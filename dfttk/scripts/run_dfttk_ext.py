# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen import MPRester, Structure
from pymatgen.io.vasp.inputs import Potcar
#from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import dfttk.pythelec as pythelec
from dfttk.pythelec import class_thelecMDB
import warnings
import copy
import os
import sys
import shutil

def ext_thelec(args):
    print ("Postprocess for thermodynamic properties, Seebeck, Lorenz number etc. Yi Wang\n")
    """
    Postprocess for thermodynamic properties, Seebeck, Lorenz number etc

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
        LAUNCH = args.LAUNCH
            Launch to lpad or not
        MAX_JOB = args.MAX_JOB
            Max job to submit
        SETTINGS = args.SETTINGS
            Settings file
        WRITE_OUT_WF = args.WRITE_OUT_WF
            Write out wf file or not
    """
    t0 = args.t0
    t1 = args.t1
    td = args.td
    xdn = args.xdn
    xup = args.xup
    ndosmx = args.ndosmx
    natom = args.natom
    gauss = args.gauss
    dope = args.dope
    doscar = args.doscar
    outf = args.outf
    qhamode = args.qhamode
    eqmode = args.eqmode
    elmode = args.elmode
    metadata = args.metadata
    everyT = args.everyT
    plot = args.plot
    expt = args.expt
    xlim = args.xlim
    if abs(dope)<5.e-9:
        ndosmx = 100001
        gaussian = 10000.

    #call API
    if metadata != None:
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        if 1==0:
            volumes, energies, thermofile = pythelec.thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, 
                outf, db_file, 
                metadata=metadata, qhamode=qhamode, eqmode=eqmode, elmode=elmode, everyT=everyT, plot=plot)
        else:
            proc = class_thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom,           
                outf, db_file, 
                metadata=metadata, qhamode=qhamode, eqmode=eqmode, elmode=elmode, everyT=everyT, plot=plot)
            volumes, energies, thermofile = proc.run_console()

        print("\nFull thermodynamic properties have outputed into:", thermofile) 
        if plot: 
            from dfttk.analysis.ywplot import plotAPI
            plotAPI(thermofile, volumes, energies, expt=expt, xlim=xlim)
    else:
        pythelec.thelecAPI(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, doscar)

def run_ext_thelec(subparsers):
    # begin process by Yi Wang, July 23, 2020
    #SUB-PROCESS: thelec
    pthelec = subparsers.add_parser("thelec", help="Calculate Seebeck and Lorenz number.")
    pthelec.add_argument("-T0", "-t0", dest="t0", nargs="?", type=float, default=0.0,
                      help="Low temperature limit. \n"
                           "Default: 0")
    pthelec.add_argument("-T1", "-t1", dest="t1", nargs="?", type=float, default=2000,
                      help="High temperature limit. \n"
                           "Default: 1300")
    pthelec.add_argument("-dT", "-td", dest="td", nargs="?", type=float, default=10,
                      help="Temperature increment. \n"
                           "Default: 10")
    pthelec.add_argument("-xdn", "--xdn", dest="xdn", nargs="?", type=float, default=-100,
                      help="Low band energy limit. \n"
                           "Default: -100 (eV)")
    pthelec.add_argument("-xup", "--xup", dest="xup", nargs="?", type=float, default=100,
                      help="High band energy limit. \n"
                           "Default: 100")
    pthelec.add_argument("-dope", "--dope", dest="dope", nargs="?", type=float, default=0.0,
                      help="dope level (electrons). \n"
                           "Default: 0")
    pthelec.add_argument("-ne", "--ndosmx", dest="ndosmx", nargs="?", type=int, default=10001,
                      help="new DOS mesh. \n"
                           "Default: 10001")
    pthelec.add_argument("-natom", "--natom", dest="natom", nargs="?", type=int, default=1,
                      help="number of atoms in the DOSCAR. \n"
                           "Default: 1")
    pthelec.add_argument("-e", "--everyT", dest="everyT", nargs="?", type=int, default=1,
                      help="number of temperature points skipped from QHA analysis. \n"
                           "Default: 1")
    pthelec.add_argument("-gauss", "--gauss", dest="gauss", nargs="?", type=float, default=1000.,
                      help="densing number near the Fermi energy. \n"
                           "Default: 1000")
    pthelec.add_argument("-i", "--doscar", dest="doscar", nargs="?", type=str, default="DOSCAR",
                      help="DOSCAR filename. \n"
                           "Default: DOSCAR")
    pthelec.add_argument("-o", "-outf", dest="outf", nargs="?", type=str, default="fvib_ele",
                      help="output filename for calculated thermoelectric properties. \n"
                           "Default: fvib_ele")
    pthelec.add_argument("-metadata", "-metadata", dest="metadata", nargs="?", type=str, default=None,
                      help="metadata: MongoDB data id. \n"
                           "Default: None")
    pthelec.add_argument("-qhamode", "-qhamode", dest="qhamode", nargs="?", type=str, default=None,
                      help="quasiharmonic mode: debye, phonon, or yphon. \n"
                           "Default: debye")
    pthelec.add_argument("-eq", "--eqmode", dest="eqmode", nargs="?", type=int, default=0,
                      help="Mode to calculate LTC. 0: Symmetrical Central differential;  \n"
                           "                       4: 4-parameter BM fitting.  \n"
                           "                       5: 5-parameter BM fitting.  \n"
                           "Default: 0")
    pthelec.add_argument("-el", "--elmode", dest="elmode", nargs="?", type=int, default=1,
                      help="Mode to interpolate thermal electronic contribution:"
                           "                       0: interp1d;  \n"
                           "                       1: UnivariateSpline.  \n"
                           "Default: 1")
    pthelec.add_argument("-plot", "-plot", dest="plot", action='store_true', default=False,
                      help="plot the figure. \n"
                           "Default: False")
    pthelec.add_argument("-expt", "-expt", dest="expt", nargs="?", type=str, default=None,
                      help="json file path for experimental thermodynamic properties for plot. \n"
                           "Default: None")
    pthelec.add_argument("-xlim", "-xlim", dest="xlim", nargs="?", type=float, default=None,
                      help="Up temperature limit for plot. \n"
                           "Default: None")
    pthelec.set_defaults(func=ext_thelec)
    # end process by Yi Wang, July 23, 2020
