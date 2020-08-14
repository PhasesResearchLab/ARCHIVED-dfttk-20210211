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
from dfttk.pythelec import thelecMDB
from dfttk.pythfind import thfindMDB
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
    gaussian = args.gaussian
    dope = args.dope
    doscar = args.doscar
    outf = args.outf
    qhamode = args.qhamode
    eqmode = args.eqmode
    elmode = args.elmode
    metatag = args.metatag
    everyT = args.everyT
    noel = args.noel
    plot = args.plot
    smooth = args.smooth
    expt = args.expt
    xlim = args.xlim
    if abs(dope)<5.e-9:
        ndosmx = 100001
        gaussian = 10000.

    #call API
    if metatag != None:
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        proc = thelecMDB(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, db_file, 
            noel=noel, metatag=metatag, qhamode=qhamode, eqmode=eqmode, elmode=elmode, everyT=everyT, 
            smooth=smooth, plot=plot)
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
    pthelec.add_argument("-gauss", "--gauss", dest="gaussian", nargs="?", type=float, default=1000.,
                      help="densing number near the Fermi energy. \n"
                           "Default: 1000")
    pthelec.add_argument("-i", "--doscar", dest="doscar", nargs="?", type=str, default="DOSCAR",
                      help="DOSCAR filename. \n"
                           "Default: DOSCAR")
    pthelec.add_argument("-o", "-outf", dest="outf", nargs="?", type=str, default="fvib_ele",
                      help="output filename for calculated thermoelectric properties. \n"
                           "Default: fvib_ele")
    pthelec.add_argument("-noel", "-noel", dest="noel", action='store_true', default=False,
                      help="do not consider the thermal electron contribution. \n"
                           "Default: False")
    pthelec.add_argument("-metatag", "-metatag", dest="metatag", nargs="?", type=str, default=None,
                      help="metatag: MongoDB metadata tag field. \n"
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
    pthelec.add_argument("-s", "-smooth", dest="smooth", action='store_true', default=False,
                      help="smooth the LTC. \n"
                           "Default: False")
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
 
    #further extension for finding phonon calculation
    run_ext_thfind(subparsers)


def run_ext_thfind(subparsers):
    #SUB-PROCESS: thfind
    pthfind = subparsers.add_parser("thfind", help="find the metadata tag that has finished.")
    pthfind.add_argument("-q", "--qhamode", dest="qhamode", nargs="?", type=str, default='phonon',
                      help="Collection. 'phonon', 'qha'.\n"
                           "Default: 'phonon'")
    pthfind.add_argument("-w", "--within", dest="within", nargs="?", type=str, default=None,
                      help="find calculations within element list\n"
                           "Default: None")
    pthfind.add_argument("-all", "--containall", dest="containall", nargs="?", type=str, default=None,
                      help="find calculations must contain all elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-any", "--containany", dest="containany", nargs="?", type=str, default=None,
                      help="find calculations contain any elements in the list\n"
                           "Default: None")
    pthfind.add_argument("-v", "--nV", dest="nV", nargs="?", type=int, default=6,
                      help="Return phonon calculations finished for number of volumes larger or equals to. \n"
                           "Default: 6")
    pthfind.add_argument("-ss", "--supercellsize", dest="supercellN", nargs="?", type=int, default=0,
                      help="only return phonon calculation with supercell size larger than. \n"
                           "Default: 0")
    pthfind.add_argument("-g", "--get", dest="get", action='store_true', default=False,
                      help="get the thermodyamic data for all found entries. \n"
                           "Default: False")
    pthfind.add_argument("-ty", "--toyphon", dest="toyphon", action='store_true', default=False,
                      help="extract the superfij.out used by Yphon for all found entries. \n"
                           "Default: False")
    pthfind.add_argument("-T0", "-t0", dest="t0", nargs="?", type=float, default=0.0,
                      help="Low temperature limit. \n"
                           "Default: 0")
    pthfind.add_argument("-T1", "-t1", dest="t1", nargs="?", type=float, default=2000,
                      help="High temperature limit. \n"
                           "Default: 1300")
    pthfind.add_argument("-dT", "-td", dest="td", nargs="?", type=float, default=10,
                      help="Temperature increment. \n"
                           "Default: 10")
    pthfind.add_argument("-xdn", "--xdn", dest="xdn", nargs="?", type=float, default=-100,
                      help="Low band energy limit. \n"
                           "Default: -100 (eV)")
    pthfind.add_argument("-xup", "--xup", dest="xup", nargs="?", type=float, default=100,
                      help="High band energy limit. \n"
                           "Default: 100")
    pthfind.add_argument("-dope", "--dope", dest="dope", nargs="?", type=float, default=0.0,
                      help="dope level (electrons). \n"
                           "Default: 0")
    pthfind.add_argument("-ne", "--ndosmx", dest="ndosmx", nargs="?", type=int, default=10001,
                      help="new DOS mesh. \n"
                           "Default: 10001")
    pthfind.add_argument("-natom", "--natom", dest="natom", nargs="?", type=int, default=1,
                      help="number of atoms in the DOSCAR. \n"
                           "Default: 1")
    pthfind.add_argument("-e", "--everyT", dest="everyT", nargs="?", type=int, default=1,
                      help="number of temperature points skipped from QHA analysis. \n"
                           "Default: 1")
    pthfind.add_argument("-gauss", "--gauss", dest="gaussian", nargs="?", type=float, default=1000.,
                      help="densing number near the Fermi energy. \n"
                           "Default: 1000")
    pthfind.add_argument("-i", "--doscar", dest="doscar", nargs="?", type=str, default="DOSCAR",
                      help="DOSCAR filename. \n"
                           "Default: DOSCAR")
    pthfind.add_argument("-o", "-outf", dest="outf", nargs="?", type=str, default="fvib_ele",
                      help="output filename for calculated thermoelectric properties. \n"
                           "Default: fvib_ele")
    pthfind.add_argument("-noel", "-noel", dest="noel", action='store_true', default=False,
                      help="do not consider the thermal electron contribution. \n"
                           "Default: False")
    pthfind.add_argument("-metatag", "-metatag", dest="metatag", nargs="?", type=str, default=None,
                      help="metatag: MongoDB metadata tag field. \n"
                           "Default: None")
    pthfind.add_argument("-qhamode", "-qhamode", dest="qhamode", nargs="?", type=str, default=None,
                      help="quasiharmonic mode: debye, phonon, or yphon. \n"
                           "Default: debye")
    pthfind.add_argument("-eq", "--eqmode", dest="eqmode", nargs="?", type=int, default=0,
                      help="Mode to calculate LTC. 0: Symmetrical Central differential;  \n"
                           "                       4: 4-parameter BM fitting.  \n"
                           "                       5: 5-parameter BM fitting.  \n"
                           "Default: 0")
    pthfind.add_argument("-el", "--elmode", dest="elmode", nargs="?", type=int, default=1,
                      help="Mode to interpolate thermal electronic contribution:"
                           "                       0: interp1d;  \n"
                           "                       1: UnivariateSpline.  \n"
                           "Default: 1")
    pthfind.add_argument("-s", "-smooth", dest="smooth", action='store_true', default=False,
                      help="smooth the LTC. \n"
                           "Default: False")
    pthfind.add_argument("-plot", "-plot", dest="plot", action='store_true', default=False,
                      help="plot the figure. \n"
                           "Default: False")
    pthfind.add_argument("-expt", "-expt", dest="expt", nargs="?", type=str, default=None,
                      help="json file path for experimental thermodynamic properties for plot. \n"
                           "Default: None")
    pthfind.add_argument("-xlim", "-xlim", dest="xlim", nargs="?", type=float, default=None,
                      help="Up temperature limit for plot. \n"
                           "Default: None")
    pthfind.set_defaults(func=ext_thfind)


def ext_thfind(args):
    """
    find the metadata tag that has finished.

    Parameters
        STR_FOLDER = args.STRUCTURE_FOLDER
            folder/file containing structures
        MATCH_PATTERN = args.MATCH_PATTERN
            Match patterns for structure file, e.g. *POSCAR
        RECURSIVE = args.RECURSIVE
            recursive or not
        WORKFLOW = args.WORKFLOW
            workflow, current only get_wf_gibbs
    """
    proc=thfindMDB(args)
    tags = proc.run_console()
    if args.get:
        for tag in tags:
            print("\nDownloading data by metadata tag:", tag, "\n")
            args.metatag = tag
            ext_thelec(args)
