# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen import MPRester, Structure
from pymatgen.io.vasp.inputs import Potcar
from dfttk.wflows import get_wf_gibbs, get_wf_EV_bjb, get_wf_gibbs_robust
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import dfttk.pythelec as pythelec
import warnings
import copy
import os
import sys
import shutil

def thelec(args):
    print ("\n  Calculate Seebeck and Lorenz number. Yi Wang\n")
    """
    calculate Seebeck and Lorenz number

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
    metadata = args.metadata
    everyT = args.everyT
    if abs(dope)<5.e-9:
        ndosmx = 100001
        gaussian = 10000.

    #call API
    if metadata != None:
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
        pythelec.thelecMDB(xdn, xup, dope, ndosmx, gaussian, natom, outf, db_file, 
            metadata=metadata, qhamode=qhamode, everyT=everyT)
    else:
        pythelec.thelecAPI(t0, t1, td, xdn, xup, dope, ndosmx, gaussian, natom, outf, doscar)

def ywrun(subparsers):
    # begin process by Yi Wang, July 23, 2020
    #SUB-PROCESS: thelec
    pthelec = subparsers.add_parser("thelec", help="Calculate Seebeck and Lorenz number.")
    pthelec.add_argument("-T0", "--t0", dest="t0", nargs="?", type=float, default=0.0,
                      help="Low temperature limit. \n"
                           "Default: 0")
    pthelec.add_argument("-T1", "--t1", dest="t1", nargs="?", type=float, default=1300,
                      help="High temperature limit. \n"
                           "Default: 1300")
    pthelec.add_argument("-dT", "--td", dest="td", nargs="?", type=float, default=10,
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
    pthelec.add_argument("-qhamode", "-qhamode", dest="qhamode", nargs="?", type=str, default="debye",
                      help="quasiharmonic mode: debye, phonon, or yphon. \n"
                           "Default: debye")
    pthelec.set_defaults(func=thelec)
    # end process by Yi Wang, July 23, 2020
