# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen import MPRester, Structure
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import os
import sys

def global_settings():
    ################ PARAMETERS FOR WF #############################
    #str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
    #  If None, it will use the configuration in fireworks 
    db_file=None
    #list, the MAGMOM of the structure, e.g. [4.0, 4.0, -4.0, -4.0]
    magmom = None
    #int, the number of initial deformations, e.g. 7
    num_deformations = 7
    #list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.05, 0.1) means (0.95, 1.1) in volume
    deformation_fraction = (-0.1, 0.1)
    #float, minimum ratio of Volumes spacing, e.g. 0.03
    volume_spacing_min = 0.03
    #bool, run phonon(True) or not(False)
    phonon=False
    #list(3x3), the supercell matrix for phonon, e.g. [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
    phonon_supercell_matrix=None
    #float, the mimimum of temperature in QHA process, e.g. 5
    t_min=5
    #float, the maximum of temperature in QHA process, e.g. 2000
    t_max=2000
    #float, the step of temperature in QHA process, e.g. 5
    t_step=5
    #float, acceptable value for average RMS, recommend >= 0.005
    tolerance = 0.01
    #str, the vasp command, if None then find in the FWorker configuration
    vasp_cmd=None
    #dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
    metadata=None
    #float, the tolerannce for symmetry, e.g. 0.05
    symmetry_tolerance = 0.05
    #bool, set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    passinitrun=False
    #bool, Whether run isif=2 calculation before isif=4 running
    run_isif2=False
    #bool, Whether pass isif=4 calculation.
    pass_isif4=False
    #Set the path already exists for new static calculations; if set as '', will try to get the path from db_file
    relax_path=''
    #dict, dict of class ModifyIncar with keywords in Workflow name. e.g.
    """
    modify_incar_params = { 'Full relax': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                            'PreStatic': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                            'PS2': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}}, 
                            'static': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
    """
    modify_incar_params={}
    #dict, dict of class ModifyKpoints with keywords in Workflow name, similar with modify_incar_params
    modify_kpoints_params={}
    #bool, print(True) or not(False) some informations, used for debug
    verbose=False
    GLOBAL_SETTING = {"db_file": db_file, "magmom": magmom, "num_deformations": num_deformations,
                  "deformation_fraction": deformation_fraction, "volume_spacing_min": volume_spacing_min,
                  "phonon": phonon, "phonon_supercell_matrix": phonon_supercell_matrix,
                  "t_min": t_min, "t_max": t_max, "t_step": t_step, "tolerance": tolerance,
                  "vasp_cmd": vasp_cmd, "metadata": metadata, "symmetry_tolerance": symmetry_tolerance,
                  "passinitrun": passinitrun, "run_isif2": run_isif2, "pass_isif4": pass_isif4,
                  "relax_path": relax_path, "modify_incar_params": modify_incar_params, 
                  "modify_kpoints_params": modify_kpoints_params, "verbose": verbose}
    return GLOBAL_SETTING

def get_abspath(path):
    # Get the abs path
    if path.startswith("~"):
        path = path.replace("~", os.environ["HOME"], 1)
    else:
        path = os.path.abspath(path)
    return path

def get_structure_file(STR_FOLDER=".", RECURSIVE=False, MATCH_PATTERN="*"):
    ## Get the file name of structure file
    STR_FOLDER = get_abspath(STR_FOLDER)
    if not os.path.exists(STR_FOLDER):
        raise IOError("ERROR: The specified folder/file does not exist.")
    if os.path.isfile(STR_FOLDER):
        STR_FILES = [STR_FOLDER]
    else:
        if RECURSIVE:
            STR_FILES = recursive_glob(STR_FOLDER, MATCH_PATTERN)
        else:
            STR_FILES = os.listdir(STR_FOLDER)
            for file_i in STR_FILES:
                if os.path.isdir(file_i):
                    STR_FILES.remove(file_i)
    return STR_FILES

def get_setting_file(STR_FILENAME, NEW_SETTING="SETTINGS"):
    SETTING_FILENAMES = ["SETTINGS", "settings"
                         "SETTINGS-" + STR_FILENAME, "settings-" + STR_FILENAME,
                         STR_FILENAME + "-SETTINGS", STR_FILENAME + "-settings"]
    replace_dict = {"SETTINGS": NEW_SETTING.upper(), "settings": NEW_SETTING.lower()}
    SETTING_FILENAMES = [multi_replace(item, replace_dict) for item in SETTING_FILENAMES]
    return SETTING_FILENAMES

def dfttk_run(args):
    STR_FOLDER = args.STRUCTURE_FOLDER  # folder/file containing structures
    MATCH_PATTERN = args.MATCH_PATTERN  # Match patterns for structure file, e.g. *POSCAR
    RECURSIVE = args.RECURSIVE          # recursive or not
    WORKFLOW = args.WORKFLOW            # workflow, current only get_wf_gibbs
    LAUNCH = args.LAUNCH               # Launch to lpad or not
    MAX_JOB = args.MAX_JOB              # Max job to submit
    SETTINGS = args.SETTINGS            # Settings file    
    WRITE_OUT_WF = args.WRITE_OUT_WF    # Write out wf file or not

    ## Load global settings
    GLOBAL_SETTING = global_settings()

    ## Get the file names of files
    STR_FILES = get_structure_file(STR_FOLDER=STR_FOLDER, RECURSIVE=RECURSIVE, MATCH_PATTERN=MATCH_PATTERN)

    ## Initial wfs and metadatas
    wfs = []
    metadatas = {}

    ## generat the wf
    for STR_FILE in STR_FILES:
        (STR_PATH, STR_FILENAME_WITH_EXT) = os.path.split(STR_FILE)
        (STR_FILENAME, STR_EXT) = os.path.splitext(STR_FILENAME_WITH_EXT)
        str_filename = STR_FILENAME.lower()
        if (str_filename.endswith("-" + SETTINGS.lower()) or 
           str_filename.startswith( SETTINGS.lower() + "-") or 
           (str_filename == SETTINGS.lower())):
            print(STR_FILE + " is a setting file, not structure file, and skipped when reading the structure.")
        elif STR_FILE == os.path.abspath(__file__):
            #This is current file
            pass
        else:
            try:
                structure = Structure.from_file(STR_FILE)

                locals().update(GLOBAL_SETTING)
                SETTING_FILENAMES = get_setting_file(STR_FILENAME, SETTINGS)
                for SETTING_FILENAME in SETTING_FILENAMES:
                    SETTING_FULL_FILENAME = os.path.join(STR_PATH, SETTING_FILENAME + STR_EXT)
                    if os.path.exists(SETTING_FULL_FILENAME):
                        locals().update(loadfn(SETTING_FULL_FILENAME))
                if magmom:
                    structure.add_site_property('magmom', magmom)
                if not db_file:
                    from fireworks.fw_config import config_to_dict
                    db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]    

                if WORKFLOW == "get_wf_gibbs":
                    #Currently, only this workflow is supported
                    wf = get_wf_gibbs(structure, num_deformations=num_deformations, deformation_fraction=deformation_fraction, 
                                phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix,  t_min=t_min, t_max=t_max, 
                                t_step=t_step, tolerance=tolerance, volume_spacing_min=volume_spacing_min,vasp_cmd=vasp_cmd, 
                                db_file=db_file, metadata=metadata, name='EV_QHA', symmetry_tolerance=symmetry_tolerance, 
                                run_isif2=run_isif2, pass_isif4=pass_isif4, passinitrun=passinitrun, relax_path=relax_path, 
                                modify_incar_params=modify_incar_params, modify_kpoints_params=modify_kpoints_params, 
                                verbose=verbose)
                metadatas[STR_FILE] = wf.as_dict()["metadata"]
                wfs.append(wf)

                if WRITE_OUT_WF:
                    dfttk_wf_filename = os.path.join(STR_PATH, "dfttk_wf-" + STR_FILENAME + ".yaml")
                    wf.to_file(dfttk_wf_filename)
            except Exception as e:
                print("The name or the contant of " + STR_FILE + " is not supported by dfttk, and skipped.")
                print("Ref. https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure.from_file")
    #Write Out the metadata for POST purpose
    dumpfn(metadatas, "METADATAS.yaml")

    if LAUNCH:
        from fireworks import LaunchPad
        lpad = LaunchPad.auto_load()

        for wflow in wfs:
            lpad.add_wf(wflow)
        if MAX_JOB:
            # Not False or Empty
            if MAX_JOB == 1:
                os.system("qlaunch singleshot")
            else:
                os.system("qlaunch -m 2")


if __name__ == '__main__':
    from dfttk._version import get_versions
    print("DFTTK version: " + get_versions()["version"])
    print("Copyright \u00a9 Phases Research Lab (https://www.phaseslab.com/)")

    parser = argparse.ArgumentParser(description='Run DFTTK jobs.')

    subparsers = parser.add_subparsers()
    prun = subparsers.add_parser("run", help="Run dfttk.")
    prun.add_argument("-f", "--structure_folder", dest="STRUCTURE_FOLDER", type=str, default=".",
        help="The folder/file containing the structure, default is current folder('.')")
    prun.add_argument("-mh", "--match_pattern", dest="MATCH_PATTERN", type=str, default="*",
        help="The match pattern for the file of structure, e.g. *POSCAR* for any file contain 'POSCAR',"
             " default is everything (*) except those yaml/json file:"
             " settings.{yaml/json} or start with settings- or end with -settings.{yaml/json} or specified by -s parameter")
    prun.add_argument("-s", "--setting", dest="SETTINGS", type=str, default="SETTINGS",
        help="Specify the name of settings file (yaml or json file)"
             "The file name(without ext) will start with SETTINGS- or end with -SETTINGS or SETTINGS"
             "SETTINGS.yaml is global settings in the folder, and *-SETTINGS.yaml or SETTINGS-*.yaml is individual settings"
             "e.g. -s set, then the file SET.yaml, *-set.json, set-*.yaml will all be taken as setting file")
    prun.add_argument("-r", "--recursive", dest="RECURSIVE", action="store_true", 
        help="Recursive the path or not")
    prun.add_argument("-wf", "--workflow", dest="WORKFLOW", type=str, default="get_wf_gibbs",
        help="Specify the workflow to run. (NOTE: currently, only get_wf_gibbs is supported.)")
    prun.add_argument("-l", "--launch", dest="LAUNCH", action="store_true",
        help="Launch the wf to launchpad")
    prun.add_argument("-m", "--max_job", dest="MAX_JOB", nargs="?", default=False,
        help="Run the job, only works with -l parameter."
             "False(default): Not submit to queue,"
             "1: qlaunch singleshot"
             "N(N>1): qlaunch -m N")
    prun.add_argument("-o", "--write_out_wf", dest="WRITE_OUT_WF", action="store_true",
        help="Write out the workflow")
    prun.set_defaults(func=dfttk_run)

    args = parser.parse_args()
    #print(args)
    #args = prun.parse_args()
    #print(args)
    #print(args.LAUNCH)

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)
