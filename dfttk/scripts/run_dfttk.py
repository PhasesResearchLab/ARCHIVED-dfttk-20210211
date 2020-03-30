# -*- coding: utf-8 -*-
# The template for batch run of DFTTK
import argparse
from pymatgen import MPRester, Structure
from pymatgen.io.vasp.inputs import Potcar
from dfttk.wflows import get_wf_gibbs
from dfttk.utils import recursive_glob
from dfttk.structure_builders.parse_anrl_prototype import multi_replace
from monty.serialization import loadfn, dumpfn
import os
import sys
import shutil

def get_abspath(path):
    """
    Get the absolute path.

    Parameter
        path: str (path-like)
            An absolute or relative path. It can start from home (~)
    Return
        path: str (path-like)
            The absolute path(start from root "/").
    """
    if path.startswith("~"):
        path = path.replace("~", os.environ["HOME"], 1)
    else:
        path = os.path.abspath(path)
    return path

def creat_folders(folder):
    """
    Create folders if not exist, leave a warning if exists
    """
    if os.path.exists(folder):
        print("WARNING: " + folder + " exists!")
    else:
        os.makedirs(folder)

def get_structure_file(STR_FOLDER=".", RECURSIVE=False, MATCH_PATTERN="*"):
    """
    Get all file names in STR_FOLDER path [and its sub-folder (RECURSIVE)] with condition (MATCH_PATTEN)

    Parameter
        STR_FOLDER: str (path-like)
            The path containing files, Default: .
        RECURSIVE: bool
            Including the sub-folder (True) or not (Default: False)
        MATCH_PATTERN: str
            The match pattern for the file name
    Return
        STR_FILES: list[str (filename-like)]
            The list of all the filenames matching the conditions
    """
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
            STR_FILES = [os.path.join(STR_FOLDER, file_i) for file_i in STR_FILES]
    return STR_FILES

def get_setting_file(STR_FILENAME, NEW_SETTING="SETTINGS"):
    """
    Get the filename (without ext) of setting file
    (By default: The setting file should be SETTINGS or start with SETTINGS- or end with -SETTINGS (case insensitive))

    Parameter
        STR_FILENAME: str
            The individual tags for the setting file
        NEW_SETTING: str
            The str to replace "SETTINGS"
    Return
        SETTING_FILENAMES: list
            The list of setting files (without ext)
    """
    SETTING_FILENAMES = ["SETTINGS", "settings"
                         "SETTINGS-" + STR_FILENAME, "settings-" + STR_FILENAME,
                         STR_FILENAME + "-SETTINGS", STR_FILENAME + "-settings"]
    replace_dict = {"SETTINGS": NEW_SETTING.upper(), "settings": NEW_SETTING.lower()}
    SETTING_FILENAMES = [multi_replace(item, replace_dict) for item in SETTING_FILENAMES]
    return SETTING_FILENAMES

def run(args):
    """
    Run dfttk
    Currently, only support get_wf_gibbs

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
    STR_FOLDER = args.STRUCTURE_FOLDER  # folder/file containing structures
    MATCH_PATTERN = args.MATCH_PATTERN  # Match patterns for structure file, e.g. *POSCAR
    RECURSIVE = args.RECURSIVE          # recursive or not
    WORKFLOW = args.WORKFLOW            # workflow, current only get_wf_gibbs
    LAUNCH = args.LAUNCH               # Launch to lpad or not
    MAX_JOB = args.MAX_JOB              # Max job to submit
    SETTINGS = args.SETTINGS            # Settings file    
    WRITE_OUT_WF = args.WRITE_OUT_WF    # Write out wf file or not

    ## Load global settings
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
                os.system("qlaunch rapidfire -m " + str(MAX_JOB))

def config(args):
    """
    Config dfttk
    It can be used to config atomate and pymatgen
    """
    import dfttk.scripts.config_dfttk as dfttkconfig
    TEST_CONFIG = args.TEST_CONFIG

    ALL = args.ALL
    PATH_TO_STORE_CONFIG = get_abspath(args.PATH_TO_STORE_CONFIG)

    ATOMATE = args.ATOMATE
    VASP_CMD_FLAG = args.VASP_CMD_FLAG
    CONFIG_FOLDER = args.CONFIG_FOLDER
    QUEUE_SCRIPT = args.QUEUE_SCRIPT
    QUEUE_TYPE = args.QUEUE_TYPE

    PYMATGEN = args.PYMATGEN
    VASP_PSP_DIR = args.VASP_PSP_DIR
    MAPI_KEY = args.MAPI_KEY
    DEFAULT_FUNCTIONAL = args.DEFAULT_FUNCTIONAL
    
    if ALL:
        ATOMATE = True
        PYMATGEN = True

    TEST_CONFIG_MAP = {"all": [True, True], "atomate": [True, False], 
                       "pymatgen": [False, True], "none": [False, False]}
    [test_atomate, test_pymagen] = TEST_CONFIG_MAP[TEST_CONFIG.lower()]
    if test_atomate or test_pymagen:
        # -t parameter exists in the input
        dfttkconfig.test_config(test_pymagen=True, test_atomate=True)
        exit()

    if PATH_TO_STORE_CONFIG is None:
        PATH_TO_STORE_CONFIG = dfttkconfig.default_path()
    PATH_TO_STORE_CONFIG = get_abspath(PATH_TO_STORE_CONFIG)

    if ATOMATE:
        dfttkconfig.config_atomate(path_to_store_config=PATH_TO_STORE_CONFIG, config_folder=CONFIG_FOLDER, 
            queue_script=QUEUE_SCRIPT, queue_type=QUEUE_TYPE, vasp_cmd_flag=VASP_CMD_FLAG)

    if PYMATGEN:
        dfttkconfig.config_pymatgen(psp_dir=VASP_PSP_DIR, def_fun=DEFAULT_FUNCTIONAL, 
            mapi=MAPI_KEY, path_to_store_psp=os.path.join(PATH_TO_STORE_CONFIG, "vasp_psp"))

def run_dfttk():
    """
    dfttk command
    """
    from dfttk._version import get_versions
    print("DFTTK version: " + get_versions()["version"])
    print("Copyright \u00a9 Phases Research Lab (https://www.phaseslab.com/)\n")

    parser = argparse.ArgumentParser(description='Run DFTTK jobs.')
    
    subparsers = parser.add_subparsers()

    #SUB-PROCESS: run
    prun = subparsers.add_parser("run", help="Run dfttk.")
    prun.add_argument("-f", "--structure_folder", dest="STRUCTURE_FOLDER", type=str, default=".",
                      help="The folder/file containing the structure,\n"
                           "Default: '.' ")
    prun.add_argument("-mh", "--match_pattern", dest="MATCH_PATTERN", type=str, default="*",
                      help="The match pattern for structure file, and it should be place in quotes."
                           " e.g. '*POSCAR*'. Default: * -- everything except SETTING files, ref. -s")
    prun.add_argument("-s", "--setting", dest="SETTINGS", type=str, default="SETTINGS",
                      help="Specify the name of SETTINGS files (yaml or json file)\n"
                           "Default: SETTINGS (case insensitive and without ext) \n"
                           "The following filename will be treat as SETTINGS file \n"
                           "\t SETTINGS (global settings in the folder)\n"
                           "\t Start with SETTINGS- (individual settings for struct)\n"
                           "\t End with -SETTINGS (individual settings)")
    prun.add_argument("-r", "--recursive", dest="RECURSIVE", action="store_true", 
                      help="Recursive the path.")
    prun.add_argument("-wf", "--workflow", dest="WORKFLOW", type=str, default="get_wf_gibbs",
                      help="""Specify the workflow to run.\n
                           Default: get_wf_gibbs \n
                           (NOTE: currently, only get_wf_gibbs is supported.)""")
    prun.add_argument("-l", "--launch", dest="LAUNCH", action="store_true",
                      help="Launch the wf to launchpad")
    prun.add_argument("-m", "--max_job", dest="MAX_JOB", nargs="?", type=int, default=0,
                      help="Run the job, only works when -l is specified.\n"
                           "Default: 0 (Not submit to queue) \n"
                           "1: qlaunch singleshot (single job) \n"
                           "N(N>1): qlaunch rapidfire -m N")
    prun.add_argument("-o", "--write_out_wf", dest="WRITE_OUT_WF", action="store_true",
                      help="Write out the workflow")
    prun.set_defaults(func=run)

    #SUB-PROCESS: config
    pconfig = subparsers.add_parser("config", help="Config dfttk.")
    pconfig.add_argument("-all", "--all", dest="ALL", action="store_true",
                         help="Configure atomate and pymatgen.")
    pconfig.add_argument("-p", "--prefix", dest="PATH_TO_STORE_CONFIG", default=".",
                         help="The folder to store the config files.\n"
                              "Default: . (current folder)")
    pconfig.add_argument("-a", "--atomate", dest="ATOMATE", action="store_true",
                         help="Configure atomate.")
    pconfig.add_argument("-c", "--config_folder", dest="CONFIG_FOLDER", default=".",
                         help="The folder containing config files, "
                              "at least contain db.json and my_launchpad.yaml. Default: '.'")
    pconfig.add_argument("-q", "--queue_script", dest="QUEUE_SCRIPT", default="vaspjob.pbs",
                         help="The filename of the script for sumitting vasp job. "
                              "It will search in current folder and sub-folders. Default: vaspjob.pbs")
    pconfig.add_argument("-qt", "--queue_type", dest="QUEUE_TYPE", type=str, default="pbs",
                         help="The type of queue system. Note, only pbs is supported now. Default: pbs")
    pconfig.add_argument("-v", "--vasp_cmd_flg", dest="VASP_CMD_FLAG", type=str, default="vasp_std",
                         help="The flag to distinguish vasp_cmd to othe commands in queue_script. Default: vasp_std")
    pconfig.add_argument("-mp", "--pymatgen", dest="PYMATGEN", action="store_true",
                         help="Configure pymatgen.")
    pconfig.add_argument("-psp", "--vasp_psp_dir", dest="VASP_PSP_DIR", type=str, default="psp",
                         help="The path of pseudopotentials. Default: psp")
    pconfig.add_argument("-mapi", "--mapi_key", dest="MAPI_KEY", type=str,
                         help="The API key of Materials Projects")
    pconfig.add_argument("-df", "--default_functional", dest="DEFAULT_FUNCTIONAL", type=str, default="PBE",
                         choices=sorted(Potcar.FUNCTIONAL_CHOICES),
                         help="The default functional. Default: PBE")
    pconfig.add_argument("-t", "--test_config", dest="TEST_CONFIG", nargs="?", const="all", default="none"
                         choices=["all", "pymatgen", "atomate"],
                         help="Test for configurations. Note: currently only support for pymatgen.")
    pconfig.set_defaults(func=config)


    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == '__main__':
    run_dfttk()
