# -*- coding: utf-8 -*-
#
#Configure script for dfttk

from dfttk.scripts.run_dfttk import get_abspath, creat_folders
from pymatgen.io.vasp.inputs import PotcarSingle, Potcar
from pymatgen import MPRester
from dfttk.utils import recursive_glob
from monty.serialization import loadfn, dumpfn
import getpass
import shutil
import math
import os
import re

def replace_file(filename, old_str, new_str):
    """
    Replace the old_str with new_str in filename

    Parameter
        filename: str (filename)
            The file name of the file to be replace
        old_str: str
            The string need to be replaced
        new_str: str
            The new string
    Return
        None, but the file (name with filename) is updated
    """
    with open(filename, "r") as fid:
        lines = fid.readlines()
    with open(filename, "w+") as fid:
        for line in lines:
            a = re.sub(old_str, new_str, line)
            fid.writelines(a)

def default_path():
    #Get the install path of dfttk
    import dfttk
    return os.path.dirname(dfttk.__file__)

def add_path_var(**kwarg):
    """
    Add path var to ~/.bashrc
    
    If the key of kwarg exists in ~/.bashrc, it will update it, 
    if not exists, it will add "export key=kwarg[key]"
    Note: it will backup the original ~/.bashrc to ~/.bashrc.dfttk.bak
    """
    homepath = os.environ["HOME"]
    bashrc = homepath + "/.bashrc"
    shutil.copyfile(bashrc, bashrc + ".dfttk.bak")  #backup the .bashrc file
    fid = open(bashrc, "r+")
    lines = fid.readlines()
    fid.seek(0,0)
    for line in lines:
        if line.startswith("export"):
            path_var = line.strip().split()[1].split("=")[0]
            if path_var in kwarg:
                path_val = kwarg.pop(path_var)
                line = "export " + path_var + "=" + path_val + "\n"
        fid.write(line)
    for key in kwarg:
        line = "export " + key + "=" + kwarg[key] + "\n"
        fid.write(line)
    fid.close()                        

def get_shortest_path(path_list):
    """
    Return the shortest path (the uppest path)
    (Used for find the uppest path for the files with the same filename)

    Parameter
        path_list: list[path-like str], e.g. ["./db.json", "./test/db.json", "./config/db.json"]
    Return
        return the shortest path, e.g. "./db.json"
    """
    path_len = [len(path_i) for path_i in path_list]
    return path_list[path_len.index(min(path_len))]

def get_config_file(config_folder=".", queue_script="vaspjob.pbs"):
    """
    Get the file name of the config file in config_folder and its sub-folders

    Parameter
        config_folder: str (path like)
            The folder containing the config file or containing config file in its sub-folder
        queue_script: str (filename like)
            The filename of the queue script
    Return
        config_file: dict
            The dict of the filename(key) and corresponding path(value).
            If the config file is not in config_folder, store None.
            If more than one file (the same filename), it will store the uppest file [Ref. get_shortest_path]
    """
    config_file = {}
    required_file = ["db.json", "my_launchpad.yaml"]
    option_file = ["FW_config.yaml", "my_fworker.yaml", "my_qadapter.yaml", queue_script]
    files = required_file + option_file
    for file in files:
        file_list = recursive_glob(config_folder, file)
        if len(file_list) == 0:
            if file in required_file:
                raise IOError("ERROR: " + file + " file is required for configuration of dfttk.")
            else:
                print("WARNING: " + file + " file does not exist, the default setting will be used")
                config_file[file] = None
        else:
            config_file[file] = get_abspath(get_shortest_path(file_list))
    return config_file

def parase_pbs_script(filename = "vaspjob.pbs", vasp_cmd_flag="vasp_std"):
    """
    Parse pbs script for config usage

    Parameter
        filename: str (filename-like)
            The filename of the pbs script
        vasp_cmd_flag: str
            The flag to distinguish vasp_cmd with other commands
    Return
        param_dict: dict
            The dict of parameters. It support vasp_cmd, pre_rocket, post_rocket
    """
    s = {"-q": "queue", "-A": "account", "-N": "job_name", "-V": "env",
         "-G": "group_name"}
    param_dict = {"post_rocket": [], "pre_rocket": []}
    flag_post = False
    with open(filename, "r") as fid:
        for eachline in fid:
            eachline = eachline.strip()
            if eachline.startswith("#PBS"):
                line_list = re.split("\s+", eachline)
                if line_list[1] == "-l":
                    if line_list[2].startswith("walltime"):
                        # walltime
                        param_dict["walltime"] = line_list[2].split("=")[1]
                    else:
                        for item in line_list[2].split(":"):
                            key = item.split("=")[0]
                            # nodes, ppn, pmem
                            value = item.split("=")[1]
                            param_dict[key] = value
                else:
                    if line_list[1] in s:
                        param_dict[s[line_list[1]]] = line_list[2]
            elif vasp_cmd_flag in eachline:
                param_dict["vasp_cmd"] = eachline
                flag_post = True
            elif eachline.startswith("cd $") or eachline.startswith("#") or (not eachline):
                #The cd $PBS_O_WORKDIR, or comments(#) or empty
                pass
            else:
                if flag_post:
                    param_dict["post_rocket"].append(eachline)
                else:
                    param_dict["pre_rocket"].append(eachline)
    for param in param_dict:
        try:
            len_param = len(param_dict[param])
            if len_param == 0:
                param_dict[param] = None
            elif len_param == 1:
                param_dict[param] = param_dict[param][0]
        except Exception as e:
            pass
    return param_dict

def parse_queue_script(template="vaspjob.pbs", type="pbs", vasp_cmd_flag="vasp_std"):
    """
    Parse the queue script. Currently only pbs is supported

    Parameter
        template: str (filename-like)
            The filename of the queue script. Default: vaspjob.pbs
        type: str
            The type of queue system. Default: pbs
        vasp_cmd_flag: str
            The flag to distinguish vasp_cmd to other commands in the queue script. Default: vasp_std
    Return
    """
    param_dict = {}
    if type == "pbs":
        param_dict = parase_pbs_script(filename=template, vasp_cmd_flag=vasp_cmd_flag)
    else:
        raise ValueError("Only PBS is supported now. Other system will coming soon...")
    return param_dict

def parse_psp_name(psp_name):
    """
    Parse the name of vasp's psp

    Parameter
        psp_name: str
            The name of vasp's psp, e.g. GGA, LDA, potpaw_LDA
    Return
        psp_name_norm: str
            The normalized psp name
    """
    psp_name = psp_name.upper()
    psp_name_list = re.split(r'\.|\_|\-|\=|\+|\*|\s',psp_name)
    flag_us = False
    for psp_name_i in psp_name_list:
        if "US" in psp_name_i:
            flag_us = True
            break
    if "LDA" in psp_name_list:
        if "52" in psp_name_list:
            psp_name_norm = "POT_LDA_PAW_52"
        elif "54" in psp_name_list:
            psp_name_norm = "POT_LDA_PAW_54"
        elif flag_us:
            psp_name_norm = "POT_LDA_US"
        else:
            psp_name_norm = "POT_LDA_PAW"
    elif "PBE" in psp_name_list:
        if "52" in psp_name_list:
            psp_name_norm = "POT_GGA_PAW_PBE_52"
        elif "54" in psp_name_list:
            psp_name_norm = "POT_GGA_PAW_PBE_54"
        else:
            psp_name_norm = "POT_GGA_PAW_PBE"
    elif "GGA" in psp_name_list:
        if flag_us:
            psp_name_norm = "POT_GGA_US_PW91"
        else:
            psp_name_norm = "POT_GGA_PAW_PW91"
    else:
        print("WARNING: This is not a proper name of vasp's pseudopotential. It will return None.")
        psp_name_norm = None
    return psp_name_norm


def handle_potcar_gz(psp_dir="psp", path_to_store_psp="psp_pymatgen", aci=True):
    """
    Compress and move the pseudopotential to a specified path(path_to_store_psp)
    (The compress is done by running "pmg config -p psp_dir path_to_store_psp" command)

    Parameter
        psp_dir: str (path-like)
            The origial path containing psp. Both original and uncompressed are ok.
            The name of the compressed file or the sub-folder containing psps must be in the following list
            ["potpaw_PBE", "POT_GGA_PAW_PBE","potpaw_PBE_52", "POT_GGA_PAW_PBE_52","potpaw_PBE_54", "POT_GGA_PAW_PBE_54",
             "potpaw_PBE.52", "POT_GGA_PAW_PBE_52","potpaw_PBE.54", "POT_GGA_PAW_PBE_54","potpaw_LDA", "POT_LDA_PAW",
             "potpaw_LDA.52", "POT_LDA_PAW_52","potpaw_LDA.54", "POT_LDA_PAW_54","potpaw_LDA_52", "POT_LDA_PAW_52",
             "potpaw_LDA_54", "POT_LDA_PAW_54","potUSPP_LDA", "POT_LDA_US","potpaw_GGA", "POT_GGA_PAW_PW91",
             "potUSPP_GGA", "POT_GGA_US_PW91"]
            For more details, Ref:https://pymatgen.org/installation.html#potcar-setup
            The example of the structure of the psp_dir:
            e.g. psp_dir
                 ├── potpaw_LDA.54.tar.gz
                 └── potpaw_PBE.54.tar.gz
              or: psp_dir
                  ├── potpaw_LDA_54
                  │   ├── Ac
                  │   ├── Ag
                  │   └── ...
                  ├── potpaw_PBE_54
                  │   ├── Ac
                  │   ├── Ag
                  │   └── ...
                  └── ...
        path_to_store_psp: str (path-like)
            The destination to store the compressed psp. Default: psp_pymatgen
    Return
        None
    """
    #aci* For ACI at PSU only
    aci_pp_path = "/opt/aci/sw/vasp/5.4.1.05Feb16_intel-16.0.3_impi-5.1.3/pp"
    aci_name_map = {"USPP_GAA": "POT_GGA_US_PW91"}  #A typo in ACI cluster
    # file_str is not abspath, is relative path
    file_str = os.listdir(psp_dir)
    psp_uncompress = os.path.join(psp_dir, "psp_uncompress")
    creat_folders(psp_uncompress)
    if aci:
        #For ACI at PSU only
        if not os.path.exists(aci_pp_path):
            raise Exception("ERROR: The -aci parameters only work for ACI cluster at PSU.")
        pp_paths = os.listdir(aci_pp_path)
        for pp_path in pp_paths:
            dst_path_name = parse_psp_name(pp_path)
            if not dst_path_name:
                if pp_path in aci_name_map:
                    dst_path_name = aci_name_map[pp_path]
            shutil.copytree(os.path.join(aci_pp_path, pp_path), os.path.join(psp_uncompress, dst_path_name))
    for file_i in file_str:
        dst_path_name = parse_psp_name(file_i)
        if not dst_path_name:
            continue
        psp_old = os.path.join(psp_dir, file_i)
        if os.path.isdir(psp_old):
            psp_new = os.path.join(psp_uncompress, dst_path_name)
            if os.path.exists(psp_new):
                print("WARNING: Potential exists, and current potential will over write it.")
                shutil.rmtree(psp_new)
            shutil.copytree(psp_old, psp_new)
        elif os.path.isfile(psp_old):
            psp_new = os.path.join(psp_uncompress, dst_path_name)
            creat_folders(psp_new)
            if file_i.endswith(".tar.gz") or file_i.endswith(".tgz"):
                os.system("tar -zxvf " + psp_old + " -C " + psp_new)
            elif file_i.endswith(".tar"):
                os.system("tar -xvf " + psp_old + " -C " + psp_new)
    # config the POTCAR
    os.system("pmg config -p " + psp_uncompress + " " + path_to_store_psp)
    # Remove the uncompress folder
    try:
        shutil.rmtree(psp_uncompress)
    except:
        os.system("chmod +w " + os.path.join(psp_uncompress, "*/*"))
        shutil.rmtree(psp_uncompress)

def config_pymatgen(psp_dir="psp", def_fun="PBE", mapi=None, path_to_store_psp="psp_pymatgen", aci=False):
    """
    Config pymatgen. 
    If the key is exists in ~/.pmgrc.yaml and not empty, skip

    Parameter
        psp_dir: str (path-like)
            Ref: handle_potcar_gz
        def_fun: str
            The default functional. Default: PBE
        mapi: str
            The API of Materials Project. Default: None. Ref. https://materialsproject.org/open 
        path_to_store_psp: str (path-like)
            The destination to store the compressed psp. default: psp_pymatgen
    Return
    """
    keys_required = ["PMG_DEFAULT_FUNCTIONAL", "PMG_MAPI_KEY", "PMG_VASP_PSP_DIR"]
    keys_dict = {"PMG_DEFAULT_FUNCTIONAL": def_fun, "PMG_VASP_PSP_DIR": path_to_store_psp, "PMG_MAPI_KEY": mapi}
    pmg_config_file = os.path.join(os.environ["HOME"], ".pmgrc.yaml")
    keys_exist = []
    params = {}
    if os.path.exists(pmg_config_file):
        pmg_config = loadfn(pmg_config_file)
        for key in keys_required:
            if key in pmg_config:
                if pmg_config[key]:
                    # Not empty or None
                    params[key] = pmg_config[key]
                    keys_exist.append(key)
        keys_required = list(set(keys_required).difference(set(keys_exist)))
        if len(keys_required) == 0:
            print("The pymatgen has been configurated before.")
            return
        else:
            #Backup the .pmgrc.yaml file
            shutil.copyfile(pmg_config_file, pmg_config_file + ".dfttk.bak")
    for key in keys_required:
        params[key] = keys_dict[key]
    dumpfn(params, pmg_config_file, default_flow_style=False)
    if "PMG_VASP_PSP_DIR" in keys_required:
        #No configuration for psp path
        handle_potcar_gz(psp_dir=psp_dir, path_to_store_psp=path_to_store_psp, aci=aci)

def update_configfile(filename, base_file):
    """
    Update the filename base on base_file.
    Update all except the path/file which exists in filename but not in base_file
    
    Parameter
        filename: str (filename-like)
            The filename of config file
        base_file: str (filename-like)
            The base file
    Return
        None
    """
    ori_file = loadfn(filename)
    base_file = loadfn(base_file)
    for item in base_file:
        flag_update = True
        if item in ori_file:
            if isinstance(base_file[item], str) and isinstance(ori_file[item], str):
                if os.path.exists(ori_file[item]) and (not os.path.exists(base_file[item])):
                    flag_update = False
        if flag_update:
            ori_file[item] = base_file[item]
    with open(filename, 'w+') as f:
        if filename.endswith(".json"):
            from json import dump
            dump(ori_file, f, indent=4)
        elif filename.endswith(".yaml"):
            from yaml import dump
            dump(ori_file, f, default_flow_style=False, sort_keys=False, indent=4)

def config_atomate(path_to_store_config=".", config_folder="config", queue_script="vaspjob.pbs", 
    queue_type="pbs", vasp_cmd_flag="vasp_std"):
    """
    Configuration for atomate

    Parameter
        path_to_store_config: str (path-like)
            The path to store the configuration files. Default: .
        config_folder: str (path-like)
            The path of config files. Default: config
        queue_script: str (filename-like)
            The submitting script for the queue system. Default: vaspjob.pbs
        queue_type: str
            Th type of queue system. Now, only support pbs. Default: pbs
        vasp_cmd_flag: str
            The flag of vasp_cmd to distinguish the vasp_cmd with other commands in the queue script
    Return
        None
    """
    config_file = get_config_file(config_folder=config_folder, queue_script=queue_script)
        
    creat_folders(path_to_store_config + "/config")
    creat_folders(path_to_store_config + "/logs")
   
    if config_file[queue_script]:
        #If the pbs file exists
        param_dict = parse_queue_script(template=config_file[queue_script], 
                                        type=queue_type, vasp_cmd_flag=vasp_cmd_flag)
    else:
        param_dict = {}
    param_dict["path_to_store_config"] = path_to_store_config

    required_file = ["db.json", "my_launchpad.yaml"]
    option_file = ["FW_config.yaml", "my_fworker.yaml", "my_qadapter.yaml"]
    FileModule = {"FW_config.yaml": "ConfigFW", "my_fworker.yaml": "ConfigFworker", "my_qadapter.yaml": "ConfigQadapter",
                  "db.json": "ConfigDb", "my_launchpad.yaml": "ConfigLaunchFile"}
    files = required_file + option_file
    for file in files:
        eval(FileModule[file] + "(**param_dict).write_file()")
        if config_file[file]:
            update_configfile(os.path.join(param_dict["path_to_store_config"], "config/" + file), config_file[file])
        #else:
        #    eval(FileModule[file] + "(**param_dict).write_file()")
    #Add environment var
    FW_CONFIG_FILE_VAL = os.path.join(path_to_store_config, "config/FW_config.yaml")
    add_path_var(FW_CONFIG_FILE=FW_CONFIG_FILE_VAL)

class ConfigTemplate(object):
    """docstring for ConfigTemplate"""
    def __init__(self, **kwargs):
        super(ConfigTemplate, self).__init__()
        #The input should be a dict
        PATH_TO_STORE_CONFIG = kwargs.get("path_to_store_config", ".")
        PATH_TO_STORE_CONFIG = get_abspath(PATH_TO_STORE_CONFIG)
        self.PATH_TO_STORE_CONFIG = PATH_TO_STORE_CONFIG
        self.VASP_CMD = kwargs.get("vasp_cmd", "mpirun vasp_std")
        self.NNODES = kwargs.get("nodes", 1)
        self.PPNODE = kwargs.get("ppn", 24)
        self.WALLTIME = kwargs.get("walltime", "48:00:00")
        self.QUEUE = kwargs.get("queue", "open")
        self.PMEM = kwargs.get("pmem", "8gb")
        #For ACI, please provide the submit script
        #maybe you need set PRE_ROCKET = "module load intel/16.0.3 impi/5.1.3 vasp/5.4.1.05Feb16"
        self.PRE_ROCKET = kwargs.get("pre_rocket", "null")
        self.POST_ROCKET = kwargs.get("post_rocket", "null")

    def write_file(self):
        filename = os.path.join(self.PATH_TO_STORE_CONFIG, "config/" + self.FILENAME)
        with open(filename, 'w') as f:
            if filename.endswith(".json"):
                from json import dump
                dump(self.DATA, f, indent=4)
            elif filename.endswith(".yaml"):
                from yaml import dump
                dump(self.DATA, f, default_flow_style=False, sort_keys=False, indent=4)

class ConfigDb(ConfigTemplate):
    """docstring for ConfigDb"""
    def __init__(self, **kwargs):
        super(ConfigDb, self).__init__(**kwargs)
        self.FILENAME = "db.json"
        self.DATA = {
            "database": "dfttk_tests",
            "collection": "tasks",
            "host": "localhost",
            "port": 27017,
            "aliases": {}}

class ConfigLaunchFile(ConfigTemplate):
    """docstring for MyLaunchFile"""
    def __init__(self, **kwargs):
        super(ConfigLaunchFile, self).__init__(**kwargs)
        self.FILENAME = "my_launchpad.yaml"
        self.DATA = {"host": "localhost",
            "port": 27017,
            "name": "dfttk-fws",
            "ssl_ca_file": "null",
            "strm_lvl": "INFO",
            "user_indices": "[]",
            "wf_user_indices": "[]"}

class ConfigFW(ConfigTemplate):
    """docstring for ConfigFW"""
    def __init__(self, **kwargs):
        super(ConfigFW, self).__init__(**kwargs)
        self.FILENAME = "FW_config.yaml"
        self.DATA = {"CONFIG_FILE_DIR": os.path.join(self.PATH_TO_STORE_CONFIG, "config"),
            "LAUNCHPAD_LOC": os.path.join(self.PATH_TO_STORE_CONFIG, "config/my_launchpad.yaml"),
            "FWORKER_LOC": os.path.join(self.PATH_TO_STORE_CONFIG, "config/my_fworker.yaml"),
            "QUEUEADAPTER_LOC": os.path.join(self.PATH_TO_STORE_CONFIG, "config/my_qadapter.yaml"),
            "QUEUE_JOBNAME_MAXLEN": 15,
            "ADD_USER_PACKAGES": ["atomate.vasp.firetasks", "atomate.feff.firetasks"]
        }

class ConfigQadapter(ConfigTemplate):
    """docstring for ConfigQadapter"""
    def __init__(self, **kwargs):
        super(ConfigQadapter, self).__init__(**kwargs)
        self.FILENAME = "my_qadapter.yaml"
        self.DATA = {"_fw_name": "CommonAdapter",
            "_fw_q_type": "PBS",
            "rocket_launch": "rlaunch -c " + os.path.join(self.PATH_TO_STORE_CONFIG, "config") + " rapidfire",
            "nnodes": self.NNODES,
            "ppnode": self.PPNODE,
            "pmem": self.PMEM,
            "walltime": self.WALLTIME,
            "queue": self.QUEUE,
            "account": "open",
            "job_name": "dfttk",
            "pre_rocket": self.PRE_ROCKET,
            "post_rocket": self.POST_ROCKET,
            "logdir": os.path.join(self.PATH_TO_STORE_CONFIG, "logs")
        }

class ConfigFworker(ConfigTemplate):
    """docstring for ConfigFworker"""
    def __init__(self, **kwargs):
        super(ConfigFworker, self).__init__(**kwargs)
        self.FILENAME = "my_fworker.yaml"
        self.DATA = {"name": "ACI",
            "category": '',
            "query": '{}',
            "env":
                {"db_file": os.path.join(self.PATH_TO_STORE_CONFIG, "config/db.json"),
                 "vasp_cmd": self.VASP_CMD,
                 "scratch_dir": "null",
                 "incar_update":{"NCORE": math.floor(math.sqrt(int(self.PPNODE)))}}
            }

def test_config(test_pymagen=True, test_atomate=True):
    if test_pymagen:
        test_config_pymatgen()
    if test_atomate:
        test_config_atomate()

def test_config_pymatgen():
    path_mp_config = os.path.join(os.environ["HOME"], ".pmgrc.yaml")
    if os.path.exists(path_mp_config):
        mp_config = loadfn(path_mp_config)

        print("##########Start to test the PMG_MAPI_KEY paramter##########")
        help_cmd = "-mapi YOUR_MP_API_KEY"
        param = "PMG_MAPI_KEY"
        if param in mp_config:
            try:
                with MPRester(mp_config[param]) as mpr:
                    mpr.get_structure_by_material_id("mp-66")
                Tips().set_properly(MAPI_KEY=mp_config[param])
            except Exception as e:
                Tips().set_improper(param, help_cmd)
        else:
            Tips().set_not_exist(param, help_cmd)

        print("\n##########Start to test the PMG_DEFAULT_FUNCTIONAL paramter##########")
        help_cmd = "-df DEFAULT_FUNCTIONAL"
        param = "PMG_DEFAULT_FUNCTIONAL"
        flag_functional = False
        FUNCTIONAL_CHOICES = Potcar.FUNCTIONAL_CHOICES
        if param in mp_config:
            DF = mp_config[param]
            if DF in FUNCTIONAL_CHOICES:
                Tips().set_properly(PMG_DEFAULT_FUNCTIONAL=DF)
                flag_functional = True
            else:
                Tips().set_improper(param, help_cmd)
        else:
            Tips().set_not_exist(param, help_cmd)

        print("\n##########Start to test the PMG_VASP_PSP_DIR paramter##########")
        help_cmd = "-psp VASP_PSP_DIR"
        param = "PMG_VASP_PSP_DIR"
        flag_psp = False
        FUNCTIONAL_DIR = list(PotcarSingle.functional_dir.values())
        functional_name = {v: k for k, v in PotcarSingle.functional_dir.items()}
        if param in mp_config:
            psp_dir = mp_config[param]
            if os.path.isdir(psp_dir):
                flag_psp = False
                functional_available = []
                sub_dirs = os.listdir(psp_dir)
                for item in sub_dirs:
                    if item in FUNCTIONAL_DIR:
                        functional = functional_name[item]
                        try:
                            p = PotcarSingle.from_symbol_and_functional("Fe", functional)
                            functional_available.append(functional)
                            flag_psp = True
                        except Exception as e:
                            print("There are some problems of " + functional + " in " + os.path.join(psp_dir, item))
                if flag_psp:
                    Tips().functional_info(param)
                    print("\t The supported functional is/are: " + ", ".join(functional_available))
                    flag_psp = True
                else:
                    print("There is no available functional in " + psp_dir)
                    Tips().set_improper(param, help_cmd)
            else:
                print("The path " + psp_dir + " not exists.")
                Tips().set_improper(param, help_cmd)
        else:
            Tips().set_not_exist(param, help_cmd)
    else:
        Tips().set_not_exist("~/.pmgrc.yaml", "-mp")
    if flag_psp and flag_functional:
        if not (mp_config["PMG_DEFAULT_FUNCTIONAL"] in functional_available):
            print("\nWARNING: The default functional is not supported.\n")

def test_config_atomate():
    #TODO
    pass

class Tips(Exception):
    """docstring for DfttkConfigError"""
    def __init__(self):
        super(Tips, self).__init__()

    def functional_info(self, param):
        print("SUCCESSFUL: The " + param + " is set properly.")

    def set_properly(self, **kwargs):
        for item in kwargs:
            print("SUCCESSFUL: Your " + item + " is " + kwargs[item])

    def set_improper(self, param, command):
        print("ERROR: The setting of " + param + " is inappropriate. (This is not a valid " + param + ")")
        self.ref(command)

    def set_not_exist(self, param, command):
        print("ERROR: The " + param + " not exists.")
        self.ref(command)

    def ref(self, command):
        print("\t Please config it using 'dfttk config " + command + "'")
        print("\t Ref. https://github.com/hitliaomq/dfttk_tutorial/tree/master/config")