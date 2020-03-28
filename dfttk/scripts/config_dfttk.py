# -*- coding: utf-8 -*-
#
#Configure script for dfttk

from dfttk.scripts.run_dfttk import get_abspath, creat_folders
from dfttk.utils import recursive_glob
from monty.serialization import loadfn, dumpfn
import getpass
import shutil
import math
import os
import re

def replace_file(filename, old_str, new_str):
    with open(filename, "r") as fid:
        lines = fid.readlines()
    with open(filename, "w+") as fid:
        for line in lines:
            a = re.sub(old_str, new_str, line)
            fid.writelines(a)

def default_path():
    import dfttk
    return os.path.dirname(dfttk.__file__)

def add_path_var(**kwarg):
    """
    Add path var to .bashrc
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
    """
    path_len = [len(path_i) for path_i in path_list]
    return path_list[path_len.index(min(path_len))]

def get_config_file(config_folder=".", queue_script="vaspjob.pbs"):
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

def parase_pbs_script(filename = "vaspjob.pbs"):
    s = {"-q": "queue", "-A": "account", "-N": "job_name", "-V": "env"}
    param_dict = {}
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
            elif eachline.startswith("module"):
                if "pre_rocket" in param_dict:
                    param_dict["pre_rocket"].append(eachline)
                else:
                    param_dict["pre_rocket"] = [eachline]
            elif eachline.startswith("mpirun"):
                # start with mpirun
                param_dict["vasp_cmd"] = eachline
    if "pre_rocket" in param_dict:
        if len(param_dict["pre_rocket"]) == 1:
            param_dict["pre_rocket"] = param_dict["pre_rocket"][0]
    return param_dict

def parse_submit_template(template="vaspjob.pbs", type="pbs"):
    param_dict = {}
    if type == "pbs":
        param_dict = parase_pbs_script(filename=template)
    else:
        print("Only PBS is supported now. Other system will coming soon...")
    return param_dict

def handle_potcar_gz(psp_dir="psp", path_to_store_psp="psp_pymatgen"):
    #The name_mappings is copied from pymatgen.cli.pmg_config
    name_mappings = {
        "potpaw_PBE": "POT_GGA_PAW_PBE",
        "potpaw_PBE_52": "POT_GGA_PAW_PBE_52",
        "potpaw_PBE_54": "POT_GGA_PAW_PBE_54",
        "potpaw_PBE.52": "POT_GGA_PAW_PBE_52",
        "potpaw_PBE.54": "POT_GGA_PAW_PBE_54",
        "potpaw_LDA": "POT_LDA_PAW",
        "potpaw_LDA.52": "POT_LDA_PAW_52",
        "potpaw_LDA.54": "POT_LDA_PAW_54",
        "potpaw_LDA_52": "POT_LDA_PAW_52",
        "potpaw_LDA_54": "POT_LDA_PAW_54",
        "potUSPP_LDA": "POT_LDA_US",
        "potpaw_GGA": "POT_GGA_PAW_PW91",
        "potUSPP_GGA": "POT_GGA_US_PW91"
    }
    name_list = list(name_mappings.keys()) + list(name_mappings.values())
    # file_str is not abspath, is relative path
    file_str = os.listdir(psp_dir)
    psp_uncompress = os.path.join(psp_dir, "psp_uncompress")
    creat_folders(psp_uncompress)
    for file_i in file_str:
        psp_old = os.path.join(psp_dir, file_i)
        if os.path.isdir(psp_old):
            if file_i in name_list:
                psp_new = os.path.join(psp_uncompress, file_i)
                shutil.copytree(psp_old, psp_new)
        elif os.path.isfile(psp_old):
            if file_i.endswith(".tar.gz"):
                file_without_gz = ".".join(file_i.split(".")[0:-2])
                if file_without_gz in name_list:
                    psp_new = os.path.join(psp_uncompress, file_without_gz)
                    creat_folders(psp_new)
                    # uncompress
                    os.system("tar -zxvf " + psp_old + " -C " + psp_new)
    # config the POTCAR
    os.system("pmg config -p " + psp_uncompress + " " + path_to_store_psp)
    # Remove the uncompress folder
    try:
        shutil.rmtree(psp_uncompress)
    except:
        os.system("chmod +w " + os.path.join(psp_uncompress, "*/*"))
        shutil.rmtree(psp_uncompress)

def config_pymatgen(psp_dir="psp", def_fun="PBE", mapi=None, path_to_store_psp="psp_pymatgen"):
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
        handle_potcar_gz(psp_dir=psp_dir, path_to_store_psp=path_to_store_psp)

def config_atomate(path_to_store_config=".", config_folder="config", queue_script="vaspjob.pbs", queue_type="pbs"):

    config_file = get_config_file(config_folder=config_folder, queue_script=queue_script)
        
    creat_folders(path_to_store_config + "/config")
    creat_folders(path_to_store_config + "/logs")
   
    if config_file[queue_script]:
        #If the pbs file exists
        param_dict = parse_submit_template(template=config_file[queue_script], type=queue_type)
    else:
        param_dict = {}
    param_dict["path_to_store_config"] = path_to_store_config

    required_file = ["db.json", "my_launchpad.yaml"]
    option_file = ["FW_config.yaml", "my_fworker.yaml", "my_qadapter.yaml"]
    FileModule = {"FW_config.yaml": "ConfigFW", "my_fworker.yaml": "ConfigFworker", "my_qadapter.yaml": "ConfigQadapter"}
    files = required_file + option_file
    for file in files:
        if config_file[file]:
            shutil.copyfile(config_file[file], os.path.join(param_dict["path_to_store_config"], "config/" + file))
        else:
            eval(FileModule[file] + "(**param_dict).write_file()")
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


