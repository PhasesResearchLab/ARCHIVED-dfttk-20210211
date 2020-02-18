#!python
#

#Temp test

import os
import re
import json
from pathlib import Path

atat_sqsdb_path = "/storage/home/mjl6505/package/atat/atat/data/sqsdb"

def SQSDatabaseATAT(atat_sqsdb_path):
    """
    Generate the SQS database using the build-in sqsdb in ATAT
    """
    sqsfilename = "bestsqs.out"
    for diri in os.listdir(atat_sqsdb_path):
        sqsgen_path = os.path.join(atat_sqsdb_path, diri)
        if os.path.isdir(sqsgen_path):
            # prototype = [name_prototype, Strukturbericht_mark]
            prototype = diri.split("_")
            for atatsqs_path in os.listdir(sqsgen_path):
                sqs_path = os.path.join(sqsgen_path, atatsqs_path)
                if os.path.isdir(sqs_path):
                    sqs_config = parse_atatsqs_path(atatsqs_path)
                    sqsfile_fullpath = os.path.join(sqs_path, sqsfilename)
                    if os.path.exists(sqsfile_fullpath):
                        with open(sqsfile_fullpath, "r") as fid:
                            #sqs_str = fid.read()
                            sqs_dict = lat_in_to_sqs(fid.read()).as_dict()
                            print(sqs_dict)
            '''
            sqs_folders, sqs_config = read_sqsgen_in(sqsgen_path)
            for sqs_folder in sqs_folders:
                sqsfile_fullpath = os.path.join(sqsgen_path, sqs_folder, sqsfilename)
                print(sqsfile_fullpath)
                if os.path.exists(sqsfile_fullpath):
                    with open(sqsfile_fullpath, "r") as fid:
                        #sqs_str = fid.read()
                        sqs_dict = lat_in_to_sqs(fid.read()).as_dict()
                        print(sqs_dict)
            '''

def parse_atatsqs_path(atatsqs_path):
    """
    Parse the path of atat sqsdb

    Parameters
    ----------
        atatsqs_path: str
            The path of atat sqsdb, e.g. sqsdb_lev=2_a=0.5,0.5_f=0.5,0.5
    Return
    ------
        sqs_config: dict
            The dict of the configuration of sqs
            It contains the following keys:
                level: The level of sqs, usually 0 for pure elements, 1 for 50%-50%, ...
                configuration: The configuration of sqs, e.g. ['a', 'c']
                occupancies: The occupancies of sqs, e.g. [0.5, 0.5]
    """
    sqs_config = {}
    configuration = []
    occupancies = []
    path_list = atatsqs_path.split("_")
    for path_param in path_list:
        path_val = path_param.split("=")
        if path_val[0] == "sqsdb":
            pass
        elif path_val[0] == "lev":
            level = path_val[1]
        else:
            configuration.append(path_val[0])
            occupancies.append(path_val[1].split(","))
    sqs_config["level"] = level
    sqs_config["configuration"] = configuration
    sqs_config["occupancies"] = occupancies
    return sqs_config

SQSDatabaseATAT(atat_sqsdb_path)
#



