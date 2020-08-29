#!python
# This script is used to query the mongodb
# 
from monty.serialization import loadfn
from fireworks.fw_config import config_to_dict

DB_FILE = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
PRO_COLLECTION_MAP = {'static': 'tasks', 'relax': 'relax_scheme', 'check_symmetry': 'relaxations',
                      'phonon': 'phonon', 'qha': 'qha', 'borncharge': 'borncharge'}

def is_property_exist_in_db(metadata, property, db_file=None):
    pass