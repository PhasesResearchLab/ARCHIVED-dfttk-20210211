#!python
# This script is used to query the mongodb
# 
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from monty.serialization import loadfn
from fireworks.fw_config import config_to_dict
from dfttk.utils import sort_x_by_y

def get_static_structure_by_metadata(metadata, db_file=None):
    '''
    Get the static structure by metadata

    Parameters
    ----------
        metadata: dict
            The metadata use for searching the database
        db_file: filepath-like
            The file path of db.json(the settings for mongodb file)
            if it is None or >>db_file<<, then using the db_file of fireworks configurations
    Returns
    -------
        structure_list: list
            The list of different structures.
            The structures are sorted by energy, the first one is the equilibrium structure
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    static_items = list(vasp_db.db['tasks'].find({'metadata': metadata}))
    structure_list = [Structure.from_dict(itemi['output']['structure']) for itemi in static_items]
    energies = [itemi['output']['energy_per_atom'] for itemi in static_items]
    band_gap = [itemi['output']['direct_gap'] for itemi in static_items]
    structure_list = sort_x_by_y(structure_list, energies)
    band_gap = sort_x_by_y(band_gap, energies)
    energies = sorted(energies)
    return (structure_list, energies, band_gap)


def is_property_exist_in_db(metadata, db_file=None, property='static'):
    '''
    Search the MongoDB for specific property by metadata
    '''
    PRO_COLLECTION_MAP = {'static': 'tasks', 'relax': 'relax_scheme', 'check_symmetry': 'relaxations',
                          'phonon': 'phonon', 'qha': 'qha', 'borncharge': 'borncharge'}
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    if property == 'static':
        return get_static_structure_by_metadata(metadata=metadata, db_file=db_file)
    else:
        vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
        collection = PRO_COLLECTION_MAP[property]
        search_items = vasp_db.db[collection].find({'metadata': self.metadata})
        if search_items:
            #Not empty
            return get_static_structure_by_metadata(metadata=metadata, db_file=db_file)
        else:
            return False