#!python
# This script is used to query the mongodb
# 
from pymatgen import Structure
from atomate.vasp.database import VaspCalcDb
from monty.serialization import loadfn
from fireworks.fw_config import config_to_dict
from dfttk.utils import sort_x_by_y
from warnings import warn

PRO_COLLECTION_MAP = {'static': 'tasks', 'relax': 'relax_scheme', 'check_symmetry': 'relaxations',
                      'phonon': 'phonon', 'qha': 'qha', 'qha_phonon': 'qha_phonon', 'borncharge': 'borncharge'}

def get_eq_structure_by_metadata(metadata, db_file=None):
    '''
    Get the equilibrium structure by metadata

    Parameters
    ----------
        metadata: dict
            The metadata use for searching the database
        db_file: filepath-like
            The file path of db.json(the settings for mongodb file)
            if it is None or >>db_file<<, then using the db_file of fireworks configurations
    Returns
    -------
        eq_structure: pymatgen.Strucutre object
            the equilibrium structure
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    static_items = list(vasp_db.db['tasks'].find({'metadata': metadata}).sort('_id', 1))
    structure_list = [itemi['output']['structure'] for itemi in static_items]
    if structure_list:
        #not empty
        eq_structure = Structure.from_dict(structure_list[0])
        return eq_structure

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
            The structures are sorted by volume
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    static_items = list(vasp_db.db['tasks'].find({'metadata': metadata}))
    structure_list = [Structure.from_dict(itemi['output']['structure']) for itemi in static_items]
    volumes = [itemi['output']['structure']['lattice']['volume'] for itemi in static_items]
    energies = [itemi['output']['energy_per_atom'] for itemi in static_items]
    band_gap = [itemi['output']['direct_gap'] for itemi in static_items]
    structure_list = sort_x_by_y(structure_list, volumes)
    band_gap = sort_x_by_y(band_gap, volumes)
    energies = sort_x_by_y(energies, volumes)
    volumes = sorted(volumes)
    return (structure_list, energies, band_gap)


def is_property_exist_in_db(metadata, db_file=None, property='static'):
    '''
    Search the MongoDB for specific property by metadata
    '''    
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

def remove_data_by_metadata(tag, db_file=None, rem_mode='vol', forcedelete=False):
    '''
    rem_mode: str/list
        allvol: remove all volume collection
        vol: remove volume collection except bandstructure and dos
        property: other collections except volume
        all: all
        aeccar: remove all aeccar related
        chgcar: remove chgcar
        locpot: remove locpot
        
    '''
    if (db_file is None) or (db_file == '>>db_file<<'):
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    metadata = {'tag': tag}

    VOL1_COLLECTION = ['aeccar0', 'aeccar1', 'aeccar2', 'chgcar', 'locpot']
    VOL2_COLLECTION = ['bandstructure', 'dos']
    VOL_COLLECTION = VOL1_COLLECTION + VOL2_COLLECTION
    OTHER_COLLECTION = ['borncharge', 'phonon', 'qha', 'qha_phonon', 'relax',
                        'relax_scheme', 'relaxations', 'tasks']
    flag_remove = False
    if forcedelete:
        flag_remove = True
    else:
        if input('Are you sure? This will delete all data related with metadata.tag={}. (Y/N)'.format(tag))[0].upper() == 'Y':
            flag_remove = True
    if flag_remove:
        if isinstance(rem_mode, str):
            rem_mode = rem_mode.lower()
            if rem_mode == 'all':
                collections = VOL_COLLECTION + OTHER_COLLECTION
            elif rem_mode == 'allvol':
                collections = VOL_COLLECTION
            elif rem_mode == 'vol':
                collections = VOL1_COLLECTION
            elif rem_mode == 'property':
                collections = OTHER_COLLECTION
            elif rem_mode == 'aeccar':
                collections = ['aeccar0', 'aeccar1', 'aeccar2']
            else:
                collections = [rem_mode]
        elif isinstance(rem_mode, list):
            collections = rem_mode
        else:
            raise ValueError('Unsupported remove mode, please provide a str or list')

        for collectioni in collections:
            if collectioni in VOL_COLLECTION:
                collectioni_file = collectioni + '_fs.files'
                collectioni_chunk = collectioni + '_fs.chunks'
                #It has files and chunks
                static_items = list(vasp_db.db['tasks'].find({'metadata': metadata}))
                for itemi in static_items:
                    task_id = itemi['task_id']
                    files_id = list(vasp_db.db[collectioni_file].find({'metadata.task_id': task_id}))
                    if files_id:
                        vasp_db.db[collectioni_chunk].remove({'files_id': files_id[0]['_id']})
                        vasp_db.db[collectioni_file].remove({'metadata.task_id': task_id})
                        print('The volume data with metadata.tag={} in {} collection is removed'.format(tag, collectioni))
            else:
                vasp_db.db[collectioni].remove({'metadata': metadata})
                print('The data with metadata.tag={} in {} collection is removed'.format(tag, collectioni))
