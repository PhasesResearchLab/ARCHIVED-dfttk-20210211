"""Custom utilities for interacting with the Materials Project.

Mostly for getting and manipulating structures. With all of the function definitions and docstrings,
these are more verbose """

import fnmatch
import os

from pymatgen import MPRester
from fireworks import LaunchPad
import numpy as np
import scipy

# TODO: wrap MPRester calls in a try-except block to catch errors and retry automatically

eV_per_atom_to_J_per_mol = scipy.constants.eV*scipy.constants.Avogadro
J_per_mol_to_eV_per_atom = 1/(scipy.constants.eV*scipy.constants.Avogadro)

def mp_structures_from_ids(mp_ids, API_KEY=None):
    """Returns a list of structures from MP ids

    Args:
        mp_ids ([str]): list of Materials Project ids in the form of 'mp-###'
        API_KEY (str): your Materials Project API_KEY. Will try to use environment key if None.

    Returns:
        List of Structure objects
    """
    structs = []
    with MPRester(API_KEY) as mpr:
        for mp_id in mp_ids:
            structs.append(mpr.get_structure_by_material_id(mp_id))
    return structs


def mp_structures_from_system(system, API_KEY=None):
    """Supply a chemical system (e.g. Fe-Cr) and get all of the structures back

    Args:
        system (str): system name (e.g. Fe-Cr)
        API_KEY (str): your Materials Project API_KEY. Will try to use environment key if None.

    Returns:
        List of Structure objects
    """
    with MPRester(API_KEY) as mpr:
            structs = mpr.get_structures(system)
    return structs


def mp_structures_and_energies_from_system(system, API_KEY=None):
    """Supply a chemical system (e.g. Fe-Cr) and get dicts of the structures and properties back

    Args:
        system (str): system name (e.g. Fe-Cr), but could also be mp-ids or formula
        API_KEY (str): your Materials Project API_KEY. Will try to use environment key if None.

    Returns:
        List of {"material_id": id, "pretty_formula": formula, "energy_per_atom"}
    """
    with MPRester(API_KEY) as mpr:
            entries = mpr.get_data(system)
    return entries


def mp_sorted_structures_from_system(system, filter_energy=0.2, API_KEY=None):
    """Supply a chemical system (e.g. Fe-Cr) and get back Structures sorted by energy above hull

    Energies too far above the hull can be removed with the filter_energy

    Args:
        system (str): system name (e.g. Fe-Cr), but could also be mp-ids or formula
        filter_energy (float): Maximum energy above hull allowed in eV
        API_KEY (str): your Materials Project API_KEY. Will try to use environment key if None.

    Returns:
        List of Structure objects sorted by energy above hull
    """
    entries = mp_structures_and_energies_from_system(system, API_KEY=API_KEY)
    # if Structure objects cannot be created from the entries
    mp_ids = [entry["material_id"] for entry in entries]
    energies_above_hull = [entry["e_above_hull"] for entry in entries]
    sorted_mp_ids = [mp_id for energy, mp_id in sorted(zip(energies_above_hull, mp_ids)) if energy <= filter_energy]
    sorted_structs = mp_structures_from_ids(sorted_mp_ids)

    return sorted_structs


def get_launchpad(launchpad_file=None):
    """
    Returns a LaunchPad object. If the launchpad_file is None, then try to auto load from environment

    Args:
        launchpad_file (File-like): A file-like or file path to the LaunchPad file.

    Returns:
        LaunchPad
    """
    if launchpad_file:
        if isinstance(launchpad_file, file):
            # a file object was found
            ext = launchpad_file.name.split('.')[-1]
            if ext == 'yaml':
                launchpad = LaunchPad.from_format(launchpad_file.read(), f_format='yaml')
            else:
                # assume json
                launchpad = LaunchPad.from_format(launchpad_file.read())
        else:
            # assume launchpad_file is a path
            launchpad = LaunchPad.from_file(launchpad_file)
    else:
        launchpad = LaunchPad.auto_load()
    return  launchpad


def update_fws_spec(wf, spec_dict, fw_name_constraint=None):
    """
    Update the fireworks matching the name constraint with the passed spec_dict. Can be used for
    generically updating the spec as long as update can be expressed as a dictionary.

    Args:
        wf (Workflow): The original workflow object
        spec_dict (dict): the keys and values to update in the spec, e.g. {'_queueadapter': {'walltime': '24:00:00'}}
        fw_name_constraint (str): a constraint on the FW name

    Returns:
        Workflow
    """
    for fw in wf.fws:
        if fw_name_constraint is None or fw_name_constraint in fw.name:
            fw.spec.update(spec_dict)
    return wf


def recursive_glob(start, pattern):
    """
    Recursively glob for the given pattern from the start directory.

    Taken from ESPEI.

    Args:
        start (str): Path of the directory to walk while for file globbing
        pattern (str): Filename pattern to match in the glob

    Returns:
        [str]: List of matched filenames
    """
    matches = []
    for root, dirnames, filenames in os.walk(start):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return sorted(matches)


def sort_x_by_y(x, y):
    """Sort a list of x in the order of sorting y"""
    return [xx for _, xx in sorted(zip(y, x), key=lambda pair: pair[0])]


def supercell_scaling_by_target_atoms(structure, min_atoms=60, max_atoms=120,
                                      target_shape='sc', lower_search_limit=-2, upper_search_limit=2,
                                      verbose=False):
    """
    Find a the supercell scaling matrix that gives the most cubic supercell for a
    structure, where the supercell has between the minimum and maximum nubmer of atoms.

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell of a structure
    min_atoms : target number of atoms in the supercell, defaults to 5
    max_atoms : int
        Maximum number of atoms allowed in the supercell
    target_shape : str
        Target shape of supercell. Could choose 'sc' for simple cubic or 'fcc' for face centered
        cubic. Default is 'sc'.
    lower_search_limit : int
        How far to search below the 'ideal' cubic scaling. Default is -2.
    upper_search_limit : int
        How far to search below the 'ideal' cubic scaling. Default is 2.
    verbose : bool
        Whether to print extra details on the cell shapes and scores. Useful for debugging.

    Returns
    -------
    numpy.ndarray
        2d array of a scaling matrix, e.g. [[3,0,0],[0,3,0],[0,0,3]]

    Notes
    -----
    The motiviation for this is for use in phonon calculations and defect calculations.
    It is important that defect atoms are far enough apart that they do not interact.
    Scaling unit cells that are not cubic by even dimensions might result in interacting
    defects. An example would be a tetragonal cell with 2x8x8 Ang lattice vectors being
    made into a 2x2x2 supercell. Atoms along the first dimension would not be very far
    apart.

    We are using a pure Python implementation from ASE, which is not very fast for a given
    supercell size. This allows for a variable supercell size, so it's going to be slow
    for a large range of atoms.

    The search limits are passed directloy to ``find_optimal_cell_shape``.
    They define the search space for each individual supercell based on the "ideal" scaling.
    For example, a cell with 4 atoms and a target size of 110 atoms might have an ideal scaling
    of 3x3x3. The search space for a lower and upper limit of -2/+2 would be 1-5. Since the
    calculations are based on the cartesian product of 3x3 matrices, large search ranges are
    very expensive.
    """
    from ase.build import get_deviation_from_optimal_cell_shape, find_optimal_cell_shape

    # range of supercell sizes in number of unitcells
    supercell_sizes = range(min_atoms//len(structure), max_atoms//len(structure) + 1)

    optimal_supercell_shapes = []  # numpy arrays of optimal shapes
    optimal_supercell_scores = []  # will correspond to supercell size

    # find the target shapes
    for sc_size in supercell_sizes:
        optimal_shape = find_optimal_cell_shape(structure.lattice.matrix, sc_size, target_shape, upper_limit=upper_search_limit, lower_limit=lower_search_limit, verbose = True)
        optimal_supercell_shapes.append(optimal_shape)
        optimal_supercell_scores.append(get_deviation_from_optimal_cell_shape(optimal_shape, target_shape))

    if verbose:
        for i in range(len(supercell_sizes)):
            print('{} {:0.4f} {}'.format(supercell_sizes[i], optimal_supercell_scores[i], optimal_supercell_shapes[i].tolist()))

    # find the most optimal cell shape along the range of sizes
    optimal_sc_shape = optimal_supercell_shapes[np.argmin(optimal_supercell_scores)]

    return optimal_sc_shape

def recursive_flatten(l):
    if l == []:
        return l
    if isinstance(l[0], list):
        return recursive_flatten(l[0]) + recursive_flatten(l[1:])
    return l[:1] + recursive_flatten(l[1:])


def mget(d, path):
    """Get from a dict using dot notation

    Parameters
    ----------
    d : dict
        Nested dictionary structure
    path : str
        Dot separated property, e.g. output.structure.lattice.volume

    Returns
    -------
    Object
        Value of the dictionary

    Examples
    --------
    >>> nested_dict = {'top_level': {'second_level': 'my_value'}}
    >>> mget(nested_dict, 'top_level.second_level')
    'my_value'

    """

    keys = path.split('.')
    current_path = ""
    curr_dict = d
    for k in keys:
        if not isinstance(curr_dict, dict):
            raise ValueError("Object at path \"{}\" (type: {}) is not a dictionary and \"{}\" in path \"{}\" cannot be looked up.".format(current_path[:-1], type(curr_dict), k, path))
        current_path += k
        try:
            curr_dict = curr_dict[k]
        except KeyError:
            raise KeyError("Cannot access key \"{}\" in path \"{}\". Possible keys are [{}].".format(k, current_path, ", ".join(curr_dict.keys())))

        current_path += "."
    return curr_dict


def get_mat_info(struct):
    """
    Get some basic information of the structure, e.g. name, configuration

    Parameters
    ----------
        struct: pymatgen.structure

    Returns
    -------
        name: string
            The name of the structure
        configuration: list
            The configuration of the structure
        occupancy:  list
            The occupancy of the structure
        site_ratio: list
            the site-ratio of the structure
    """
    name = struct.formula
    configuration = []
    occupancy = []
    site_ratio = []
    for e, a in struct.composition.items():
        configuration.append([str(e)])
        occupancy.append([1.0])
        site_ratio.append([a])
    return name, configuration, occupancy, site_ratio
 

def mark_adopted_TF(tag, db_file, adpoted):
    from atomate.vasp.database import VaspCalcDb
    vasp_db = VaspCalcDb.from_db_file(db_file, admin = True)
    if vasp_db:
        vasp_db.collection.update({'metadata.tag': tag}, {'$set': {'adopted': adpoted}}, upsert = True, multi = True)
        vasp_db.db['phonon'].update({'metadata.tag': tag}, {'$set': {'adopted': adpoted}}, upsert = True, multi = True)


def mark_adopted(tag, db_file, volumes):
    mark_adopted_TF(tag, db_file, False)             # Mark all as adpoted
    from atomate.vasp.database import VaspCalcDb
    vasp_db = VaspCalcDb.from_db_file(db_file, admin = True)
    for volume in volumes:
        vasp_db.collection.update({'$and':[ {'metadata.tag': tag}, {'output.structure.lattice.volume': volume} ]},
                                  {'$set': {'adopted': True}}, upsert = True, multi = False)            # Mark only one
        vasp_db.db['phonon'].update({'$and':[ {'metadata.tag': tag}, {'volume': volume} ]},
                                    {'$set': {'adopted': True}}, upsert = True, multi = False)


def consistent_check_db(db_file, tag):
    '''
    In the subsequent running(run DFTTK again with the same tag exists in Mongo DB), 
    if phonon method is committed, it'd better to check the lengths of 
    "task" and "phonon" collections.
    '''
    from atomate.vasp.database import VaspCalcDb
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    num_task = vasp_db.collection.count_documents({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]})
    num_phonon = vasp_db.db['phonon'].count_documents({'$and':[ {'metadata.tag': tag}, {'adopted': True} ]})
    if num_task == num_phonon:
        return(True)
    else:
        print('The records length of "task"(%s) differs to the length of "phonon"(%s) in mongodb.' 
              %(num_task, num_phonon))
        return(False)


def check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4):
    from atomate.vasp.database import VaspCalcDb
    vasp_db = VaspCalcDb.from_db_file(db_file, admin=True)
    if relax_path != '':
        if os.path.exists(relax_path):
            return(relax_path, run_isif2, pass_isif4)
    
    if vasp_db.db["relax"].count_documents({'metadata.tag': tag}) > 0:
        items = vasp_db.db["relax"].find({'metadata.tag': tag}).sort([('_id', -1)]).limit(1)
        if os.path.exists(items[0]['path']):
            print('Relax result "%s" with "run_isif2 = %s" and "run_isif4 = %s" has been found, and will be used for new static calculations.' 
                  %(items[0]['path'], items[0]['run_isif2'], items[0]['pass_isif4']))
            return(items[0]['path'], items[0]['run_isif2'], items[0]['pass_isif4'])
        else:
            print('Relax result "%s" has been found but NOT exists. Change tag and try again!' %relax_path)
            return('', run_isif2, pass_isif4)
    else:
        print('No relax result found.')
        return('', run_isif2, pass_isif4)
    

def add_modify_incar_by_FWname(wf, modify_incar_params):
    from atomate.vasp.powerups import add_modify_incar
    for keyword in modify_incar_params.keys():
        add_modify_incar(wf, modify_incar_params = modify_incar_params[keyword], fw_name_constraint = keyword)


def add_modify_kpoints(original_wf, modify_kpoints_params, fw_name_constraint=None):
    """
    Every FireWork that runs VASP has a ModifyIncar task just beforehand. For example, allows
    you to modify the INCAR based on the Worker using env_chk or using hard-coded changes.

    Args:
        original_wf (Workflow)
        modify_incar_params (dict) - dict of parameters for ModifyIncar.
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    from atomate.utils.utils import get_fws_and_tasks
    from dfttk.ftasks import ModifyKpoints
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp"):
        original_wf.fws[idx_fw].tasks.insert(idx_t, ModifyKpoints(modify_kpoints_params = modify_kpoints_params))
    return original_wf


def add_modify_kpoints_by_FWname(wf, modify_kpoints_params):
    for keyword in modify_kpoints_params.keys():
        add_modify_kpoints(wf, modify_kpoints_params = modify_kpoints_params[keyword], fw_name_constraint = keyword)


import re
from pymatgen import Structure
class metadata_in_POSCAR():
    '''
    First line in POSCAR is like: SIGMA1;[0.5,0.5]16[0.25,0.75]32...;SQS;    Occupancies of 0 in [,] could not omitted
    meta = metadata_in_POSCAR('POSCAR')
    for config in configs:                          #configs writen like [['V', 'Ni'], ['Cr']]
        metadata = meta.get_metadata(config) 
    '''
    def __init__(self, filename='POSCAR'):
        self.poscarfile = filename
        ss = self.parse_poscar()
        if len(ss) <= 0:
            self.phase_name = ''
        else:
            self.phase_name = ss[0]
        if len(ss) < 3:
            return
        self.occupancies = []
        self.site_ratios = []
        digiarrs = ss[1].split('[')
        for digis in digiarrs:
            occupancy = []
            if digis == '':
                continue
            digis = re.split(',|]', digis)
            for i in range(len(digis) - 1):
                occupancy.append(float(digis[i])) 
            self.site_ratios.append(int(digis[-1]))
            self.occupancies.append(occupancy)
        self.method = ss[2]
    
    def parse_poscar(self):
        '''
        To parse the first line in POSCAR
        Each tag word segmented by a semicolon(";")
        Tag word could be as followings:
            
        '''
        if not os.path.exists(self.poscarfile):
            print('''
#####################################################################################################
#                                                                                                   #
     The file "%s" does NOT exist, please biuld it with some instructions in first line!         
#                                                                                                   #
#####################################################################################################
                  ''' %(self.poscarfile))
            ss = list()
        else:
            file = open(self.poscarfile)
            firstline = file.readline()
            file.close
            firstline = firstline.strip('\n')
            firstline = firstline.replace(' ', '') 
            firstline = firstline.upper()
            ss = re.split('[:;~]', firstline)
            i = len(ss) - 1
            while i >= 0:
                if ss[i] == '':
                    ss.pop(i)
                i -= 1
        return(ss)  
    
    def get_metadata(self, config):
        '''
        configs writen like [['V', 'Ni'], ['Cr']]
        '''
        m = len(config)
        n = len(self.occupancies)
        if m != n:
            print('Material configuration number(%s) is NOT equal to the occupancy number(%s), please check!' 
                  %(m, n))
            return()
        for i in range(m):
            if len(config[i]) > len(self.occupancies[i]):
                print('Wrong configuration in %s, please check!' %config)
                return()
        if not self.check_POSCAR(config):
            return()
        metadata = {
    		'phase': self.phase_name,
    		'sublattice_model': {
    			'configuration': config,
    			'occupancies': self.occupancies,
    			'site_ratios': self.site_ratios
    			},
            'method': self.method
            }
        return(metadata)
    
    def check_POSCAR(self, config):      
        '''             
        First line must like [1,0]32 to match the elements in POSCAR, 0 could not ignored.
        '''
        # To check the sum of occupancies
        for m in range(len(self.site_ratios)):
            sum_occupancy = 0
            for n in range(len(self.occupancies[m])):
                sum_occupancy += self.occupancies[m][n]
            if abs(sum_occupancy - 1) > 1e-10:
                print('The sum of occupancies in %s is NOT equal to 1, please check!' %self.occupancies[m])
                return(False)

        # To check config and occupancy
        
        temp_struct = Structure.from_file(self.poscarfile)
        namelist_elements = []
        numlist_elements = []
        for e, a in temp_struct.composition.items():
            namelist_elements.append(e)
            numlist_elements.append(a)                        # [8.0, 24.0]
        num_element_firstline = 0
        for ocs in self.occupancies:
            num_element_firstline += len(ocs) 
        if len(numlist_elements) != num_element_firstline:
            print('The number of element kind(%s) in first line of POASCAR is NOT same one in the structure(%s), maybe "0" occupancy should be added.'
                  %(num_element_firstline, len(numlist_elements)))
            return(False)
        
        index = 0
        for m in range(len(self.site_ratios)):
            if len(config[m]) < len(self.occupancies[m]):
                num_occupancy = 0
                for n in range(len(self.occupancies[m])):
                    num_occupancy += self.occupancies[m][n] * self.site_ratios[m]
                if abs(num_occupancy - numlist_elements[index]) > 1e-10:
                    print('The sum of sites in %s%s is NOT equal to %s(Element: %s), please check!' 
                          %(self.occupancies[m], self.site_ratios[m], 
                            numlist_elements[index], namelist_elements[index]))
                    return(False)
                index += len(self.occupancies[m])  
            else:
                for n in range(len(self.occupancies[m])):
                    if abs(numlist_elements[index] - self.occupancies[m][n] * self.site_ratios[m]) > 1e-10:
                        print('The sites in %s * %s is NOT equal to %s(Element: %s), please check!' 
                              %(self.occupancies[m][n], self.site_ratios[m], 
                                numlist_elements[index], namelist_elements[index]))
                        return(False)
                    index += 1   
        return(True)
        