"""Custom utilities for interacting with the Materials Project.

Mostly for getting and manipulating structures. With all of the function definitions and docstrings,
these are more verbose """

import fnmatch
import os

from pymatgen import MPRester
from fireworks import LaunchPad
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
