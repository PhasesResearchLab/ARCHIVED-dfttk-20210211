"""Custom utilities for interacting with the Materials Project.

Mostly for getting and manipulating structures. With all of the function definitions and docstrings,
these are more verbose """

from pymatgen import MPRester

# TODO: wrap MPRester calls in a try-except block to catch errors and retry automatically

def mp_structures_from_ids(mp_ids, API_KEY=None):
    """Returns a list of structures from MP ids

    Args:
        mp_ids ([str]): list of Materials Project ids in the form of 'mp-###'

    Kwargs:
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

    Kwargs:
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

    Kwargs:
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

    Kwargs:
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