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

    The search limits are passed directloy to ``find_optimal_cell_shape_pure_python``.
    They define the search space for each individual supercell based on the "ideal" scaling.
    For example, a cell with 4 atoms and a target size of 110 atoms might have an ideal scaling
    of 3x3x3. The search space for a lower and upper limit of -2/+2 would be 1-5. Since the
    calculations are based on the cartesian product of 3x3 matrices, large search ranges are
    very expensive.
    """
    from ase.build import get_deviation_from_optimal_cell_shape, find_optimal_cell_shape_pure_python

    # range of supercell sizes in number of unitcells
    supercell_sizes = range(min_atoms//len(structure), max_atoms//len(structure) + 1)

    optimal_supercell_shapes = []  # numpy arrays of optimal shapes
    optimal_supercell_scores = []  # will correspond to supercell size

    # find the target shapes
    for sc_size in supercell_sizes:
        optimal_shape = find_optimal_cell_shape_pure_python(structure.lattice.matrix, sc_size, target_shape, upper_limit=upper_search_limit, lower_limit=lower_search_limit)
        optimal_supercell_shapes.append(optimal_shape)
        optimal_supercell_scores.append(get_deviation_from_optimal_cell_shape(optimal_shape, target_shape))

    if verbose:
        for i in range(len(supercell_sizes)):
            print('{} {:0.4f} {}'.format(supercell_sizes[i], optimal_supercell_scores[i], optimal_supercell_shapes[i].tolist()))

    # find the most optimal cell shape along the range of sizes
    optimal_sc_shape = optimal_supercell_shapes[np.argmin(optimal_supercell_scores)]

    return optimal_sc_shape

