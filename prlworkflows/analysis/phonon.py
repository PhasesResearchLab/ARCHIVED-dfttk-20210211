"""
Phonon analysis using phonopy
"""

import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import Vasprun as PhonopyVasprun
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from prlworkflows.utils import J_per_mol_to_eV_per_atom


def get_all_force_sets(displacement_vasprun_files):
    """Will be replaced by reading the outputs from Fireworks"""
    force_sets = []
    for fp in displacement_vasprun_files:
        vasprun = PhonopyVasprun(fp)
        force_sets.append(vasprun.read_forces())
    return force_sets


def create_supercell_displacements(structure, supercell_matrix, displacement_distance=0.01):
    """

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    displacement_distance : float
        Distance of each displacement. Defaults to 0.01, consistent with phonopy.


    Returns
    -------
    list
        List of 2-tuples containing (supercell structure, displacement_dict)
    """
    ph_unitcell = get_phonopy_structure(structure)
    ph = Phonopy(ph_unitcell, supercell_matrix)
    ph.generate_displacements(distance=displacement_distance)
    supercells = [get_pmg_structure(sc) for sc in ph.get_supercells_with_displacements()]
    disp_dataset = ph.get_displacement_dataset()
    disp_dicts = disp_dataset['first_atoms']
    # phonopy does some weird stuff with the displacement dataset, so we will sanitize it
    for disp_dict in disp_dicts:
        disp_dict['number'] = int(disp_dict['number'])
        disp_dict['displacement'] = disp_dict['displacement'].tolist()  # numpy array to list
        disp_dict['direction'] = [int(d) for d in disp_dict['direction']]  # directions are always -1, 0, 1
    return supercells, disp_dicts


def get_f_vib_phonopy(structure, supercell_matrix, displacement_dicts, force_sets=None,
                     qpoint_mesh=(50, 50, 50), t_min=5, t_step=5, t_max=2000.0,):
    """
    Return F_vib(T)

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    force_sets : list
        List of forces in the supercell. Each element in the list should be a 2d array of
        shape (number of atoms, 3), for 3 cartesian coordinates. Note that the force sets
        must correspond to the displacement dicts.
    displacement_dicts : list
        List of displacment dictionaries. A displacement dictionary is a dictionary of
        ``{'direction': [1, 0, 0], 'displacement': np.array([0.01, 0, 0]), 'number': 0}``.
        Note that number is atom number, NOT sort order!
    qpoint_mesh : list
        Mesh of q-points to calculate thermal properties on.
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)

    Returns
    -------
    tuple
        Tuple of (temperature, F_vib)

    """
    # reconstruct the displacement dataset:
    disp_dataset = {'first_atoms': [ds for ds in displacement_dicts], 'natom': len(structure)}

    ph_unitcell = get_phonopy_structure(structure)
    # I don't think I need a primitive matrix here, but it needs to be tested.
    # If we do need a primitive matrix, is this the matrix of the primitive cell even if the unit cell is conventional? Or is it just the unit cell matrix?
    ph = Phonopy(ph_unitcell, supercell_matrix)
    # set the forces and displacements
    ph.set_displacement_dataset(disp_dataset)
    # force_sets can also be passed as 'forces' in displacement dicts, however they are updated here if passed
    if force_sets is not None:
        ph.set_forces(force_sets)
    # make the force constants from the forces and displacements
    ph.produce_force_constants()
    # calculate the thermal properties
    ph.set_mesh(qpoint_mesh)
    ph.set_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    # the thermal properties are for the unit cell
    temperatures, f_vib, s_vib, cv_vib = ph.get_thermal_properties()
    # convert the units into our expected eV/atom-form (and per K)
    f_vib *= J_per_mol_to_eV_per_atom*1000
    s_vib *= J_per_mol_to_eV_per_atom
    cv_vib *= J_per_mol_to_eV_per_atom
    return temperatures, f_vib, s_vib, cv_vib

test = False
# sample code for how to run:
if test:
    from pymatgen import Structure
    struct = Structure.from_file('/Users/brandon/Projects/phonopy-runs/Al-QHA/run_on_cluster/QHA-01/POSCAR')

    # sorting really important here. Displacement dict order and vasprun order must agree.
    displacement_vasprun_files = [
        '/Users/brandon/Projects/phonopy-runs/Al-QHA/run_on_cluster/QHA-01/disp-001/vasprun.xml'
    ]

    displacement_dicts = [
        {'number': 0, 'displacement': np.array([ 0.01,  0.  ,  0.  ]), 'direction': [1, 0, 0]},
    ]

    fs = get_all_force_sets(displacement_vasprun_files)
    t, f, cv = get_f_vib_phonopy(struct, np.eye(3, dtype=np.int), fs, displacement_dicts)


