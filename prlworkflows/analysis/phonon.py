"""
Phonon analysis using phonopy
"""

# from the web for preprocessing:
# import numpy as np
# from phonopy import Phonopy

# phonon = Phonopy(unitcell,
#                  [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
#                  primitive_matrix=[[0, 0.5, 0.5],
#                                    [0.5, 0, 0.5],
#                                    [0.5, 0.5, 0]])
# phonon.generate_displacements(distance=0.03) # default is 0.1
# supercells = phonon.get_supercells_with_displacements()
# # save teh displacement dicts with
# ds = phonon.get_displacement_dataset()
# # attach each displacement in ds['first_atoms'] to the individual supercell generated

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.vasp import Vasprun as PhonopyVasprun

def get_all_force_sets(displacement_vasprun_files):
    """Will be replaced by reading the outputs from Fireworks"""
    force_sets = []
    for fp in displacement_vasprun_files:
        vasprun = PhonopyVasprun(fp)
        force_sets.append(vasprun.read_forces())
    return force_sets

def pmg_structure_to_phonopy_atoms(structure):
    """Convert a pymatgen Structure to PhonopyAtoms"""
    return PhonopyAtoms(symbols=[str(s.specie) for s in structure],
                 scaled_positions=structure.frac_coords,
                 cell=structure.lattice.matrix)

def get_f_vib_phonopy(structure, supercell_matrix, force_sets, displacement_dicts,
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

    """
    # reconstruct the displacement dataset:
    disp_dataset = {'first_atoms': [ds for ds in displacement_dicts], 'natom': len(structure)}

    # TODO: Do the displacement dicts need to be ordered or just correspond? Verify with a test
    ph_unitcell = pmg_structure_to_phonopy_atoms(structure)
    # I don't think I need a primitive matrix here, but it needs to be tested.
    # If we do need a primitive matrix, is this the matrix of the primitive cell even if the unit cell is conventional? Or is it just the unit cell matrix?
    ph = Phonopy(ph_unitcell, supercell_matrix)
    # set the forces and displacements
    ph.set_displacement_dataset(disp_dataset)
    ph.set_forces(force_sets)
    # make the force constants from the forces and displacements
    ph.produce_force_constants()
    # calculate the thermal properties
    ph.set_mesh(qpoint_mesh)
    ph.set_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    # the thermal properties are for the unit cell
    temperatures, f_vib, s_vib, cv_vib = ph.get_thermal_properties()
    return f_vib

test = False
# sample code for how to run:
if test:
    from pymatgen import Structure
    struct = Structure.from_file('/Users/brandon/Projects/phonopy-runs/Al-QHA/run_on_cluster/QHA-01/POSCAR')

    # TODO: sorting probably really important here. Displacement dict order and vasprun order must agree.
    displacement_vasprun_files = [
        '/Users/brandon/Projects/phonopy-runs/Al-QHA/run_on_cluster/QHA-01/disp-001/vasprun.xml'
    ]

    displacement_dicts = [
        {'number': 0, 'displacement': np.array([ 0.01,  0.  ,  0.  ]), 'direction': [1, 0, 0]},
    ]

    fs = get_all_force_sets(displacement_vasprun_files)
    t, f, cv = get_f_vib_phonopy(struct, np.eye(3, dtype=np.int), fs, displacement_dicts)


