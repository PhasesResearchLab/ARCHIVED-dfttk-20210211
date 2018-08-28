"""
Phonon analysis using phonopy
"""

from phonopy import Phonopy
from phonopy.interface.vasp import Vasprun as PhonopyVasprun
from pymatgen.io.phonopy import get_phonopy_structure
from dfttk.utils import J_per_mol_to_eV_per_atom


def get_f_vib_phonopy(structure, supercell_matrix, vasprun_path,
                     qpoint_mesh=(50, 50, 50), t_min=5, t_step=5, t_max=2000.0,):
    """
    Return F_vib(T) for the unitcell in eV/atom

    Parameters
    ----------
    structure : pymatgen.Structure
        Unitcell (not supercell) of interest.
    supercell_matrix : numpy.ndarray
        3x3 matrix of the supercell deformation, e.g. [[3, 0, 0], [0, 3, 0], [0, 0, 3]].
    vasprun_path : str
        String pointing to a vasprun.xml file from a force constants run
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
        Tuple of (temperature, F_vib, S_vib, Cv_vib, force_constants)

    """
    # get the force constants from a vasprun.xml file
    vasprun = PhonopyVasprun(vasprun_path)
    force_constants, elements = vasprun.read_force_constants()

    ph_unitcell = get_phonopy_structure(structure)
    ph = Phonopy(ph_unitcell, supercell_matrix)
    # set the force constants we found
    ph.set_force_constants(force_constants)
    # calculate the thermal properties
    ph.set_mesh(qpoint_mesh)
    ph.set_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    # the thermal properties are for the unit cell
    temperatures, f_vib, s_vib, cv_vib = ph.get_thermal_properties()
    # convert the units into our expected eV/atom-form (and per K)
    f_vib *= J_per_mol_to_eV_per_atom*1000
    s_vib *= J_per_mol_to_eV_per_atom
    cv_vib *= J_per_mol_to_eV_per_atom
    return temperatures, f_vib, s_vib, cv_vib, ph.force_constants

