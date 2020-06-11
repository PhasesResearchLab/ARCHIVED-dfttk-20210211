"""
Phonon analysis using phonopy
"""

from phonopy import Phonopy
from phonopy.interface.vasp import Vasprun as PhonopyVasprun
from pymatgen.io.phonopy import get_phonopy_structure
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
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
    ph.run_mesh(qpoint_mesh)
    ph.run_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    # the thermal properties are for the unit cell
    tp_dict = ph.get_thermal_properties_dict()
    temperature = tp_dict['temperatures']
    # convert the units into our expected eV/atom-form (and per K)
    f_vib = tp_dict['free_energy'] * J_per_mol_to_eV_per_atom*1000
    s_vib = tp_dict['entropy'] * J_per_mol_to_eV_per_atom
    cv_vib = tp_dict['heat_capacity'] * J_per_mol_to_eV_per_atom
    return temperatures, f_vib, s_vib, cv_vib, ph.force_constants


def get_phonon_band_dos(structure, supercell_matrix, force_constants, qpoint_mesh=(50, 50, 50), band_paths=None, 
                        npoints=51, labels=None, phonon_dos=True, phonon_band=True, phonon_pdos=False, 
                        save_data=False, save_fig=False):
    '''
    Return the phonon dos and band

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
    band_paths :  list, multi dimention
        Sets of end points of paths, e.g. [[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0.5, 0.5]], [[0.5, 0.25, 0.75], [0, 0, 0]]]
    Returns
    -------
    '''
    volume = structure.volume
    formula = structure.composition.reduced_formula
    pre_filename = "{}-phonon-Vol{:.2f}".format(formula, volume)

    unitcell = get_phonopy_structure(structure)
    ph = Phonopy(unitcell, supercell_matrix)
    ph.set_force_constants(force_constants)

    phonon_dos = np.vstack((ph._total_dos._frequency_points, ph._total_dos._dos))

    if phonon_band:
        if band_paths:
            qpoints, connections = get_band_qpoints_and_path_connections(band_paths, npoints=npoints)
            ph.run_band_structure(qpoints, path_connections=connections, labels=labels)
        else:
            ph.auto_band_structure()
        if save_fig:
            fig_band = ph.plot_band_structure()
            fig_band.savefig(fname = '{}-band.png'.format(pre_filename))
            fig_band.close()
        if save_data:
            ph.write_yaml_band_structure(filename = '{}-band.yaml'.format(pre_filename))

    #for dos
    phonon_total_dos = None
    if phonon_dos:
        ph.run_mesh(qpoint_mesh)
        ph.run_total_dos()
        phonon_total_dos = np.vstack((ph._total_dos._frequency_points, ph._total_dos._dos))
        if save_fig:
            fig_dos = ph.plot_total_dos()
            fig_dos.savefig(fname = '{}-dos.png'.format(pre_filename))
            fig_dos.close()
        if save_data:
            ph.write_total_dos(filename = '{}-dos.dat'.format(pre_filename))
    #for pdos.
    if phonon_pdos:
        ph.run_mesh(qpoint_mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        ph.run_projected_dos()
        if save_fig:
            ph.plot_projected_dos().savefig(fname = '{}-pdos.png'.format(filename))
        if save_data:
            ph.write_projected_dos(filename = '{}-pdos.dat'.format(filename))

    return phonon_total_dos