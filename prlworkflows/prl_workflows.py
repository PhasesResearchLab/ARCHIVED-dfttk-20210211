"""Custom Phases Research Lab Workflows
"""

import numpy as np
from uuid import uuid4

from fireworks import Workflow, Firework
from atomate.vasp.config import VASP_CMD, DB_FILE

from prlworkflows.analysis.phonon import create_supercell_displacements
from prlworkflows.prl_firetasks import CalculatePhononThermalProperties, QHAAnalysis
from prlworkflows.prl_fireworks import PRLOptimizeFW, PRLStaticFW, PRLPhononDisplacementFW
from prlworkflows.input_sets import PRLRelaxSet, PRLStaticSet


def get_wf_ev_curve(structure, num_deformations=7, deformation_fraction=0.05, vasp_cmd=None, db_file=None, metadata=None, name='EV_Curve'):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
    num_deformations: int
    deformation_fraction: float
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    metadata = metadata or {}
    if 'tag' not in metadata.keys():
        metadata['tag'] = '{}'.format(str(uuid4()))

    deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)

    # follow a scheme of
    # 1. ISIF 2
    # 2. ISIF 4
    # 3. Static
    fws = []
    for i, deformation in enumerate(deformations):
        struct = structure.copy()
        struct.scale_lattice(struct.volume*deformation)
        vis = PRLRelaxSet(struct, user_incar_settings={'ISIF': 2})
        isif_2_fw = PRLOptimizeFW(structure, job_type='normal', name='structure_{}-relax-isif_2'.format(i), prev_calc_loc=True, vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata)
        fws.append(isif_2_fw)

        vis = PRLRelaxSet(struct, user_incar_settings={'ISIF': 4})
        isif_4_fw = PRLOptimizeFW(structure, job_type='normal', name='structure_{}-relax-isif_4'.format(i), prev_calc_loc=True, vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata, parents=isif_2_fw)
        fws.append(isif_4_fw)

        vis = PRLStaticSet(struct)
        static = PRLStaticFW(structure, name='structure_{}-static'.format(i), vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata, parents=isif_4_fw)
        fws.append(static)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_phonon_single_volume(structure, supercell_matrix, smearing_type='methfessel-paxton', displacement_distance=0.01, vasp_cmd=None, name='Phonon'):
    """
    A workflow to calculate the phonon vibration properties for a single volume.

    Parameters
    ----------
    structure : pymatgen.Structure
        Input structure. Should be a unitcell, not a supercell.
    supercell_matrix : numpy.ndarray
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]].
    smearing_type : str
        It must be one of 'gaussian', 'methfessel-paxton', or 'tetrahedron'. The default
        is 'methfessel-paxton', which uses a SIGMA of 0.2 and is well suited for metals,
        but it should not be used for semiconductors or insulators. Using 'tetrahedron'
        or 'gaussian' gives a SIGMA of 0.05. Any further customizations should use a custom workflow.
    displacement_distance : float
        Distance of each displacement. Defaults to 0.01, consistent with phonopy.
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    name : str
        Name of the workflow

    Returns
    -------
    Workflow
        A workflow for a single volume phonon calculation.

    """
    vasp_cmd = vasp_cmd or VASP_CMD

    supercells, displacement_dicts = create_supercell_displacements(structure, supercell_matrix, displacement_distance=displacement_distance)

    fws = []
    for sc, disp_dict in zip(supercells, displacement_dicts):
        fws.append(PRLPhononDisplacementFW(sc, disp_dict, name="phonon_displacement", smearing_type=smearing_type, vasp_cmd=vasp_cmd))

    thermal_props = Firework([CalculatePhononThermalProperties()], parents=fws, name='CalculateThermalProperties',
                               spec={'unitcell': structure, 'supercell_matrix': supercell_matrix},)
    fws = fws + [thermal_props]  # have to add here to create a new object, otherwise the thermal_props parents includes self.
    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    return Workflow(fws, name=wfname)


def get_wf_gibbs(structure, num_deformations=7, deformation_fraction=(-0.05, 0.1), phonon=False, phonon_supercell_matrix=None, phonon_kwargs=None, t_min=5, t_max=2000, t_step=5, vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA'):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
    num_deformations: int
    deformation_fraction: float
        Can be a float (a single value) or a 2-type of a min,max deformation fraction.
        Default is (-0.05, 0.1) leading to volumes of 0.95-1.10. A single value gives plus/minus
        by default.
    phonon : bool
        Whether to do a phonon calculation. Defaults to False, meaning the Debye model.
    phonon_supercell_matrix : list
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]]. Must be specified if phonon is specified.
    phonon_kwargs : dict
        Keyword arguments to send to the GeneratePhononDisplacements Firetask.
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE
    phonon_kwargs = phonon_kwargs or {}
    # update these to be consistent with the QHA.
    phonon_kwargs.update({
        't_min': t_min,
        't_max': t_max,
        't_step': t_step,
        'supercell_matrix': phonon_supercell_matrix,
    })

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    if 'tag' not in metadata.keys():
        metadata['tag'] = tag

    if isinstance(deformation_fraction, (list, tuple)):
        deformations = np.linspace(1+deformation_fraction[0], 1+deformation_fraction[1], num_deformations)
    else:
        deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)

    # follow a scheme of
    # 1. ISIF 2
    # 2. ISIF 4
    # 3. Static
    fws = []
    static_calcs = []
    for i, deformation in enumerate(deformations):
        struct = structure.copy()
        struct.scale_lattice(struct.volume*deformation)
        vis = PRLRelaxSet(struct, user_incar_settings={'ISIF': 2})
        isif_2_fw = PRLOptimizeFW(structure, job_type='normal', name='structure_{}-relax-isif_2'.format(i), prev_calc_loc=True, vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata)
        fws.append(isif_2_fw)

        vis = PRLRelaxSet(struct, user_incar_settings={'ISIF': 4})
        isif_4_fw = PRLOptimizeFW(structure, job_type='normal', name='structure_{}-relax-isif_4'.format(i), prev_calc_loc=True, vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata, parents=isif_2_fw)
        fws.append(isif_4_fw)

        vis = PRLStaticSet(struct)
        static = PRLStaticFW(structure, name='structure_{}-static'.format(i), vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata, parents=isif_4_fw, phonon_detour=phonon, phonon_kwargs=phonon_kwargs)
        fws.append(static)
        static_calcs.append(static)

    qha_fw = Firework(QHAAnalysis(phonon=phonon, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag), parents=static_calcs, name="{}-qha_analysis".format(structure.composition.reduced_formula))
    fws.append(qha_fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)
