"""Custom Phases Research Lab Workflows
"""

import numpy as np
from uuid import uuid4

from fireworks import Workflow, Firework
from atomate.vasp.config import VASP_CMD, DB_FILE

from prlworkflows.prl_firetasks import QHAAnalysis
from prlworkflows.prl_fireworks import PRLOptimizeFW, PRLStaticFW, PRLPhononFW
from prlworkflows.input_sets import PRLRelaxSet, PRLStaticSet, PRLForceConstantsSet


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


def get_wf_gibbs(structure, num_deformations=7, deformation_fraction=(-0.05, 0.1),
                 phonon=False, phonon_supercell_matrix=None,
                 t_min=5, t_max=2000, t_step=5,
                 vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA'):
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
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
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
    qha_calcs = []
    # for each FW, we set the structure to the original structure to verify to ourselves that the
    # volume deformed structure is set by input set.
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

        if phonon:
            vis = PRLForceConstantsSet(struct)
            phonon_fw = PRLPhononFW(structure, phonon_supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step,
                     name='structure_{}-phonon'.format(i), vasp_input_set=vis,
                     vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata,
                     prev_calc_loc=True, parents=static)
            fws.append(phonon_fw)
            qha_calcs.append(phonon_fw)
        else:
            qha_calcs.append(static)

    qha_fw = Firework(QHAAnalysis(phonon=phonon, t_min=t_min, t_max=t_max, t_step=t_step, db_file=db_file, tag=tag), parents=qha_calcs, name="{}-qha_analysis".format(structure.composition.reduced_formula))
    fws.append(qha_fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)
