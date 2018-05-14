"""Custom Phases Research Lab Workflows
"""

import numpy as np
from uuid import uuid4

from fireworks import Workflow, Firework
from atomate.vasp.config import VASP_CMD, DB_FILE

from prlworkflows.analysis.phonon import create_supercell_displacements
from prlworkflows.prl_firetasks import CalculatePhononThermalProperties
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


def get_wf_phonon(structure, supercell_matrix, smearing_type='methfessel-paxton', displacement_distance=0.01, vasp_cmd=None, name='Phonon'):
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

