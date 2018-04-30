"""Custom Phases Research Lab Workflows
"""
from pymatgen import Structure, Lattice
import numpy as np
from uuid import uuid4

from fireworks import Workflow
from atomate.vasp.config import VASP_CMD, DB_FILE

from prlworkflows.prl_fireworks import PRLOptimizeFW, PRLStaticFW
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
    db_file : str
        Points to the database JSON file
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

