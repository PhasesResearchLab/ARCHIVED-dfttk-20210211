"""Custom Phases Research Lab Workflows
"""

import numpy as np
from uuid import uuid4

from fireworks import Workflow
from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA
from atomate.vasp.powerups import add_common_powerups, add_modify_incar, add_wf_metadata
from atomate.vasp.workflows.base.gibbs import get_wf_gibbs_free_energy

from prlworkflows.prl_fireworks import OptimizeFW
from prlworkflows.input_sets import PRLRelaxSet, PRLStaticSet

def get_wf_robust_optimization(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None,
                               tag="", metadata=None, name='structure optimization'):
    """
    Return an optimization workflow tailored to calculations for unstable structures.

    Calculation steps:
    1. Relax volume until the volume change is less than 5%
    2. Relax ions
    3. Relax cell shape and ions

    Parameters
    ----------
    structure : Structure
        A pymatgen Structure to be optimized
    vasp_input_set : DictVaspInputSet
        VASP input set for the relaxation
    vasp_cmd : str
        Command to run
    db_file : str
        Path to file containing the database credentials.
    tag : str
        Some unique string that will be appended to the names of the fireworks so that the data from
        those tagged fireworks can be queried later during the analysis.
    metadata : dict
        meta data
    name : str
        Name for the workflow

    Returns
    -------
    Workflow

    Notes
    -----
    In the future, this may have some advanced symmetry checking. For now the workflow runs through
    the entire optimization procedure and it is up to the user to pick which ones to continue
    their calculations with.

    Changelog:
    0.1: 7-2-4 scheme with custodian full opt ISIF 7, custodian normal ISIF 2 and 4.
    """
    metadata = metadata or {}
    metadata.update({'robust_optimization_version': 0.1 }) # SEMVER naming scheme
    vasp_input_set = vasp_input_set or PRLRelaxSet(structure)

    vol_relax_fw = OptimizeFW(structure, isif=7, job_type='full_opt_run', name=name, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata)
    ion_relax_fw = OptimizeFW(structure, isif=2, name=name, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, parents=vol_relax_fw, metadata=metadata)
    ion_shape_relax_fw = OptimizeFW(structure, isif=4, name=name, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, parents=ion_relax_fw, metadata=metadata)

    fws = [vol_relax_fw, ion_relax_fw, ion_shape_relax_fw]

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_optimization(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None,
                        tag="", metadata=None, name='structure optimization'):
    """
    Return an optimization workflow for stable structures using a standard double relaxation.

    Parameters
    ----------
    structure : Structure
        A pymatgen Structure to be optimized
    vasp_input_set : DictVaspInputSet
        VASP input set for the relaxation
    vasp_cmd : str
        Command to run
    db_file : str
        Path to file containing the database credentials.
    tag : str
        Some unique string that will be appended to the names of the fireworks so that the data from
        those tagged fireworks can be queried later during the analysis.
    metadata : dict
        meta data
    name : str
        Name for the workflow

    Returns
    -------
    Workflow
    """
    metadata = metadata or {}
    vasp_input_set = vasp_input_set or PRLRelaxSet(structure)

    relax_fw = OptimizeFW(structure, isif=3, job_type='double_relaxation_run', name=name, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata)
    fws = [relax_fw]

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    return Workflow(fws, name=wfname, metadata=metadata)


def wf_gibbs_free_energy(structure, c=None):
    """The same as workflow Gibbs free energy preset in atomate proper, but allows for choosing
    between standard and robust optimization with the ``ROBUST`` boolean (default: False) in the config dict.

    Also uses a much more reasonable 7 deformations with +/- 10% volume

    Parameters
    ----------
    structure : Structure
        input structure
    c : dict
        workflow config dict

    Returns
    -------
        Workflow
    """
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
    robust_optimization = c.get("ROBUST", False)
    optimize = c.get("OPTIMIZE", False)
    eos = c.get("EOS", "vinet")
    qha_type = c.get("QHA_TYPE", "debye_model")
    # min and max temp in K, default setting is to compute the properties at 300K only
    t_min = c.get("T_MIN", 300.0)
    t_max = c.get("T_MAX", 300.0)
    t_step = c.get("T_STEP", 100.0)
    pressure = c.get("PRESSURE", 0.0)
    poisson = c.get("POISSON", 0.25)
    e_diff = c.get("EDIFF", 1e-6)
    anharmonic_contribution = c.get("ANHARMONIC_CONTRIBUTION", False)
    metadata = c.get("METADATA", {})
    # 7 deformed structures: from -10% to +10% volume
    defos = [(np.identity(3)*(1+x)**(1.0/3.0)).tolist() for x in np.linspace(-0.1, 0.1, 7)]
    deformations = c.get("DEFORMATIONS", defos)
    user_kpoints_settings = {"grid_density": 7000}

    tag = "gibbs group: >>{}<<".format(str(uuid4()))

    # static input set for the transmute firework
    uis_static = {
        "ISIF": 2,
        "ISTART": 1,
    }
    lepsilon = False
    if qha_type not in ["debye_model"]:
        lepsilon = True
        try:
            from phonopy import Phonopy
        except ImportError:
            raise RuntimeError("'phonopy' package is NOT installed but is required for the final "
                               "analysis step; you can alternatively switch to the qha_type to "
                               "'debye_model' which does not require 'phonopy'.")
    vis_static = PRLStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                             user_kpoints_settings=user_kpoints_settings,
                             user_incar_settings=uis_static)
    # get gibbs workflow and chain it to the optimization workflow
    wf_gibbs = get_wf_gibbs_free_energy(structure, user_kpoints_settings=user_kpoints_settings,
                                        deformations=deformations, vasp_cmd=vasp_cmd, db_file=db_file,
                                        eos=eos, qha_type=qha_type, pressure=pressure, poisson=poisson,
                                        t_min=t_min, t_max=t_max, t_step=t_step, metadata=metadata,
                                        anharmonic_contribution=anharmonic_contribution,
                                        tag=tag, vasp_input_set=vis_static)

    # chaining
    if optimize:
        # input set for structure optimization
        vis_relax = PRLRelaxSet(structure, force_gamma=True)
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": user_kpoints_settings})
        vis_relax = vis_relax.__class__.from_dict(v)
        name = "{} structure optimization".format(tag)
        if robust_optimization:
            wf = get_wf_robust_optimization(structure, vasp_cmd=vasp_cmd, db_file=db_file, name=name, vasp_input_set=vis_relax)
        else:
            # optimization only workflow
            wf = get_wf_optimization(structure, vasp_cmd=vasp_cmd, db_file=db_file, name=name, vasp_input_set=vis_relax)
        wf.append_wf(wf_gibbs, wf.leaf_fw_ids)
    else:
        wf = wf_gibbs

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"EDIFF": e_diff,
                                                                    "LAECHG": False}})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf
