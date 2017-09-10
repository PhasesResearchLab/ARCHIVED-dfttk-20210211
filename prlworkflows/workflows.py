"""Custom Phases Research Lab Workflows

For now, mostly contains some vanilla (non-custodian) optimizations.
"""

import numpy as np
from uuid import uuid4

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from FireWorks import Workflow, FWAction
from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA
from atomate.vasp.powerups import add_common_powerups, add_modify_incar, add_wf_metadata
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.workflows.base.gibbs import get_wf_gibbs_free_energy

from prlworkflows.firetasks import CheckSymmetry
from prlworkflows.fireworks import OptimizeFW

def get_wf_robust_optimization(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None,
                               tag="", metadata=None, name='structure optimization'):
    """
    Return an optimization workflow tailored to calculations for unstable structures.

    Calculation steps:
    1. Relax volume
    2. Relax ions, if symmetry breaks then return the structure from #1
    3. Relax cell shape and ions, if symmetry breaks then return the structure from #2

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

    vasp_input_set = vasp_input_set or MPRelaxSet(structure)

    # volume relax
    vol_relax_fw = OptimizeFW(structure, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, isif=7)
    # ion relax
    ion_relax_fw = OptimizeFW(structure, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, isif=2, parents=vol_relax_fw)
    # ion and shape relax, will be added as a detour to #2 on sucess
    ion_shape_relax_fw = OptimizeFW(structure, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd, db_file=db_file, isif=4)

    # add the ion shape relax as a detour to the ion wf
    fail_action = FWAction()  # continue on failing
    ion_relax_fw.tasks.append(CheckSymmetry(rename_on_fail=True, pass_action=FWAction(detours=[ion_shape_relax_fw]), fail_action=fail_action))
    ion_shape_relax_fw.tasks.append(CheckSymmetry(rename_on_fail=True, fail_action=fail_action))

    fws = [vol_relax_fw, ion_relax_fw]

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

    Returns:
        Workflow
    """
    c = c or {}

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
    robust_optimization = c.get("ROBUST", False)
    eos = c.get("EOS", "vinet")
    qha_type = c.get("QHA_TYPE", "debye_model")
    # min and max temp in K, default setting is to compute the properties at 300K only
    t_min = c.get("T_MIN", 300.0)
    t_max = c.get("T_MAX", 300.0)
    t_step = c.get("T_STEP", 100.0)
    pressure = c.get("PRESSURE", 0.0)
    poisson = c.get("POISSON", 0.25)
    anharmonic_contribution = c.get("ANHARMONIC_CONTRIBUTION", False)
    metadata = c.get("METADATA", None)

    # 7 deformed structures: from -10% to +10% volume
    defos = [(np.identity(3)*(1+x)**(1.0/3.0)).tolist() for x in np.linspace(-0.1, 0.1, 7)]
    deformations = c.get("DEFORMATIONS", defos)
    user_kpoints_settings = {"grid_density": 7000}

    tag = "gibbs group: >>{}<<".format(str(uuid4()))

    # input set for structure optimization
    vis_relax = MPRelaxSet(structure, force_gamma=True)
    v = vis_relax.as_dict()
    v.update({"user_kpoints_settings": user_kpoints_settings})
    vis_relax = vis_relax.__class__.from_dict(v)

    if robust_optimization:
        wf = get_wf_robust_optimization(structure, vasp_cmd=vasp_cmd, db_file=db_file, )
    else:
        # optimization only workflow
        wf = get_wf(structure, "optimize_only.yaml",
                    params=[{"vasp_cmd": vasp_cmd,  "db_file": db_file,
                             "name": "{} structure optimization".format(tag)}],
                    vis=vis_relax)

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
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
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
    wf.append_wf(wf_gibbs, wf.leaf_fw_ids)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600,
                                                                    "EDIFF": 1e-6,
                                                                    "LAECHG": False}})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf
