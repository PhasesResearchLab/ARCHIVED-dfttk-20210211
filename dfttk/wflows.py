"""
Custom DFTTK Workflows
"""

import numpy as np
from uuid import uuid4
from copy import deepcopy
from fireworks import Workflow, Firework
from atomate.vasp.config import VASP_CMD, DB_FILE
from dfttk.fworks import OptimizeFW, StaticFW, RobustOptimizeFW
from dfttk.ftasks import CheckRelaxScheme
from dfttk.input_sets import PreStaticSet, RelaxSet
from dfttk.EVcheck_QHA import EVcheck_QHA, PreEV_check
from dfttk.utils import check_relax_path, add_modify_incar_by_FWname, add_modify_kpoints_by_FWname


def _get_deformations(def_frac, num_def):
    if isinstance(def_frac, (list, tuple)):
        return np.linspace(1 + def_frac[0], 1 + def_frac[1], num_def)
    else:
        return np.linspace(1 - def_frac, 1 + def_frac, num_def)


def get_wf_EV_bjb(structure, deformation_fraction=(-0.08, 0.12),
                  num_deformations=11, override_symmetry_tolerances=None, metadata=None):
    """
    Perform an E-V curve, robustly relaxating all structures on the curve.

    Parameters
    ----------
    structure : pymatgen.Structure

    deformation_fraction : tuple, optional
        Min and max volume scaling factors
    num_deformations : int, optional
        Number of volumes on the E-V curve
    override_symmetry_tolerances : Dict, optional
        See ``tol_X`` below.
    tol_energy : float, optional
        Symmetry is considered broken if the energy decrease for a relaxation
        step exceeds this value.
    tol_strain : float, optional
        Symmetry is considered broken if the the non-isotropic strain for a
        relaxation step exceeds this value.
    tol_bond : float, optional
        Symmetry is considered broken if the bond length per atom change for a
        relaxation step exceeds this value.

    """
    deformations = _get_deformations(deformation_fraction, num_deformations)*structure.volume

    fws = []
    for defo in deformations:
        struct = deepcopy(structure)
        struct.scale_lattice(defo)
        full_relax_fw = RobustOptimizeFW(struct, isif=5, vasp_cmd=VASP_CMD, db_file=DB_FILE)
        fws.append(full_relax_fw)
    if metadata is not None and all(x in metadata for x in ('phase_name', 'sublattice_configuration')):
        # create a nicer name for the workflow
        subl_config = ':'.join(','.join(subl_comps) for subl_comps in metadata['sublattice_configuration'])
        wfname = f"{metadata['phase_name']}:{subl_config}:{structure.composition.reduced_formula}"
    else:
        wfname = f"unknown:{structure.composition.reduced_formula}:unknown"
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf

def get_wf_gibbs_robust(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1), phonon=False, isif4=False,
                        phonon_supercell_matrix=None, override_symmetry_tolerances=None, t_min=5, t_max=2000, 
                        t_step=5, tolerance=0.01, volume_spacing_min=0.03, vasp_cmd=None, db_file=None, 
                        metadata=None, name='EV_QHA', override_default_vasp_params=None, modify_incar_params={},
                        modify_kpoints_params={}, verbose=False):
    """
    E - V
    curve

    workflow
    Parameters
    ------
    structure: pymatgen.Structure
        The initial structure
    num_deformations: int
        The number of deformation
    deformation_fraction: float
        Can be a float (a single value) or a 2-type of a min,max deformation fraction.
        Default is (-0.1, 0.1) leading to volumes of 0.90-1.10. A single value gives plus/minus by default.
    phonon : bool
        Whether to do a phonon calculation. Defaults to False, meaning the Debye model.
    phonon_supercell_matrix : list
        3x3 array of the supercell matrix, e.g. [[2,0,0],[0,2,0],[0,0,2]]. Must be specified if phonon is specified.
    override_symmetry_tolerances : dict
        The symmetry tolerance. It contains three keys, default: {'tol_energy':0.025, 'tol_strain':0.05, 'tol_bond':0.1}
    t_min : float
        Minimum temperature
    t_step : float
        Temperature step size
    t_max : float
        Maximum temperature (inclusive)
    tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    override_default_vasp_params: dict
        Override vasp parameters for all vasp jobs. e.g override_default_vasp_params={'user_incar_settings': {'ISIF': 4}}
    modify_incar_params : dict
        Override vasp settings in firework level
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
   """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    override_symmetry_tolerances = override_symmetry_tolerances or {'tol_energy':0.025, 'tol_strain':0.05, 'tol_bond':0.10}
    override_default_vasp_params = override_default_vasp_params or {}

    #TODO: Delete this paragraph, and ensure that any sub-function which call db_file is a Firetask
    if db_file == ">>db_file<<":
        #In PengGao's version, some function used the absolute db_file
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

    site_properties = deepcopy(structure).site_properties

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))
    metadata.update({'tag': tag})

    deformations = _get_deformations(deformation_fraction, num_deformations)
    vol_spacing = max((deformations[-1] - deformations[0]) / (num_deformations - 1), volume_spacing_min)

    fws = []
    
    vis_relax = RelaxSet(structure, isif=7, **override_default_vasp_params)
    full_relax_fw = RobustOptimizeFW(structure, isif=7, override_symmetry_tolerances=override_symmetry_tolerances,
                       prev_calc_loc=False, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file, isif4=isif4,
                       override_default_vasp_params=override_default_vasp_params, name='Full relax',
                       metadata=metadata, tag=tag)
    fws.append(full_relax_fw)

    #TODO: add a phonon after the relax
    #      How to pass the static structure to phonon
    '''
    if phonon:
        visphonon = ForceConstantsSet(struct)
        phonon_fw = PhononFW(struct, phonon_supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step,
                 name='structure_%.3f-phonon' %(vol_add), vasp_input_set=visphonon,
                 vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata,
                 prev_calc_loc='static', parents=full_relax_fw)
        fws.append(phonon_fw)
        calcs.append(phonon_fw)
    '''

    check_relax_fw = Firework(CheckRelaxScheme(db_file=db_file, tag=tag), parents=full_relax_fw)
    fws.append(check_relax_fw)

    check_result = Firework(EVcheck_QHA(db_file=db_file, tag=tag, deformations=deformations, site_properties=site_properties,
                                        tolerance=tolerance, threshold=14, vol_spacing=vol_spacing, vasp_cmd=vasp_cmd, 
                                        metadata=metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon=phonon, 
                                        phonon_supercell_matrix=phonon_supercell_matrix, verbose = verbose, 
                                        override_symmetry_tolerances=override_symmetry_tolerances,
                                        override_default_vasp_params=override_default_vasp_params,
                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params),
                            parents=check_relax_fw, name='%s-EVcheck_QHA' %structure.composition.reduced_formula)
    fws.append(check_result)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf


def get_wf_gibbs(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1), run_isif2=False,
                 phonon=False, phonon_supercell_matrix=None, pass_isif4=False,
                 t_min=5, t_max=2000, t_step=5, tolerance = 0.01, volume_spacing_min = 0.03,
                 vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA', symmetry_tolerance = 0.05,
                 passinitrun=False, relax_path='', modify_incar_params={},
                 modify_kpoints_params={}, verbose=False):
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
    tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    passinitrun : bool
        Set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    relax_path : str
        Set the path already exists for new static calculations; if set as '', will try to get the path from db_file.
    modify_incar_params : dict
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
    run_isif2: bool
        Whether run isif=2 calculation before isif=4 running.
    pass_isif4: bool
        Whether pass isif=4 calculation.
   """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    if db_file == ">>db_file<<":
        #In PengGao's version, some function used the absolute db_file
        from fireworks.fw_config import config_to_dict
        from monty.serialization import loadfn
        db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

    site_properties = deepcopy(structure).site_properties

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))

    if isinstance(deformation_fraction, (list, tuple)):
        deformations = np.linspace(1+deformation_fraction[0], 1+deformation_fraction[1], num_deformations)
        vol_spacing = max((deformation_fraction[1] - deformation_fraction[0]) / (num_deformations - 0.999999) + 0.001,
                          volume_spacing_min)
    else:
        deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)
        vol_spacing = max(deformation_fraction / (num_deformations - 0.999999) * 2 + 0.001,
                          volume_spacing_min)

    fws = []
    if 'tag' not in metadata.keys():
        metadata['tag'] = tag
    relax_path, run_isif2, pass_isif4 = check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4)

    if (relax_path == ''):
        # follow a scheme of
        # 1. Full relax + symmetry check
        # 2. If symmetry check fails, detour to 1. Volume relax, 2. inflection detection
        # 3. Inflection detection
        # 4. Static EV
        # 5. Phonon EV
        # for each FW, we set the structure to the original structure to verify to ourselves that the
        # volume deformed structure is set by input set.

        vis_relax = RelaxSet(structure)
        print('Full relax will be running ...')
        full_relax_fw = OptimizeFW(structure, symmetry_tolerance=symmetry_tolerance, job_type='normal', name='Full relax',
                                   prev_calc_loc=False, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file,
                                   metadata=metadata, record_path = True, run_isif2=run_isif2, pass_isif4=pass_isif4,
                                   modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params,
                                   spec={'_preserve_fworker': True})
        fws.append(full_relax_fw)
    else:
        full_relax_fw = None

    check_result = Firework(EVcheck_QHA(db_file = db_file, tag = tag, relax_path = relax_path, deformations = deformations, site_properties=site_properties,
                                        tolerance = tolerance, threshold = 14, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd, 
                                        metadata = metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon = phonon, symmetry_tolerance = symmetry_tolerance,
                                        phonon_supercell_matrix = phonon_supercell_matrix, verbose = verbose, run_isif2=run_isif2, pass_isif4=pass_isif4,
                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params),
                            parents=full_relax_fw, name='%s-EVcheck_QHA' %structure.composition.reduced_formula)
    fws.append(check_result)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf


def get_wf_gibbs_SQS(structure, num_deformations=7, deformation_fraction=(-0.1, 0.1),
                 phonon=False, phonon_supercell_matrix=None, run_isif2=False, pass_isif4=False,
                 t_min=5, t_max=2000, t_step=5, tolerance = 0.01, volume_spacing_min = 0.03,
                 vasp_cmd=None, db_file=None, metadata=None, name='EV_QHA', symmetry_tolerance = 0.05,
                 passinitrun=False, relax_path='', modify_incar_params={},
                 modify_kpoints_params={}, verbose=False):
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
    tolerance: float
        Acceptable value for average RMS, recommend >= 0.005.
    volume_spacing_min: float
        Minimum ratio of Volumes spacing
    vasp_cmd : str
        Command to run VASP. If None (the default) is passed, the command will be looked up in the FWorker.
    db_file : str
        Points to the database JSON file. If None (the default) is passed, the path will be looked up in the FWorker.
    name : str
        Name of the workflow
    metadata : dict
        Metadata to include
    passinitrun : bool
        Set True to pass initial VASP running if the results exist in DB, use carefully to keep data consistent.
    relax_path : str
        Set the path already exists for new static calculations; if set as '', will try to get the path from db_file.
    modify_incar_params : dict
        User can use these params to modify the INCAR set. It is a dict of class ModifyIncar with keywords in Workflow name.
    modify_kpoints_params : dict
        User can use these params to modify the KPOINTS set. It is a dict of class ModifyKpoints with keywords in Workflow name.
        Only 'kpts' supported now.
    run_isif2: bool
        Whether run isif=2 calculation before isif=4 running, pass to get_gibbs.
    pass_isif4: bool
        Whether pass isif=4 calculation, pass to get_gibbs.
    """
    vasp_cmd = vasp_cmd or VASP_CMD
    db_file = db_file or DB_FILE

    metadata = metadata or {}
    tag = metadata.get('tag', '{}'.format(str(uuid4())))

    if isinstance(deformation_fraction, (list, tuple)):
        deformations = np.linspace(1+deformation_fraction[0], 1+deformation_fraction[1], num_deformations)
        vol_spacing = max((deformation_fraction[1] - deformation_fraction[0]) / (num_deformations - 0.999999) + 0.001,
                          volume_spacing_min)
    else:
        deformations = np.linspace(1-deformation_fraction, 1+deformation_fraction, num_deformations)
        vol_spacing = max(deformation_fraction / (num_deformations - 0.999999) * 2 + 0.001,
                          volume_spacing_min)

    if 'tag' not in metadata.keys():
        metadata['tag'] = tag
    relax_path, run_isif2, pass_isif4 = check_relax_path(relax_path, db_file, tag, run_isif2, pass_isif4)

    fws = []
    prestatic_calcs = []
    if relax_path != '':
        pass
    else:
        for i, deformation in enumerate(np.linspace(0.85, 1.15, 7)):
            structure1 = deepcopy(structure)
            structure1.scale_lattice(deformation * structure.volume)
            vis_PreStatic = PreStaticSet(structure1)
            prestatic = StaticFW(structure=structure1, scale_lattice=deformation, name='VR_%.3f-PreStatic' %deformation,
                               prev_calc_loc=False, vasp_input_set=vis_PreStatic, vasp_cmd=vasp_cmd, db_file=db_file,
                               metadata=metadata, Prestatic=True)

            fws.append(prestatic)
            prestatic_calcs.append(prestatic)

    check_result = Firework(PreEV_check(db_file = db_file, tag = tag, relax_path = relax_path, deformations =deformations, structure = structure,
                                        tolerance = tolerance, threshold = 14, vol_spacing = vol_spacing, vasp_cmd = vasp_cmd, run_isif2=run_isif2,
                                        metadata = metadata, t_min=t_min, t_max=t_max, t_step=t_step, phonon = phonon, symmetry_tolerance = symmetry_tolerance,
                                        phonon_supercell_matrix = phonon_supercell_matrix, verbose = verbose, pass_isif4=pass_isif4,
                                        modify_incar_params=modify_incar_params, modify_kpoints_params = modify_kpoints_params),
                            parents=prestatic_calcs, name='%s-PreEV_check' %structure.composition.reduced_formula)
    fws.append(check_result)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    wf = Workflow(fws, name=wfname, metadata=metadata)
    add_modify_incar_by_FWname(wf, modify_incar_params = modify_incar_params)
    add_modify_kpoints_by_FWname(wf, modify_kpoints_params = modify_kpoints_params)

    return wf
