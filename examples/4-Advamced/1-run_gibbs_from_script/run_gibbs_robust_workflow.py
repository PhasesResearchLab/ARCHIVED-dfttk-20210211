# coding:utf-8
# The template for single run get_wf_gibbs_robust

################## COMMOM SETTING #############################
TEMPLATE_STRUCTURE_FILENAME = "POSCAR.Si"
#bool, submit(False) or not(True) the workflows to LaunchPad
DRY_RUN = False
#bool, write out(True) or not(False) the workflow (single workflow) to yaml file
WRITE_OUT_WF = False

################ PARAMETERS FOR WF #############################
#str, the absolute path of db.json file, e.g. /storage/home/mjl6505/atomate/config/db.json
#  If None, it will use the configuration in fireworks 
db_file = None
#list, the MAGMOM of the structure, e.g. [4.0, 4.0, -4.0, -4.0]
magmom = None
#int, the number of initial deformations, e.g. 7
num_deformations = 7
#list/tuple(min, max) or float(-max, max), the maximum amplitude of deformation, e.g. (-0.05, 0.1) means (0.95, 1.1) in volume
deformation_fraction = (-0.1, 0.1)
#float, minimum ratio of Volumes spacing, e.g. 0.03
volume_spacing_min = 0.03
#bool, run phonon(True) or not(False)
phonon = False
#list(3x3), the supercell matrix for phonon, e.g. [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
phonon_supercell_matrix = None
phonon_supercell_matrix_min = None
phonon_supercell_matrix_max = None
optimize_sc = False
#run phonon always, no matter ISIF=4 passed or not
force_phonon  = False
#The tolerance for phonon stable
stable_tor = 0.01
#float, the mimimum of temperature in QHA process, e.g. 5
t_min = 5
#float, the maximum of temperature in QHA process, e.g. 2000
t_max = 2000
#float, the step of temperature in QHA process, e.g. 5
t_step = 5
#float, acceptable value for average RMS, recommend >= 0.005
eos_tolerance = 0.01
#str, the vasp command, if None then find in the FWorker configuration
vasp_cmd = None
#dict, metadata to be included, this parameter is useful for filter the data, e.g. metadata={"phase": "BCC_A2", "tag": "AFM"}
metadata = None
#It is for RobustOptimizeFW, if run ISIF=4 followed ISIF=7
isif4 = False
#The level for robust optimization
level = 1
#float, the tolerannce for symmetry, e.g. 0.05
symmetry_tolerance = 0.05
#The symmetry tolerance, including three keys, 
#e.g. override_symmetry_tolerances={'tol_strain': 0.05, 'tol_energy': 0.025, 'tol_bond': 0.10}
override_symmetry_tolerances = None
#Global settings for all vasp job, e.g.
#override_default_vasp_params = {'user_incar_settings': {}, 'user_kpoints_settings': {}, 'user_potcar_functional': str}
#If some value in 'user_incar_settings' is set to None, it will use vasp's default value
override_default_vasp_params = {}
#dict, dict of class ModifyIncar with keywords in Workflow name. e.g.
"""
modify_incar_params = { 'Full relax': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                        'PreStatic': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
                        'PS2': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}}, 
                        'static': {'incar_update': {"LAECHG":False,"LCHARG":False,"LWAVE":False}},
"""
modify_incar_params = {}
#dict, dict of class ModifyKpoints with keywords in Workflow name, similar with modify_incar_params
modify_kpoints_params = {}
#bool, print(True) or not(False) some informations, used for debug
verbose = False

###################### DO NOT CHANGE THE FOLLOWING LINES ##############################
from pymatgen import Structure
from dfttk.wflows import get_wf_gibbs_robust
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn, dumpfn
from fireworks import LaunchPad

structure = Structure.from_file(TEMPLATE_STRUCTURE_FILENAME)

if magmom:
    structure.add_site_property('magmom', magmom)

if not db_file:    
    db_file = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]

wf = get_wf_gibbs_robust(structure, num_deformations=num_deformations, deformation_fraction=deformation_fraction,
         phonon=phonon, phonon_supercell_matrix=phonon_supercell_matrix, t_min=t_min, t_max=t_max, t_step=t_step, 
         eos_tolerance=eos_tolerance, volume_spacing_min=volume_spacing_min, vasp_cmd=vasp_cmd, db_file=db_file, 
         isif4=isif4, metadata=metadata, name='EV_QHA', override_symmetry_tolerances=override_symmetry_tolerances,
         override_default_vasp_params=override_default_vasp_params, modify_incar_params=modify_incar_params,
         modify_kpoints_params=modify_kpoints_params, verbose=verbose, phonon_supercell_matrix_min=phonon_supercell_matrix_min,
         phonon_supercell_matrix_max=phonon_supercell_matrix_max, optimize_sc=optimize_sc, level=level,
         force_phonon=force_phonon, stable_tor=stable_tor)

if not DRY_RUN:
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
if WRITE_OUT_WF:
    dumpfn(wf.to_dict(), 'dfttk_wf.yaml')
