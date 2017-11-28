from pymatgen import Structure
import pymatgen
from pymatgen.io.vasp.inputs import Incar
from fireworks import FWorker, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket
from prlworkflows.prl_fireworks import OptimizeFW
from prlworkflows.input_sets import PRLRelaxSet
from prlworkflows.prl_workflows import wf_gibbs_free_energy, get_wf_robust_optimization
from prlworkflows.utils import update_fws_spec
import pytest
import shutil
import os
import time

MODULE_DIR = os.path.dirname(__file__)

# TODO: does not support use_fake_vasp yet. VASP will fail.

DEBUG_MODE = False

POSCAR_STR = """Si2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si
2
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 Si"""
STRUCT = Structure.from_str(POSCAR_STR, fmt='POSCAR')
TEST_DIR = os.path.join(MODULE_DIR, 'tmp_fw_test_dir')
LPAD = LaunchPad.from_dict({'host': 'localhost', 'logdir': None, 'name': 'prlworkflows_unittest', 'password': None, 'port': 27017, 'ssl_ca_file': None, 'strm_lvl': 'DEBUG', 'user_indices': [], 'username': None, 'wf_user_indices': []})

# TODO: enable debug mode by having a launchpad that does not reset
# Can this be done by still having other tests pass?
# Should we only run one test?
# Stop on failure?
@pytest.fixture
def lpad():
    """A LaunchPad object for test instances to use. Always gives a clean (reset) LaunchPad. """
    LPAD.reset(None, require_password=False, max_reset_wo_password=5)
    yield LPAD
    LPAD.connection.close()
    return

@pytest.fixture
def launch_dir():
    test_dir = TEST_DIR + '-' + str(time.time()).split('.')[0]
    os.mkdir(test_dir)
    os.chdir(test_dir)
    yield test_dir
    os.chdir('..')
    shutil.rmtree(test_dir)
    return


@pytest.fixture
def patch_pmg_psp_dir():
    current_psp_dir = pymatgen.SETTINGS.get('PMG_VASP_PSP_DIR')
    if current_psp_dir is None:
        pymatgen.SETTINGS['PMG_VASP_PSP_DIR'] = os.path.join(MODULE_DIR, 'test_potcars')
    yield
    pymatgen.SETTINGS['PMG_VASP_PSP_DIR'] = current_psp_dir


@pytest.fixture(scope='module')
def launch_dir_debug():
    test_dir = TEST_DIR + '-' + str(time.time()).split('.')[0]
    os.mkdir(test_dir)
    os.chdir(test_dir)
    yield test_dir
    os.chdir('..')
    return

if DEBUG_MODE:
    launch_dir = launch_dir_debug


@pytest.fixture
def fworker():
    scratch_dir = os.path.join(MODULE_DIR, 'scratch_dir')
    os.mkdir(scratch_dir)
    yield FWorker(env={"db_file": os.path.join(MODULE_DIR, "db.json"), 'scratch_dir': scratch_dir})
    shutil.rmtree(scratch_dir)


def test_full_opt_fw_writes_correct_fw_for_UIS_in_set_constructor(patch_pmg_psp_dir, launch_dir, lpad, fworker):
    s = PRLRelaxSet(STRUCT, user_incar_settings={'ISIF': 4})
    fw = OptimizeFW(STRUCT, vasp_input_set=s, job_type='full_opt_run', vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=fworker)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 4}
    assert all([incar[k] == v for k, v in desired_parameters.items()])


def test_full_opt_fw_writes_isif_setting_takes_effect(patch_pmg_psp_dir, launch_dir, lpad, fworker):
    fw = OptimizeFW(STRUCT, isif=7, job_type='full_opt_run', vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=fworker)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 7}
    assert all([incar[k] == v for k, v in desired_parameters.items()])


def test_full_opt_fw_writes_isif_setting_does_take_effects_with_VIS(patch_pmg_psp_dir, launch_dir, lpad, fworker):
    s = PRLRelaxSet(STRUCT)
    fw = OptimizeFW(STRUCT, vasp_input_set=s, isif=5, job_type='full_opt_run', vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=fworker)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 5}
    assert all([incar[k] == v for k, v in desired_parameters.items()])


def test_fw_spec_modified_by_powerup():
    wf = get_wf_robust_optimization(STRUCT)
    wf = update_fws_spec(wf, {'_preserve_fworker': True})
    assert all([fw.spec['_preserve_fworker'] == True for fw in wf.fws])


def test_prl_gibbs_wf(patch_pmg_psp_dir, launch_dir, lpad, fworker):
    wf = wf_gibbs_free_energy(STRUCT, {'VASP_CMD': None})
    lpad.add_wf(wf)
    os.mkdir('scratch')
    # TODO: make this actually run by using run_no_vasp
    launch_rocket(lpad, fworker=fworker)

def test_prl_gibbs_optimization():
    without_optimize = wf_gibbs_free_energy(STRUCT, {'OPTIMIZE': False,'ROBUST': False})
    with_optimize = wf_gibbs_free_energy(STRUCT, {'OPTIMIZE': True,'ROBUST': False})
    with_robust_optimize = wf_gibbs_free_energy(STRUCT, {'OPTIMIZE': True,'ROBUST': True})
    assert len(with_optimize.fws) == len(without_optimize.fws)+1
    assert len(with_robust_optimize.fws) == len(without_optimize.fws)+3
