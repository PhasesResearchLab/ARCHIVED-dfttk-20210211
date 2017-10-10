from pymatgen import Structure
from pymatgen.io.vasp.inputs import Incar
from fireworks import FWorker, Workflow, LaunchPad
from fireworks.core.rocket_launcher import launch_rocket
from prlworkflows.prl_fireworks import FullOptFW
from prlworkflows.input_sets import PRLRelaxSet
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
# TODO: make FWORKER a fixture and clean the db with every function call
FWORKER = FWorker(env={"db_file": os.path.join(MODULE_DIR, "db.json")})

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


def test_full_opt_fw_writes_correct_fw_for_UIS_in_set_constructor(launch_dir, lpad):
    s = PRLRelaxSet(STRUCT, user_incar_settings={'ISIF': 4})
    fw = FullOptFW(STRUCT, vasp_input_set=s, vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=FWORKER)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 4}
    assert all([incar[k] == v for k, v in desired_parameters.items()])


def test_full_opt_fw_writes_isif_setting_takes_effect(launch_dir, lpad):
    fw = FullOptFW(STRUCT, isif=7, vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=FWORKER)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 7}
    assert all([incar[k] == v for k, v in desired_parameters.items()])


def test_full_opt_fw_writes_isif_setting_does_not_take_effect_with_VIS(launch_dir, lpad):
    s = PRLRelaxSet(STRUCT)
    fw = FullOptFW(STRUCT, vasp_input_set=s, isif=7, vasp_cmd=None)
    wf = Workflow([fw])
    lpad.add_wf(wf)
    launch_rocket(lpad, fworker=FWORKER)
    incar = Incar.from_file(os.path.join(launch_dir, 'INCAR.gz'))
    desired_parameters = {'ISIF': 3}
    assert all([incar[k] == v for k, v in desired_parameters.items()])

