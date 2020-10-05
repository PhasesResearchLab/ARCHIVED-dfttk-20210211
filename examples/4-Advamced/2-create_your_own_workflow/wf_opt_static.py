#!python
#A work flow run OptimizeFW followed by StaticFW
#

from dfttk.fworks import OptimizeFW, StaticFW
from fireworks import Workflow, LaunchPad
from pymatgen import Structure
from uuid import uuid4
from fireworks.fw_config import config_to_dict
from monty.serialization import loadfn, dumpfn

DB_FILE = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["db_file"]
VASP_CMD = loadfn(config_to_dict()["FWORKER_LOC"])["env"]["vasp_cmd"]

def wf_opt_static(structure, vasp_cmd=None, db_file=None, metadata=None, tag=None):
    #Set the default value of input parameters
    if db_file is None:
        db_file = DB_FILE
    if vasp_cmd is None:
        vasp_cmd = VASP_CMD

    if not metadata:
        metadata = {}
    tag = metadata.get('tag', None)
    if tag is None:
        tag = str(uuid4())
        metadata['tag'] = tag

    fws = []

    #Create Optimize fireworks
    optfw = OptimizeFW(structure, isif=3, vasp_cmd=vasp_cmd, db_file=db_file, metadata=metadata, tag=tag)
    fws.append(optfw)

    #Create Static fireworks, note: we need the optimize finished and then run static, set parents=optfw
    #       In addition, we need to set the result of optimize as the input of static, set prev_calc_loc=True, parents=optfw
    staticfw = StaticFW(structure, isif=2, prev_calc_loc=True, parents=optfw, vasp_cmd=vasp_cmd, 
        db_file=db_file, metadata=metadata, tag=tag)
    fws.append(staticfw)

    #Create workflow
    wfname = "{}:{}".format(structure.composition.reduced_formula, "opt-static")
    wf = Workflow(fws, name=wfname, metadata=metadata)
    return wf

####Use of above workflow############
TEMPLATE_STRUCTURE_FILENAME = "POSCAR.Si"
#bool, submit(False) or not(True) the workflows to LaunchPad
DRY_RUN = False
#bool, write out(True) or not(False) the workflow (single workflow) to yaml file
WRITE_OUT_WF = True

structure = Structure.from_file(TEMPLATE_STRUCTURE_FILENAME)
wf = wf_opt_static(structure)

if not DRY_RUN:
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
if WRITE_OUT_WF:
    dumpfn(wf.to_dict(), 'dfttk_wf.yaml')