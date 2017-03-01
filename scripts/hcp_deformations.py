import numpy as np
from pymatgen import MPRester
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from Fireworks import LaunchPad

API_KEY = None # If None, will use the one from your ~/.pmgrc.yaml. Otherwise string of API_KEY
LAUNCHPAD_FILE = None # If None, will try to autoload from FW_CONFIG_FILE environment variable. Otherwise string of file path
MP_ID = 'mp-72' # materials project id of the structure. mp-72 is the P6/mmm Ti
# mp-46 is P6_3/mmc 
# mp-6985 is Fm-3m
# mp-73 is Im-3m

number_of_deformations = 5 #number of deformations in a and c directions
# for cubic systems, you'd have to do the matricies slightly differently.
# see comment below

if API_KEY:
    mpr = MPRester(API_KEY)
else:
    mpr = MPRester()

if LAUNCHPAD_FILE:
    lpad = LaunchPad.from_file(LAUNCHPAD_FILE)
else:
    lpad = LaunchPad.auto_load()

# transformation matricies
# here the values in the linspace are the scaling factors for the a and c directions. E.g. -0.1 = 90% of original a
a_transformations = [np.array([[1+x, 0, 0],[0,1,0],[0,0,1]]) for x in np.linspace(-0.1, 0.1, number_of_deformations)]
c_transformations = [np.array([[1,0,0],[0,1,0],[0,0,x+1]]) for x in np.linspace(-0.1, 0.1, number_of_deformations)]

# here the values in the linspace are the scaling factor for cell volume. E.g. -0.1 = 90%
# cubic_transformations = [np.eye(3)*((x+1)**(1.0/3.0)) for x in np.linspace(-0.1, 0.1, number_of_deformations)]

# construct the workflow
deformations = np.append(a_transformations, c_transformations, axis=0) # if cubic, don't need this. Just pass cubic_transformations
name_formula = struct.composition.reduced_formula
name_spg = SpacegroupAnalyzer(struct).get_space_group_symbol()

wf = get_wf_deformations(struct, deformations, name="{}-{}-deformation".format(name_formula, name_spg))

lpad.add_wf(wf)
