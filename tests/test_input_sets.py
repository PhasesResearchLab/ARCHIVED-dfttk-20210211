#!python
#
#
#from dfttk.input_sets import PreStaticSet, RelaxSet, StaticSet, ForceConstantsSet
from pymatgen import Structure
from dfttk.ftasks import WriteVaspFromIOSetPrevStructure
from dfttk.input_sets import RelaxSet
from copy import deepcopy

struct = Structure.from_file('POSCAR')

"""
print(struct)

struct.to(fmt="poscar", filename="poscar-test")
"""

force_gamma = True
override_default_vasp_params = None
vasp_input_set = None

t = []
site_properties = deepcopy(struct).site_properties
print(site_properties)
vasp_input_set = vasp_input_set or RelaxSet(struct, force_gamma=force_gamma)

t.append(WriteVaspFromIOSetPrevStructure(vasp_input_set=vasp_input_set, site_properties=site_properties))
print(t)

print(struct)
scal_f = struct.volume * 1.01 ** 3
print(scal_f)
print(type(scal_f))
struct.scale_lattice(struct.volume * 1.01 ** 3)
print(struct)