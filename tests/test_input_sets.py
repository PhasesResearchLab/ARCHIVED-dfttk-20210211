#!python
#
#
#from dfttk.input_sets import PreStaticSet, RelaxSet, StaticSet, ForceConstantsSet
from pymatgen import Structure
import sys
import os
dfttkhome = os.path.abspath(os.path.join('..'))
sys.path.append(dfttkhome)
from dfttk.structure_builders.substitutions import get_density_from_pt


POSCAR_STR = """FeNi3
1.0
3.5391385555 0.0000000000 0.0000000000
0.0000000000 3.5391385555 0.0000000000
0.0000000000 0.0000000000 3.5391385555
Fe Ni
1 3
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.500000000 0.000000000
0.500000000 0.000000000 0.500000000"""
struct = Structure.from_str(POSCAR_STR, fmt='POSCAR')
print(struct.species)
for i in struct.species:
    print(str(i))

ele_dict = struct.composition.get_el_amt_dict()
print(ele_dict)
for ele in ele_dict:
    print(ele)


def scale_struct(struct):
    """Scale the structure according to the weighted average density of each element.

    Parameters
    ----------
    struct : pymatgen.Structure
    density_dict : dict
        Dictionary of {element: density}, e.g. {'Fe': 9, 'Ti': 4}. The units do
        not matter as long as the densities in the dict are internally consistent.

    Returns
    -------
    pymatgen.Structure
        Modifies the structure in place, but also returns for convenience.

    """

    species_amnt_dict = struct.composition.get_el_amt_dict()  # dict of {'V': 10.0, 'Ni': 30.0}
    density_dict = get_density_from_pt(species_amnt_dict)
    # densities is dict of densities, {'V': 6.313, 'Ni': 9.03}
    expected_density = float(sum([density_dict[species]*amnt for species, amnt in species_amnt_dict.items()]))/struct.composition.num_atoms
    current_density = struct.density
    current_volume = struct.volume
    expected_volume = current_volume/expected_density*current_density
    struct.scale_lattice(float(expected_volume))
    return struct

density_dict = {'Fe': 7.874, 'Ni': 8.908}
print(struct)
struct = scale_struct(struct)
print(struct.lattice.a)

"""
print(struct)

struct.to(fmt="poscar", filename="poscar-test")
"""


'''
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
'''