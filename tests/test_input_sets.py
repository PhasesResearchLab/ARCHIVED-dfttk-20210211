#!python
#
#
from dfttk.input_sets import PreStaticSet, RelaxSet, StaticSet, ForceConstantsSet
from pymatgen import Structure


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
structure = Structure.from_str(POSCAR_STR, fmt='POSCAR')

vis_relax = RelaxSet(structure)
print(vis_relax)