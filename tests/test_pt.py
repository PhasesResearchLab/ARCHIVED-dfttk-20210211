#!python
#

from pymatgen.core.periodic_table import Element
import sys
import os
dfttkhome = os.path.abspath(os.path.join('..'))
sys.path.append(dfttkhome)
from dfttk.structure_builders.substitutions import get_density_from_pt


print(get_density_from_pt(['Fe', 'Ni']))

Fe = Element("Ti").density_of_solid
print(Fe)
print(type(Fe))
print(float(Fe))
print(type(float(Fe)))