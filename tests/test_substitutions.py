#!python
#

#import pytest
import sys
import os
dfttkhome = os.path.abspath(os.path.join('..'))
sys.path.append(dfttkhome)
import dfttk.structure_builders.substitutions as substitutions
from pymatgen import Structure


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
stuct = Structure.from_str(POSCAR_STR, fmt='POSCAR')


def test_canonicalize_config():
    configuration = [['Ni', 'Fe'], ['Fe']]
    occupancies = [[0.75, 0.25], [1.0]]
    (new_configuration, new_occupancies) = substitutions.canonicalize_config(configuration, occupancies)
    assert new_configuration == [['Fe', 'Ni'], ['Cr', 'Fe', 'Ni']]
    assert new_occupancies == [[0.25, 0.75], [0.2, 0.1, 0.7]]

def test_get_ele_list_from_struct():
    ele_list = get_ele_list_from_struct(struct)
    assert ele_list == ['Fe', 'Ni', 'Ni', 'Ni']

def test_get_density_from_pt():
    den_dict = {'Nb': 8.57, 'Ti': 4.507}
    #test for the list of elements
    ele_list = ['Nb', 'Ti']
    density_dict = get_density_from_pt(ele_list)
    for ele in ele_list:
        assert den_dict[ele] == density_dict[ele]
    #test for the dict of tlements
    ele_dict = {'Nb': 3, 'Ti': 1}
    density_dict = get_density_from_pt(ele_dict)
    for ele in ele_dict:
        assert den_dict[ele] == density_dict[ele]

def test_scale_struct():
    struct = scale_struct(struct)
    assert struct.lattice.a == 3.5443397446212437

def test_gen_replacement_dict():
    old_config = [['Fe', 'Cr'], 'Ni']
    new_config = [['Ti', 'V'], 'Zr']
    replacement_dict = gen_replacement_dict(old_config, new_config)
    assert replacement_dict['Fe'] == 'Ti'
    assert replacement_dict['Cr'] == 'V'
    assert replacement_dict['Ni'] == 'Zr'