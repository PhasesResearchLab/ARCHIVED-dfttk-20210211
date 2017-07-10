"""
This module tests the functionality of the sqs module.

The tests are focused on conversion of structures from the ATAT format to
pymatgen Structures and matching those Structure objects to ESPEI-style
sublattice models and symmetry from user-input.
"""

import numpy as np
import pytest
from pymatgen import Lattice

from prlworkflows.sqs import SQS, enumerate_sqs
from prlworkflows.sqs_db import lat_in_to_sqs


ATAT_FCC_L12_LATTICE_IN = """1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
-1.000000 1.000000 -1.000000
1.000000 -1.000000 -1.000000
-2.000000 -2.000000 0.000000
-1.000000 -2.000000 -1.000000 a_A
-2.000000 -1.000000 -1.000000 a_A
-1.000000 -1.000000 -1.000000 a_B
-0.000000 -1.000000 -1.000000 a_B
-1.000000 -0.000000 -1.000000 a_B
-2.000000 -2.000000 -1.000000 a_B
-1.000000 -1.000000 -2.000000 a_A
-2.000000 -2.000000 -2.000000 a_A
-1.000000 -1.500000 -0.500000 c_A
-1.000000 -1.500000 -1.500000 c_A
-1.000000 -0.500000 -0.500000 c_A
-0.000000 -0.500000 -0.500000 c_A
-0.000000 -0.500000 -1.500000 c_A
-2.000000 -1.500000 -0.500000 c_A
-1.000000 -0.500000 -1.500000 c_A
-2.000000 -1.500000 -1.500000 c_A
-1.500000 -1.000000 -1.500000 c_A
-1.500000 -1.000000 -0.500000 c_A
-0.500000 -1.000000 -0.500000 c_A
-0.500000 0.000000 -1.500000 c_A
-0.500000 -0.000000 -0.500000 c_A
-1.500000 -2.000000 -0.500000 c_A
-0.500000 -1.000000 -1.500000 c_A
-1.500000 -2.000000 -1.500000 c_A
-0.500000 -1.500000 -1.000000 c_A
-1.500000 -0.500000 -1.000000 c_A
-0.500000 -0.500000 -1.000000 c_A
-1.500000 -2.500000 -1.000000 c_A
-2.500000 -1.500000 -1.000000 c_A
-1.500000 -1.500000 -1.000000 c_A
-0.500000 -0.500000 -2.000000 c_A
-1.500000 -1.500000 -2.000000 c_A"""

ATAT_ROCKSALT_B1_LATTICE_IN = """1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
2.000000 1.000000 1.000000
-2.000000 1.000000 1.000000
-0.000000 -0.500000 0.500000
-0.500000 1.500000 2.000000 a_A
-0.000000 1.000000 2.000000 a_A
-0.000000 1.500000 2.500000 a_A
0.500000 1.000000 1.500000 a_B
0.500000 1.500000 2.000000 a_B
-1.500000 1.000000 1.500000 a_A
1.000000 0.500000 1.500000 a_A
1.000000 1.000000 2.000000 a_A
-1.000000 0.500000 1.500000 a_B
-1.000000 1.000000 2.000000 a_B
1.500000 1.000000 1.500000 a_B
-0.500000 0.500000 1.000000 a_B
-0.500000 1.000000 1.500000 a_A
-0.000000 -0.000000 1.000000 a_B
-0.000000 0.500000 1.500000 a_B
0.500000 0.500000 1.000000 a_A
-1.000000 1.000000 1.500000 b_B
-0.500000 0.500000 1.500000 b_A
-0.500000 1.000000 2.000000 b_A
0.000000 0.500000 1.000000 b_A
0.000000 1.000000 1.500000 b_B
0.000000 1.500000 2.000000 b_B
0.500000 -0.000000 1.000000 b_B
0.500000 0.500000 1.500000 b_A
0.500000 1.000000 2.000000 b_A
-1.500000 0.500000 1.500000 b_B
1.000000 0.500000 1.000000 b_A
1.000000 1.000000 1.500000 b_A
-1.000000 0.500000 1.000000 b_B
1.500000 0.500000 1.500000 b_B
-0.500000 -0.000000 1.000000 b_B
-0.000000 -0.000000 0.500000 b_A"""


# noinspection PyProtectedMember,PyProtectedMember
def test_atat_bestsqs_is_correctly_parsed_to_sqs():
    """lattice.in files in the ATAT format should be converted to SQS correctly."""
    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    specie_types = {specie.symbol for specie in structure.types_of_specie}
    assert specie_types == {'Xaa', 'Xab', 'Xca'}
    assert np.all(structure.sublattice_model == [['a', 'b'], ['a']])
    assert np.all(structure.sublattice_site_ratios == [8, 24])
    assert np.all(structure._sublattice_names == ['a', 'c'])
    assert structure.is_abstract

    structure = lat_in_to_sqs(ATAT_ROCKSALT_B1_LATTICE_IN)
    specie_types = {specie.symbol for specie in structure.types_of_specie}
    assert specie_types == {'Xaa', 'Xab', 'Xba', 'Xbb'}
    assert np.all(structure.sublattice_model == [['a', 'b'], ['a', 'b']])
    assert np.all(structure.sublattice_site_ratios == [16, 16])
    assert np.all(structure._sublattice_names == ['a', 'b'])
    assert structure.is_abstract


def test_sqs_obj_correctly_serialized():
    """Tests that the as_dict method of the SQS object correctly includes metadata and is able to be seralized/unserialized."""
    sqs = SQS(Lattice.cubic(5), ['Xaa', 'Xab'], [[0,0,0],[0.5,0.5,0.5]],
              sublattice_model=[['a', 'b']],
              sublattice_names=['a'],
              sublattice_site_ratios=[1])
    assert sqs.is_abstract

    # first seralization
    s1 = SQS.from_dict(sqs.as_dict())
    assert sqs == s1
    assert s1.is_abstract
    assert s1.sublattice_model == [['a', 'b']]
    assert s1._sublattice_names == ['a']
    assert s1.sublattice_site_ratios == [1]

    # second serialization
    s2 = SQS.from_dict(sqs.as_dict())
    assert sqs == s2
    assert s2.is_abstract
    assert s2.sublattice_model == [['a', 'b']]
    assert s2._sublattice_names == ['a']
    assert s2.sublattice_site_ratios == [1]

    # test that we can make it concrete
    s2.make_concrete([['Fe', 'Ni']])
    assert {s.symbol for s in s2.types_of_specie} == {'Fe', 'Ni'}



def test_higher_order_sqs_list_from_database():
    """List of SQS objects that match the criteria should be extracted from the database.

    This tests that phases with multicomponent (3+) sublattices are matched with SQS for the lower order subsystems.
    """
    raise NotImplementedError


def test_multiple_sqs_list_from_database():
    """List of SQS objects that match the criteria should be extracted from the database.

    This tests that phases with multiple solution sublattices can match different SQS that describe
    each sublattice.
    """
    raise NotImplementedError


def test_abstract_sqs_is_properly_substituted_with_sublattice_model():
    """Test that an abstract SQS can correctly be make concrete."""
    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)

    structure.make_concrete([['Fe', 'Ni'], ['Al']])
    assert {s.symbol for s in structure.types_of_specie} == {'Al', 'Fe', 'Ni'}
    assert structure.is_abstract is False

    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    structure.make_concrete([['Al', 'Al'], ['Al']])
    assert {s.symbol for s in structure.types_of_specie} == {'Al'}
    assert structure.is_abstract is False
    with pytest.raises(ValueError):
        structure.make_concrete([['Fe', 'Fe'], ['Fe']])


def test_abstract_sqs_scales_volume_when_made_concrete():
    """SQS should scale in volume by default, but optionally not when made concrete"""

    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    structure.make_concrete([['Fe', 'Ni'], ['Al']])
    assert np.isclose(structure.volume, 421.08774505083824)

    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    structure.make_concrete([['Fe', 'Ni'], ['Al']], scale_volume=False)
    assert np.isclose(structure.volume, 8.0)


def test_sqs_is_properly_enumerated_for_a_higher_order_sublattice_model():
    """Tests that a sublattice model of higher order than an SQS properly enumerated"""
    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    structures = enumerate_sqs(structure, [['Al', 'Ni'], ['Fe', 'Cr']], endmembers=False)
    assert len(structures) == 4

    structure = lat_in_to_sqs(ATAT_ROCKSALT_B1_LATTICE_IN)
    structures = enumerate_sqs(structure, [['Al', 'Ni'], ['Fe', 'Cr']], endmembers=True)
    assert len(structures) == 16

def test_sqs_is_properly_enumerated_for_a_multiple_solution_sublattice_model():
    """Tests that a sublattice model with multiple solution sublattices is properly enumerated"""
    structure = lat_in_to_sqs(ATAT_ROCKSALT_B1_LATTICE_IN)
    structures = enumerate_sqs(structure, [['Al', 'Ni'], ['Fe', 'Cr']], endmembers=False)
    assert len(structures) == 4

    structure = lat_in_to_sqs(ATAT_ROCKSALT_B1_LATTICE_IN)
    structures = enumerate_sqs(structure, [['Al', 'Ni'], ['Fe', 'Cr']], endmembers=True)
    assert len(structures) == 16
    assert all([isinstance(s, SQS) for s in structures])
    assert all([not s.is_abstract for s in structures])


def test_enumerating_sqs_with_lower_order_subl_raises():
    """If a lower order sublattice model is passed be enumerated in an SQS, it should raise."""
    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    with pytest.raises(ValueError):
        enumerate_sqs(structure, [['Fe'], ['Al']])
