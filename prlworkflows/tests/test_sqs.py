"""
This module tests the functionality of the sqs module.

The tests are focused on conversion of structures from the ATAT format to
pymatgen Structures and matching those Structure objects to ESPEI-style
sublattice models and symmetry from user-input.
"""

import numpy as np
import pytest

from prlworkflows.sqs import SQS
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
-1.500000 -1.500000 -2.000000 c_A
"""


def test_atat_bestsqs_is_correctly_parsed_to_sqs():
    """lattice.in files in the ATAT format should be converted to SQS correctly."""
    structure = lat_in_to_sqs(ATAT_FCC_L12_LATTICE_IN)
    specie_types = {specie.symbol for specie in structure.types_of_specie}
    assert specie_types == {'Xaa', 'Xab', 'Xca'}
    assert np.all(structure.sublattice_model == [['a', 'b'], ['a']])
    assert np.all(structure.sublattice_site_ratios == [8, 24])
    assert np.all(structure._sublattice_names == ['a', 'c'])
    assert structure.is_abstract


def test_sqs_obj_correctly_serialized():
    """Tests that the as_dict method of the SQS object correctly includes metadata"""
    raise NotImplementedError


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


def test_sqs_is_properly_enumerated_for_a_higher_order_sublattice_model():
    """Tests that a sublattice model of higher order than an SQS properly enumerated"""
    raise NotImplementedError


def test_sqs_is_properly_enumerated_for_a_multiple_solution_sublattice_model():
    """Tests that a sublattice model with multiple solution sublattices is properly enumerated"""
    raise NotImplementedError
