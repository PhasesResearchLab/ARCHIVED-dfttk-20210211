"""
This module tests the functionality of the sqs module.

The tests are focused on conversion of structures from the ATAT format to
pymatgen Structures and matching those Structure objects to ESPEI-style
sublattice models and symmetry from user-input.
"""

from prlworkflows.sqs import SQS
from prlworkflows.sqs_db import latt_in_to_cif


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

ATAT_FCC_L12_CIF = """
"""

ATAT_FCC_L12_CIF_RENAMED = """
"""


def test_atat_sqs_is_correctly_serialized_to_cif():
    """lattice.in files in the ATAT format should be converted to CIF format, optionally with element renaming."""
    assert ATAT_FCC_L12_CIF == latt_in_to_cif(ATAT_FCC_L12_LATTICE_IN)
    assert ATAT_FCC_L12_CIF_RENAMED == latt_in_to_cif(ATAT_FCC_L12_LATTICE_IN, rename=True)


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
    raise NotImplementedError
    # TODO: make sure to test that the is_abstract is set to False after


def test_sqs_is_properly_enumerated_for_a_higher_order_sublattice_model():
    """Tests that a sublattice model of higher order than an SQS properly enumerated"""
    raise NotImplementedError


def test_sqs_is_properly_enumerated_for_a_multiple_solution_sublattice_model():
    """Tests that a sublattice model with multiple solution sublattices is properly enumerated"""
    raise NotImplementedError
