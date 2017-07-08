"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.

The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

from pymatgen import Structure


def substitute_sqs(structure, subl_model):
    """Takes a single abstract SQS and makes a concrete Structure.

    Parameters
    ----------
    structure : Structure
        Abstract SQS.
    subl_model : [[str]]
        List of strings of species names. **Note that order does matter!**
        [["Al", "Fe"]] and [["Fe", "Al"]] will produce different results!

    Returns
    -------
    Structure
        A concrete SQS structure with the
    """
    pass


def enumerate_sqs(structure, subl_model, endmembers=True):
    """
    Return a list of all of the concrete Structure objects from an abstract Structure and concrete sublattice model.

    Parameters
    ----------
    structure : Structure
        Abstract SQS.
    subl_model : [[str]]
        List of strings of species names, in the style of ESPEI `input.json`. This sublattice model
        can be of higher dimension than the SQS, e.g. a [["Al", "Fe", "Ni"]] for a fcc 75/25 binary SQS
        will generate the following Structures:
        Al0.75Fe0.25, Al0.75Ni0.25
        Fe0.75Al0.25, Fe0.75Ni0.25
        Ni0.75Al0.25, Ni0.75Fe0.25
        *Note that the ordering of species the sublattice model does not matter!*
    endmembers: bool
        Include endmembers in the enumerated structures if True. Defaults to True.

    Returns
    -------
    [Structure]
        List of all Structure objects that can be created from the sublattice model
    """
    pass
