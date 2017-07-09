"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.

The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

from pymatgen import Structure


# TODO: override the upstream constructors to call super and then add the sublattice configuration
# TODO: override the as_dict method to add the extra metadata, if necessary.
# TODO: implement making the concrete structure abstract again
class SQS(Structure):
    """A pymatgen Structure with special features for SQS.
    """

    def __init__(self, *args, **kwargs):
        # TODO: add support for adding sublattice model metadata (model, site ratios, symmetry, version info)
        # TODO: check for any DummySpecies and set is_abstract based on the result
        super(SQS, self).__init__(*args, **kwargs)
        self.sublattice_model = None
        self._sublattice_names = None
        self.sublattice_site_ratios = None

    @property
    def is_abstract(self):
        return all([specie.symbol.startswith('X') for specie in self.types_of_specie])

    def make_concrete(self, subl_model):
        """Modify self to be a concrete SQS based on the sublattice model.

        Parameters
        ----------
        subl_model : [[str]]
            List of strings of species names. Must exactly match the shape of self.sublattice_model.
            **Note that order does matter!** [["Al", "Fe"]] and [["Fe", "Al"]] will produce different results!

        """
        def _subl_error():
            raise ValueError('Concrete sublattice model {} does not match size of abstract sublattice model {}'.format(subl_model, self.sublattice_model))
        if not self.is_abstract:
            raise ValueError('SQS cannot be made concrete because it already is concrete with species.'.format({s.symbol for s in self.types_of_specie}))
        if len(subl_model) != len(self.sublattice_model):
            _subl_error()
        # build the replacement dictionary
        # we have to look up the sublattice names to build the replacement species names
        replacement_dict = {}
        for abstract_subl, concrete_subl, subl_name in zip(self.sublattice_model, subl_model, self._sublattice_names):
            if len(abstract_subl) != len(concrete_subl):
                _subl_error()
            for abstract_specie, concrete_specie in zip(abstract_subl, concrete_subl):
                specie = 'X' + subl_name + abstract_specie
                replacement_dict[specie] = concrete_specie
        self.replace_species(replacement_dict)


def enumerate_sqs(structure, subl_model, endmembers=True):
    """
    Return a list of all of the concrete Structure objects from an abstract Structure and concrete sublattice model.

    Parameters
    ----------
    structure : SQS
        SQS object. Must be abstract.
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
    [SQS]
        List of all comcrete SQS objects that can be created from the sublattice model.
    """
    pass
