"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.

The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

import itertools
import copy

import numpy as np
from pymatgen import Structure


# TODO: implement making the concrete structure abstract again
class SQS(Structure):
    """A pymatgen Structure with special features for SQS.
    """

    def __init__(self, *args, **kwargs):
        """Create a SQS object

        Parameters
        ----------
        args :
            args to pass to Structure
        sublattice_model : np.ndarray
            Abstract sublattice model in the ESPEI style, e.g. `[['a', 'b'], ['a']]`.
        sublattice_names : np.ndarray
            Names of the sublattices, or the second character in the species names, e.g. `['a', 'c']`.
        sublattice_site_ratios : np.ndarray
            Site ratios of the sublattices, e.g. `[8.0, 24.0]`
        kwargs :
            kwargs to pass to Structure
        """
        self.sublattice_model = kwargs.pop('sublattice_model', None)
        self._sublattice_names = kwargs.pop('sublattice_names', None)
        self.sublattice_site_ratios = kwargs.pop('sublattice_site_ratios', None)
        super(SQS, self).__init__(*args, **kwargs)

    @property
    def is_abstract(self):
        return all([specie.symbol.startswith('X') for specie in self.types_of_specie])

    def make_concrete(self, subl_model, scale_volume=True):
        """Modify self to be a concrete SQS based on the sublattice model.

        Parameters
        ----------
        subl_model : [[str]]
            List of strings of species names. Must exactly match the shape of self.sublattice_model.
            **Note that order does matter!** [["Al", "Fe"]] and [["Fe", "Al"]] will produce different results!
        scale_volume : bool
            If True, scales the volume of the cell so the ions have at least their minimum atomic radii between them.
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
        if scale_volume:
            idx = np.nonzero(self.distance_matrix)  # exclude distance from self
            # construct a minimal distance matrix based on average ionic radii
            radius_matrix = np.zeros((len(self.sites), len(self.sites)))
            for i in range(len(self.sites)):
                for j in range(len(self.sites)):
                    radius_matrix[i, j] = (
                        self.sites[i].specie.data['Atomic radius'] + self.sites[j].specie.data['Atomic radius'])
            # find the scale factor and scale the lattice
            sf = np.max(radius_matrix[idx] / self.distance_matrix[idx]) ** 3
            self.scale_lattice(self.volume * sf)

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(SQS, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_model'] = self.sublattice_model
        d['sublattice_names'] = self._sublattice_names
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        sqs = super(SQS, cls).from_dict(d, fmt=fmt)
        sqs.sublattice_model = d.get('sublattice_model')
        sqs._sublattice_names = d.get('sublattice_names')
        sqs.sublattice_site_ratios = d.get('sublattice_site_ratios')
        return sqs


def enumerate_sqs(structure, subl_model, endmembers=True, scale_volume=True):
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
    endmembers : bool
        Include endmembers in the enumerated structures if True. Defaults to True.
    scale_volume : bool
        If True, scales the volume of the cell so the ions have at least their minimum atomic radii between them.

    Returns
    -------
    [SQS]
        List of all concrete SQS objects that can be created from the sublattice model.
    """
    # error checking
    if len(subl_model) != len(structure.sublattice_model):
        raise ValueError('Passed sublattice model ({}) does not agree with the passed structure ({})'.format(subl_model, structure.sublattice_model))
    if np.any([len(subl) for subl in subl_model] < [len(subl) for subl in structure.sublattice_model]):
        raise ValueError('The passed sublattice model ({}) is of lower order than the passed structure supports ({})'.format(subl_model, structure.sublattice_model))
    possible_subls = []
    for subl, abstract_subl in zip(subl_model, structure.sublattice_model):
        if endmembers:
            subls = itertools.product(subl, repeat=len(abstract_subl))
        else:
            subls = itertools.permutations(subl, r=len(abstract_subl))
        possible_subls.append(subls)
    unique_subl_models = itertools.product(*possible_subls)
    # create a list of concrete structures with the generated sublattice models
    structs = []
    for model in unique_subl_models:
        s = copy.deepcopy(structure)
        s.make_concrete(model, scale_volume)
        structs.append(s)
    return structs
