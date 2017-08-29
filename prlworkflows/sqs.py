"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.

The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

from __future__ import division

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
        sublattice_model : [[str]]
            Abstract sublattice model in the ESPEI style, e.g. `[['a', 'b'], ['a']]`.
        sublattice_names : [[str]]
            Names of the sublattices, or the second character in the species names, e.g. `['a', 'c']`.
        kwargs :
            kwargs to pass to Structure
        """
        self.sublattice_model = kwargs.pop('sublattice_model', None)
        self._sublattice_names = kwargs.pop('sublattice_names', None)
        super(SQS, self).__init__(*args, **kwargs)

    @property
    def espei_sublattice_model(self):
        """
        Return ESPEI-formatted sublattice model [['a', 'b'], 'a']
        """
        # short function to convert [['A', 'B'], ['A']] to [['A', 'B'], 'A'] as in ESPEI format
        canonicalize_sublattice = lambda sl: sl[0] if len(sl) == 1 else sl
        return [canonicalize_sublattice(sl) for sl in self.sublattice_model]

    @property
    def is_abstract(self):
        return all([specie.symbol.startswith('X') for specie in self.types_of_specie])

    @property
    def normalized_sublattice_site_ratios(self):
        """Return normalized sublattice site ratio. E.g. [[0.25, 0.25], [0.1666, 0.1666, 0.1666]]
        """
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = self.composition.as_dict()
        site_ratios = [[comp_dict['X'+name+e+'0+']/self.num_sites for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

    @property
    def sublattice_site_ratios(self):
        """Return normalized sublattice site ratio. E.g. [[0.25, 0.25], [0.1666, 0.1666, 0.1666]]
        """
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = {k: int(v) for k, v in self.composition.reduced_composition.as_dict().items()}
        site_ratios = [[comp_dict['X'+name+e+'0+'] for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

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
            raise ValueError('SQS cannot be made concrete because it already is concrete with species {}.'.format({s.symbol for s in self.types_of_specie}))
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
        if scale_volume and not self.is_abstract:
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

    def get_endmember_space_group_info(self, symprec=1e-2, angle_tolerance=5.0):
        """
        Return endmember space group info..

        Args:
            symprec (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 1e-2.
            angle_tolerance (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 5 degrees.

        Returns:
            spacegroup_symbol, international_number
        """
        endmember_subl = [['X' + subl_name for _ in subl] for subl, subl_name in
                     zip(self.sublattice_model, self._sublattice_names)]
        self_copy = copy.deepcopy(self)
        self_copy.make_concrete(endmember_subl)
        endmember_space_group_info = self_copy.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)
        del self_copy
        return endmember_space_group_info

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(SQS, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_model'] = self.sublattice_model
        d['sublattice_names'] = self._sublattice_names
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        endmember_symmetry = self.get_endmember_space_group_info()
        d['symmetry'] = {'symbol': endmember_symmetry[0], 'number': endmember_symmetry[1]}
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        sqs = super(SQS, cls).from_dict(d, fmt=fmt)
        sqs.sublattice_model = d.get('sublattice_model')
        sqs._sublattice_names = d.get('sublattice_names')
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
