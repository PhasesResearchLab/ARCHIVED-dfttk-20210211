from pymatgen.io.vasp.sets import DictSet, _load_yaml_config

class PRLRelaxSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 6000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].update({
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
    })

    def __init__(self, structure, **kwargs):
        self.kwargs = kwargs
        isif = self.kwargs.pop('isif', 3)
        super(PRLRelaxSet, self).__init__(
            structure, PRLRelaxSet.CONFIG, **kwargs)
        self._config_dict['INCAR']['ISIF'] = isif
