from pymatgen.io.vasp.sets import DictSet, _load_yaml_config

class PRLRelaxSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['INCAR'].pop('ENCUT')
    CONFIG['KPOINTS'].update({
        'grid_density': 6000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].update({
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'PREC': 'HIGH',
        'ALGO': 'NORMAL'
    })

    def __init__(self, structure, **kwargs):
        self.kwargs = kwargs
        super(PRLRelaxSet, self).__init__(
            structure, PRLRelaxSet.CONFIG, **kwargs)

class PRLStaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['INCAR'].pop('ENCUT')
    CONFIG['KPOINTS'].update({
        'grid_density': 6000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'ISMEAR': -5,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ENCUT': 520,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": True,
        "LORBIT": 11,
        "LVHAR": True,
        "LWAVE": False,
        "ICHARG": 0,
    })

    def __init__(self, structure, **kwargs):
        self.kwargs = kwargs
        super(PRLStaticSet, self).__init__(
            structure, PRLStaticSet.CONFIG, **kwargs)
