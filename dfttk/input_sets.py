import six
import os
import subprocess

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.sets import DictSet, get_vasprun_outcar, get_structure_from_prev_run, _load_yaml_config


# We set potentials back to the VASP recommended, since we primarily don't work on oxides.
POTCAR_UPDATES = {
        'Be': 'Be',  # 2 electrons, default was Be_sv (4 electrons)
        'Cu': 'Cu',  # 11 electrons, default was Cu_pv (17 electrons)
        'Fe': 'Fe',  # 8 electrons, default was Fe_pv (14 electrons)
        'Mg': 'Mg',  # 2 electrons, default was Mg_pv (8 electrons)
        'Ni': 'Ni',  # 10 electrons, default was Ni_pv (16 electrons)
        'Mo': 'Mo_sv',  # 14 electrons, default was Mo_pv (12 electrons)
        'Nb': 'Nb_sv',  # 13 electrons, default was Nb_pv (11 electrons)
        'Os': 'Os',  # 8 electrons, default was Os_pv (14 electrons)
        'Re': 'Re',  # 7 electrons, default was Re_pv (13 electrons)
        'Ti': 'Ti_sv',  # 12 electrons, default Ti_pv (10 electrons)
        'V': 'V_sv',  # 13 electrons, default V_pv (11 electrons)
    }

class RelaxSet(DictSet):
    """
    Set for performing relaxations.

    The smearing must be set to give the correct forces.
    The default is tuned for metal relaxations.
    Kpoints have a 8000 kpoints per reciprocal atom default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    # we never are comparing relaxations, only using them for optimizing structures.
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-5,
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'PREC': 'Accurate',
        'ALGO': 'NORMAL',
        'LWAVE': False,
        'LCHARG': False,
        'ISIF': 3,
        "ICHARG": 2,
        'ENCUT': 520,
    })
    # now we reset the potentials
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, volume_relax=False, **kwargs):
        """If volume relax is True, will do volume only, ISIF 7"""
        self.kwargs = kwargs
        self.volume_relax = volume_relax
        super(RelaxSet, self).__init__(structure, RelaxSet.CONFIG, **kwargs)
        if self.volume_relax:
            self.incar.update({'ISIF': 7})


class ForceConstantsSet(DictSet):
    """
    Set for calculating force constants calculations.

    Force constants are calculated by the finite difference method with symmetry considered.

    The smearing must be set to give the correct forces.
    The default is tuned for metals.

    Kpoints have a 8000 kpoints per reciprocal atom default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    # we never are comparing relaxations, only using them for optimizing structures.
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-6,
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': False,
        'ISIF': 0,  # only calculate the forces, stress tensor is not needed
        'IBRION': 6,  # calculate force constants by finite differences with symmetry
        'POTIM': 0.015,  # displacement distance
        'NFREE': 2,  # how many displacments to do. 2 gives +POTIM and -POTIM
        'NSW': 1,  # backwards compatibility setting
        'PREC': 'Accurate',
        'ALGO': 'NORMAL',
        'SYMPREC': 1e-4,  # some supercells seem to have issues with primcel VASP algorithm
        "ICHARG": 2,
    })
    # now we reset the potentials
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        self.kwargs = kwargs
        super(ForceConstantsSet, self).__init__(
            structure, ForceConstantsSet.CONFIG, **kwargs)


class StaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-6,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': -5,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": True,
        "LWAVE": False,
        "LORBIT": 11,
        "LVHAR": True,
        "ICHARG": 2,
        "NEDOS": 5001,
    })
    # now we reset the potentials
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = kwargs
        super(StaticSet, self).__init__(structure, StaticSet.CONFIG, **kwargs)


class ATATIDSet():
    """Set tuned for Inflection Detection runs using ATAT with correct smearing for metals.
    Kpoints have a 8000 reciprocal density default.

    Overrides write_input to write the INCAR, KPPRA, USEPOT and DOSTATIC to the vasp.wrap instead.
    """

    def __init__(self, structure, grid_density=8000):
        self.structure = structure
        self.grid_density = grid_density


    def write_input(self, output_dir):
        """Write vasp.wrap and generate the other commands with robustrelax_vasp -mk"""
        # TODO: handle magmoms
        EDIFF_PER_ATOM = 1e-6
        EDIFF = len(self.structure)*EDIFF_PER_ATOM
        vasp_wrap = """[INCAR]
        EDIFF = {0}
        PREC = Accurate
        ALGO = Fast
        ENCUT = 520
        ISMEAR = 1
        SIGMA = 0.2
        IBRION = -1
        NSW = 1
        ISPIN = 2
        NELMIN = 4
        ISIF = 2
        LREAL = False
        ISYM = 0
        ICHARG = 1
        ISTART = 2
        USEPOT = PAWPBE
        KPPRA = {1}
        """.format(EDIFF, self.grid_density)
        with open(os.path.join(output_dir, 'vaspid.wrap'), 'w') as fp:
            fp.write(vasp_wrap)

class ForcesSet(DictSet):
    """Set tuned for generic force calculations (Gaussian smearing).
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].pop('ENCUT',None) 
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-8,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LCHARG": False,
        "LORBIT": 11,
        "LVHAR": True,
        "LWAVE": False,
        "ICHARG": 2,
        "NEDOS": 5001,
    })
    # now we reset the potentials
    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = kwargs
        super(ForcesSet, self).__init__(structure, ForcesSet.CONFIG, **kwargs)


