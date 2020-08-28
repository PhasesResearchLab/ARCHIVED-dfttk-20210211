import six
import os
import subprocess

from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.sets import DictSet, get_vasprun_outcar, get_structure_from_prev_run, _load_yaml_config


PeriodicTable = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Uut": 113, "Fl": 114, "Uup": 115, "Lv": 116, "Uus": 117, "Uuo": 118} 

magnetic_elements = [e for e in PeriodicTable.keys() if PeriodicTable[e]>=23 and PeriodicTable[e]<=28]
magnetic_elements.extend ([e for e in PeriodicTable.keys() if PeriodicTable[e]>=44 and PeriodicTable[e]<=46])
magnetic_elements.extend ([e for e in PeriodicTable.keys() if PeriodicTable[e]>=58 and PeriodicTable[e]<=71])
magnetic_elements.extend ([e for e in PeriodicTable.keys() if PeriodicTable[e]>=91 and PeriodicTable[e]<=118])

def magnetic_check(els):
    return any(item in magnetic_elements for item in els)

"""
#print (magnetic_elements)
t1 = ['Pb', 'Ti', 'O']
t2 = ['Pb', 'Ti', 'O', 'Fe']
print (magnetic_check(t1))
print (magnetic_check(t2))
"""

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

# the potential reset by Yi Wang, 08/24/20 
POTCAR_UPDATES = {
        'Mo': 'Mo_sv',  # 14 electrons, default was Mo_pv (12 electrons)
        'Nb': 'Nb_sv',  # 13 electrons, default was Nb_pv (11 electrons)
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
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['KPOINTS'].update({
        'grid_density': 1000,
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

    def __init__(self, structure, volume_relax=False, isif=None, **kwargs):
        """If volume relax is True, will do volume only, ISIF 7"""
        self.kwargs = kwargs
        self.volume_relax = volume_relax
        self.isif = isif
        uis = kwargs.get('user_incar_settings', {})
        if self.volume_relax and self.isif is not None:
            raise ValueError("isif cannot have a value while volume_relax is True.")
        if self.volume_relax:
            uis['ISIF'] = 7
        if self.isif is not None:
            uis['ISIF'] = self.isif
        if kwargs.get('user_incar_settings', None):
            kwargs.pop('user_incar_settings')
        super(RelaxSet, self).__init__(structure, RelaxSet.CONFIG, sort_structure=False, user_incar_settings=uis, **kwargs)


class PreStaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['KPOINTS'].update({
        'grid_density': 4000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-5,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': 1,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LAECHG": True,
        "LCHARG": False,
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
        super(PreStaticSet, self).__init__(structure, PreStaticSet.CONFIG, sort_structure=False, **kwargs)



class ForceConstantsSet(DictSet):
    """
    Set for calculating force constants calculations.

    Force constants are calculated by the finite difference method with symmetry considered.

    The smearing must be set to give the correct forces.
    The default is tuned for metals.

    Kpoints have a 8000 kpoints per reciprocal atom default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    # we never are comparing relaxations, only using them for optimizing structures.
    CONFIG['KPOINTS'].update({
        'grid_density': 4000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density') # to be explicit
    CONFIG['INCAR'].pop('ENCUT')  # use the ENCUT set by PREC
    CONFIG['INCAR'].update({
        'EDIFF_PER_ATOM': 1e-6,
        'ISMEAR': 1,
        'SIGMA': 0.2,
        'LREAL': 'Auto',
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
        uis = kwargs.get('user_incar_settings', {})
        super(ForceConstantsSet, self).__init__(
            structure, ForceConstantsSet.CONFIG, sort_structure=False, user_incar_settings=uis, **kwargs)


class StaticSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
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

    def __init__(self, structure, isif=2, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        self.isif = isif
        uis = kwargs.get('user_incar_settings', {})
        uis['ISIF'] = isif
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = kwargs
        super(StaticSet, self).__init__(structure, StaticSet.CONFIG, sort_structure=False, user_incar_settings=uis, **kwargs)


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
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
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
        super(ForcesSet, self).__init__(structure, ForcesSet.CONFIG, sort_structure=False, **kwargs)




class BornChargeSet(DictSet):
    """Set tuned for metal relaxations (correct smearing).
    Add `isif` parameter to the set to easily allow for overriding ISIF setting.
    Kpoints have a 6000 reciprocal density default.
    """
    CONFIG = _load_yaml_config("MPRelaxSet")
    CONFIG['POTCAR_FUNCTIONAL'] = 'PBE'
    CONFIG['KPOINTS'].update({
        'grid_density': 8000,
    })
    CONFIG['KPOINTS'].pop('reciprocal_density')  # to be explicit

    default = {
        'EDIFF_PER_ATOM': 1e-6,
        'ENCUT': 520,  # MP compatibility
        'ISMEAR': 0,
        "NSW": 0,
        "IBRION": -1,
        'LREAL': False,
        'ALGO': 'NORMAL',
        # other settings from MPStaticSet
        "LCHARG": False,
        "LWAVE": False,
        "ICHARG": 2,
        "LEPSILON": True,
    }


    CONFIG['POTCAR'].update(POTCAR_UPDATES)

    def __init__(self, structure, isif=2, **kwargs):
        # pop the old kwargs, backwards compatibility from the complex StaticSet
        self.isif = isif
        uis = kwargs.get('user_incar_settings', {})
        uis['ISIF'] = isif
        old_kwargs = ['prev_incar', 'prev_kpoints', 'grid_density', 'lepsilon', 'lcalcpol']
        for k in old_kwargs:
            try:
                kwargs.pop(k)
            except KeyError:
                pass
        self.kwargs = kwargs

        poscar = structure.to(fmt="poscar")
        unitcell_l = str(poscar).split('\n')
        els = [ e for e in unitcell_l[5].split(' ') if e!=""]
        if not magnetic_check(els): self.default['ISPIN':1]
        self.CONFIG['INCAR'].update(self.default)

        super(BornChargeSet, self).__init__(structure, BornChargeSet.CONFIG, sort_structure=False, user_incar_settings=uis, **kwargs)

