#!python

import dfttk.utils as dfttkutils
from dfttk.input_sets import RelaxSet
from pymatgen import Structure

POSCAR_STR_check_symbol = """FCC_Fe_WithDiffMagMom
1.0
3.5391385555 0.0000000000 0.0000000000
0.0000000000 3.5391385555 0.0000000000
0.0000000000 0.0000000000 3.5391385555
Fe Fe
2 2
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.500000000 0.000000000
0.500000000 0.000000000 0.500000000"""
def test_check_symbol():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    symbols, natom = dfttkutils.check_symbol(InputSet)
    assert(symbols == ["Fe", "Fe"])
    assert(natom == ["2", "2"])

def test_update_pos_by_symbols():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    poscar_str = dfttkutils.update_pos_by_symbols(InputSet, write_file=False)
    syms = poscar_str.split("\n")[5]
    natom = poscar_str.split("\n")[6]
    assert(syms == "Fe Fe")
    assert(natom == "2 2")

def test_update_pot_by_symbols():
    struc = Structure.from_str(POSCAR_STR_check_symbol, fmt="POSCAR")
    magmoms = [4.0, 4.0, -4.0, -4.0]
    struc.add_site_property("magmom", magmoms)
    InputSet = RelaxSet(struc)
    potcar = dfttkutils.update_pot_by_symbols(InputSet, write_file=False)
    syms = potcar.symbols
    assert(syms == ["Fe", "Fe"])
