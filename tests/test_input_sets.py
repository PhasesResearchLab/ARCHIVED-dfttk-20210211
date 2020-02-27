#!python
#
import pytest
import dfttk.input_sets as input_sets
from pymatgen import Structure

struct = Structure.from_file("POSCAR")
vis = input_sets.RelaxSet(struct)
vis.write_input("test_write")