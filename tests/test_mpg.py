#!python
#
#Tests for some function of pymatgen

from pymatgen import Structure
from pymatgen.analysis.defects.generators import VacancyGenerator


POSCAR_STR = """FeNi3
1.0
3.5391385555 0.0000000000 0.0000000000
0.0000000000 3.5391385555 0.0000000000
0.0000000000 0.0000000000 3.5391385555
Fe Ni Cr Ti
1 1 1 1
Direct
0.000000000 0.000000000 0.000000000
0.000000000 0.500000000 0.500000000
0.500000000 0.500000000 0.000000000
0.500000000 0.000000000 0.500000000"""
struc = Structure.from_str(POSCAR_STR, fmt='POSCAR')

def test_vacancy_gen():
    #struc = Structure.from_file("VO2")
    vac_gen = VacancyGenerator(struc)

    vacs = list(vac_gen)
    for vac in vacs:
        print(vac.site)

    multiplicities = {str(v.site.specie): v.multiplicity for v in vacs}
    print(multiplicities)

test_vacancy_gen()