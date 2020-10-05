# Run dfttk for multi-structures

**This shows the `-mh` and `-r` parameters in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -mh MATCH_PATTERN, --match_pattern MATCH_PATTERN
                        The match pattern for structure file, and it should be
                        place in quotes. e.g. '*POSCAR*'. Default: * --
                        everything except SETTING files, ref. -s
  -r, --recursive       Recursive the path.

```

`-f` parameter, specify the structure or folder containing structure

- When the parameter value of `-f` is a file  it run this file; if it is a folder, it will search all supported files in the specified folder.



## Multi structure with `mh`

```bash
dfttk run -mh 'POSCAR*' -o
```

This command will write out the workflow of  `POSCAR.Si` and `POSCAR.mp-13`

## Run dfttk recursively

```bash
dfttk run -r -o
```

This command will write out the workflow of all of the supported structures (including it's sub- or sub-sub- ... folders)



## Supported file type

- It support **POSCAR**, **CONTCAR**, **CHGCAR**, **LOCPOT**, **cif**, **vasprun.xml**, **pymatgen’s structures**, **mcsqs's structure**, more details, ref [pymatgen.core.structure.IStructure.from_file](https://pymatgen.org/pymatgen.core.structure.html#pymatgen.core.structure.IStructure.from_file)
  - For **POSCAR** or **CONTCAR**, the name should be **\*POSCAR\*** or **\*CONTCAR\*** or **\*.vasp**
  - For **CHGCAR** or **LOCPOT**, the name should be **\*CHGCAR\*** or **\*LOCPOT\***
  - For **cif**, the name should be **\*.cif*** or **\*.mcif\***
  - For **vasprun.xml**, the name should be **vasprun\*.xml\***
  - For **pymatgen's structure**, the name shold be **\*.json** or **\*.yaml**
  - For **atat's structure**, the name should be **\*rndstr.in\*** or **\*lat.in\*** or **\*bestsqs\***

[Return](../)