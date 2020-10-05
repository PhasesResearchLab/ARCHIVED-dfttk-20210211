# Run dfttk with phonon

**This shows the `-ph`  parameters in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -ph, --phonon         Run phonon. This is equivalent with set phonon=True in
                        SETTINGS file

```



## Run `phonon` with `SETTINGS.yaml` file

```yaml
phonon : True
phonon_supercell_matrix : [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
```

For different structures, we need to set different `SETTINGS-STRUCUTRS.yaml` files to set different `phonon_supercell_matrix`



## Run `phonon` with `-ph` parameter

```bash
dfttk run -r -ph -o
```

This command will find the proper supercell matrix for different structures

1. For BCC-Nb

```yaml
#For bcc-Nb
 lattice_matrix:
 - [3.32052, 0.0, 0.0]
 - [0.0, 3.32052, 0.0]
 - [0.0, 0.0, 3.32052]
phonon_supercell_matrix:
- [3, 0, 0]
- [0, 3, 0]
- [0, 0, 3]
```

2. For HCP-Zr

```yaml
#For HCP-Zr
lattice_matrix:
- [3.239232, 0.0, 0.0]
- [-1.619616, 2.805257, 0.0]
- [0.0, 0.0, 5.17222]
phonon_supercell_matrix:
- [3, 0, 0]
- [2, 4, 0]
- [0, 0, 2]
```

3. For Orth-CaTiO3

```yaml
 lattice_matrix:
 - [5.40776, 0.0, 0.0]
 - [0.0, 5.507429, 0.0]
 - [0.0, 0.0, 7.694012]
phonon_supercell_matrix:
- [2, 0, 0]
- [0, 2, 0]
- [0, 0, 1]
```

4. For P1-MnP4

```yaml
 lattice_matrix:
 - [5.090496, 0.0, -0.460003]
 - [-2.532782, 5.2766, 0.075695]
 - [0.0, 0.0, 16.351185]
phonon_supercell_matrix:
- [2, 0, 0]
- [1, 2, 0]
- [0, 0, 1]
```

[Return](../)