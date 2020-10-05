# Run dfttk with settings



**This shows the `-s` parameter and `SETTINGS.yaml` file in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -s SETTINGS, --setting SETTINGS
                        Specify the name of SETTINGS files (yaml or json file)
                        Default: SETTINGS (case insensitive and without ext)
                        The following filename will be treat as SETTINGS file
                        SETTINGS (global settings in the folder) Start with
                        SETTINGS- (individual settings for struct) End with
                        -SETTINGS (individual settings)

```

`-s` parameter, specify the name for settings, default: SETTINGS

- `SETTINGS.yaml` or `SETTINGS.json`: The global settings for current folder

- `SETTINGS-*.yaml(json)` or `*-SETTINGS.yaml(json)`: The individual settings, `*` is the name of the structure

- Case insensitive. (both SETTINGS and settings are OK)

- The value of `-s` parameter will replace SETTINGS.



## Change `INCAR` parameters and `POTCAR FUNCTIONAL`

Prepare the `SETTINGS.yaml` file, as follows:

```yaml
override_default_vasp_params : {'user_incar_settings': {'LREAL': 'Auto'}, 'user_potcar_functional': 'LDA', 'user_potcar_settings': {'Fe': 'Fe_sv'}}
```

Then run the following command:

```bash
dfttk run -o
```

This command will write out the workflow of  `POSCAR.mp-13_Fe` and changed the vasp settings



## Run phonon

Prepare the `SETTINGS.yaml` file as follows:

```yaml
phonon : True
phonon_supercell_matrix : [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
```

Then run the following command:

```bash
dfttk run -o
```

This workflow will run the `phonon` automatically



## Individual settings for multi structures

1. The `SETTINGS.yaml` is set for global

   ```yaml
   t_step : 10
   ```

2. For individual settings, use `SETTINGS-FILENAME_OF_STRUCTURE.yaml`, e.g.:

   `SETTINGS-POSCAR.mp-13_Fe.yaml`

   ```yaml
   override_default_vasp_params : {'user_incar_settings': {'LREAL': 'Auto'}, 'user_potcar_functional': 'LDA', 'user_potcar_settings': {'Fe': 'Fe_sv'}}
   ```

   `SETTINGS-POSCAR.Si.yaml`

   ```yaml
   phonon : True
   phonon_supercell_matrix : [[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]
   ```

3. Then run the following command

   ```bash
   dfttk run -o
   ```

   

## `-s` parameter

​	`-s` parameter allows the user change the filename of `SETTINGS.yaml` file, e.g. change into `SET.yaml`

​	Then run the following command:

```bash
dfttk run -s SET -o
```

[Return](../)

