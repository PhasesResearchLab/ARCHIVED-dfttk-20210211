# Configuration of DFTTK

## Content

- [Preparation](#Preparation)
- [Configuration for atomate](#Configuration-for-atomate)
- [Configuration for pymatgen](#Configuration-for-pymatgen)
- [Configure all with one command](#Configure-all-with-one-command)
- [Example](#Example)
- [Validation for configuration](#Validation-for-configuration)
- [Help for `dfttk config` command](#Help-for-dfttk-config-command)

---

## Preparation

For configure dfttk, the following file you need to prepare.

```
current_folder
├── psp                          [specified by -psp]
│   ├── pseudopotential_content  [required if you didnot configurate pymatgen]
│   └── ...
├── config                       [specified by -c parameter]
│   ├── db.json                  [required]
│   ├── my_launchpad.yaml        [required]
│   ├── FW_config.yaml           [optional]
│   ├── my_fworker.yaml          [optional]
│   └── my_qadapter.yaml         [optional]
└── vaspjob.pbs                  [optional, specified by -q parameter]
```

- **pseudopotential_content:** The vasp pseudopotential.

  - It can be compressed file (`*.tar.gz`) or uncompressed folder. The name should be in the following list. e.g `potpaw_PBE.tar.gz` (compressed file)

    ```python
    ["potpaw_PBE", "POT_GGA_PAW_PBE", "potpaw_PBE_52", "POT_GGA_PAW_PBE_52", "potpaw_PBE_54", "POT_GGA_PAW_PBE_54", "potpaw_PBE.52", "POT_GGA_PAW_PBE_52", "potpaw_PBE.54", "POT_GGA_PAW_PBE_54", "potpaw_LDA", "POT_LDA_PAW", "potpaw_LDA.52", "POT_LDA_PAW_52", "potpaw_LDA.54", "POT_LDA_PAW_54", "potpaw_LDA_52", "POT_LDA_PAW_52", "potpaw_LDA_54", "POT_LDA_PAW_54", "potUSPP_LDA", "POT_LDA_US", "potpaw_GGA", "POT_GGA_PAW_PW91", "potUSPP_GGA": "POT_GGA_US_PW91"]
    ```

    The file structure should be as follows:

    ```shell
    e.g. psp
         ├── potpaw_LDA.54.tar.gz
         ├── potpaw_PBE.54.tar.gz
         └── ...
     or: psp
         ├── potpaw_LDA_54
         │   ├── Ac
         │   └── ...
         ├── potpaw_PBE_54
         │   ├── Ac
         │   └── ...
         └── ...
    ```

    The original pseudopotential file, please ask those who are in charge of vasp in your group.

- **MAPI_KEY:** The API key of [Materials Project](https://materialsproject.org), ref. [API KEY](https://materialsproject.org/open)

- **vaspjob.pbs:** The submitting script of your queue system. **Currently, only pbs system is supported**

- **other config files:** At least, you need prepare two files, `db.json` and `my_launchpad.yaml`. The template is shown in the `config` folder.

  - For more details, please ref. [Configure database connections and computing center parameters](https://atomate.org/installation.html#configure-fireworks)

- **PRL GROUP NOTES:**

  - For **db.json** and **my_launchpad.yaml** file, please ask **Brandon** for help.
  - **vaspjob.pbs** for ACI can be download from github
  - If you have the pseudopotential, please download it from group's box(Group Documents/VASP_potential).

[TO TOP](#Content)

---

## Configuration for atomate

- Config manual, ref [Configure database connections and computing center parameters](https://atomate.org/installation.html#configure-fireworks)

- `dfttk config -a`

  ```shell
  usage: dfttk config -a [-p PATH_TO_STORE_CONFIG] [-c CONFIG_FOLDER] [-q QUEUE_SCRIPT] 
                         [-qt QUEUE_TYPE] [-v VASP_CMD_FLAG] 
  ```

  The meaning of the parameters, please ref [Help for `dfttk config` command](#Help-for-dfttk-config-command) or just run `dfttk config -h`

[TO TOP](#Content)

---

## Configuration for pymatgen

- Config manual, ref [POTCAR setup in pymatgen](https://pymatgen.org/installation.html#potcar-setup)

- `dfttk config -mp`

  ```shell
  usage: dfttk config -mp [-p PATH_TO_STORE_CONFIG] [-psp VASP_PSP_DIR] [-mapi MAPI_KEY] 
                      [-df DEFAULT_FUNCTIONAL]
  ```

  The meaning of the parameters, please ref [Help for `dfttk config` command](#Help-for-dfttk-config-command) or just run `dfttk config -h`

[TO TOP](#Content)

---

## Configure all with one command

- `dfttk config -all`

  ```shell
  usage: dfttk config [-h] [-all] [-p PATH_TO_STORE_CONFIG] [-a]
                      [-c CONFIG_FOLDER] [-q QUEUE_SCRIPT] [-qt QUEUE_TYPE]
                      [-v VASP_CMD_FLAG] [-mp] [-psp VASP_PSP_DIR]
                      [-mapi MAPI_KEY] [-df DEFAULT_FUNCTIONAL]
  ```

  The meaning of the parameters, please ref [Help for `dfttk config` command](#Help-for-dfttk-config-command) or just run `dfttk config -h`

[TO TOP](#Content)

---

## Example

```shel
dfttk config -all -p test_config -mapi test_api
```

In `test_config` folder

```shell
test_config
├── config
│   ├── db.json
│   ├── FW_config.yaml
│   ├── my_fworker.yaml
│   ├── my_launchpad.yaml
│   └── my_qadapter.yaml
├── logs
└── vasp_psp
    ├── POT_GGA_PAW_PBE_54
    │   ├──Ac
    │   └── ...
    └── POT_LDA_PAW_54
        ├──Ac
        └── ...
```

In `~/.bashrc`

```shell
export FW_CONFIG_FILE=/gpfs/scratch/mjl6505/test/dfttk/tutorial/config/test_config/config/FW_config.yaml
```

in `~/.pmgrc.yaml`

```yaml
PMG_DEFAULT_FUNCTIONAL: PBE
PMG_MAPI_KEY: test_api
PMG_VASP_PSP_DIR: /gpfs/scratch/mjl6505/test/dfttk/tutorial/config/test_config/vasp_psp
```

[TO TOP](#Content)

---

## Validation for configuration

- `dfttk config -t`

  This command will validate the configuration, and give some error or tips for the incorrect configuration

  ```shell
  dfttk config -t [{all,pymatgen,atomate,none}]
  ```

  - Currently only support validation for the configuration of pymatgen

[TO TOP](#Content)

---

## Help for `dfttk config` command

```shell
dfttk config -h
```

```shell
DFTTK version: 0.1+99.ga2da70f.dirty
Copyright © Phases Research Lab (https://www.phaseslab.com/)

usage: dfttk config [-h] [-all] [-p PATH_TO_STORE_CONFIG] [-a]
                    [-c CONFIG_FOLDER] [-q QUEUE_SCRIPT] [-qt QUEUE_TYPE]
                    [-v VASP_CMD_FLAG] [-mp] [-psp VASP_PSP_DIR]
                    [-mapi MAPI_KEY]
                    [-df {LDA,LDA_52,LDA_54,LDA_US,PBE,PBE_52,PBE_54,PW91,PW91_US,Perdew-Zunger81}]
                    [-t [{all,pymatgen,atomate,none}]]

optional arguments:
  -h, --help            show this help message and exit
  -all, --all           Configure atomate and pymatgen.
  -p PATH_TO_STORE_CONFIG, --prefix PATH_TO_STORE_CONFIG
                        The folder to store the config files. Default: .
                        (current folder)
  -a, --atomate         Configure atomate.
  -c CONFIG_FOLDER, --config_folder CONFIG_FOLDER
                        The folder containing config files, at least contain
                        db.json and my_launchpad.yaml. Default: '.'
  -q QUEUE_SCRIPT, --queue_script QUEUE_SCRIPT
                        The filename of the script for sumitting vasp job. It
                        will search in current folder and sub-folders.
                        Default: vaspjob.pbs
  -qt QUEUE_TYPE, --queue_type QUEUE_TYPE
                        The type of queue system. Note, only pbs is supported
                        now. Default: pbs
  -v VASP_CMD_FLAG, --vasp_cmd_flg VASP_CMD_FLAG
                        The flag to distinguish vasp_cmd to othe commands in
                        queue_script. Default: vasp_std
  -mp, --pymatgen       Configure pymatgen.
  -psp VASP_PSP_DIR, --vasp_psp_dir VASP_PSP_DIR
                        The path of pseudopotentials. Default: psp
  -mapi MAPI_KEY, --mapi_key MAPI_KEY
                        The API key of Materials Projects
  -df {LDA,LDA_52,LDA_54,LDA_US,PBE,PBE_52,PBE_54,PW91,PW91_US,Perdew-Zunger81}, --default_functional {LDA,LDA_52,LDA_54,LDA_US,PBE,PBE_52,PBE_54,PW91,PW91_US,Perdew-Zunger81}
                        The default functional. Default: PBE
  -t [{all,pymatgen,atomate,none}], --test_config [{all,pymatgen,atomate,none}]
                        Test for configurations. Note: currently only support
                        for pymatgen.
```

