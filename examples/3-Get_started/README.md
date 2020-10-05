# Get Started

## Content

- [1. Run dfttk for a single structure or a single folder](./1-Single_Structure/)
- [2. Run dfttk for multi structures](./2-Multi_Structures/)
- [3. Run dfttk with different settings](./3-Run_With_Settings)
- [4. Run dfttk with phonon](./4-Run_With_Phonon)
- [5. Run different workflow](./5-Different_Workflow)
- [6. Write out workflow](./6-Wirte_Out_Workflow)
- [7. Launch to launchpad and submit to queue](./7-Launch_and_Submit)
- [Help for `dfttk run` command](#Help-for-dfttk-run-command)

## Help for `dfttk run` command

```shell
dfttk run -h
```

```shell
DFTTK version: 0.1+198.ga48e239
Copyright © Phases Research Lab (https://www.phaseslab.com/)

usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -h, --help            show this help message and exit
  -f STRUCTURE_FOLDER, --structure_folder STRUCTURE_FOLDER
                        The folder/file containing the structure, Default: '.'
  -mh MATCH_PATTERN, --match_pattern MATCH_PATTERN
                        The match pattern for structure file, and it should be
                        place in quotes. e.g. '*POSCAR*'. Default: * --
                        everything except SETTING files, ref. -s
  -s SETTINGS, --setting SETTINGS
                        Specify the name of SETTINGS files (yaml or json file)
                        Default: SETTINGS (case insensitive and without ext)
                        The following filename will be treat as SETTINGS file
                        SETTINGS (global settings in the folder) Start with
                        SETTINGS- (individual settings for struct) End with
                        -SETTINGS (individual settings)
  -r, --recursive       Recursive the path.
  -wf WORKFLOW, --workflow WORKFLOW
                        Specify the workflow to run. Default: robust (run
                        get_wf_gibbs_robust workflow) (NOTE: currently, only
                        robust and born are supported.)
  -ph, --phonon         Run phonon. This is equivalent with set phonon=True in
                        SETTINGS file
  -l, --launch          Launch the wf to launchpad
  -m [MAX_JOB], --max_job [MAX_JOB]
                        Run the job, only works when -l is specified. Default:
                        0 (Not submit to queue) 1: qlaunch singleshot (single
                        job) N(N>1): qlaunch rapidfire -m N
  -o, --write_out_wf    Write out the workflow
```

[TO TOP](#Content)