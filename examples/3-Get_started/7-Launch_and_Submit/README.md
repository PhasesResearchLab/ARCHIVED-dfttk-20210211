# Write out workflow



**This shows the `-l` and `-m` parameters in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -l, --launch          Launch the wf to launchpad
  -m [MAX_JOB], --max_job [MAX_JOB]
                        Run the job, only works when -l is specified. Default:
                        0 (Not submit to queue) 1: qlaunch singleshot (single
                        job) N(N>1): qlaunch rapidfire -m N

```



`-l` parameter, launch to launchpad

`-m` parameter, launch to queue system and determine the number of jobs

- Add `-l` parameter, it will submit the job to launchpad

- Add `-m N` parameter, it will submit the job to queue system. `N` specifies the number of jobs run at the same time. (Note: This only work when `-l` parameter is added.) 

  ```shell
  dfttk run -l
  dfttk run -l -m 1
  dfttk run -l -m 2
  ```

  - When `-m 1`, it will run `qlaunch singleshot` (fireworks command)
  - When `-m N` (N>1), it will run `qlaunch rapidfire -m N` (fireworks command)
  - For more details, please ref. [Fireworks: Launch Rockets through a queue](https://materialsproject.github.io/fireworks/queue_tutorial.html)



e.g.

```bash
dfttk run -r -l -m 2
```

This will submit all structures to launchpad and queue

[Return](..)