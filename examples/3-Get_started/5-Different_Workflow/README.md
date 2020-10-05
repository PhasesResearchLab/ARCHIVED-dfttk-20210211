# Write out workflow



**This shows the `-wf` parameter in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -wf WORKFLOW, --workflow WORKFLOW
                        Specify the workflow to run. Default: robust (run
                        get_wf_gibbs_robust workflow) (NOTE: currently, only
                        robust and born are supported.)

```



This command will allow you to run born charge workflow

```bash
dfttk run -wf born -o
```

[Return](../)