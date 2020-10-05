# Write out workflow



**This shows the `-o` parameter in `dfttk run`**

```bash
usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-l] [-m [MAX_JOB]] [-o]

optional arguments:
  -o, --write_out_wf    Write out the workflow

```



`-o` parameter, write out work flow

- When add `-o` parameter, the work flow for every structure will be write out

- The file name of the work flow is `dfttk_wf-` + filename of the structure + `.yaml` (e.g. `dfttk_wf-POSCAR.yaml`)

  ```shell
  dfttk run -f POSCAR -o
  ```

  It will write out the work flow in current folder as `dfttk_wf-POSCAR.yaml`

- The file will write out in the same folder with corresponding structure.

  ```shell
  dfttk run -f ./test_folder/POSCAR -o
  ```

  It will write out the work flow as `./test_folder/dfttk_wf-POSCAR.yaml`

[Return](../)