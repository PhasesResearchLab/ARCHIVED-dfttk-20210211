# dfttk_tutorial
 Tutorials for dfttk (https://github.com/PhasesResearchLab/dfttk)

## Content

- [The simplest way to install dfttk](#The-simplest-way-to-install-dfttk)
  - [Installation](#Installation)
  - [Trouble Shooting](#Trouble-Shooting)
- [Installation of dfttk](./1-Installation)
- [Configuration for DFTTK](./2-Configuration/)
- [Get Started](./3-Get_started)



## The simplest way to install dfttk

### Installation

```bash
#Take ACI as an example
#1. Install dfttk
module load python
pip install --user dfttk

#2. Add the following lines to ~/.bashrc
export PATH=$HOME/.local/bin:$PATH

#3. Update ~/.bashrc
source ~/.bashrc

#4. Configuration of dfttk 
# (Prepare db.json and my_launchpad.yaml [do not put them in scratch folder]
#  Run the following command at the location)
dfttk config -all -aci
```

### Trouble shooting

1. If the code can't find vasp, please load vasp previously or provide `vaspjob.pbs`

   ```bash
   #For ACI
   module load intel
   module load impi
   module load vasp
   ```

2. Please ensure that the version of automate>0.9.4.

3. In ACI, if the following error occured, 

   ```bash
   RuntimeError: Running cythonize failed!
   ```

   Please update cython by following command

   ```bash
   pip install --user --upgrade cython
   ```

   

For details, please ref. the following sections

## [Installation of dfttk](./1-Installation)

## [Configuration for dfttk](./2-Configuration/)

## [Get Started](./3-Get_started/)