# Installation of DFTTK

## Content

- [1. Create a Virtual Environment (Optional)](#Create-a-Virtual-Environment-(Optional))
  - [a. Using Anaconda or Miniconda](#Using-Anaconda-or-Miniconda)
  - [b. Using virtualenv](#Using-virtualenv)
- [2. Install DFTTK](#Install-DFTTK)
  - [a. pip](#pip)
  - [b. develop version](#develop-version)
  - [c. No virtual environment](#No-virtual-environment)
- [3. Notice](#Notice)

## Create a Virtual Environment (Optional)

**Note: python version requirements: 3.6+ ([Atomate requirements](https://atomate.org/installation.html#create-a-python-3-virtual-environment))**

### Using [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

(1)  Using the anaconda on ACI 

```bash
#1.1 Load anconda environment (by default it will load python/3.6.3-anaconda5.0.1)
module load python
#1.2 Check the version of python
python --version
# It should shows: Python 3.6.3 :: Anaconda, Inc.

#2.1 Creat virtual environment
# conda create -n ENV_NAME python=VERSION 
conda create -n dfttk python=3.6 
#2.2 Activate 
source activate dfttk 
#2.3 Deactivate 
source deactivate 

#3.1 If you want to load the environment auto, please add the following sentence into ~/.basrc
module load python
source activate dfttk
#3.2 Update ~/.bashrc
source ~/.bashrc

#4. After installed dfttk, total space is 782Mb
```

Or

(2) Install your own conda

```bash
#Install of Miniconda
#1.1 Download Miniconda
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#1.2 Install Miniconda
bash Miniconda3-latest-Linux-x86_64.sh
#1.3 Check the version of python
python --version

#2.1 Creat virtual environment
# conda create -n ENV_NAME python=VERSION 
conda create -n dfttk python=3.8 
#2.2 Activate 
conda activate dfttk 
#2.3 Deactivate 
conda deactivate 

#3.1 Update ~/.bashrc
source ~/.bashrc
#3.2 It will add the environment variables in ~/.basrc, if you don't want the auto load, please delete/comment the following lines (There are some difference with yours, different usernames: mjl6505)
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/storage/home/mjl6505/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/storage/home/mjl6505/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/storage/home/mjl6505/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/storage/home/mjl6505/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

#4. After installed dfttk, total space is 994Mb
```

[To Top](#Content)

### Using [virtualenv](https://github.com/pypa/virtualenv)

```bash
#Here, we assumpt that all operations is in you home dir

#1.1 Load the python
module load python
#1.2 Install virtualenv
python -m pip install --user virtualenv

#2.1 Create virtual environment
#virtualenv --python=PYTHON_VERSION ENV_NAME 
virtualenv --python=python3.6 dfttk 
#2.2 Activate 
source dfttk/bin/activate 
#2.3 Deactivate 
deactivate 

#3.1 If you want to load the environment auto, please add the following sentence into ~/.basrc
source $HOME/dfttk/bin/activate
#3.2 Update ~/.bashrc
source ~/.bashrc

#4. After installed dfttk, total space is 458Mb
```

[To Top](#Content)

## Install DFTTK

Firstly, load your virtual environment (as described in above Section)

###  pip

```bash
pip install dfttk
```

**Note**: currently, only **version<=0.2.1** is supported by pip, the new version will be updated after more tests and release.

### develop version

```bash
#1. Download dfttk code
git clone https://github.com/phasesresearchlab/dfttk 
cd dfttk 
#2. Install dfttk
pip install -e . 
```

**Note**: Currently, there is some bugs in the development version in the group's branch, if you want to install it now, please download my branch of dfttk. (https://github.com/hitliaomq/dfttk)

### No virtual environment

```bash
#1. Load python
module load python
#2. Install dfttk (use pip or install develop version)
pip install --user dfttk
#3. After install, add following lines in your ~/.bashrc
module load python
export PATH=$HOME/.local/bin:$PATH
#Then update the ~/.bashrc
source ~/.bashrc

#3. Note. In ACI, if the following error occured, 
 RuntimeError: Running cythonize failed!
# please update cython by following command
 pip install --user --upgrade cython
```

[To Top](#Content)

## Notice

**Note: ACI will update linux version(from RedHat 6 to RedHat 7) recently (Sep. 14), and after the update, maybe there will exist some difference.**

```tex
Users who rely on the software already installed on the ICDS-ACI system likely won’t experience any problems and their transition should be seamless.
But if you:
	Compile your own software or build your own software from source
	Install your own Python or R modules
…then it’s likely that you’ll need to recompile or reinstall these elements for them to work properly in the new operating system.
```

[To TOP](#Content)

