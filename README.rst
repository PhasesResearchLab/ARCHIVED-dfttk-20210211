Usage of thfind module
usage: dfttk thfind [-h] [-q [QHAMODE]] [-w [WITHIN]] [-all [CONTAINALL]]
                    [-any [CONTAINANY]] [-v [NV]] [-ss [SUPERCELLN]] [-get]
                    [-py] [-T0 [T0]] [-T1 [T1]] [-dT [TD]] [-xdn [XDN]]
                    [-xup [XUP]] [-dope [DOPE]] [-ne [NDOSMX]]
                    [-natom [NATOM]] [-e [EVERYT]] [-gauss [GAUSSIAN]]
                    [-i [DOSCAR]] [-o [OUTF]] [-noel] [-metatag [METATAG]]
                    [-qhamode [QHAMODE]] [-eq [EQMODE]] [-el [ELMODE]] [-s]
                    [-plot] [-g] [-expt [EXPT]] [-xlim [XLIM]]

optional arguments:
  -h, --help            show this help message and exit
  -q [QHAMODE], --qhamode [QHAMODE]
                        Collection. 'phonon', 'qha'. Default: 'phonon'
  -w [WITHIN], --within [WITHIN]
                        find calculations within element list Default: None
  -all [CONTAINALL], --containall [CONTAINALL]
                        find calculations must contain all elements in the
                        list Default: None
  -any [CONTAINANY], --containany [CONTAINANY]
                        find calculations contain any elements in the list
                        Default: None
  -v [NV], --nV [NV]    Return phonon calculations finished for number of
                        volumes larger or equals to. Default: 6
  -ss [SUPERCELLN], --supercellsize [SUPERCELLN]
                        only return phonon calculation with supercell size
                        larger than. Default: 0
  -get, --get           get the thermodyamic data for all found entries.
                        Default: False
  -py, --pyphon         use Yphon to recalculate vibrational properties.
                        Default: False
  -T0 [T0], -t0 [T0]    Low temperature limit. Default: 0
  -T1 [T1], -t1 [T1]    High temperature limit. Default: 1300
  -dT [TD], -td [TD]    Temperature increment. Default: 10
  -xdn [XDN], --xdn [XDN]
                        Low band energy limit. Default: -100 (eV)
  -xup [XUP], --xup [XUP]
                        High band energy limit. Default: 100
  -dope [DOPE], --dope [DOPE]
                        dope level (electrons). Default: -1.e-8 for numerical
                        stability
  -ne [NDOSMX], --ndosmx [NDOSMX]
                        new DOS mesh. Default: 10001
  -natom [NATOM], --natom [NATOM]
                        number of atoms in the DOSCAR. Default: 1
  -e [EVERYT], --everyT [EVERYT]
                        number of temperature points skipped from QHA
                        analysis. Default: 1
  -gauss [GAUSSIAN], --gauss [GAUSSIAN]
                        densing number near the Fermi energy. Default: 1000
  -i [DOSCAR], --doscar [DOSCAR]
                        DOSCAR filename. Default: DOSCAR
  -o [OUTF], -outf [OUTF]
                        output filename for calculated thermoelectric
                        properties. Default: fvib_ele
  -noel, -noel          do not consider the thermal electron contribution.
                        Default: False
  -metatag [METATAG], -metatag [METATAG]
                        metatag: MongoDB metadata tag field. Default: None
  -qhamode [QHAMODE], -qhamode [QHAMODE]
                        quasiharmonic mode: debye, phonon, or yphon. Default:
                        debye
  -eq [EQMODE], --eqmode [EQMODE]
                        Mode to calculate LTC. 0: Symmetrical Central
                        differential; 4: 4-parameter BM fitting. 5:
                        5-parameter BM fitting. Default: 0
  -el [ELMODE], --elmode [ELMODE]
                        Mode to interpolate thermal electronic contribution:
                        0: interp1d; 1: UnivariateSpline. Default: 0
  -s, -smooth           smooth the LTC. Default: False
  -plot, -plot          plot the figure. Default: False
  -g, --debug           turn on debug mode by reducing the mesh. Default:
                        False
  -expt [EXPT], -expt [EXPT]
                        json file path for experimental thermodynamic
                        properties for plot. Default: None
  -xlim [XLIM], -xlim [XLIM]


=========================================
DFTTK: Density Functional Theory Tool Kit
=========================================

**Ultimate goals:** For a given structure and elements, calculate the free energy with respect to possible internal degree of freedoms.

- Features

 - High-throughput. It can run plenty of structures with one simple command.
 - Simple. Only the structure file is required.

- The following workflows are currently implemented:

 - Gibbs energy workflow for stable structures
 - Minimum volume finding workflow

**Note:** This repo contains the custom workflows developed by the Phases Research Lab that do not fit into the scope of the public atomate repository.


Installation
============

DFTTK requires Python 3. Python 2 support for NumPy ends 2019-01-01.

Create virtual environment (optional)
-------------------------------------

Anaconda_ or Miniconda_ is required. (Another option is using virtualenv_)

.. code-block:: bash

    #conda create -n ENV_NAME python=VERSION
    conda create -n dfttk python=3.6
    #Activate
    conda activate dfttk
    #Deactivate
    conda deactivate

.. _virtualenv: https://github.com/pypa/virtualenv
.. _Anaconda: https://www.anaconda.com/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html


Install dfttk
-------------

To install dfttk, there are several ways.

- pip

From the command line, run ``pip install dfttk``

- conda

Anaconda packages of DFTTK are currently not supported. If you are using Anaconda, you should be able to install with pip.

- development versions

.. code-block:: bash

    git clone https://github.com/phasesresearchlab/dfttk
    cd dfttk
    pip install -e .

Configuration
=============

Preparation
-----------

Prepare following files.

.. code-block:: bash

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


Config
------

After prepared above files, simply run

.. code-block:: bash

    dfttk config -all

**PRL GROUP NOTE:** If you use ACI cluster at PRL group, you can get the pseudopotentials from ACI

.. code-block:: bash

    dfttk config -all -aci

For more details, ref. Configuration_

.. _Configuration: docs/Configuration.md

Use
===

1. Prepare structure file(s), e.g. POSCAR
2. Simply run ``dfttk run -f POSCAR``

For more details, ref. Get_Started_

.. _Get_Started: docs/Get_started.md

Contributing
============

See CONTRIBUTING.rst_

.. _CONTRIBUTING.rst: CONTRIBUTING.rst

License
-------

DFTTK is MIT licensed. See LICENSE_

.. _LICENSE: LICENSE
