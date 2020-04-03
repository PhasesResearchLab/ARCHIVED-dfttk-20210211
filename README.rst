==========================================
DFTTK: Density Functional Theory Tool Kits
==========================================

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

.. code-block:: bass

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
