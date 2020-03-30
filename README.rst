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

- virtualenv

 - virtualenv_ is required

    .. code-block:: bash
    
        virtualenv --python=python3.6 dfttk
        #Activate
        source dfttk/bin/activate
        #Deactivate
        deactivate

- conda

 - Anaconda_ or Miniconda_ is required

    .. code-block:: bash

        conda create -n dfttk python=3.6
        #Activate
        conda active dfttk
        #Deactivate
        conda deactivate

- FAQ

 Q: The path pythonX (from --python=pythonX) does not exist (for virtualenv)

 A: Please check the python version of the system ``python --version``. If its python 2.x, please update it to python3.6+. (**PRL Group Notes:** For **ACI** account, please you just need to load the correct python by ``module load python``)


.. _virtualenv: https://github.com/pypa/virtualenv
.. _Anaconda: https://www.anaconda.com/
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html


Install dfttk
-------------

- pip

From the command line, run ``pip install dfttk``

- conda

Anaconda packages of DFTTK are currently not supported. If you are using Anaconda, you should be able to install with pip.

- development versions

.. code-block:: bash

    git clone https://github.com/phasesresearchlab/dfttk
    cd dfttk
    pip install -e .

Use
===

``from dfttk import get_wf_gibbs``. Examples forthcoming.

Contributing
============

See CONTRIBUTING.rst_

.. _CONTRIBUTING.rst: CONTRIBUTING.rst

License
-------

DFTTK is MIT licensed. See LICENSE_

.. _LICENSE: LICENSE
