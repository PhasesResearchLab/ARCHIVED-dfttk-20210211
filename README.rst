=====
DFTTK
=====

This repo contains the custom workflows developed by the Phases Research Lab that do not fit into the scope of the public atomate repository.

The following workflows are currently implemented:
- Gibbs energy workflow for stable structures

or under development:
- Minimum volume finding workflow

Installation
============

DFTTK requires Python 3. Python 2 support for NumPy ends 2019-01-01.

pip
---

From the command line, run ``pip install dfttk``

conda
-----

Anaconda packages of DFTTK are currently not supported. If you are using Anaconda, you should be able to install with pip.

development versions
--------------------

1. ``git clone https://github.com/phasesresearchlab/dfttk``
2. ``cd dfttk``
3. ``pip install -e .``

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
