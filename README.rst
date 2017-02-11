=============
PRL Workflows
=============

This repo contains the custom workflows developed by the Phases Research Lab that cannot be public or have not been submitted to the public atomate repository.

The following workflows are currently implemented:

or under development:
- A basic (half) NEB calculation

You can also find some scripts for basic atomate calculations in the `scripts` folder. Simply configure the settings as you want them run the script using Python on a system that is configured to connect to your LaunchPad.

Installation
------------

1. You should ideally download this repo to your codes directory in your atomate installation 
2. Install it as editable using pip (``pip install -e .``)
3. Add the folder containing the *Firetasks* you are interested in to the ``ADD_USER_PACKAGES`` entry in ``FW_config.yaml`` formatted like a Python import. For example:

::

    ADD_USER_PACKAGES:
      - PRLWorkflows.neb

Contributing
------------

See CONTRIBUTING.rst_

.. _CONTRIBUTING.rst: CONTRIBUTING.rst
