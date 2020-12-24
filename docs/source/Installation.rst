Installation
============

- Requirements

  DFTTK requires YPHON, MongodB, Python 3. Python 2 support for NumPy ends 2019-01-01. 

- Release version

.. code-block:: bash

    pip install dfttk
    mkdir dfttk #if dfttk folder not exist. 
    cd dfttk
    dfttk config -mp -aci #a folder named "config" will be created where running environmental info saved

- Development version

.. code-block:: bash

    git clone https://github.com/PhasesResearchLab/dfttk.git
    cd dfttk
    pip install -e .
    cd dfttk
    dfttk config -mp -aci #a folder named "config" will be created where running environmental info saved

- YPHON

To postprocess the finite properties, the Yphon package is required. Yphon can be installed by run

.. code-block:: bash

    cd ~
    git clone https://github.com/yiwang62/YphonPackage
    #Note: Usually the precompiled binaries should be executable in the common Linux/Unix environment. If not, do the following:

.. code-block:: bash

    cd YphonPackage/YPHON/YPHON 
    make
    #Note: If errors reported in the compiling stage, insert one line #define R_OK 1 after #include

For csh user: the command search path should be changed by inserting line below into the .cshrc  (.tcshrc) file

.. code-block:: bash

    set path = (. ~/YphonPackage/YPHON/YPHON $BIN_PATH $path)

For bsh user: the command search path should be changed by inserting the lines below into the .bash_profile (.bashrc) file

.. code-block:: bash

    PATH=.:~/YphonPackage/YPHON/YPHON:$BIN_PATH:$PATH
    export PATH


- MongoDB 

  Ask the MongoDB system manager to set up the access credential information by downlond a python code named
  ``mongodb_user.py`` from https://github.com/PhasesResearchLab/dfttk/tree/master/dfttk/scripts followed run it by

.. code-block:: bash

    python mongodb_user.py

  The run will prompt the MongoDB manager to input an userid for your user. After you input userid 
  and hit enter, the following lines 

.. code-block:: python

    use psuid-fws
    db.createUser({user: "psuid", pwd: "B5nRcUvoCZ92", roles: [{role: "dbOwner", db: "psuid-fws"}]})
    use psuid-results
    db.createUser({user: "psuid", pwd: "BeFihJ2mrKGm", roles: [{role: "dbOwner", db: "psuid-results"}]})
    db.createUser({user: "psuid-ro", pwd: "QIvaUT9ca6H8", roles: [{role: "read", db: "psuid-results"}]})

  These lines will be used by the MongoDB manager the MongoDB user setup, see the ``VM setup`` section in this document.

  Meanwhile, a file named ``db.json`` containing in the JSON format similiar to the following lines which should be sent to your user.

.. _JSONLint: https://jsonlint.com

.. code-block:: JSON

    {
        "database": "psuid-results",
        "collection": "tasks",
        "admin_user": "psuid",
        "admin_password": "BeFihJ2mrKGm",
        "readonly_user": "psuid-ro",
        "readonly_password": "QIvaUT9ca6H8",
        "host": "146.186.149.69",
        "port": 27018,
        "aliases": {}
    }
    #note: Sometimes, when using "copy/paste" with Windows, some invisible characters may be hidden by linux "vi". Make sure show/delete the invisible characters by vi command ":set list". 

  Save this as a json file named "db.json" under the "dfttk/config" folder that created by "dfttk config -mp -aci" command mentioned above. 
