Installation
============

- Requirements

  DFTTK requires MongodB, Python 3. Python 2 support for NumPy ends 2019-01-01. 

- Release version

  pip install dfttk

- Development version

  git clone https://github.com/PhasesResearchLab/dfttk.git

  cd dfttk

  pip install -e .

- Run environmental setup

  dfttk config -mp -aci 

- MongoDB 

  Ask the MongoDB system manager to set up the access credential file named "db.json" containing: 

{

    "database": "psuid-results",

    "collection": "tasks",

    "admin_user": "psuid",

    "admin_password": "mz3vWHwKGvEU",

    "readonly_user": "psuid-ro",

    "readonly_password": "tv4BZQmHBinD",

    "host": "146.186.149.69",

    "port": 27018,

    "aliases": {}
}

