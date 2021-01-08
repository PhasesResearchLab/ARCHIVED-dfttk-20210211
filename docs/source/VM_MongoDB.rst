Virtual Machines (VM) hosting
=============================

This section is for the dfttk users who want to host their MongoDB by themselves.

In computing, a virtual machine (VM) is an emulation of a computer system. Virtual machines are based on computer architectures and provide functionality of a physical computer. Their implementations may involve specialized hardware, software, or a combination. See https://en.wikipedia.org/wiki/Virtual_machine

In our case, we have changed the default tcp port from 27017 into 27018 due to historical reason

VM operation
------------

  Connect your VM by ssh (ssh youruserid@146.186.149.69 in our case). One should use VPN if you have firewall for your system

- To add user to your VM linux system

  Run the follwing command under Linux

.. code-block:: bash

    sudo adduser newuser
    usermod -aG sudo newuser #add admin user to your VM

Note on adding other users to access the VM. In order to add other users to access the VM from the Morpheus portal you would need to add them to the VM. For the VM itself you would need to add their user ID to the /etc/security/access.conf file and if they need sudo access you would need to add their ID to the /etc/sudoers.d/sudo-users file as well. 

MongoDB operation
-----------------

- MongoDB installation

  Run the follwing command

.. code-block:: bash

    apt update
    apt install mongodb

More information about installing and configuring MongoDB 3.6 is available here:
https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/

If you would like the most recent version of MongoDB (4.4), please reference this guide for installation and configuration instructions for Ubuntu 18.04:
https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/.
#Note: In step 2 of the link listed above, please reference the command for Ubuntu 18.04 (Bionic)

Additional configuration info for MongoDB 4.4 can be found here:
https://www.digitalocean.com/community/tutorials/how-to-install-mongodb-on-ubuntu-18-04-source

- Start mongdb service

.. code-block:: bash

    mongod --bind_ip_all -port 27018 &

- Shut down mongdb service

.. code-block:: bash

    db.adminCommand( { shutdown: 1 } )

For more details on MongodB user management, see https://docs.mongodb.com/manual/tutorial/enable-authentication/

- Connect to MongoDB by port 27018 for management

.. code-block:: bash

    mongo -port 27018

- Quit from MongoDB

  hit ``Ctrl+d``

- Create admin user for mongdb

  After connected to your mongoDB by ``mongo -port 27018``, input the following lines

.. code-block:: python

    use admin
    db.createUser(
      {
        user: "admin",
        pwd: "xxxxxxxxx", // xxxxxxxx is the admin password of your choice
        roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]
      }
    )

- Create general user

  Connect to your mongoDB as admin user by ``mongo --port 27018 --authenticationDatabase "admin" -u "admin" -p``, 
  followed by inputting the following lines

.. code-block:: python

    use userid-fws
    db.createUser({user: "userid", pwd: "B5nRcUvoCZ92", roles: [{role: "dbOwner", db: "userid-fws"}]})
    use userid-results
    db.createUser({user: "userid", pwd: "BeFihJ2mrKGm", roles: [{role: "dbOwner", db: "userid-results"}]})
    db.createUser({user: "userid-ro", pwd: "QIvaUT9ca6H8", roles: [{role: "read", db: "userid-results"}]})

  These lines can be produced by dfttk by run a python code named ``mongodb_user.py`` which 
  can be downlonded from
  https://github.com/PhasesResearchLab/dfttk/tree/master/dfttk/scripts
  After download the code, one can run it by 

.. code-block:: bash

    python mongodb_user.py

  The run will prompt the MongoDB system manager to input an userid for the user. After you input 
  userid and hit enter, one gets the above outputs in the screen. 

  Meanwhile, a file named ``db.json`` in the JSON format containing something similiar to 
  the following lines which should be sent to the MongoDB user.

.. _JSONLint: https://jsonlint.com

.. code-block:: JSON

    {
        "database": "userid-results",
        "collection": "tasks",
        "admin_user": "userid",
        "admin_password": "BeFihJ2mrKGm",
        "readonly_user": "userid-ro",
        "readonly_password": "QIvaUT9ca6H8",
        "host": "146.186.149.69",
        "port": 27018,
        "aliases": {}
    }

  The MongoDB user should save this data in a json file named "db.json" under the path 
  "dfttk/config" that created by "dfttk config -mp -aci" command.

- Remove user

.. code-block:: python

    db.removeUser(username)

- Check if mongodb is running, use

.. code-block:: python

    ps -ef | grep mongo

