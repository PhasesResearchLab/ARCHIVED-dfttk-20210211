MongoDB VM hosting
==================

This section is for the dfttk users who want to host their MongoDB by themselves.

MongoDB is one of the most popular document-oriented databases under the banner of NoSQL database. The schema-free implementation of MongoDB eliminates the prerequisites of defining a fixed structure required by the SQL database. 

In computing, a virtual machine (VM) is an emulation of a computer system. Virtual machines are based on computer architectures and provide functionality of a physical computer. Their implementations may involve specialized hardware, software, or a combination. See https://en.wikipedia.org/wiki/Virtual_machine

Our MongoDB databases are currently hosted by `Penn State's VM hosting service <https://cyberinfrastructure.psu.edu/?q=node/161>`_ that provides cost-effective, reliable VM for departments, colleges, and research units at Penn State University.

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

<<<<<<< HEAD
More broadly, the `MongoDB security checklist <https://docs.mongodb.com/manual/administration/security-checklist/>`_ should be followed. Critically, you should `enable access control <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_ and set up an authentication database and credentials for at least one administrator.

Some key things to configure in ``/etc/mongodb.conf``:

- Set the path of the database to the correct location. You may need to ``sudo chown mongodb /path/to/database``  for permissions to work
- Add the public IP of the server to ``bind_ip`` setting so it reads like ``127.0.0.1,XXX.XXX.XXX.XX`` for your server's IP ``XXX.XXX.XXX.XX``.
- Set a non-standard port if you want one
- Enable authentication (assuming you `enabled authentication <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_ and set up an administrator database and user)


MongoDB operation
-----------------

Managing the MongoDB service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**This section assumes that you already have your VM set up and your are managing it under linux environment.**

The MongoDB server is managed by a systemd service that is managed through ``systemctl``. Common commands are:

- Check status of mongodb service::

   sudo systemctl status mongodb

- Start mongodb service::

   sudo systemctl start mongodb
=======
.. code-block:: bash
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

    mongod --bind_ip_all -port 27018 &

- Shut down mongdb service

.. code-block:: bash

    db.adminCommand( { shutdown: 1 } )

For more details on MongodB user management, see https://docs.mongodb.com/manual/tutorial/enable-authentication/

- Connect to MongoDB by port 27018 for management

<<<<<<< HEAD
**This section deals with your MongoDB database management to locally or remotely operate on it.** This is to say you are going to manage your database from your local computer by ``mongo`` or ``mongosh``. 

#note mongosh by `MongoDB Shell <https://www.mongodb.com/try/download/shell?jmp=docs>`_ is the quickest way to connect, configure, query, and work with your MongoDB database 


- Create admin user for mongodb

These are the **one time setup** steps for MongoDB:

With access control enabled, ensure you have a user with userAdmin or userAdminAnyDatabase role in the admin database. This user can administrate user and roles such as: create users, grant or revoke roles from users, and create or modify customs roles.
=======
.. code-block:: bash

    mongod -port 27018

- Quit from MongoDB

  hit ``Ctrl+d``
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

- Create admin user for mongdb

<<<<<<< HEAD
=======
After connected to your mongoDB by ``mongo -port 27018``, input the following lines
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

1. Connect to the instance by open another terminal in your VM and connect a mongo shell to the instance::

    mongod --port 27018

after the prompt ">" input::

    use admin
    db.createUser(
      {
        user: "admin",
        pwd: "xxxxxxxxx", // xxxxxxxx is the admin password of your choice
        roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]
      }
    )

2. Re-start the MongoDB instance with access control

    a. Shut down the mongod instance
    b. Exit the mongo shell by run the command ``exit`` or give an EOF (``Ctrl+D``)

    c. Start the mongod with access control enabled by

      - adding the security.authorization configuration file setting

        .. code-block:: bash

          security:
              authorization: enabled

      - or If you start the mongod from the command line

        .. code-block:: bash

          mongod --auth --port 27018 


- Create general user

<<<<<<< HEAD
Assuming the service is running and configured with authentication (see above), Connect to your mongoDB as admin user locally by::
=======
Connect to your mongoDB as admin user by 
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

   mongo --port 27018 --authenticationDatabase "admin" -u "admin" -p

or remotelly by::

   mongo 146.186.149.69:27018 --authenticationDatabase admin -u <admin username> -p <admin password>
 
or remotelly use ``mongosh`` by::

   mongosh --username <admin username> --password --authenticationDatabase admin --host 146.186.149.69 --port 27018

followed by inputting the following lines after the prompt ">"::

    use userid-fws
    db.createUser({user: "userid", pwd: "B5nRcUvoCZ92", roles: [{role: "dbOwner", db: "userid-fws"}]})
    use userid-results
    db.createUser({user: "userid", pwd: "BeFihJ2mrKGm", roles: [{role: "dbOwner", db: "userid-results"}]})
    db.createUser({user: "userid-ro", pwd: "QIvaUT9ca6H8", roles: [{role: "read", db: "userid-results"}]})

These lines can be produced by dfttk by run a python code named ``mongodb_user.py`` which 
can be downlonded from
https://github.com/PhasesResearchLab/dfttk/tree/master/dfttk/scripts
<<<<<<< HEAD
After download the code, one can run it by::
=======
After download the code, one can run it by 

.. code-block:: bash
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

    python mongodb_user.py

The run will prompt the MongoDB system manager to input an userid for the user. After you input 
userid and hit enter, one gets the above outputs in the screen. 

<<<<<<< HEAD
Meanwhile, a file named ``db.json`` in the JSON format containing something similiar to
the following lines which should be sent to the MongoDB user::
=======
Meanwhile, a file named ``db.json`` in the JSON format containing something similiar to 
the following lines which should be sent to the MongoDB user.

.. _JSONLint: https://jsonlint.com

.. code-block:: bash
>>>>>>> parent of ceaf9b3f... Merge branch 'master' of github.com:PhasesResearchLab/dfttk into 20200807

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

The MongoDB user should save this data in a json file named ``db.json`` under the path 
``dfttk/config`` that created by ``dfttk config -mp -aci`` command.

- Remove user::

    db.removeUser(username)


