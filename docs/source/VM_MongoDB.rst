MongoDB VM hosting
==================

This section is for the dfttk users who want to host their MongoDB by themselves.

MongoDB is one of the most popular document-oriented databases under the banner of NoSQL database. The schema-free implementation of MongoDB eliminates the prerequisites of defining a fixed structure required by the SQL database.

In computing, a virtual machine (VM) is an emulation of a computer system. Virtual machines are based on computer architectures and provide functionality of a physical computer. Their implementations may involve specialized hardware, software, or a combination. See https://en.wikipedia.org/wiki/Virtual_machine

Our MongoDB databases are currently by Penn State's VM hosting service that provides cost-effective, reliable VM for departments, colleges, and research units at Penn State University.

VM operation
------------

  Connect your VM by ssh (ssh youruserid@146.186.149.69 in our case). One should use VPN if you have firewall for your system

- To add user to your VM linux system

  Run the follwing command under Linux

.. code-block:: bash

    sudo adduser newuser
    usermod -aG sudo newuser #add admin user to your VM

Note on adding other users to access the VM. In order to add other users to access the VM from the Morpheus portal you would need to add them to the VM. For the VM itself you would need to add their user ID to the /etc/security/access.conf file and if they need sudo access you would need to add their ID to the /etc/sudoers.d/sudo-users file as well.

MongoDB installation/configuration
----------------------------------

These are the **one time setup** steps for MongoDB:

Assuming you are on Ubuntu and using the ``apt`` package manager::

   apt update
   apt install mongodb

This may be an outdated MongoDB. To follow the best security practices and get the latest features, you should install something more recent. More information about installing and configuring MongoDB 3.6 is available here:
https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/


Additional configuration info for MongoDB 4.4 can be found here:
https://www.digitalocean.com/community/tutorials/how-to-install-mongodb-on-ubuntu-18-04-source


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

The MongoDB server is managed by a systemd service that is managed through ``systemctl``. Common commands are:

- Check status of mongodb service::

   sudo systemctl status mongodb

- Start mongodb service::

   sudo systemctl start mongodb

- Shut down mongodb service::

   sudo systemctl stop mongodb

- Restart mongodb service (if it's already running)::

   sudo systemctl restart mongodb


Connecting to the MongoDB server console
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Assuming the service is running and configured with authentication (see above), connect to the MongoDB console by::

   mongo 146.186.149.69:27018 --authenticationDatabase admin -u <admin username> -p <admin password>

- To exit, run the command ``exit`` or give an EOF (``Ctrl+D``)

User management
~~~~~~~~~~~~~~~

For more details on MongoDB user management, see https://docs.mongodb.com/manual/tutorial/enable-authentication/

- Create admin user for mongodb

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

Connect to your mongoDB as admin user by

.. code-block:: bash

    mongo --port 27018 --authenticationDatabase "admin" -u "admin" -p

followed by inputting the following lines

.. code-block:: bash

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

.. code-block:: bash

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

- Remove user

.. code-block:: python

    db.removeUser(username)

- Check if mongodb is running, use

.. code-block:: python

    ps -ef | grep mongo

