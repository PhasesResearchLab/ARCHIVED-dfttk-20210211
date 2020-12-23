Virtual Machines (VM) hosting
=============================

This section is for the dfttk users who want to host their MongoDB by themselves.

In computing, a virtual machine (VM) is an emulation of a computer system. Virtual machines are based on computer architectures and provide functionality of a physical computer. Their implementations may involve specialized hardware, software, or a combination. See https://en.wikipedia.org/wiki/Virtual_machine

VM setup
--------

  Connect your VM by ssh (ssh yourpsuid@146.186.149.69 in our case). One should use VPN if you have firewall for your system

  - To add user to your VM linux system

.. code-block:: bash

    sudo adduser newuser
    usermod -aG sudo newuser #add admin user to your VM

Note on adding other users to access the VM. In order to add other users to access the VM from the Morpheus portal you would need to add them to the VM. For the VM itself you would need to add their user ID to the /etc/security/access.conf file and if they need sudo access you would need to add their ID to the /etc/sudoers.d/sudo-users file as well. 

MongoDB
-------

MongoDB installation

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

Start mongdb service

.. code-block:: bash

    mongod --bind_ip_all -port 27018 &

Connect to MongoDB by port 27018

.. code-block:: bash

    mongo -port 27018

Shut down

.. code-block:: bash

    db.adminCommand( { shutdown: 1 } )

MongodB user management, see https://docs.mongodb.com/manual/tutorial/enable-authentication/

Create admin user

.. code-block:: python

    use admin
    db.createUser(
      {
        user: "admin",
        pwd: "xxxxxxxxx", // or cleartext password
        roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]
      }
    )

Create general user

.. code-block:: python

    mongo --port 27018 --authenticationDatabase "admin" -u "admin" -p
    use psuid-fws
    db.createUser({user: "psuid", pwd: "B5nRcUvoCZ92", roles: [{role: "dbOwner", db: "psuid-fws"}]})
    use psuid-results
    db.createUser({user: "psuid", pwd: "BeFihJ2mrKGm", roles: [{role: "dbOwner", db: "psuid-results"}]})
    db.createUser({user: "psuid-ro", pwd: "QIvaUT9ca6H8", roles: [{role: "read", db: "psuid-results"}]})

Remove user

.. code-block:: python

    db.removeUser(username)

Check if mongodb is running

.. code-block:: python

    ps -ef | grep mongo

