==================
DFT job monitoring
==================

- Check run results

To get a brief summary for data in the database, use

.. code-block:: bash

  dfttk thfind -v 1

this will report the metadata, number of finished phonon calculations, number of finished static calculations, supercell size used for the phonon calculations, brief phase name; and the following command reports the location where the DFT calculations were submitted

.. code-block:: bash

  dfttk thfind -v 1 -jobpath parentfolder 

where parentfolder is the parent folder whose subfolders hold the individual DFT calculation job submission data.

- Batch postprocessing results
 
.. code-block:: bash

  dfttk thfind -get -py -td -50 -plot find_or_DFT -eq 4 -renew -ss 30 -w Pb-Ti-O 

where ``-get`` instruct thfind to call ``thelec`` module to postprocess the data, ``-py`` to use Yphon to recalculate the phonon density of states based on the force constants saved in the ``phonon`` collection, ``-td -50`` to use the self-adapted temperature mesh, ``-eq 4`` to use Birch-Murnaghan fitting to get equilibrium volume, ``-ss 30`` means only handle results with supercell size larger thatn 30, -w elist means only extract for compounds with all elements in the list of elist. Other enhanced screening conditions can be combinations of the follows

.. code-block:: bash

  -all elist # compounds must contain all elements in the list of elist;
  -any elist # compounds contain any elements in the list of elist;
  -xall elist # compounds must not contain all elements in the list of elist;
  -xany elist # compounds must not contain an y elements in the list of elist.

Batch postprocess may take longer time to finish. If it is the case, it is recommended to submit batch job. When submit batch, make sure compatibilities of non-ascii character by including the following in the job script:

.. code-block:: bash

  export LC_ALL='en_US.utf8' #for bsh;
  setenv LC_ALL en_US.utf8 #for csh

- Monitor workflows

To check running status of all submitted dfttk jobs, use

.. code-block:: bash

  lpad get_wflows

To find the running calculations

.. code-block:: bash

  lpad get_fws -s RUNNING 

To find the jobdir for a job with specific fw_id (The number portion of job id reported by ``lpad get_wflows``), use 
lpad get_launchdir  fw_id

- FIZZLED jobs

``FIZZLED`` means failed. The following command can summarize all FIZZLED jobs:

.. code-block:: bash

  lpad get_fws -s FIZZLED; or
  lpad get_wflows -s FIZZLED

One can rerun ``FIZZLED`` by 

.. code-block:: bash

  lpad rerun_fws -s FIZZLED

For more details, see ref. FIZZLED.

To get help for the subcommands, try ``lpad get_fws -h`` or ``lpad get_wflows -s FIZZLED -h``

