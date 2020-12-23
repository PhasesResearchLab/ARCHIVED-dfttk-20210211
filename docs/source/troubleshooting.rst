***************
Troubleshooting
***************


This document covers how to handle jobs that have fizzled. There are two main sections: general troubleshooting workflow to find the root cause of issue and specific issues and their fixes.


Troubleshooting Workflow
========================


**My job has fizzled!** The following steps can help you get information about how your job. You can imagine it as a decision tree. Check one thing before moving on to the next one.

1. Check that the job ran and has raised an exception with a traceback.

   Run ``lpad get_fws -i <ID> -d more``, replacing ``<ID>`` with the integer id of the Firework.
   Search the output for ``_exception``.
   What you see is the Python exception that was raised when running the Firework.

   *TIP:*: Searching works well when you pipe the output to ``less`` with ``lpad get_fws -i <ID> -d more | less`` and search using ``/``.

   .. todo:: If you don't see a traceback, that means... (this is the first step, but does this actually happen?)

   
2. Check the traceback is not a common error.

   See the `Common Errors section <CommonErrors>`_


.. _CommonErrors:

Common Errors
=============

Custodian VasprunXMLValidator failed
------------------------------------

In this error, you get a traceback that looks something like:

.. code-block:: python

   Traceback (most recent call last):
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 320, in run
       self._run_job(job_n, job)
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 428, in _run_job
       raise CustodianError(s, True, v)
   custodian.custodian.CustodianError: (CustodianError(...), 'Validation failed: <custodian.vasp.validators.VasprunXMLValidator object at 0x2af45b1d3908>')

   During handling of the above exception, another exception occurred:

   Traceback (most recent call last):
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/fireworks/core/rocket.py", line 262, in run
       m_action = t.run_task(my_spec)
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/atomate/vasp/firetasks/run_calc.py", line 204, in run_task
       c.run()
     File "/storage/home/bjb54/.conda/envs/wfs/lib/python3.7/site-packages/custodian/custodian.py", line 330, in run
       .format(self.total_errors, ex))
   RuntimeError: 0 errors reached: (CustodianError(...), 'Validation failed: <custodian.vasp.validators.VasprunXMLValidator object at 0x2af45b1d3908>'). Exited...


With the key being that Custodian fails to validate the ``vasprun.xml``. After running VASP, Custodian will try to parse the ``vasprun.xml`` file using pymatgen.

There are usually two possible triggers for this failure:

1. VASP failed to run at all (more common) or quit in a way that custodian did not detect (less common)
2. The ``vasprun.xml`` file could not be parsed by pymatgen.

To investigate this, first check that VASP ran (e.g. the ``OUTCAR`` shows that the run completed successfully).
If VASP did not run, find out why and fix that issue.
If VASP did run successfully, it was probably an issue parsing the ``vasprun.xml`` file.
Try parsing the ``vasprun.xml`` file using the ``pymatgen.io.vasp.outputs.Vasprun`` class.
If it throws an error when you try to parse, that's what made Custodian fail and you should fix that.
