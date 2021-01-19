=========
run DFTTK
=========

To run DFTTK, the user should follow pattern:

.. code-block:: bash

    dfttk 'module' [option]

where the current available modules are listed below 

+---------------+---------------------------------------------------------------+
|   module      |       function
+---------------+---------------------------------------------------------------+
|   run         | call VASP to perform 0 K E-V, phonon, elastic moduli, Born    |
|               | effective charge (dielectric tensor), etc                     |
+---------------+---------------------------------------------------------------+
|   config      | set basic dfttk run envirnment, such as path to               |
|               | pseudopotential, batch job prototype                          |
+---------------+---------------------------------------------------------------+
|   db_remove   | Remove certain data in MongoDB to save storage cost           |
+---------------+---------------------------------------------------------------+
|   thelec      | Extract data from MongoDB database to get thermodynamic       |
|               | properties, elastic constant, etc                             |
+---------------+---------------------------------------------------------------+
|   thfind      | Check the dfttk DFT calculation results followed by           |
|               | calling the 'thelec' module to get thermodynamic              |
|               | properties when the option '-get' is given.                   |
+---------------+---------------------------------------------------------------+
|   EVfind      | Find the metadata tags that have 0 K static calculaton        |
|               | finished.                                                     |
+---------------+---------------------------------------------------------------+

Usage of run module:
--------------------

.. code-block:: bash

usage: dfttk run [-h] [-f STRUCTURE_FOLDER] [-mh MATCH_PATTERN] [-s SETTINGS]
                 [-r] [-wf WORKFLOW] [-ph] [-tag TAG] [-a] [-l] [-m [MAX_JOB]]
                 [-o]

optional arguments:

.. code-block:: bash

  -h, --help            show this help message and exit
  -f STRUCTURE_FOLDER, --structure_folder STRUCTURE_FOLDER
                        The folder/file containing the structure, Default: '.'
  -mh MATCH_PATTERN, --match_pattern MATCH_PATTERN
                        The match pattern for structure file, and it should be
                        place in quotes. e.g. '*POSCAR*'. Default: * --
                        everything except SETTING files, ref. -s
  -s SETTINGS, --setting SETTINGS
                        Specify the name of SETTINGS files (yaml or json file)
                        Default: SETTINGS (case insensitive and without ext)
                        The following filename will be treat as SETTINGS file
                        SETTINGS (global settings in the folder) Start with
                        SETTINGS- (individual settings for struct) End with
                        -SETTINGS (individual settings)
  -r, --recursive       Recursive the path.
  -wf WORKFLOW, --workflow WORKFLOW
                        Specify the workflow to run. Default: robust (run
                        get_wf_gibbs_robust workflow) (NOTE: currently, only
                        robust and born are supported.)
  -ph, --phonon         Run phonon. This is equivalent with set phonon=True in
                        SETTINGS file
  -tag TAG, --tag TAG   Specify the tag for continue mode
  -a, --append          Append calculation according to metadata, e.g.
                        appending volumes or phonon
  -l, --launch          Launch the wf to launchpad
  -m [MAX_JOB], --max_job [MAX_JOB]
                        Run the job, only works when -l is specified. Default:
                        0 (Not submit to queue) 1: qlaunch singleshot (single
                        job) N(N>1): qlaunch rapidfire -m N
  -o, --write_out_wf    Write out the workflow


Usage of thelec module:
-----------------------

 - ``thelec`` : process the downloaded from the database, covering phonon/debye approach

  * try ``dfttk thelec -h`` for command available line options

.. code-block:: bash

       dfttk thelec [-h] [-py] [-T0 [T0]] [-T1 [T1]] [-dT [TD]] [-xdn [XDN]]
                    [-xup [XUP]] [-dope [DOPE]] [-ne [NDOSMX]]
                    [-natom [NATOM]] [-e [EVERYT]] [-gauss [GAUSSIAN]]
                    [-i [DOSCAR]] [-o [OUTF]] [-noel] [-metatag [METATAG]]
                    [-qhamode [QHAMODE]] [-pn [PHASENAME]] [-eq [EQMODE]]
                    [-el [ELMODE]] [-s] [-plot] [-g] [-expt [EXPT]]
                    [-xlim [XLIM]]

optional arguments:

.. code-block:: bash

  -h, --help            show this help message and exit
  -py, --pyphon         use Yphon to recalculate vibrational properties.
                        Default: False
  -T0 [T0], -t0 [T0]    Low temperature limit. Default: 0
  -T1 [T1], -t1 [T1]    High temperature limit. Default: 1300
  -dT [TD], -td [TD]    Temperature increment. Default: 10
  -xdn [XDN], --xdn [XDN]
                        Low band energy limit. Default: -100 (eV)
  -xup [XUP], --xup [XUP]
                        High band energy limit. Default: 100
  -dope [DOPE], --dope [DOPE]
                        dope level (electrons). Default: -1.e-8 for numerical
                        stability
  -ne [NDOSMX], --ndosmx [NDOSMX]
                        new DOS mesh. Default: 10001
  -natom [NATOM], --natom [NATOM]
                        number of atoms in the DOSCAR. Default: 1
  -e [EVERYT], --everyT [EVERYT]
                        number of temperature points skipped from QHA
                        analysis. Default: 1
  -gauss [GAUSSIAN], --gauss [GAUSSIAN]
                        densing number near the Fermi energy. Default: 1000
  -i [DOSCAR], --doscar [DOSCAR]
                        DOSCAR filename. Default: DOSCAR
  -o [OUTF], -outf [OUTF]
                        output filename for calculated thermoelectric
                        properties. Default: fvib_ele
  -noel, -noel          do not consider the thermal electron contribution.
                        Default: False
  -metatag [METATAG], -metatag [METATAG]
                        metatag: MongoDB metadata tag field. Default: None
  -qhamode [QHAMODE], -qhamode [QHAMODE]
                        quasiharmonic mode: debye, phonon, or yphon. Default:
                        debye
  -pn [PHASENAME], -phasename [PHASENAME]
                        assigan phase name. Default: None
  -eq [EQMODE], --eqmode [EQMODE]
                        Mode to calculate LTC. 0: Symmetrical Central
                        differential; 4: 4-parameter BM fitting. 5:
                        5-parameter BM fitting. Default: 0
  -el [ELMODE], --elmode [ELMODE]
                        Mode to interpolate thermal electronic contribution:
                        0: interp1d; 1: UnivariateSpline. Default: 0
  -s, -smooth           smooth the LTC. Default: False
  -plot, -plot          plot the figure. Default: False
  -g, --debug           turn on debug mode by reducing the mesh. Default:
                        False
  -expt [EXPT], -expt [EXPT]
                        json file path for experimental thermodynamic
                        properties for plot. Default: None
  -xlim [XLIM], -xlim [XLIM]
                        Up temperature limit for plot. Default: None


Usage of thfind module:
-----------------------

 - ``thfind`` : search what have in the database and batch processing the data

  * try ``dfttk thfind -h`` for command available line options

.. code-block:: bash

       dfttk thfind [-h] [-q [QHAMODE]] [-w [WITHIN]] [-all [CONTAINALL]]
                    [-any [CONTAINANY]] [-v [NV]] [-ss [SUPERCELLN]] [-get]
                    [-py] [-T0 [T0]] [-T1 [T1]] [-dT [TD]] [-xdn [XDN]]
                    [-xup [XUP]] [-dope [DOPE]] [-ne [NDOSMX]]
                    [-natom [NATOM]] [-e [EVERYT]] [-gauss [GAUSSIAN]]
                    [-i [DOSCAR]] [-o [OUTF]] [-noel] [-metatag [METATAG]]
                    [-qhamode [QHAMODE]] [-eq [EQMODE]] [-el [ELMODE]] [-s]
                    [-plot] [-g] [-expt [EXPT]] [-xlim [XLIM]]

optional arguments:

.. code-block:: bash

  -h, --help            show this help message and exit
  -q [QHAMODE], --qhamode [QHAMODE]
                        Collection. 'phonon', 'qha'. Default: 'phonon'
  -w [WITHIN], --within [WITHIN]
                        find calculations within element list Default: None
  -all [CONTAINALL], --containall [CONTAINALL]
                        find calculations must contain all elements in the
                        list Default: None
  -any [CONTAINANY], --containany [CONTAINANY]
                        find calculations contain any elements in the list
                        Default: None
  -v [NV], --nV [NV]    Return phonon calculations finished for number of
                        volumes larger or equals to. Default: 6
  -ss [SUPERCELLN], --supercellsize [SUPERCELLN]
                        only return phonon calculation with supercell size
                        larger than. Default: 0
  -get, --get           get the thermodyamic data for all found entries.
                        Default: False
  -py, --pyphon         use Yphon to recalculate vibrational properties.
                        Default: False
  -T0 [T0], -t0 [T0]    Low temperature limit. Default: 0
  -T1 [T1], -t1 [T1]    High temperature limit. Default: 1300
  -dT [TD], -td [TD]    Temperature increment. Default: 10
  -xdn [XDN], --xdn [XDN]
                        Low band energy limit. Default: -100 (eV)
  -xup [XUP], --xup [XUP]
                        High band energy limit. Default: 100
  -dope [DOPE], --dope [DOPE]
                        dope level (electrons). Default: -1.e-8 for numerical
                        stability
  -ne [NDOSMX], --ndosmx [NDOSMX]
                        new DOS mesh. Default: 10001
  -natom [NATOM], --natom [NATOM]
                        number of atoms in the DOSCAR. Default: 1
  -e [EVERYT], --everyT [EVERYT]
                        number of temperature points skipped from QHA
                        analysis. Default: 1
  -gauss [GAUSSIAN], --gauss [GAUSSIAN]
                        densing number near the Fermi energy. Default: 1000
  -i [DOSCAR], --doscar [DOSCAR]
                        DOSCAR filename. Default: DOSCAR
  -o [OUTF], -outf [OUTF]
                        output filename for calculated thermoelectric
                        properties. Default: fvib_ele
  -noel, -noel          do not consider the thermal electron contribution.
                        Default: False
  -metatag [METATAG], -metatag [METATAG]
                        metatag: MongoDB metadata tag field. Default: None
  -qhamode [QHAMODE], -qhamode [QHAMODE]
                        quasiharmonic mode: debye, phonon, or yphon. Default:
                        debye
  -eq [EQMODE], --eqmode [EQMODE]
                        Mode to calculate LTC. 0: Symmetrical Central
                        differential; 4: 4-parameter BM fitting. 5:
                        5-parameter BM fitting. Default: 0
  -el [ELMODE], --elmode [ELMODE]
                        Mode to interpolate thermal electronic contribution:
                        0: interp1d; 1: UnivariateSpline. Default: 0
  -s, -smooth           smooth the LTC. Default: False
  -plot, -plot          plot the figure. Default: False
  -g, --debug           turn on debug mode by reducing the mesh. Default:
                        False
  -expt [EXPT], -expt [EXPT]
                        json file path for experimental thermodynamic
                        properties for plot. Default: None
  -xlim [XLIM], -xlim [XLIM]

