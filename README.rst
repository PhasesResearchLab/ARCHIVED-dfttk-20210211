=========================================
DFTTK: Density Functional Theory Tool Kit
=========================================

**Ultimate goals:** High-throughput postprocessing the result produce by the DFTTK package
Copyright © Phases Research Lab (https://www.phaseslab.com/)

This branch is based the robust_relax branch

- The following extensions are currently implemented:
 - thelec : process the downloaded from the database, covering phonon/debye approach
  * try "dfttk thelec -h" for command available line options
 - thfind : search what have in the database and batch processing the data
  * try "dfttk thfind -h" for command available line options

- Features
 - High-throughput DFT calculation and postprosess.
 - Can postprocess plenty of data stored in MongoDB with one simple command.
 - Can look through the MongoDB database to see what phases have been calculated by various conditions
 - Compatible with Yphon package (https://www.sciencedirect.com/science/article/pii/S0010465514002288) and
   phonopy (https://phonopy.github.io/phonopy/)
 - Included more than 10 examples in the Example folder
 - Can recover corrupted data that are not completely finished by the workflow
 - Can handle the thermal electron contribution to the thermodynamic properties
 - Can plot most of common thermodynamic properties, phonon dispersion/dos  in the publishable accuracy
 - Avoid numerical differential as much as possible, made data smoothness greatly enhanced
 - Thermal expansion coefficients are calculated by thermodynamic relation with respect to volume
 - Can use the parabolic description of phonon DOS in low frequency region
 - Can handle phonons for temperature at 0 K and a few tenth K for which the original codes mostly failed
 - Can picked the results from correpted calculations due to reasons such as tow wide temperature range
 - temperature can be made parabolic with option "-td -5" for better accuracy and efficiency
 - Can call the thelec module from the thfind module for batch processing
 - Dope under rigid band approximation made it fast of thermoelectric properties
 - Can report the effect of thermal expansion/temperature on Seebeck coefficient, Lorenz number, thermal carrier concentrations
 - Can produce/plot more than 20 thermodynamic properties, as a function of temperature, including:
   atomic volume
   free energy
   entropy
   enthalpy
   lineat thermal expansion coefficient
   isbartic specific heat
   constant volume specific heat
   lattice only isbartic specific heat
   bulk modulus
   Debye temperature
   Debye temperature by Debye model (if it can find)
   Heat capacity diveded temperature vs square of temperature
   electronic contribution to the free energy
   electronic contribution to the entropy
   electronic contribution to the specific heat
   Seebeck coefficient
   Lorenz number
   Absolute thermal electric force



.. code-block:: bash

   Usage of thelec module:
   usage: dfttk thelec [-h] [-py] [-T0 [T0]] [-T1 [T1]] [-dT [TD]] [-xdn [XDN]]
                       [-xup [XUP]] [-dope [DOPE]] [-ne [NDOSMX]]
                       [-natom [NATOM]] [-e [EVERYT]] [-gauss [GAUSSIAN]]
                       [-i [DOSCAR]] [-o [OUTF]] [-noel] [-metatag [METATAG]]
                       [-qhamode [QHAMODE]] [-pn [PHASENAME]] [-eq [EQMODE]]
                       [-el [ELMODE]] [-s] [-plot] [-g] [-expt [EXPT]]
                       [-xlim [XLIM]]

   optional arguments:
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
   usage: dfttk thfind [-h] [-q [QHAMODE]] [-w [WITHIN]] [-all [CONTAINALL]]
                       [-any [CONTAINANY]] [-v [NV]] [-ss [SUPERCELLN]] [-get]
                       [-py] [-T0 [T0]] [-T1 [T1]] [-dT [TD]] [-xdn [XDN]]
                       [-xup [XUP]] [-dope [DOPE]] [-ne [NDOSMX]]
                       [-natom [NATOM]] [-e [EVERYT]] [-gauss [GAUSSIAN]]
                       [-i [DOSCAR]] [-o [OUTF]] [-noel] [-metatag [METATAG]]
                       [-qhamode [QHAMODE]] [-eq [EQMODE]] [-el [ELMODE]] [-s]
                       [-plot] [-g] [-expt [EXPT]] [-xlim [XLIM]]

   optional arguments:
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

