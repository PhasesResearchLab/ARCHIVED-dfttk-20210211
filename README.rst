=========================================
DFTTK: Density Functional Theory Tool Kit
=========================================

**Ultimate goals:** High-throughput postprocessing the result produce by the DFTTK package
Copyright © Phases Research Lab (https://www.phaseslab.com/)

This branch is based the robust_relax branch

- The following extensions are currently implemented:

    - thelec process the downloaded from the database, covering phonon/debye approach

        * try "dfttk thelec -h" for command available line options

    - thfind search what have in the database and batch processing the data

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
 - Can produce/plot more than 20 thermodynamic properties, as a function of temperature, including: atomic volume, free energy, entropy, enthalpy, lineat thermal expansion coefficient, isbartic specific heat, constant volume specific heat, lattice only isbartic specific heat, bulk modulus, Debye temperature, Debye temperature by Debye model (if it can find), Heat capacity diveded temperature vs square of temperature, electronic contribution to the free energy, electronic contribution to the entropy, electronic contribution to the specific heat, Seebeck coefficient, Lorenz number, Absolute thermal electric force
