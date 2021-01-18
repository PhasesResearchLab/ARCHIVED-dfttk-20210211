=========================================
DFTTK: Density Functional Theory Tool Kit
=========================================

High-throughput calculations
----------------------------

By its definition, the term of “first-principles” represents a philosophy that the prediction is to be based on a basic, fundamental proposition or assumption that cannot be deduced from any other proposition or assumption.  This implies that the computational formulations are based on the most fundamental theory of quantum mechanics - Schrödinger equation or density functional theory (DFT) and the inputs to the calculations must be based on well-defined physical constants – the nuclear and electronic charges.  In another word, once the atomic species of an assigned material are known, the theory should predict the energy of all possible crystalline structures, without invoking any phenomenological fitting parameters.  

However, to perform DFT calculations in reality, it still needs the user to have extensive experiences on a variety of parameter choices and a lot of human handling on numerical or system exceptions. In the last decade, we have been working on solving the problem by integrating our experiences accumulated on high-throughout DFT calculations into a software package named as DFTTK (DFT based toolkits) and opened to the community (https://www.dfttk.org). 

The exampled workflows of DFTTK
-------------------------------

1.      Structure maker by protype and elemental substitution;
2.      Robust 0 K equilibrium volume optimization;
3.      0 K energy-volume curve;
4.      Quasiharmonic phonon calculation; 
5.      Born effective charge calculation;
6.      Elastic constant calculations.

DFTTK features
--------------

To perform DFT calculation using DFTTK, the user only needs to name the structure file called POSCAR by VASP, either prepared by user or produced by DFTTK  by elemental substation on given prototype. DFTTK is developed on atomate from the Materials Project4 which is built on three open-source Python libraries. The main benefits of atomate are its flexibility and data management platform, in particular the numerical convergence control and computational exception handling. DFTTK is able to predict properties at finite temperatures by phonon or Debye model for both stoichiometric and solution phases, featured by:

1.      High-throughput DFT calculation and postprocess;
2.      Postprocess plenty of data stored in MongoDB with one simple command;
3.      Compatible with Yphon package and phonopy;
4.      Can recover data from certain fizzled calculations;
5.      Can account thermal electron contribution to thermodynamic properties;
6.      Can calculate thermodynamic properties at 0 K and a few tenth K;
7.      Can perform doping calculations for semiconductors or thermoelectric materials under rigid band approximation;
8.      Can account the effect of thermal expansion/temperature on Seebeck coefficient, Lorenz number, thermal carrier concentrations;
9.      Automatic plot figures for more than 20 thermodynamic properties in the publishable resolution, including atomic volume, free energy, entropy, enthalpy, linear thermal expansion coefficient, isbartic specific heat, constant volume specific heat, lattice only specific heat, bulk modulus, Debye temperature, Seebeck coefficient, Lorenz number, absolute thermal electric force, etc.

