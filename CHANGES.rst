=========
Changelog
=========

0.2 (2020-03-30)
================

New features

(Contributor: @bocklund_ , @Peng_Gao, @hitliaomq_ )

* The relax scheme is optimized. (from ``ISIF=3`` to ``ISIF=2`` followed by ``ISIF=4``) (@Peng_Gao)
* Change the static workflow to dynamic workflow. (``EVcheck_QHA.py`` increase the data points atomately if the fitting of initial points is incorrect) (@Peng_Gao)
* Support run dfttk by command. (Add ``dfttk run [options]``) (@hitliaomq_)
* Support configrate dfttk automately. (Add ``dfttk config [options]``) (@hitliaomq_)
* Documents' enhance. (@hitliaomq_)
* Bug fix. (@bocklund_, @Peng_Gao, @hitliaomq_)

0.1 (2018-08-28)
================

Initial release. Includes

(Contributor: @bocklund_, @mxf469_)

* Gibbs workflow for stable structures
* Analysis code and libraries for calculation quasiharmonic Gibbs energies with 0K, vibrational and thermal electronic contributions
* Useful utilities for interfacing with structure, calculations and the Materials Project

.. _bocklund: https://github.com/bocklund
.. _mxf469: https://github.com/mxf469
.. _hitliaomq: https://github.com/hitliaomq