Physical Constants
==================
    
:ref:`O2scl <o2scl>`

The constants from GSL are reworked with the type ``const double`` and
placed in namespaces called :ref:`o2scl_mks <Namespace o2scl_mks>`,
:ref:`o2scl_cgs <Namespace o2scl_cgs>`, and :ref:`o2scl_mksa
<Namespace o2scl_cgsm>` The GSL MKSA constants are identical to the
MKS constants and thus are not duplicated here. The numerical
constants from ``gsl_num`` and some other additional constants
are given in the namespace :ref:`o2scl_const <Namespace o2scl_const>`,

Some of the numerical values have been updated from recently
released data. Electron, neutron, proton, and atomic mass have
been updated with CODATA 2018 values. Also electron charge,
gravitational constant, plancks_constant_hbar, are updated. The
astronomical unit has been updated with the result from \ref
Luzum11 (and possibly other values need updating as well).

These physical constants are often used to create unit conversion
factors, described in :ref:`Unit conversions`.

Namespace o2scl_const
---------------------

:ref:`Top <Physical Constants>`

.. doxygennamespace:: o2scl_const

Namespace o2scl_mks
-------------------

:ref:`Top <Physical Constants>`

.. doxygennamespace:: o2scl_mks
   
Namespace o2scl_cgs
-------------------

:ref:`Top <Physical Constants>`

.. doxygennamespace:: o2scl_cgs
   
Namespace o2scl_cgsm
--------------------

:ref:`Top <Physical Constants>`

.. doxygennamespace:: o2scl_cgsm
   
