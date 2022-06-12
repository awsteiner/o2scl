Physical Constants
==================
    
:ref:`O2scl <o2scl>`

Constant contents
-----------------

- :ref:`Constant introduction`
- :ref:`Namespace o2scl_const`
- :ref:`Namespace o2scl_mks`
- :ref:`Namespace o2scl_cgs`
- :ref:`Namespace o2scl_cgsm`

Constant introduction
---------------------
     
The constants from GSL are reworked with the type ``const double`` and
placed in namespaces called :ref:`o2scl_mks <Namespace o2scl_mks>`,
:ref:`o2scl_cgs <Namespace o2scl_cgs>`, and :ref:`o2scl_mksa
<Namespace o2scl_cgsm>`. The GSL MKSA constants are identical to the
MKS constants and thus are not duplicated here. The numerical
constants from ``gsl_num`` and some other additional constants are
given in the namespace :ref:`o2scl_const <Namespace o2scl_const>`,
Some of the numerical values have been updated with CODATA 2018
values.

The :ref:`find_constants <find_constants>` class contains a
simple constant database which can be searched at compiled time
and also provides the constant database to ``acol -constants``
(see :ref:`The acol Command Line Utility`).

These physical constants can also be used to create unit conversion
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
   
