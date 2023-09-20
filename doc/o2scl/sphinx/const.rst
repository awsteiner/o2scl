Physical Constants
==================
    
:ref:`O2scl <o2scl>`

Constant contents
-----------------

- :ref:`Constant introduction`
- :ref:`Namespace o2scl_const`

Constant introduction
---------------------

In order to avoid confusing numerical differences when using
multiprecision arithmetic, physical constants are template functions
which return a value given a user-specified floating-point type.
Physical constants are promoted to higher precision by adding zeros in
the base-10 representation. The constants are in the namespace
:ref:`o2scl_const <Namespace o2scl_const>`. The numerical values are
periodically updated with CODATA, the Particle Data Book, and other
databases. A handful of constants have a short name for their
values in double precision, see below. 

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

   
