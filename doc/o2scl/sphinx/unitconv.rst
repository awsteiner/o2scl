:ref:`O2scl <o2scl>`

Unit Conversions
================

There is a class which performs conversion between units specified
with strings in :ref:`convert_units <convert_units>`. If ``popen()``
is supported then the class can use a ``system()`` call to GNU units
to obtain the proper conversion factors. Unit conversions are cached
so that future requests of the same conversion factor do not require
additional ``system()`` calls.
