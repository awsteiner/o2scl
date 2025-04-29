Library Settings
================

:ref:`O2scl <o2scl>`

There are a couple library settings which are handled by a global
object of type :ref:`lib_settings_class <lib_settings_class>` (see
below).

There are several data files that are used by various classes in the
library. The installation procedure should ensure that these files are
automatically found. However, if these data files are moved after
installation, then a call to
:cpp:func:`o2scl::lib_settings_class::set_data_dir()` can adjust the
library to use the new directory. It is assumed that the directory
structure within the data directory has not changed.

If you have changed the directory structure of the O₂scl
data files, then many of the I/O classes allow you to specify the full
pathname of the file you wish to load.

Global settings object
----------------------

.. _o2scl_settings:

.. doxygenvariable:: o2scl::o2scl_settings
