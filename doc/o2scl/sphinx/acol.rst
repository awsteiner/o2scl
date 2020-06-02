:ref:`O2scl <o2scl>`

The acol Command-line Utility
=============================

O\ :sub:`2`\ scl contains a command-line utility, \c acol, designed to
facilitate the manipulation of various objects stored in HDF5 files.
It can handle integers, characters, double-precision floating point
numbers, size_t objects, arrays of any of these types, or O\ :sub:`2`\
scl objects of type :ref:`table <table>`, :ref:`hist <hist>`,
:ref:`table3d <table3d>`, :ref:`hist_2d <hist_2d>`, :ref:`tensor_grid
<tensor_grid>` :ref:`contour_line <contour_line>` .

``acol`` can only operate with one object at a time. The
basic workflow is:

- create an object from scratch or read it from an HDF5 file
- perform operations on that object
- output the object to the screen or write it to an HDF5 file

The available command list can be obtained using ``'help'`` or
``'commands'`` and changes depending on what type of object is
currently in memory. In order to list the commands which would be
available given a particular type, give ``'commands'`` the type as an
argument, i.e. ``acol -commands table`` . In order to
get detailed help on how a command operates on a particular type,
give the type and the command as arguments to help, e.g.
``acol -help table interp``. There are some commands
which are available for all types, and obtaining the help
information for these commands does not require a type argument,
i.e. ``acol -commands`` or ``acol -help read``.

``acol`` can sometimes, but not always read and write HDF5
files generated outside of O\ :sub:`2`\ scl.

``acol`` has a command, ``run``, which allows you to run
a set of commands which are given in a separate file. An example
script in the ``extras`` directory of the documentation is 
named ``acol.scr``. The associated output is a useful demonstration
of the capabilities of ``acol``.

.. literalinclude:: static/acol.out
