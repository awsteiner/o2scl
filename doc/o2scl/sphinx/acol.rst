The acol Command-line Utility
=============================

:ref:`O2scl <o2scl>`

acol contents
-------------

- :ref:`acol introduction`
- :ref:`acol functions`
- :ref:`acol types`
- :ref:`Value specifications`
- :ref:`Vector specifications`
- :ref:`String list specifications`
- :ref:`Multiple vector specifications`
- :ref:`acol example`
  
acol introduction
-----------------
  
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

acol functions
--------------

Some ``acol`` commands can work with mathematical function arguments.
Functions can be created using the operators and functions listed below.
Examples are ``x==5 && y<1``, ``acos(-1)``, and ``sin(x>5)``.
Comparison operators result in either 1.0 (true) or 0.0 (false).

Operators: ``() ^ * / % + - == != < > && || << >> >= <=``

Functions: ``exp(x) log(x) log10(x) sin(x) cos(x) tan(x) sqrt(x) abs(x)
asin(x) acos(x) atan(x) sinh(x) cosh(x) tanh(x) asinh(x) acosh(x)
atanh(x) floor(x)``

There is also a command called ``function`` which works with several
different types to generate data based on functions. Use ``acol -help
function`` to get more information on these type-specific commands.

acol types
----------

The types which can be handled by ``acol`` are either related to C++
internal types: ``char, double, double[], int, int[], size_t,
size_t[], std::string``, and ``std::string[]``, or O\ :sub:`2`\ scl
types :ref:`hist <hist>`, :ref:`hist_2d <hist_2d>`,
:ref:`prob_dens_mdim_amr <prob_dens_mdim_amr>`, :ref:`table <table>`,
:ref:`table3d <table3d>`, :ref:`tensor <tensor>` (including
``double``, ``int``, and ``size_t`` versions), :ref:`tensor_grid
<tensor_grid>` :ref:`uniform_grid <uniform_grid>`, and
vector<:ref:`contour_line <contour_line>`>. 

Value specifications
--------------------

Some acol commands take "value specifications" as arguments, which are
string expressions which represent a fixed numerical value. These
specifications are handled by :cpp:func:`o2scl_hdf::value_spec()`. The
first part of the specification is a "type of specification" followed
by a colon, followed by arguments which depend on the type. If no
colon is present, then a "func:" prefix is assumed. The different
types for a value specification are:

1. ``<numeric value or function>`` - Value equal to the result of
<function>, e.g. "7.6" or "sin(0.5)". See :ref:`acol Functions` for a
list of functions that can be used.

2. ``hdf5:<object name>:[addl. spec.]`` - Read an HDF5 value and obtain
the value from object named <object name>. For some object types,
additional specifications are required to specify which value should
be used. A list of object types and additional specifications and more
detail is given below.

====================== ======================== ========================
Type                   Additional specification Description
====================== ======================== ========================
               double  (none)
                  int  (none)
               size_t  (none)
             double[]  index
                int[]  index
             size_t[]  index
 uniform_grid<double>  index
                table  column name, row index
====================== ======================== ========================

3. ``shell:<shell command>`` - Set the value equal to the first result
obtained using the specified shell command.

Vector specifications
---------------------

Some acol commands take "value specifications" as arguments, which are
string expressions which represent a list of numerical values. These
specifications are handled by :cpp:func:`o2scl_hdf::vector_spec()`.
The different parts of the string are separated by a colon, and the
first part specifes the type of vector specification. The different
types are:

1. ``val:<value>`` - Create a vector with one element equal to <value>,
which may be a number or a simple function, e.g. ``val:sin(0.5)``.

2. ``list:<entry 0>,<entry 1>, ..., <entry n-1>`` - Create a vector with
a simple list of numbers or functions, e.g.
'list:3.0,1.0e-3,sqrt(2.0)'.

3. ``func:<N>:<function of i>`` - Create a vector by specifying the
length of the vector and a function used to fill the elements. For
example: 'func:41:sin(i/20.0*acos(-1))'.

4. ``grid:<begin>,<end>,<width>,["log"]`` - Create a vector equal to a
uniform grid, e.g. use 'grid:1.0,10.0,1.0' for a 10-element vector
filled with the numbers 1 to 10.

5. ``text:<filename>:<column index>`` - Read a text file and extract a
vector of numbers from a column of the text file (starting with zero
for the first column), ignoring any header rows which contain
non-numeric values. For example 'text:~/temp.dat:2' will construct a
vector from the third column of the file 'temp.dat' in the user's home
directory (using :cpp:func:`wordexp_single_file()` which calls the
system ``wordexp()`` function to expand the tilde).

6. ``hdf5:<file name>:<object name>:[addtional spec.]`` - Read an HDF5
file and obtain a vector from the object with the specified name. The
remaining parts of the string contain additional information which may
be needed depending on the type of object stored in the HDF5 file. A
list of object types and additional specifications and more detail is
given below.

==================== ======================== =======================
Type                 Additional specification Description
==================== ======================== =======================
              double (none)                   Implies vector of size 1
            double[] (none)
                hist (none)                   Vector of histogram weights
                 int (none)                   Implies vector of size 1
               int[] (none)
              size_t (none)                   Implies vector of size 1
            size_t[] (none)
               table <column>                 Selected column from table
               table <row> <col pat>          Selected row and columns  
uniform_grid<double> (none)
==================== ======================== =======================

For table <row> <col pat>, the first additional specification is a row
number, which can be negative to refer to counting from the end of the
table. The second additional specification is a pattern of column
names using either '*' or '?'.

String list specifications
--------------------------

Some acol commands take "string list specifications" as arguments,
which are string expressions which represent a list of strings. These
specifications are handled by :cpp:func:`o2scl_hdf::strings_spec()`.
The different parts of the string are separated by a colon, and the
first part specifes the type of vector specification. The different
types are:

1. ``list:<comma-separated list>`` - A list of strings

2. ``shell:<command>`` - The lines obtained from the result of a shell
command, with a maximum of 256 characters per line.

3. ``pattern:N:x[0][a][A]`` - The N strings obtained from a pattern.
Occurrences of [0] are replaced with the integer 'i' where i runs from
0 to N-1. Occurrences of [a] are replaced with 'a' through 'z' from 0
through 25, and 'aa' through 'zz' for i from 26 to 701. Occurrences of
[A] are replaced with 'A' through 'Z' from 0 through 25, and 'AA'
through 'ZZ' for i from 26 to 701.

4. hdf5: - Unfinished.

Multiple vector specifications
------------------------------

Some acol commands take "multiple vector specifications" as arguments,
which are string expressions which represent a list of vectors (which
need not have the same length). These specifications are handled by
:cpp:func:`o2scl_hdf::mult_vector_spec()`. The different parts of the
string are separated by a colon, and the first part specifes the type
of multiple vector specification. The different types are:

1. ``func:<N>:<function of i>:<function of i and j>`` - Specify the
number of vectors, a function of "i" which determines the length of
the ith vector, and a function of "i" and "j" which specifies the jth
element of the ith vector.

2. ``text:<filename pattern>:<numeric column list>`` - Read one or
more text files and extract vectors of numbers from columns of the
text file, ignoring any header rows which contain non-numeric values.
For example 'text:~/temp.dat:2-4' will construct vectors from the
third, fourth, and fifth columns of the file 'temp.dat' in the user's
home directory.

3. ``hdf5:<filename pattern>:<object name>:[additional spec.]`` - Read
one or more HDF5 files and obtain a vector from the object with the
specified name. The remaining parts of the string contain additional
information which may be needed depending on the type of object stored
in the HDF5 file.

==================== ======================== =======================
Type                 Additional specification Description
==================== ======================== =======================
              table  <column pattern>
==================== ======================== =======================

Also, many normal vector specifications (from 'acol -help
vector-spec') also work as multiple vector specifications. These
include specifications which begin with 'val:', 'list:', 'grid:', and
'table-row:'. Also included are 'hdf5:' specifications which refer to
objects of type double, double[], hist, int, int[], size_t, size_t[],
and uniform_grid<double>.
     
acol Example
------------

.. literalinclude:: static/acol.out
