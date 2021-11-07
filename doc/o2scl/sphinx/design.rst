Design Considerations
=====================

:ref:`O2scl <o2scl>`

Design contents
---------------

- :ref:`Design introduction`
- :ref:`Header file dependencies`
- :ref:`The use of templates`
- :ref:`Error handling`
- :ref:`Define constants and macros`
- :ref:`Parameter ordering`
- :ref:`Global objects`
- :ref:`Thread safety`
- :ref:`Copyright notices`

Design introduction
-------------------

The design goal is to create an object-oriented computing library
with classes that perform common numerical tasks. The most
important principle is that the library should add functionality
to the user while at the same time retaining as much freedom for
the user as possible and allowing for ease of use and extensibility. 
To that end,

- The classes which utilize user-specified functions
  should be able to operate on member functions without requiring
  a particular inheritance structure,
- The interfaces ought to be generic so that the user can create new
  classes which perform related numerical tasks through inheritance.
- The classes should not use static variables or status functions.
- Const-correctness and type-safety should be respected wherever possible.
- The design should be somewhat compatible with GSL.

Header file dependencies
------------------------
    
For reference, it is useful to know how the top-level header files
depend on each other, since it can be difficult to trace everything
down. The following are the most "top-level" header files and their
associated dependencies within O\ :sub:`2`\ scl (there are other
dependencies on GSL and the C standard library not listed here). Note
that not all of the headers in the "base" directory are listed here
(because they are less likely to cause problems)

Tier 1
  - constants.h : (none)
  - err_hnd.h : (none)

Tier 2    
  - exception.h : err_hnd.h
  - misc.h : err_hnd.h
  - find_constants.h : constants.h
  - rng.h (in the mcarlo directory) : err_hnd.h
  - tridiag.h (in the linalg directory) : err_hnd.h

Tier 3
  - string_conv.h : misc.h
  
Tier 4
  - calc_utf8.h : rng.h err_hnd.h string_conv.h find_constants.h
  - columnify.h : misc.h string_conv.h
  - uniform_grid.h : string_conv.h err_hnd.h
  - format_float.h : err_hnd.h misc.h string_conv.h

Tier 5
  - funct.h : err_hnd.h calc_utf8.h
  - vector.h : misc.h uniform_grid.h (vector_special.h)
  - mm_funct.h : calc_utf8.h
  - multi_funct.h : calc_utf8.h

Tier 6    
  - convert_units.h : find_constants.h calc_utf8.h misc.h string_conv.h
    vector.h constants.h
  - search_vec.h : err_hnd.h vector.h misc.h
  - permutation.h : vector.h
  - interp.h : search_vec.h tridiag.h vector.h

Tier 7    
  - lib_settings.h : convert_units.h find_constants.h
    
The interpolation, testing, and table headers are not
as top-level as the ones above because they depend on 
tridiagonalization in the linear algebra directory::

  interp.h : search_vec.h tridiag.h vector.h
  table.h : misc.h interp.h shunting_yard.h
  table_units.h : table.h lib_settings.h
  test_mgr.h : string_conv.h misc.h table_units.h

The use of templates
--------------------
    
Templates are used extensively, and this makes for longer
compilation times so any code that can be removed conveniently
from the header files should be put into source code files
instead. 

Error handling
--------------

Thread safety for errors
^^^^^^^^^^^^^^^^^^^^^^^^

Two approaches to thread-safe error handling which are worth
comparing: the first is GSL which uses return codes and global
function for an error handler, and the second is the Math/Special
Functions section of Boost, which uses a separate policy type for
each function. One issue is thread safety: the GSL approach is
thread safe only in the sense that one can in principle use the
return codes in different threads to track errors. What one cannot
do in GSL is use different user-defined error handlers for
different threads. The Special Functions library allows one to
choose a different Policy for every special function call, and
thus allows quite a bit more flexibility in designing
multi-threaded error handling.

Memory allocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several classes have allocate() and free() functions to allocate
and deallocate memory. If an error occurs in an allocate()
function, the function should free() the partial memory that was
allocated and then call the error handler. Functions which
deallocate memory should never fail and should never be required
to call the error handler. Similarly, class destructors should
never be required to call the error handler.

Define constants and macros
---------------------------

There are a couple define constants and macros that O\ :sub:`2`\ scl
understands, they are all in upper case and begin with the prefix
``O2SCL_``.

Range-checking for arrays and matrices is turned on by default, but
can be turned off by defining ``O2SCL_NO_RANGE_CHECK`` during the
initial configuration of the library. To see how the library was
configured at runtime, use the :cpp:var:`o2scl::o2scl_settings` class.

There is a define constant O2SCL_NO_SYSTEM_FUNC which permanently
disables the shell command ``'!'`` in :ref:`cli <cli>` (when the 
constant is defined, the shell command doesn't work even if
:cpp:var:`o2scl::cli::shell_cmd_allowed` is ``true``). 

The constant O2SCL_DATA_DIR is defined internally to provide the
directory which contains the O\ :sub:`2`\ scl data files. After
installation, this can be accessed in :cpp:var:`o2scl::o2scl_settings`.

All of the header files have their own define constant of
the form ``O2SCL_HEADER_FILE_NAME`` which ensures that
the header file is only included once.

Finally, I sometimes comment out sections of code with::

  #ifdef O2SCL_NEVER_DEFINED
  ...
  #endif

This constant should not be defined by the user as it will cause
compilation to fail.

..
  These are makefile constants not source code define constants

  The two define constants O2SCL_PARTLIB and O2SCL_EOSLIB are used
  internally to control which sublibraries are compiled together
  with the main library (see \ref install_section ). The end-user
  shouldn't have to worry about these.

Parameter ordering
------------------

In functions where this makes sense, generally input parameters will
appear first, while output parameters or parameters which handle both
input and output will appear later.
    
Global objects
--------------

There are four global objects that are created in
libo2scl:
:cpp:var:`o2scl::def_err_hnd` is the default error handler
:cpp:var:`o2scl::alt_err_hnd` is the GSL-like error handler 
:cpp:var:`o2scl::err_hnd` is the pointer to the error handler (points to
:cpp:var:`o2scl::def_err_hnd` by default)
- :cpp:var:`o2scl::o2scl_settings` to control a few library settings

All other global objects are to be avoided.

Thread safety
-------------

Most of the classes are thread-safe, meaning that two instances of
the same class will not clash if their methods are called
concurrently since static variables are only used for compile-time
constants. However, two threads cannot, in general, safely
manipulate the same instance of a class. In this respect, O\
:sub:`2`\ scl is
no different from GSL.
    
.. Documentation design
   --------------------
    
   The commands \\comment and \\endcomment delineate comments about
   the documentation that are present in the header files but don't
   ever show up in the HTML or LaTeX documentation. 

Copyright notices
-----------------

For files where it is appropriate to do so, I have followed the
prescription suggested in
http://lists.gnu.org/archive/html/help-gsl/2008-11/msg00017.html
retaining the GSL copyright notices and putting the O\
:sub:`2`\ scl notices at
the top. CERNLIB has no such standard, but their licensing information
is outlined at
http://cernlib.web.cern.ch/cernlib/conditions.html .

