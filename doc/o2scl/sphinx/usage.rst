General Usage
=============

:ref:`O2scl <o2scl>`

General usage contents
----------------------

- :ref:`Namespaces`
- :ref:`Documentation conventions`
- :ref:`Basic error handling`
- :ref:`What is an error?`
- :ref:`GSL error codes and C++ exception types`
- :ref:`Objects and scope`
- :ref:`Define constants`
- :ref:`Error handler objects`

Namespaces
----------
    
Most of the classes and functions reside in the namespace ``o2scl``.
Numerical constants (many of them based on the GSL constants) are
placed in separate namespaces (:ref:`o2scl_const <Namespace
o2scl_const>`). The O₂scl functions and classes for HDF5 output are in
the ``o2scl_hdf`` namespace. There are also two namespaces which hold
integration coefficients, :ref:`o2scl_inte_gk_coeffs <Gauss-Kronrod
integration coefficients>` and :ref:`o2scl_inte_qng_coeffs
<Non-adaptive quadrature integration coefficients>`. There are also
some namespaces for the linear algebra functions, see :ref:`Linear
Algebra` for more information on these.

Documentation conventions
-------------------------

In the following documentation, function parameters are denoted by
``parameter``, except when used in mathematical formulas as in 
:math:`\mathrm{variable}`.

Basic error handling
--------------------

Error handling is a hybrid approach combining GSL with C++ exceptions.
An abstract type has been defined which operates as a GSL-like error
hander. The default error handler is a implementation of this abstract
type which throws a C++ exception when an error is encountered. The
various exceptions, and their correspondence with the GSL error codes,
are given in :ref:`GSL error codes and C++ exception types`. By
default in O₂scl, the default GSL error handler is replaced
with the O₂scl default error handler, i.e. GSL functions
will throw C++ exceptions.

Library users can use the O₂scl handler with the macro ``O2SCL_ERR`` which
calls the O₂scl error handler. The error handler,
:cpp:var:`o2scl::err_hnd` is a global pointer to an object of type
:ref:`err_hnd_type <err_hnd_type>`. There is a global default error
handler, :cpp:var:`o2scl::def_err_hnd` :ref:`err_hnd_cpp
<err_hnd_cpp>`, which throws C++ exceptions, and an alternate default
error handler, :cpp:var:`o2scl::alt_err_hnd`, of type
:ref:`err_hnd_gsl <err_hnd_gsl>`, which outputs an error message and
aborts execution. The global error handler can be replaced by simply
assigning the address of a descendant of :ref:`err_hnd_type
<err_hnd_type>` to :cpp:var:`o2scl::err_hnd`.

.. note::
   The alternate error handler is particularly useful to avoid
   throwing exceptions across DLL boundaries, e.g. when O₂scl
   is loaded as a DLL in python. 
   
..
  11/14/20: This old text regarding exceptions is now unnecessary.

  O₂scl does not support any execution beyond the point at
  which the error handler is called. Many functions which would have
  had integer return values in GSL, now return ``void`` in O\
  :sub:`2`\ scl. Internally, O₂scl does not use ``try`` blocks,
  but these can easily be effectively employed by an O₂scl
  user.

The C++ exception classes are also mapped to the list of GSL error
codes (including a few extra ones for O₂scl), which is
given in below in :ref:`GSL error codes and C++ exception types`.
One can instruct the library to use the GSL-like O₂scl
error handler :cpp:var:`o2scl::alt_err_hnd` by default, by defining
the constant ``O2SCL_USE_GSL_HANDLER``. This is also useful if one
wants to compile without C++ exceptions (e.g. to avoid the small
overhead incurred with C++ exceptions).

What is an error?
-----------------

O₂scl assumes that errors are events which should happen
infrequently. Error handling strategies are often time-consuming
and they are not a replacement for normal code flow. However, even
with this in mind, one can still distinguish a large spectrum of
posibillities from "fatal" errors, those likely to corrupt the
stack and/or cause a dreaded "segmentation fault" and "non-fatal"
errors, those errors which might cause incorrect results, but
might be somehow recoverable. One of the purposes of error
handling is to decide if and how these different types of errors
should be handled differently.

Sometimes, it is undesirable to abort execution upon a failure to
reach numerical convergence. While these failures are treated as
errors (and by default an exception is thrown), some of the classes
which attempt to reach numerical convergence have an option (e.g.
:cpp:var:`o2scl::mroot::err_nonconv`) to turn this default behavior
off for these convergence errors. To set these "convergence" errors in
code provided by the user, the macros ``O2SCL_CONV`` and
``O2SCL_CONV_RET`` can be used. Functions which may have convergence
errors sometimes return ``int``, to indicate which convergence error
was returned when the value of ``err_nonconv`` has been set to false.

Of course, the standard C++ ``try`` and  ``catch`` mechanism of error
handling may also be used for finer-grained control. 

Another related issue is that O₂scl often calls functions
which are supplied by the user, these user-designed functions may
create errors, and the library needs to decide how to deal with them,
even though it knows little about what is actually happening inside
these user-defined functions. For this reason, O₂scl does
not typically try to handle any exceptions or errors occuring in
user-specified functions.

GSL error codes and C++ exception types
---------------------------------------

See also the description of the error codes in ``err_hnd.h``

.. doxygenenumvalue:: success
.. doxygenenumvalue:: gsl_continue

Error codes associated with :ref:`exc_exception <exc_exception>`:
		      
.. doxygenenumvalue:: gsl_failure
.. doxygenenumvalue:: exc_efailed
.. doxygenenumvalue:: exc_esanity
.. doxygenenumvalue:: exc_eunsup
.. doxygenenumvalue:: exc_eunimpl

Error codes associated with :ref:`exc_range_error <exc_range_error>`:
   
.. doxygenenumvalue:: exc_edom
.. doxygenenumvalue:: exc_erange
.. doxygenenumvalue:: exc_eundrflw

Error codes associated with :ref:`exc_runtime_error <exc_runtime_error>`:
   
.. doxygenenumvalue:: exc_efault
.. doxygenenumvalue:: exc_efactor
.. doxygenenumvalue:: exc_enomem
.. doxygenenumvalue:: exc_ebadfunc
.. doxygenenumvalue:: exc_erunaway
.. doxygenenumvalue:: exc_emaxiter
.. doxygenenumvalue:: exc_etol
.. doxygenenumvalue:: exc_eloss
.. doxygenenumvalue:: exc_eround
.. doxygenenumvalue:: exc_esing
.. doxygenenumvalue:: exc_ediverge
.. doxygenenumvalue:: exc_ecache
.. doxygenenumvalue:: exc_etable
.. doxygenenumvalue:: exc_enoprog
.. doxygenenumvalue:: exc_enoprogj
.. doxygenenumvalue:: exc_etolf
.. doxygenenumvalue:: exc_etolx
.. doxygenenumvalue:: exc_etolg
.. doxygenenumvalue:: exc_enotfound
.. doxygenenumvalue:: exc_outsidecons

Error codes associated with :ref:`exc_invalid_argument <exc_invalid_argument>`:
   
.. doxygenenumvalue:: exc_einval
.. doxygenenumvalue:: exc_ebadtol
.. doxygenenumvalue:: exc_ebadlen
.. doxygenenumvalue:: exc_enotsqr
.. doxygenenumvalue:: exc_eindex
		      
Error codes associated with :ref:`exc_overflow_error <exc_overflow_error>`:

.. doxygenenumvalue:: exc_ezerodiv
.. doxygenenumvalue:: exc_eovrflw
		      
Error codes associated with :ref:`exc_ios_failure <exc_ios_failure>`:

.. doxygenenumvalue:: exc_eof
.. doxygenenumvalue:: exc_efilenotfound

Error codes associated with :ref:`exc_logic_error <exc_logic_error>`:

.. doxygenenumvalue:: exc_ememtype
		      
Objects and scope
-----------------
    
O₂scl objects frequently take inputs which are of the form
of a reference to a smaller object. This is particularly convenient
because it allows a lot of flexibility, while providing a certain
degree of safety. In many cases, the user retains the responsibility
of ensuring that input objects do not go out of scope before they are
utilized by objects which require them. This is actually no different
than the requirements on the user imposed by GSL, for example.

Member functions which store pointers to user-specified objects should
warn that they are doing so in the documentation for the class. (AWS:
If they do not, then this is a bug and you should let me know.)

For example, say that a user wants to solve several equations using a
:ref:`mroot_hybrids <mroot_hybrids>` object and use the functions
:cpp:func:`o2scl::mroot_hybrids::set()` and
:cpp:func:`o2scl::mroot_hybrids::iterate()`. Since the function is
specified by the call to ``set()``, it is important that this function
object does not go out of scope before the call to ``iterate()``
occurs.

..
  This is now moved to a design requiremet

  Reference parameters
  --------------------
   
  When a O₂scl function contains two reference parameters for
  objects, it is not typically possible to provide the same object to
  both parameters or to provide two objects which share the same memory.
  This is particularly an issue when the associated types are template
  types, since then the O₂scl library has no way of knowing
  how memory is organized in these unspecified types. Thread safety is
  also an issue, as care must be taken if two functions which are
  running simultaneously access the same instance of any class.

Define constants
----------------

Various define constants used in O₂scl are listed below. Note
that after O₂scl is installed, ``acol -v`` will report
some details about how O₂scl was compiled.

- O2SCL_CBLAS_NAMESPACE - This flag is internally used by some of the
  headers in the ``src/linalg`` directory to select between
  ``o2scl_cblas`` and ``o2scl_cblas_bracket``. The end-user should not
  need to use this.
- O2SCL_COND_FLAG - Internally used to handle the option of compiling
  with or without Armadillo or Eigen support. The end-user should not
  need to use this.
- O2SCL_DATA_DIR - Used by O₂scl classes during
  installation to refer to the directory where data is stored. After
  installation, this is accessible through the global object named
  :ref:`o2scl_settings <o2scl_settings>` and the command ``acol -v``.
  In a typical installation, the end-user should not need to use this.
- O2SCL_DOC_DIR - Used by O₂scl classes during
  installation to refer to the directory where the documentation is
  stored. After
  installation, this is accessible through the global object named
  :ref:`o2scl_settings <o2scl_settings>` and the command ``acol -v``.
  In a typical installation, the end-user should not need to use this.
- O2SCL_HDF5_COMP - Define this flag when O₂scl is compiled
  and for code which uses O₂scl to include support for HDF5
  compression. The command ``acol -v`` reports whether or not HDF5
  compression support was enabled during compilation.
- O2SCL_MPI - Flag to allow MPI functionality in O₂scl classes
  which contain MPI code (see also O2SCL_OPENMP). All current
  MPI functionality in O₂scl is header only, thus MPI support does 
  not need to be specified to the ``configure`` script.
- O2SCL_NEVER_DEFINED - Used internally to comment out large 
  sections of code. This should never be defined by the end-user. 
- O2SCL_NO_EXCEPTIONS - If this is defined, then the error handler
  calls :cpp:func:`o2scl::err_hnd_gsl::set()` instead of throwing a
  C++ exception. Used in ``src/base/exception.cpp``. This is useful,
  for example, with the ``g++`` flag ``-fno-exceptions``.
- O2SCL_NO_SYSTEM_FUNC - If this is defined, then the :ref:`cli <cli>`
  class will never allow shell commands, independent of the 
  setting of :cpp:var:`o2scl::cli::shell_cmd_allowed`.
- O2SCL_NO_RANGE_CHECK - If this is defined, range checking is turned
  off where it is used in :ref:`table <table>`, :ref:`tensor
  <tensor>`, :ref:`permutation <permutation>` and the functions in
  ``o2scl_cblas``. Some O₂scl header files use this define
  constant and so range checking can be turned off or on separately
  from the setting that was used during installation.
- O2SCL_OLDER_COMPILER - A flag which can be defined both during
  compilation of O₂scl and compilation of code using O₂scl
  for compilers which do not have
  ``std::initializer_list``, ``std::to_string``,
  ``std::uniform_distribution`` and cannot run ``mcmc_para_ts``.
- O2SCL_PLAIN_HDF5_HEADER - If true, assume HDF5 include statements 
  should be of the form ``#include <hdf5.h>`` independent
  of the automatically determined operating system type.
- O2SCL_SET_* flags - These flags indicate support for external
  libraries. The header files o2scl/set_*.h automatically define
  these constants if O₂scl was compiled with the associated
  external libraries (i.e. using a ``--enable-*`` flag to
  the ``configure`` script). 
- O2SCL_SNAP - If true, then make some modifications to enable
  compiling O₂scl with snap.
- O2SCL_UBUNTU_PKG - If true, don't use current date and time macros
  to avoid Ubuntu packaging errors (used in
  ``src/base/lib_settings.cpp`` and ``src/hdf/acolm.cpp``). The
  end-user should not need to use this macro.
- O2SCL_USE_BOOST_FILESYSTEM - Doesn't do anything (yet).
- O2SCL_USE_GSL_HANDLER - If this is defined, then an object of type
  :ref:`err_hnd_gsl <err_hnd_gsl>` is the default error handler. Used
  in ``src/base/exception.cpp``

Error handler objects
---------------------

These objects are defined in ``src/base/err_hnd.h``:

.. doxygenvariable:: o2scl::err_hnd

.. doxygenvariable:: o2scl::alt_err_hnd
		     
This object is defined in ``src/base/exception.h``:

.. doxygenvariable:: o2scl::def_err_hnd

