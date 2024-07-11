Installation
============

:ref:`O2scl <o2scl>`

Installation Contents
---------------------

- :ref:`General notes <install_general>`
- :ref:`Compiling O₂scl from a release distribution <compile_dist>`
- :ref:`Compiling O₂scl from a release on Linux <compile_release>`
- :ref:`Compiling O₂scl from the source code <compile_source>`
- :ref:`Compiling O₂scl on Docker <compile_docker>`
- :ref:`Python support <python_support>`  
- :ref:`Optional linear algebra libraries`
- :ref:`Other optional libraries`  
- :ref:`Range-checking`
- :ref:`More configure flags`
- :ref:`Generation of documentation`
- :ref:`Uninstallation`

.. _install_general:
   
General notes
-------------

O₂scl requires Boost, GSL, libquadmath, and the HDF5 libraries (the
precise procedure for installing these libraries differs from system
to system, but some common cases and useful information is given
below). O₂scl is designed to be used with the most recent release
version of all of these libraries, but is sometimes compatible with
recent older versions. The configure script attempts to add these
libraries to LIBS and LDFLAGS during the installation of O₂scl. In
order to compile your code with O₂scl, you will need to include, e.g.
``-lo2scl -lhdf5 -lgsl -lgslcblas -lquadmath -lm``, and you may need
to include ``-I`` flags for O₂scl headers and ``-L`` flags for O₂scl
libraries. The sections below describe several different ways of
installing O₂scl.

It is important to ensure that O₂scl is compiled with the same version
of the HDF5 libraries that it is linked with when compiling code based
on O₂scl. In order to help resolve these version conflicts, the
``acol`` utility (see :ref:`The acol Command Line Utility`) reports
the two different HDF5 versions (see ``acol -v``) so that it is easy
to check that they are the same. 

..
  x 
  x- :ref:`Compiling O₂scl on Mac OSX with Homebrew <compile_homebrew>`
  x.. _compile_homebrew:
    
  Compiling O₂scl on Mac OSX with Homebrew
  ----------------------------------------
  
  The easiest way to install on Mac OSX is with homebrew. Use::
  
    brew tap awsteiner/science
    brew install o2scl
  
  to install O₂scl. There are a few options for ``brew install``. The
  option ``--with-check`` performs the build-time tests and the option
  ``--with-examples`` double checks that the examples can also be
  compiled and executed. The homebrew recipe for O₂scl uses the Mac OS X
  compiler clang. Homebrew also supports the installation of the current
  version directly from the repository using the ``--HEAD`` option to
  ``brew install``. The homebrew installation includes readline support.
  The O₂scl homebrew recipes are stored at the
  xhttps://github.com/awsteiner/homebrew-science repository.
  
  By default, a homebrew installation of O₂scl uses the OSX LLVM
  compiler. However, a homebrew installation of O₂scl will also install
  ``gcc`` because O₂scl requires ``hdf5``, and the homebrew ``hdf5``
  package requires ``gcc``.
  
  Python support in the homebrew package does not yet work yet.

.. _compile_dist:

Compiling O₂scl from a release distribution
-------------------------------------------

O₂scl installation is generally similar to that for GNU-style
libraries. The file ``INSTALL`` has some details on this procedure.
Once the dependencies are installed you should be able to run
``./configure`` and then type ``make`` and ``make install``. More
information on the ``configure`` command can also be obtained from
``./configure --help``. On some systems, you may have to add
additional flag to the ``CXXFLAGS`` environment variable manually
before the ``./configure`` script. The documentation is included in
the O₂scl release distribution and automatically installed by ``make
install``.

.. note::
   If you are trying to install O₂scl with a version of
   HDF5 earlier than 1.12 you will need to compile with
   ``-DO2SCL_HDF5_PRE_1_12``.

O₂scl requires the Boost (v1.80.0 or later) and the GSL libraries
(version 2.0 or later). If the ``configure`` script cannot find Boost
or GSL, you may have to specify their location for the associated
header files in the ``CXXFLAGS`` variable and the associated libraries
in the ``LDFLAGS`` environment variable. Running ``./configure
--help`` shows some information on this. For example, in a bash shell,
you could do something like::

  CXX="g++" CXXFLAGS="-I/dir/to/gsl/include" LDFLAGS="-L/dir/to/gsl/libs" ./configure --prefix=="/dir/to/destination_directory

Along with GSL, a CBLAS library is also required, and ``./configure``
will look for ``libcblas`` first, and if not found then it will look
for ``libgslcblas``. If neither is present, then you may have to
manually specify a CBLAS library using the ``LIBS`` and ``LDFLAGS``
environment variables.

Compiling with the readline library is optional, but it is assumed to
be present by default.

After ``make install``, you may test the library with ``make check``
or ``make o2scl-test``. At the end, the phrase ``"All O2scl tests
passed"`` indicates that the testing was successful. You may also run
``make o2scl-test`` in the individual subdirectories of the src
directory to individually test the classes and functions in that part
of O₂scl. After installation, running ``acol -v`` will output several
of the installation settings.

.. _compile_release:

Compiling O₂scl from a release on Linux
---------------------------------------

For example, to install O₂scl on Ubuntu, begin by installing g++ and
make (the ``g++`` and ``make`` packages), GSL (the ``libgsl-dev``
package), Boost (the ``libboost-all-dev`` package), GNU readline (the
``libreadline-dev`` package), HDF5 (the ``libhdf5-dev`` package), and
quadmath (the ``libquadmath0`` package). You can then install O₂scl
from one of the release distributions by using the standard GNU
``./configure`` script and then invoking ``make`` and ``make install``
(which often requires ``sudo``).
 
The HDF5 package for Ubuntu and many other Linux systems is installed
in ``hdf5/serial/hdf5.h`` instead of ``hdf5.h``, so O₂scl presumes
that Linux systems are arranged that way. If HDF5 include statements
should not have the ``hdf5/serial/`` prefix, then you can use
``-DO2SCL_HDF5_PLAIN_HEADER``, i.e.::

  CXXFLAGS="-DO2SCL_PLAIN_HDF5_HEADER" ./configure

to instruct O₂scl to look for them there (for example, on bridges at
the PSC). On many systems, one can use a parallel HDF5 library using
``-DO2SCL_HDF5_PLAIN_HEADER`` and a ``-I`` option to select the proper
location for the parallel HDF5 header files. Finally, if your version
of HDF5 is earlier than 1.12, you will need to let O₂scl know, using::

  CXXFLAGS="-DO2SCL_HDF5_PRE_1_12" ./configure

Other Linux distributions are similar. For example, in OpenSUSE, you
will need to use ``zypper`` to install ``gcc-c++, make, gsl-devel,
hdf5-devel, readline-devel``, ``libquadmath0``, and ``boost-devel``.

Note that if your boost installation is earlier than 1.70, you will
need to use the -DO2SCL_OLD_BOOST flag to get all of the tests to run
successfully.

.. _compile_source:

Compiling O₂scl from the source code
------------------------------------

If you want to install from source (without generating the
documentation), then you must first install ``g++``, ``make``,
``automake``, ``autoconf``, and ``libtool`` packages. You also need to
install all the dependencies described above (see, e.g. the section
:ref:`Compiling O₂scl from a release on Linux`). Then you can use
something along the lines of::

  git clone https://github.com/awsteiner/o2scl
  cd o2scl
  autoreconf -i
  ./configure

Then, you will either need to generate the documentation from doxygen
using ``make o2scl-doc`` or use ``make blank-doc`` to create blank
documentation. Then you can proceed using ``make`` and ``make
install`` (which may require ``sudo`` depending on your
configuration). For a full installation with parallelism, I typically
also install ``libopenmpi-dev`` and then use ``./configure
--enable-openmp``

.. _compile_docker:

Docker images for O₂scl
-----------------------

There are a few docker images for recent versions of O₂scl available
at https://hub.docker.com/r/awsteiner/o2scl . These images are based
on the docker files which are stored in the ``docker`` subdirectory,
and can be found at
https://github.com/awsteiner/o2scl/tree/main/docker . 

.. _python_support:

Python support
--------------

O₂scl can be compiled with python support by providing the option
``--enable-python`` when the library is configured. This may also
require adjusting CXXFLAGS and LDFLAGS in order to ensure the Python
headers and libraries can be found. O₂scl code which uses Python also
assumes that numpy was installed, so the headers for the numpy package
may need to be specified. For example, using g++ on MacOS may need
something of the form::

  CXX="g++-13"
  CXXFLAGS="-I/usr/local/lib/python3.11/site-packages/numpy/core/include
  `python3-config --includes`" LDFLAGS="`python3-config --ldflags`"
  ./configure --enable-python

Including Python support also requires the installation of O₂sclpy
(for example, using \c pip) to ensure that the tests pass
successfully. Thus, when including Python support it is best to
install O₂scl first, install O₂sclpy second, and then test O₂scl and
O₂sclpy last. See also :ref:`Python Integration` for more details.

Some of the docker images available at
https://hub.docker.com/r/awsteiner/o2scl include an installation of
O₂sclpy and O₂scl with Python support. 

.. 
  x .. _compile_snap:
  x- :ref:`Compiling O₂scl on Ubuntu with Snap <compile_snap>`

  Compiling O₂scl on Ubuntu with Snap
  -----------------------------------
  
  .. note:: AWS, 6/23/23: The snap package needs some work and I have
            not had the time to fix it yet.
  
  The easiest way to install on Ubuntu is with snap (see
  https://snapcraft.io/o2scl). Use::
  
    sudo snap install (--edge or --beta) --devmode o2scl
  
  The snap installation includes readline support and uses the GSL CBLAS.
  
  Using the command-line utility ``acol`` may require you to set the
  environment variable ``LD_LIBRARY_PATH``. For example, on machines
  where I use snap to install in my ``.bashrc``, I use::
  
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/snap/o2scl/current/usr/lib/x86_64-linux-gnu:/snap/o2scl/current/lib/x86_64-linux-gnu

Optional linear algebra libraries
---------------------------------

O₂scl is fully functional without any additional linear algebra
libraries. However, many classes and functions which require linear
algebra are faster with either the Eigen (http://eigen.tuxfamily.org)
or Armadillo (http://arma.sourceforge.net) libraries. Support for
these can be specified in the ``configure`` command with
``--enable-armadillo`` or ``--enable-eigen``. These libraries may
require additional ``-I`` or ``-L`` flags to be defined when O₂scl is
installed, depending on how your particular system is configured. For
example, O₂scl classes which use Armadillo use matrix decompositions
so Armadillo must be compiled with LAPACK support, and you may need to
specify the location of the LAPACK libraries manually.

..
  If you are
  installing on Mac OS X with homebrew, the options ``--with-eigen`` and
  ``with-armadillo`` can be used.

Other optional libraries
------------------------

As with the linear algebra libraries, these libraries may require
additional ``-I`` or ``-L`` flags to be defined when O₂scl is
installed, depending on how your particular system is configured. The
configure script should automatically add ``-l<library name>`` to
LDFLAGS during installation, but you will need to also add this flag
to your codes which use O₂scl.

Readline support (``-lreadline``): The command-line interface class
:ref:`cli <cli>`, and ``acol`` (see :ref:`The acol Command Line
Utility`) can both take advantage of readline support. If the library
is configured with ``--disable-readline``, then the readline library
is not used.

OpenMP support (typically involves the ``-fopenmp`` compiler flag):
O₂scl contains a few functions which use multiple threads for
faster execution. This support can be included using the
``-enable-openmp`` option to the configure script. On some systems,
this will also include explicitly specifying the OpenMP libraries
in the ``LDFLAGS`` environment variable. See more information in
:ref:`Parallel Programming with O2scl`. 
  
FFTW support (``-lfftw3``): O₂scl contains a few functions which
require FFTW support, and this can be included if ``--enable-fftw`` is
passed to the configure script.

Module support, curses support, MFPR support, cubature support, and
pugixml support are all experimental.

Range-checking
--------------

Some extra range-checking for vectors and matrices is turned on by
default. You can disable range-checking by defining
-DO2SCL_NO_RANGE_CHECK, e.g.::

  CXXFLAGS="-DO2SCL_NO_RANGE_CHECK" ./configure

More configure flags
--------------------

There are several warning flags that are useful when configuring
and compiling with O₂scl. See the GSL documentation for an 
excellent discussion, and also see the generic installation
documentation in the file ``INSTALL`` in the O₂scl top-level 
directory. For running ``configure``, for example, if you do
not have privileges to write to ``/usr/local``::

  CPPFLAGS="-O3 -I/home/asteiner/install/include" \
  LDFLAGS="-L/home/asteiner/install/lib" ./configure \
  --prefix=/home/asteiner/install

In this example, specifying ``-I/home/asteiner/install/include`` and
``-L/home/asteiner/install/lib`` above ensures that the GSL libraries
can be found. The ``--prefix=/home/asteiner/install`` argument to
``./configure`` ensures that O₂scl is installed there as well.

Generation of documentation
---------------------------

The O₂scl documentation is included in every release tarball.
It can also be generated by the user, but requires ``doxygen``, ``sphinx``,
``breathe``, ``alabaster``, ``pugixml`` and other external
applications not included in the distribution.

The most recent release documentation is available at
https://awsteiner.org/code/o2scl/html/index.html and the
current development version documentation is available at
https://awsteiner.org/code/o2scl-dev/html/index.html . The
documentation for previous releases is not on the web, but is still
stored in the release ``.tar.gz`` file.

Uninstallation
--------------

While there is no explicit "uninstall" makefile target, there are only
a couple places to check. Installation creates directories named
``o2scl`` in the include, doc and shared files directory (which
default to ``/usr/local/include``, ``/usr/local/share/doc/``, and
``/usr/local/share``) which can be removed. The ``acol`` command-line
utility is installed to ``/usr/local/bin`` . Finally, all of the
libraries are named with the prefix ``libo2scl`` and are created by
default in ``/usr/local/lib``.

