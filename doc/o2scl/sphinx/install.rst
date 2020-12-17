Installation
============

:ref:`O2scl <o2scl>`

Contents
--------

- :ref:`Compiling O2scl on Ubuntu with Snap <compile_snap>`
- :ref:`Compiling O2scl on Mac OSX with Homebrew <compile_homebrew>`
- :ref:`Compiling O2scl from a release distribution <compile_dist>`
- :ref:`Compiling O2scl from a release on Linux <compile_release>`
- :ref:`Compiling O2scl from the source code <compile_source>`
- :ref:`Compiling O2scl on Docker <compile_docker>`
- :ref:`Optional linear algebra libraries`
- :ref:`Range-checking`
- :ref:`Optional physics libraries`
- :ref:`More configure flags`
- :ref:`Generation of documentation`
- :ref:`Uninstallation`

.. _compile_snap:

Compiling O\ :sub:`2`\ scl on Ubuntu with Snap
----------------------------------------------

The easiest way to install on Ubuntu is with snap (see
https://snapcraft.io/o2scl). Use::

  sudo snap install (--edge or --beta) --devmode o2scl

The snap installation includes HDF5 support, the O\ :sub:`2`\ scl_part
and O\ :sub:`2`\ scl_eos sub-libraries, readline support, and uses the
GSL CBLAS.

Using the command-line utility ``acol`` may require you to set the
environment variable ``LD_LIBRARY_PATH``. For example, on machines
where I use snap to install in my ``.bashrc``, I use::

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/snap/o2scl/current/usr/lib/x86_64-linux-gnu:/snap/o2scl/current/lib/x86_64-linux-gnu

.. _compile_homebrew:
  
Compiling O\ :sub:`2`\ scl on Mac OSX with Homebrew
---------------------------------------------------

The easiest way to install on Mac OSX is with homebrew. Use::

  brew tap awsteiner/science
  brew install o2scl

to install O\ :sub:`2`\ scl. There are a few options for ``brew
install``. The option ``--with-check`` performs the build-time tests
and the option ``--with-examples`` double checks that the examples can
also be compiled and executed. The homebrew recipe for O\ :sub:`2`\
scl uses the Mac OS X compiler clang. Homebrew also supports the
installation of the current version directly from the repository using
the ``--HEAD`` option to ``brew install``. The homebrew installation
includes the O\ :sub:`2`\ scl_part and O\ :sub:`2`\ scl_eos
sub-libraries and readline support. The O\ :sub:`2`\ scl homebrew
recipes are stored at the
https://github.com/awsteiner/homebrew-science repository.

(By default, a homebrew installation of O\ :sub:`2`\ scl uses the OSX LLVM
compiler. However, a homebrew installation of O\ :sub:`2`\ scl will also
install ``gcc`` because O\ :sub:`2`\ scl requires ``hdf5``, and the homebrew
``hdf5`` package requires ``gcc``. The homebrew installation of 
O\ :sub:`2`\ scl is tested by Travis CI.)

.. _compile_dist:

Compiling O\ :sub:`2`\ scl from a release distribution
------------------------------------------------------

O\ :sub:`2`\ scl installation is generally similar to that for
GNU-style libraries. The file ``INSTALL`` has some details on this
procedure. Once the dependencies are installed you should be able to
run ``./configure`` and then type ``make`` and ``make install``. More
information on the ``configure`` command can also be obtained from
``./configure --help``. O\ :sub:`2`\ scl assumes some C++11 support,
so compilation may be more difficult on compilers released before
about 2018. The ``./configure`` script attempts to determine the
proper compiler flags for C++11 support, e.g. ``-std=gnu++11``. If
this fails, you may have to add the proper C++11 flag to the
``CXXFLAGS`` environment variable manually before the ``./configure``
script. The documentation is included in the O\ :sub:`2`\ scl release
distribution and automatically installed by ``make install``.

.. note::
   If you are trying to install O\ :sub:`2`\ scl with a version of
   HDF5 earlier than 1.12 you will need to compile with
   ``-DO2SCL_HDF5_PRE_1_12``.

O\ :sub:`2`\ scl requires the Boost (any relatively recent version)
and the GSL libraries (version 2.0 or later). If the
``configure`` script cannot find Boost or GSL, you may have to
specify their location for the associated header files in the
``CXXFLAGS`` variable and the associated libraries in the
``LDFLAGS`` environment variable. Running ``./configure
--help`` shows some information on this. For example, in a bash
shell, you could do something like::

  CXX="g++" CXXFLAGS="-I/dir/to/gsl/include" LDFLAGS="-L/dir/to/gsl/libs" ./configure --prefix=="/dir/to/destination_directory

Along with GSL, a CBLAS library is also required, and ``./configure``
will look for ``libcblas`` first, and if not found then it will look
for ``libgslcblas``. If neither is present, then you may have to
manually specify a CBLAS library using the ``LIBS`` and ``LDFLAGS``
environment variables.

Compiling with the readline, ncurses, and HDF5 libraries is optional,
but they are assumed to be present by default. To compile without
these libraries, you will need to use the arguments
``--disable-readline``, ``--disable-ncurses`` or ``--disable-hdf`` to
``./configure``, respectively. Note that HDF5 is currently required
for the physics sub-libraries, so ``--disable-hdf`` should be
accompanied by the ``--disable-eoslib`` and ``--disable-partlib``
flags.

After ``make install``, you may test the library with ``make check``
or ``make o2scl-test``. At the end, the phrase ``"All O2scl tests
passed"`` indicates that the testing was successful. You may also run
``make o2scl-test`` in the individual subdirectories of the src
directory to individually test the classes and functions in that part
of O\ :sub:`2`\ scl. The testing code in
``src/base/lib_settings_ts.cpp`` can be useful in finding out how O\
:sub:`2`\ scl was compiled. After ``make o2scl-test``, running
``src/base/lib_settings_ts`` will output several of the installation
settings. If HDF5 is enabled, ``acol -v`` also outputs the
installation settings.

O\ :sub:`2`\ scl uses Travis CI (see
https://travis-ci.org/awsteiner/o2scl ) to ensure that compilation and
testing works on standard Ubuntu and Mac OS X environments.

.. _compile_release:

Compiling O\ :sub:`2`\ scl from a release on Linux
--------------------------------------------------

For example, to install O\ :sub:`2`\ scl on Ubuntu, begin by
installing g++ and make (the ``g++`` and ``make`` packages),
GSL (the ``libgsl-dev`` package), Boost (the
``libboost-all-dev`` package), GNU readline (the ``libreadline-dev``
package), ncurses (the ``libncurses-dev`` packages), and HDF5 the
``libhdf5-dev`` package). You can then install O\ :sub:`2`\ scl from
one of the release distributions by using the standard GNU
``./configure`` script and then invoking ``make`` and ``make install``
(which sometimes requires ``sudo``). This installation method is
tested by the Travis CI script.
 
The HDF5 package for Ubuntu and many other Linux systems is
installed in ``hdf5/serial/hdf5.h`` instead of
``hdf5.h``, so O\ :sub:`2`\ scl presumes that Linux systems are arranged
that way. If HDF5 include statements should not have the
``hdf5/serial/`` prefix, then you can use
``-DO2SCL_HDF5_PLAIN_HEADER``, i.e.::

  CXXFLAGS="-DO2SCL_PLAIN_HDF5_HEADER" ./configure

to instruct O\ :sub:`2`\ scl to look for them there (for example, on bridges at
the PSC). On many systems, one can use a parallel HDF5 library
using ``-DO2SCL_HDF5_PLAIN_HEADER`` and a ``-I`` option
to select the proper location for the parallel HDF5 header files.
Finally, if your version of HDF5 is earlier than 1.12,
you will need to let O\ :sub:`2`\ scl know, using::

  CXXFLAGS="-DO2SCL_HDF5_PRE_1_12" ./configure

Other Linux distributions are similar. For example, in OpenSUSE, you
will need to use ``zypper`` to install ``gcc-c++, make, gsl-devel,
hdf5-devel, ncurses-devel, readline-devel``, and ``boost-devel``.

.. _compile_source:

Compiling O\ :sub:`2`\ scl from the source code
-----------------------------------------------

If you want to install from source (without generating the
documentation), then you must first install ``g++``, ``make``,
``automake``, ``autoconf``, and ``libtool`` packages. Then you can use
something along the lines of::

  git clone https://github.com/awsteiner/o2scl
  cd o2scl
  mkdir m4
  autoreconf -i
  ./configure

Then, you will either need to generate the documentation from doxygen
using ``make o2scl-doc`` or use ``make blank-doc`` to create blank
documentation. Then you can proceed using ``make`` and ``make
install`` (which may require ``sudo`` depending on your
configuration). For a full installation with parallelism, I
typically also install ``libopenmpi-dev`` and then use
``./configure --enable-openmp``

.. _compile_docker:

Compiling O\ :sub:`2`\ scl on Docker
------------------------------------

There are also some experimental dockerfiles which you can use to
install O\ :sub:`2`\ scl which can be found at
https://github.com/awsteiner/o2scl/tree/master/docker . For those on
MacOS, I recommend the guide at
https://medium.com/crowdbotics/a-complete-one-by-one-guide-to-install-docker-on-your-mac-os-using-homebrew-e818eb4cfc3
to installing docker.

Optional linear algebra libraries
---------------------------------

Most classes and functions which require linear algebra can be used
with the Eigen (http://eigen.tuxfamily.org) or Armadillo
(http://arma.sourceforge.net) vector and matrix objects. This can be
specified in the ``configure`` command with ``--enable-armadillo`` or
``--enable-eigen``. Note that the O\ :sub:`2`\ scl classes which use
Armadillo use matrix decompositions so Armadillo must be compiled with
LAPACK support, and you may need to specify the location of the LAPACK
libraries manually. If you are installing on Mac OS X with homebrew,
the options ``--with-eigen`` and ``with-armadillo`` can be used.

Range-checking
--------------

Some extra range-checking for vectors and matrices is turned on by
default. You can disable range-checking by defining
-DO2SCL_NO_RANGE_CHECK, e.g.::

  CXXFLAGS="-DO2SCL_NO_RANGE_CHECK" ./configure

Optional physics libraries
--------------------------

The separate libraries O\ :sub:`2`\ scl_eos and O\ :sub:`2`\ scl_part
are installed by default. To disable the installation of these
libraries and their associated documentation, run ``./configure`` with
the flags ``--disable-eoslib`` or ``--disable-partlib``. Note that O\
:sub:`2`\ scl_eos depends on O\ :sub:`2`\ scl_part so using
``--disable-partlib`` without ``--disable-eoslib`` will not work. Note
also that both O\ :sub:`2`\ scl_part and O\ :sub:`2`\ scl_eos require
HDF5 support.

More configure flags
--------------------

There are several warning flags that are useful when configuring
and compiling with O\ :sub:`2`\ scl. See the GSL documentation for an 
excellent discussion, and also see the generic installation
documentation in the file ``INSTALL`` in the O\ :sub:`2`\ scl top-level 
directory. For running ``configure``, for example, if you do
not have privileges to write to ``/usr/local``::

  CPPFLAGS="-O3 -I/home/asteiner/install/include" \
  LDFLAGS="-L/home/asteiner/install/lib" ./configure \
  --prefix=/home/asteiner/install

In this example, specifying ``-I/home/asteiner/install/include`` and
``-L/home/asteiner/install/lib`` above ensures that the GSL libraries
can be found. The ``--prefix=/home/asteiner/install`` argument to
``./configure`` ensures that O\ :sub:`2`\ scl is installed there as
well.

Generation of documentation
---------------------------

The O\ :sub:`2`\ scl documentation is generated with ``doxygen``,
``sphinx``, ``breathe``, and ``alabaster`` and packaged in with every
release file. In principle, the documentation can be regenerated by
the end-user, but this is not supported and requires several external
applications not included in the distribution.

The most recent release documentation is available at
https://neutronstars.utk.edu/code/o2scl/html/index.html and the
current development version documentation is available at
https://neutronstars.utk.edu/code/o2scl-dev/html/index.html . The
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

