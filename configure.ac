# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([o2scl],[0.918],[awsteiner@users.sourceforge.net])
AC_CONFIG_SRCDIR(config.h.in)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for compiler and other required programs
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Newer replacement for AC_PROG_LIBTOOL
LT_INIT

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_SELECT_ARGTYPES
AC_FUNC_STRTOD

# Check for popen
AC_CHECK_FUNC([popen],[AC_DEFINE([HAVE_POPEN],[1],
[Define if popen exists])])

# ----------------------------------------------
# Take care of library version numbers
# ----------------------------------------------

OLIB_CURRENT=8
OLIB_REVISION=0
OLIB_AGE=8

OLIB_LT_VERSION="${OLIB_CURRENT}:${OLIB_REVISION}:${OLIB_AGE}"
AC_SUBST(OLIB_LT_VERSION)

# ----------------------------------------------
# Enable or disable particle and nuclei library
# ----------------------------------------------

AC_ARG_ENABLE(partlib,
[  --enable-partlib        Include particle library @<:@default=yes@:>@],
   [case "${enableval}" in
   yes) partlib=true ;;
   no)  partlib=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-partlib) ;;
   esac],[partlib=true])

AM_CONDITIONAL(O2SCL_PARTLIB, test x$partlib = xtrue)

# ----------------------------------------------
# Enable or disable EOS library
# ----------------------------------------------

AC_ARG_ENABLE(eoslib,
[  --enable-eoslib         Include equation of state library @<:@default=yes@:>@],
   [case "${enableval}" in
   yes) eoslib=true ;;
   no)  eoslib=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-eoslib) ;;
   esac],[eoslib=true])

AM_CONDITIONAL(O2SCL_EOSLIB, test x$eoslib = xtrue)

# ----------------------------------------------
# Check for Boost headers
# ----------------------------------------------

# AWS: This doesn't seem to correctly find the Boost libraries yet

#AC_CHECK_HEADER([boost/numeric/ublas/vector.hpp],[],
#	[echo ""
#        echo "Boost not found. O2scl requires Boost."
#	echo ""
#	echo -n " The present value of CPPFLAGS is: " 
#	echo $CPPFLAGS
#	echo ""
#	exit -1
#	],[])

# ----------------------------------------------
# Check for BLAS
# ----------------------------------------------

AC_CHECK_LIB([cblas],[cblas_dgemm],[],[
  echo "Library libcblas not found. Looking for GSL cblas." 
  echo -n " The present value of LDFLAGS is: " 
  echo $LDFLAGS
  AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[
    echo "Library libgslcblas not found. O2scl requires a cblas library." 
    echo "You may be required to add a cblas library to the LIBS "
    echo "environment variable. "
    echo ""
    echo -n " The present value of LDFLAGS is: " 
    echo $LDFLAGS
    echo ""
  ],[-lm])
],[-lm])

# ----------------------------------------------
# Check for GSL library (must check for BLAS first)
# ----------------------------------------------

AC_CHECK_LIB([gsl], [gsl_vector_get], [], 
	[echo ""
        echo "GSL not found."
        echo "O2scl requires the GSL library."
	echo ""
	echo -n " The present value of LDFLAGS is: " 
	echo $LDFLAGS
        exit -1
	],[-lm])

# ----------------------------------------------
# Check for Readline libraries
# ----------------------------------------------

AC_ARG_ENABLE(readline,
[  --enable-readline       Include readline support @<:@default=yes@:>@],
	[case "${enableval}" in
	yes) readline=true ;;
	no)  readline=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-readline) ;;
	esac],[readline=true])

AM_CONDITIONAL(O2SCL_READLINE, test x$readline = xtrue)

# ----------------------------------------------
# Enable or disable HDF5 support
# ----------------------------------------------

AC_ARG_ENABLE(hdf,
	[  --enable-hdf            Include HDF5 support @<:@default=yes@:>@],
	[case "${enableval}" in
	yes) hdf=true ;;
	no)  hdf=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-hdf) ;;
	esac],[hdf=true])

AM_CONDITIONAL(O2SCL_HDF, test x$hdf = xtrue)

# ----------------------------------------------
# If enabled, check for the HDF5 library
# ----------------------------------------------

AS_IF([test "x$enable_hdf" != "xno"], [
	    
AC_CHECK_LIB([hdf5], [H5Fopen], [],
	[echo ""
        echo "HDF5 not found. Either provide the location for "
	echo "the hdf5 libraries or ensure that HDF5 support is "
	echo "disabled with the --disable-hdf option to configure. "
	echo ""
	echo -n " The present value of LDFLAGS is: " 
	echo $LDFLAGS
 	],[-lm])

# ----------------------------------------------
# ... and for the HDF5 high-level library
# ----------------------------------------------

AC_CHECK_LIB([hdf5_hl], [H5TBread_table], [],
	[echo ""
        echo "High-level HDF5 interface not found. Either provide the "
	echo "location for the hdf5_hl libraries or ensure that HDF5 "
	echo "support is disabled with the --disable-hdf option to "
	echo "configure. "
	echo ""
	echo -n " The present value of LDFLAGS is: " 
	echo $LDFLAGS
 	],[-lm])

])

# ----------------------------------------------
# Enable or disable FFTW support
# ----------------------------------------------

AC_ARG_ENABLE(fftw,
	[  --enable-fftw           Include FFTW support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) fftw=true ;;
	no)  fftw=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-fftw) ;;
	esac],[fftw=false])
	
AM_CONDITIONAL(O2SCL_FFTW, test x$fftw = xtrue)

# ----------------------------------------------
# If enabled, check for the FFTW library
# ----------------------------------------------

dnl AS_IF([test "x$enable_fftw" != "xno"],[
dnl AC_CHECK_LIB([fftw3], [fftw_ifree], [],
dnl 	[echo ""
dnl         echo "FFTW not found. Either provide the location for "
dnl 	echo "the fftw libraries or ensure that FFTW support is "
dnl 	echo "disabled with the --disable-fftw option to configure. "
dnl 	echo ""
dnl 	echo -n " The present value of LDFLAGS is: " 
dnl 	echo $LDFLAGS
dnl  	],[-lm])
dnl ])

# ----------------------------------------------
# Enable or disable support for acol
# ----------------------------------------------

AC_ARG_ENABLE(acol,
[  --enable-acol           Include 'acol' utility @<:@default=yes@:>@],
   [case "${enableval}" in
   yes) acol=true ;;
   no)  acol=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-acol) ;;
   esac],[acol=true])

AM_CONDITIONAL(O2SCL_ENABLE_ACOL, test x$acol = xtrue)

# ----------------------------------------------
# Enable or disable armadillo support
# ----------------------------------------------

AC_ARG_ENABLE(armadillo,
[  --enable-armadillo      Include Armadillo support @<:@default=no@:>@],
   [case "${enableval}" in
   yes) armadillo=true ;;
   no)  armadillo=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-armadillo) ;;
   esac],[armadillo=false])

AM_CONDITIONAL(O2SCL_ARMA, test x$armadillo = xtrue)

# ----------------------------------------------
# If enabled, check for Armadillo library
# ----------------------------------------------
 
# AWS: This doesn't correctly find the Armadillo libraries yet
 
#AS_IF([test "x$enable_armadillo" != "xno"], [
#	    
#AC_CHECK_LIB([armadillo], [arma::is_finite()],[],
#	[echo ""
#        echo "Armadillo library not found."
#	echo " The present value of LDFLAGS is: " 
#	echo $LDFLAGS
#	echo ""
#	],[])
#
#])

# ----------------------------------------------
# Enable or disable Eigen support
# ----------------------------------------------

AC_ARG_ENABLE(eigen,
	[  --enable-eigen          Include Eigen support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) eigen=true ;;
	no)  eigen=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-eigen) ;;
	esac],[eigen=false])

AM_CONDITIONAL(O2SCL_EIGEN, test x$eigen = xtrue)

# ----------------------------------------------
# Enable or disable openmp support
# ----------------------------------------------

AC_ARG_ENABLE(openmp,
	[  --enable-openmp         Include OpenMP support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) openmp=true ;;
	no)  openmp=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-openmp) ;;
	esac],[openmp=false])

AM_CONDITIONAL(O2SCL_OPENMP, test x$openmp = xtrue)

# ----------------------------------------------

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/anneal/Makefile
                 src/base/Makefile
                 src/base/arma/Makefile
                 src/base/both/Makefile
                 src/base/eigen/Makefile
                 src/base/neither/Makefile
                 src/boost/Makefile
		 src/cubature/Makefile
                 src/deriv/Makefile
                 src/eos/Makefile
                 src/fit/Makefile
		 src/fparser/Makefile
                 src/hdf/Makefile
                 src/inte/Makefile
                 src/linalg/Makefile
                 src/linalg/neither/Makefile
                 src/linalg/both/Makefile
                 src/linalg/arma/Makefile
                 src/linalg/eigen/Makefile
		 src/mcarlo/Makefile
                 src/min/Makefile
                 src/nuclei/Makefile
                 src/ode/Makefile
                 src/other/Makefile
                 src/internal/Makefile
                 src/part/Makefile
                 src/root/Makefile	
                 src/root/arma/Makefile	
                 src/root/both/Makefile	
                 src/root/eigen/Makefile	
                 src/root/neither/Makefile	
		 include/Makefile
		 include/o2scl/Makefile	
                 data/Makefile
                 data/o2scl/Makefile
                 data/o2scl/nucmass/Makefile
                 data/o2scl/rmfdata/Makefile
                 data/o2scl/skdata/Makefile
		 examples/Makefile
		 examples/plot/Makefile
                 doc/Makefile
                 doc/o2scl/Makefile
                 doc/o2scl/html/Makefile
                 doc/o2scl/html/search/Makefile
		 doc/o2scl/extras/Makefile
                 doc/o2scl/eos/Makefile
                 doc/o2scl/eos/html/Makefile
                 doc/o2scl/eos/html/search/Makefile
                 doc/o2scl/eos/extras/Makefile
                 doc/o2scl/part/Makefile
                 doc/o2scl/part/html/Makefile
                 doc/o2scl/part/html/search/Makefile
		 doc/o2scl/part/extras/Makefile])

AC_OUTPUT

