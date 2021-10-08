# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([o2scl],[0.927a1],[awsteiner@mykolab.com])
AC_CONFIG_SRCDIR(config.h.in)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for compiler and other required programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# O2scl requires c++11 support, but sometimes it isn't properly
# detected when it's present (notably in travis), so we make it
# optional here and users will have to modify their flags if this
# fails.
AX_CXX_COMPILE_STDCXX_11(,[optional])

# Newer replacement for LT_INIT
LT_INIT

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_SELECT_ARGTYPES

# Check for popen
AC_CHECK_FUNC([popen],[AC_DEFINE([HAVE_POPEN],[1],
[Define if popen exists])])

# ----------------------------------------------
# Take care of library version numbers
# ----------------------------------------------

OLIB_CURRENT=17
OLIB_REVISION=0
OLIB_AGE=17

OLIB_LT_VERSION="${OLIB_CURRENT}:${OLIB_REVISION}:${OLIB_AGE}"
AC_SUBST(OLIB_LT_VERSION)

# ----------------------------------------------
# Enable or disable GSL 2.0 support
# ----------------------------------------------

AC_ARG_ENABLE(gsl2,
[  --enable-gsl2           Include GSL v2.0 support @<:@default=yes@:>@],
   [case "${enableval}" in
   yes) gsl2=true ;;
   no)  gsl2=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-gsl2) ;;
   esac],[gsl2=true])

AM_CONDITIONAL(O2SCL_GSL2, test x$gsl2 = xtrue)

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
# Check for math library
# ----------------------------------------------

AC_CHECK_LIB([m],[cos])

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
  ],[])
],[])

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
	],[])

# ----------------------------------------------
# Check for Readline library
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
# If enabled, check for the readline library
# ----------------------------------------------

AS_IF([test "x$enable_readline" != "xno"], [

AC_CHECK_LIB([readline], [add_history], [],
	[echo ""
        echo "Readline not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])

])

# ----------------------------------------------
# Check for ncurses library
# ----------------------------------------------

AC_ARG_ENABLE(ncurses,
[  --enable-ncurses        Include ncurses support @<:@default=yes@:>@],
	[case "${enableval}" in
	yes) ncurses=true ;;
	no)  ncurses=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-ncurses) ;;
	esac],[ncurses=true])

AM_CONDITIONAL(O2SCL_NCURSES, test x$ncurses = xtrue)

# ----------------------------------------------
# If enabled, check for the ncurses library
# ----------------------------------------------

AS_IF([test "x$enable_ncurses" != "xno"], [

AC_CHECK_LIB([ncurses], [initscr], [],
	[echo ""
        echo "Ncurses not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])

])

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
	    
AC_CHECK_LIB([hdf5], [H5Fopen], [], [
  AC_CHECK_LIB([hdf5_serial], [H5Fopen], [], [
    echo ""
    echo "HDF5 not found. Either provide the location for "
    echo "the hdf5 libraries or ensure that HDF5 support is "
    echo "disabled with the --disable-hdf option to configure. "
    echo ""
    echo -n " The present value of LDFLAGS is: " 
    echo $LDFLAGS
  ],[])
],[])

])

# ----------------------------------------------
# ... and for the HDF5 high-level library
# ----------------------------------------------

AC_CHECK_LIB([hdf5_hl], [H5TBread_table], [], [
  AC_CHECK_LIB([hdf5_serial_hl], [H5TBread_table], [], [
    echo ""
    echo "High-level HDF5 interface not found. Either provide the "
    echo "location for the hdf5_hl libraries or ensure that HDF5 "
    echo "support is disabled with the --disable-hdf option to "
    echo "configure. "
    echo ""
    echo -n " The present value of LDFLAGS is: " 
    echo $LDFLAGS
  ],[])
],[])  

# ----------------------------------------------
# Enable or disable FFTW support (included here
# for later use)
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

AS_IF([test "x$enable_fftw" == "xyes"],[
AC_CHECK_LIB([fftw3], [fftw_ifree], [],
	[echo ""
        echo "FFTW not found. Either provide the location for "
	echo "the fftw libraries or ensure that FFTW support is "
	echo "disabled with the --disable-fftw option to configure. "
	echo ""
	echo -n " The present value of LDFLAGS is: " 
	echo $LDFLAGS
 	],[])
])

# ----------------------------------------------
# Enable or disable python support
# ----------------------------------------------

dnl AC_ARG_ENABLE(python,
dnl 	[  --enable-python         Include Python support @<:@default=no@:>@],
dnl 	[case "${enableval}" in
dnl 	yes) python=true ;;
dnl 	no)  python=false ;;
dnl 	*) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
dnl 	esac],[python=false])

dnl AM_CONDITIONAL(O2SCL_PYTHON, test x$python = xtrue)

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
# Enable or disable snap configurations
# ----------------------------------------------

AC_ARG_ENABLE(snap,
[  --enable-snap           Enable snap features @<:@default=no@:>@],
   [case "${enableval}" in
   yes) snap=true ;;
   no)  snap=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-snap) ;;
   esac],[snap=false])

AM_CONDITIONAL(O2SCL_SNAP, test x$snap = xtrue)

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
 
AS_IF([test "x$enable_armadillo" == "xyes"], [
	    
AC_CHECK_LIB([armadillo], [wrapper_dgemm_],[],
	[echo ""
        echo "Armadillo library not found."
	echo " The present value of LDFLAGS is: " 
	echo $LDFLAGS
	echo ""
	],[])
])

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
# Define 'host_os' to determine target system
# from https://stackoverflow.com/questions/38898591/how-to-check-the-os-with-automake

AC_CANONICAL_HOST
build_linux=no
build_windows=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([O2SCL_LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([O2SCL_WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([O2SCL_OSX], [test "$build_mac" = "yes"])

# ----------------------------------------------

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/anneal/Makefile
                 src/base/Makefile
                 src/base/arma/Makefile
                 src/base/both/Makefile
                 src/base/eigen/Makefile
                 src/base/neither/Makefile
		 src/cubature/Makefile
                 src/deriv/Makefile
                 src/eos/Makefile
                 src/fit/Makefile
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
		 bin/Makefile
		 include/Makefile
		 include/o2scl/Makefile	
                 data/Makefile
                 data/o2scl/Makefile
                 data/o2scl/nucmass/Makefile
                 data/o2scl/nucmass/frib_mex/Makefile
                 data/o2scl/nucmass/ame20/Makefile
                 data/o2scl/rmfdata/Makefile
                 data/o2scl/skdata/Makefile
		 examples/Makefile
		 examples/plot/Makefile
                 doc/Makefile
                 doc/o2scl/Makefile
                 doc/o2scl/sphinx/Makefile
                 doc/o2scl/eos/Makefile
                 doc/o2scl/eos/sphinx/Makefile
                 doc/o2scl/part/Makefile
                 doc/o2scl/part/sphinx/Makefile])

AC_OUTPUT

