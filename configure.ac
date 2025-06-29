# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([o2scl],[0.931],[awsteiner0@protonmail.com])
AC_CONFIG_SRCDIR(config.h.in)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# Checks for compiler and other required programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# 7/20/22, fixing this since o2scl uses c++14 extensions
#
AX_CXX_COMPILE_STDCXX(14,[],[optional])

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

# ------------------------------------------------------------------
# Take care of library version numbers
# ------------------------------------------------------------------

# For the time being, current and age are always
# "version number minus 910"
OLIB_CURRENT=21
OLIB_REVISION=0
OLIB_AGE=21

OLIB_LT_VERSION="${OLIB_CURRENT}:${OLIB_REVISION}:${OLIB_AGE}"
AC_SUBST(OLIB_LT_VERSION)

# ------------------------------------------------------------------
# Check for Boost headers
# ------------------------------------------------------------------

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

# ------------------------------------------------------------------
# Check for math library
# ------------------------------------------------------------------

AC_CHECK_LIB([m],[cos])

# ------------------------------------------------------------------
# Check for BLAS
# ------------------------------------------------------------------

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

# ------------------------------------------------------------------
# Check for GSL library (must check for BLAS first)
# ------------------------------------------------------------------

AC_CHECK_LIB([gsl], [gsl_vector_get], [], 
	[echo ""
        echo "GSL not found."
        echo "O2scl requires the GSL library."
	echo ""
	echo -n " The present value of LDFLAGS is: " 
	echo $LDFLAGS
        exit -1
	],[])

# ------------------------------------------------------------------
# Check for libquadmath
# ------------------------------------------------------------------

dnl AC_CHECK_LIB([quadmath], [isnanq], [], 
dnl 	[echo ""
dnl         echo "Quadmath not found."
dnl         echo "O2scl typically requires libquadmath."
dnl 	echo ""
dnl 	echo -n " The present value of LDFLAGS is: " 
dnl 	echo $LDFLAGS
dnl 	],[])

# ------------------------------------------------------------------
# Check for Readline library
# ------------------------------------------------------------------

AC_ARG_ENABLE(readline,
[  --enable-readline       Include readline support @<:@default=yes@:>@],
	[case "${enableval}" in
	yes) readline=true ;;
	no)  readline=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-readline) ;;
	esac],[readline=true])

AM_CONDITIONAL(O2SCL_READLINE, test x$readline = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the readline library
# ------------------------------------------------------------------

AS_IF([test "x$enable_readline" != "xno"], [
AC_CHECK_LIB([readline], [add_history], [],
	[echo ""
        echo "Readline not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])

])

# ------------------------------------------------------------------
# Check for module support
# ------------------------------------------------------------------

dnl AC_ARG_ENABLE(modules,
dnl [  --enable-modules        Include module support @<:@default=no@:>@],
dnl 	[case "${enableval}" in
dnl 	yes) modules=true ;;
dnl 	no)  modules=false ;;
dnl 	*) AC_MSG_ERROR(bad value ${enableval} for --enable-modules) ;;
dnl 	esac],[modules=false])

dnl AM_CONDITIONAL(O2SCL_MODULES, test x$modules = xtrue)

# ------------------------------------------------------------------
# Check for cuda support
# ------------------------------------------------------------------

AC_ARG_ENABLE(cuda,
[  --enable-cuda        Include cuda support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) cuda=true ;;
	no)  cuda=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
	esac],[cuda=false])

AM_CONDITIONAL(O2SCL_CUDA, test x$cuda = xtrue)

# ------------------------------------------------------------------
# Check for Mpfr library
# ------------------------------------------------------------------

AC_ARG_ENABLE(mpfr,
[  --enable-mpfr           Include mpfr support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) mpfr=true ;;
	no)  mpfr=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-mpfr) ;;
	esac],[mpfr=false])

AM_CONDITIONAL(O2SCL_MPFR, test x$mpfr = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the mpfr library
# ------------------------------------------------------------------

AS_IF([test "x$enable_mpfr" == "xyes"], [
AC_CHECK_LIB([mpfr], [mpfr_add], [],
	[echo ""
        echo "Mpfr not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])
])

# ------------------------------------------------------------------
# Check for Pugixml library
# ------------------------------------------------------------------

AC_ARG_ENABLE(pugixml,
[  --enable-pugixml        Pugixml support (dev only) @<:@default=no@:>@],
	[case "${enableval}" in
	yes) pugixml=true ;;
	no)  pugixml=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-pugixml) ;;
	esac],[pugixml=false])

AM_CONDITIONAL(O2SCL_PUGIXML, test x$pugixml = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the pugixml library
# ------------------------------------------------------------------

AS_IF([test "x$enable_pugixml" == "xyes"], [
AC_CHECK_LIB([pugixml], [atoi], [],
	[echo ""
        echo "Pugixml not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])
])

# ------------------------------------------------------------------
# Check for Cubature library
# ------------------------------------------------------------------

dnl AC_ARG_ENABLE(cubature,
dnl [  --enable-cubature       Include cubature support @<:@default=no@:>@],
dnl 	[case "${enableval}" in
dnl 	yes) cubature=true ;;
dnl 	no)  cubature=false ;;
dnl 	*) AC_MSG_ERROR(bad value ${enableval} for --enable-cubature) ;;
dnl 	esac],[cubature=false])

dnl AM_CONDITIONAL(O2SCL_CUBATURE, test x$cubature = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the cubature library
# ------------------------------------------------------------------

dnl AS_IF([test "x$enable_cubature" == "xyes"], [
dnl AC_CHECK_LIB([cubature], [pcubature], [],
dnl 	[echo ""
dnl         echo "Cubature not found. "
dnl 	echo ""
dnl 	echo -n " The present value of LDFLAGS is: "
dnl 	echo $LDFLAGS
dnl 	],[])
dnl ])

# ------------------------------------------------------------------
# Enable Python extensions
# ------------------------------------------------------------------

AC_ARG_ENABLE(python,
[  --enable-python         Include python support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) python=true ;;
	no)  python=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
	esac],[python=false])

AM_CONDITIONAL(O2SCL_PYTHON, test x$python = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the python library
# ------------------------------------------------------------------

AS_IF([test "x$enable_python" == "xyes"], [
        AX_PYTHON
])

# ------------------------------------------------------------------
# If enabled, check for the polylogarithm library
# ------------------------------------------------------------------

dnl AS_IF([test "x$enable_polylogarithm" == "xyes"], [
dnl AC_CHECK_LIB([polylog_cpp], [log], [],
dnl 	[echo ""
dnl         echo "Polylogarithm not found. "
dnl 	echo ""
dnl 	echo -n " The present value of LDFLAGS is: "
dnl 	echo $LDFLAGS
dnl 	],[])
dnl ])

# ------------------------------------------------------------------
# Check for ncurses library
# ------------------------------------------------------------------

AC_ARG_ENABLE(ncurses,
[  --enable-ncurses        Include ncurses support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) ncurses=true ;;
	no)  ncurses=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-ncurses) ;;
	esac],[ncurses=false])

AM_CONDITIONAL(O2SCL_NCURSES, test x$ncurses = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the ncurses library
# ------------------------------------------------------------------

AS_IF([test "x$enable_ncurses" != "xno"], [
AC_CHECK_LIB([ncurses], [initscr], [],
	[echo ""
        echo "Ncurses not found. "
	echo ""
	echo -n " The present value of LDFLAGS is: "
	echo $LDFLAGS
	],[])

])

# ------------------------------------------------------------------
# If enabled, check for the HDF5 library
# ------------------------------------------------------------------

AC_CHECK_LIB([hdf5], [H5Fopen], [], [
  AC_CHECK_LIB([hdf5_serial], [H5Fopen], [], [
    echo ""
    echo "HDF5 not found. Please provide the location for "
    echo "the hdf5 libraries."
    echo ""
    echo -n " The present value of LDFLAGS is: " 
    echo $LDFLAGS
  ],[])
],[])

# ------------------------------------------------------------------
# ... and for the HDF5 high-level library
# ------------------------------------------------------------------

AC_CHECK_LIB([hdf5_hl], [H5TBread_table], [], [
  AC_CHECK_LIB([hdf5_serial_hl], [H5TBread_table], [], [
    echo ""
    echo "High-level HDF5 interface not found. Please provide the "
    echo "location for the hdf5_hl libraries."
    echo ""
    echo -n " The present value of LDFLAGS is: " 
    echo $LDFLAGS
  ],[])
],[])  

# ------------------------------------------------------------------
# Enable or disable FFTW support (included here
# for later use)
# ------------------------------------------------------------------

AC_ARG_ENABLE(fftw,
	[  --enable-fftw           Include FFTW support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) fftw=true ;;
	no)  fftw=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-fftw) ;;
	esac],[fftw=false])
	
AM_CONDITIONAL(O2SCL_FFTW, test x$fftw = xtrue)

# ------------------------------------------------------------------
# If enabled, check for the FFTW library
# ------------------------------------------------------------------

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

# ------------------------------------------------------------------
# Enable or disable python support
# ------------------------------------------------------------------

# AC_ARG_ENABLE(python,
# 	[  --enable-python         Include Python support @<:@default=no@:>@],
# 	[case "${enableval}" in
# 	yes) python=true ;;
# 	no)  python=false ;;
# 	*) AC_MSG_ERROR(bad value ${enableval} for --enable-python) ;;
# 	esac],[python=false])

# AM_CONDITIONAL(O2SCL_PYTHON, test x$python = xtrue)

# ------------------------------------------------------------------
# Enable or disable support for acol
# ------------------------------------------------------------------

AC_ARG_ENABLE(acol,
[  --enable-acol           Include 'acol' utility @<:@default=yes@:>@],
   [case "${enableval}" in
   yes) acol=true ;;
   no)  acol=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-acol) ;;
   esac],[acol=true])

AM_CONDITIONAL(O2SCL_ENABLE_ACOL, test x$acol = xtrue)

# ------------------------------------------------------------------
# Enable or disable support for multiprecision
# ------------------------------------------------------------------

AC_ARG_ENABLE(multip,
[  --enable-multip           Add multiprecision support @<:@default=no@:>@],
   [case "${enableval}" in
   yes) multip=true ;;
   no)  multip=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-multip) ;;
   esac],[multip=false])

AM_CONDITIONAL(O2SCL_MULTIP, test x$multip = xtrue)

# ------------------------------------------------------------------
# Enable or disable snap configurations
# ------------------------------------------------------------------

AC_ARG_ENABLE(snap,
[  --enable-snap           Enable snap features @<:@default=no@:>@],
   [case "${enableval}" in
   yes) snap=true ;;
   no)  snap=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-snap) ;;
   esac],[snap=false])

AM_CONDITIONAL(O2SCL_SNAP, test x$snap = xtrue)

# ------------------------------------------------------------------
# Enable or disable armadillo support
# ------------------------------------------------------------------

AC_ARG_ENABLE(armadillo,
[  --enable-armadillo      Include Armadillo support @<:@default=no@:>@],
   [case "${enableval}" in
   yes) armadillo=true ;;
   no)  armadillo=false ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-armadillo) ;;
   esac],[armadillo=false])

AM_CONDITIONAL(O2SCL_ARMA, test x$armadillo = xtrue)

# ------------------------------------------------------------------
# If enabled, check for Armadillo library
# ------------------------------------------------------------------
 
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

# ------------------------------------------------------------------
# Enable or disable Eigen support
# ------------------------------------------------------------------

AC_ARG_ENABLE(eigen,
	[  --enable-eigen          Include Eigen support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) eigen=true ;;
	no)  eigen=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-eigen) ;;
	esac],[eigen=false])

AM_CONDITIONAL(O2SCL_EIGEN, test x$eigen = xtrue)

# ------------------------------------------------------------------
# Enable or disable openmp support
# ------------------------------------------------------------------

AC_ARG_ENABLE(openmp,
	[  --enable-openmp         Include OpenMP support @<:@default=no@:>@],
	[case "${enableval}" in
	yes) openmp=true ;;
	no)  openmp=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-openmp) ;;
	esac],[openmp=false])

AM_CONDITIONAL(O2SCL_OPENMP, test x$openmp = xtrue)

# ------------------------------------------------------------------
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

# ------------------------------------------------------------------

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/anneal/Makefile
                 src/base/Makefile
                 src/base/arma/Makefile
                 src/base/both/Makefile
                 src/base/cuda/Makefile
                 src/base/eigen/Makefile
                 src/base/fftw/Makefile
                 src/base/multip/Makefile
                 src/base/mpfr/Makefile
                 src/base/neither/Makefile
                 src/base/openmp/Makefile
                 src/base/python/Makefile
                 src/base/readline/Makefile
                 src/base/pugixml/Makefile
                 src/deriv/Makefile
                 src/eos/Makefile
                 src/fit/Makefile
                 src/hdf/Makefile
                 src/inte/Makefile
                 src/internal/Makefile
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
                 data/o2scl/python/Makefile
		 examples/Makefile
		 examples/plot/Makefile
                 doc/Makefile
                 doc/o2scl/Makefile
                 doc/o2scl/sphinx/Makefile])

AC_OUTPUT

