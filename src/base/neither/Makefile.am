
if O2SCL_EIGEN
HEADER_VAR = 
else
if O2SCL_ARMA
HEADER_VAR = 
else
HEADER_VAR = vector_special.h
endif
endif

if O2SCL_MPFR
MPFR_VAR = 
else
MPFR_VAR = set_mpfr.h
endif

if O2SCL_FFTW
FFTW_VAR = 
else
FFTW_VAR = set_fftw.h
endif

if O2SCL_PYTHON
PYTHON_VAR = 
else
PYTHON_VAR = set_python.h
endif

if O2SCL_OPENMP
OPENMP_VAR = 
else
OPENMP_VAR = set_openmp.h
endif

if O2SCL_READLINE
READLINE_VAR = 
else
READLINE_VAR = set_readline.h
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = vector_special.h set_fftw.h set_mpfr.h \
	set_openmp.h set_python.h set_readline.h

pkginclude_HEADERS = $(HEADER_VAR)

