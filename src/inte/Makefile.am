# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

INTE_SRCS = inte_kronrod_gsl.cpp 

HEADER_VAR = inte.h inte_multi.h inte_gen.h inte_qng_gsl.h inte_qag_gsl.h \
	inte_cauchy_cern.h inte_adapt_cern.h \
	inte_gsl.h inte_gauss56_cern.h inte_gauss_cern.h \
	inte_kronrod_gsl.h inte_qawc_gsl.h inte_kronrod_boost.h \
	inte_qags_gsl.h inte_qagi_gsl.h inte_qagil_gsl.h inte_qagiu_gsl.h \
	inte_qawo_gsl.h inte_qawf_gsl.h inte_qaws_gsl.h inte_singular_gsl.h \
	inte_double_exp_boost.h

TEST_VAR = inte_adapt_cern.scr inte_cauchy_cern.scr inte_gauss_cern.scr \
	inte_gauss56_cern.scr inte_qawc_gsl.scr \
	inte_qng_gsl.scr inte_qag_gsl.scr inte_qags_gsl.scr \
	inte_qagi_gsl.scr inte_qagil_gsl.scr inte_qagiu_gsl.scr \
	inte_qawo_gsl.scr inte_qawf_gsl.scr inte_qaws_gsl.scr \
	inte_kronrod_boost.scr inte_double_exp_boost.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_OSX
O2SCL_OSX_MVAR = -DO2SCL_OSX
else 
O2SCL_OSX_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) -DO2SCL_COND_FLAG $(O2SCL_OSX_MVAR) \
	$(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) -DO2SCL_COND_FLAG $(O2SCL_OSX_MVAR) \
	$(O2SCL_PYTHON_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = inte_adapt_cern_ts inte_cauchy_cern_ts inte_gauss_cern_ts \
	inte_qng_gsl_ts inte_qag_gsl_ts inte_qags_gsl_ts \
	inte_qagi_gsl_ts inte_qagil_gsl_ts inte_qagiu_gsl_ts \
	inte_qawo_gsl_ts inte_qawf_gsl_ts inte_gauss56_cern_ts \
	inte_qawc_gsl_ts inte_qaws_gsl_ts inte_kronrod_boost_ts \
	inte_double_exp_boost_ts

check_SCRIPTS = o2scl-test

if O2SCL_PYTHON

ADDL_TEST_LIBS = ../libo2scl.la -l$(PYTHON_LIB)
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR) -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
endif

else

ADDL_TEST_LIBS = ../libo2scl.la
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = 
endif

endif

inte_adapt_cern_ts_LDADD = $(ADDL_TEST_LIBS)
inte_cauchy_cern_ts_LDADD = $(ADDL_TEST_LIBS)
inte_gauss_cern_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qng_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qag_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qags_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qagi_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qagil_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qagiu_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qawc_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qawo_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qawf_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_qaws_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
inte_gauss56_cern_ts_LDADD = $(ADDL_TEST_LIBS)
inte_kronrod_boost_ts_LDADD = $(ADDL_TEST_LIBS)
inte_double_exp_boost_ts_LDADD = $(ADDL_TEST_LIBS)

inte_adapt_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_cauchy_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_gauss_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qng_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qag_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qags_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qagi_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qagil_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qagiu_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qawc_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qawo_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qawf_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_qaws_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_gauss56_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_kronrod_boost_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
inte_double_exp_boost_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

inte_adapt_cern.scr: inte_adapt_cern_ts$(EXEEXT) 
	./inte_adapt_cern_ts$(EXEEXT) > inte_adapt_cern.scr
inte_cauchy_cern.scr: inte_cauchy_cern_ts$(EXEEXT) 
	./inte_cauchy_cern_ts$(EXEEXT) > inte_cauchy_cern.scr
inte_gauss_cern.scr: inte_gauss_cern_ts$(EXEEXT) 
	./inte_gauss_cern_ts$(EXEEXT) > inte_gauss_cern.scr
inte_qng_gsl.scr: inte_qng_gsl_ts$(EXEEXT) 
	./inte_qng_gsl_ts$(EXEEXT) > inte_qng_gsl.scr
inte_qag_gsl.scr: inte_qag_gsl_ts$(EXEEXT) 
	./inte_qag_gsl_ts$(EXEEXT) > inte_qag_gsl.scr
inte_qags_gsl.scr: inte_qags_gsl_ts$(EXEEXT) 
	./inte_qags_gsl_ts$(EXEEXT) > inte_qags_gsl.scr
inte_qagi_gsl.scr: inte_qagi_gsl_ts$(EXEEXT) 
	./inte_qagi_gsl_ts$(EXEEXT) > inte_qagi_gsl.scr
inte_qagil_gsl.scr: inte_qagil_gsl_ts$(EXEEXT) 
	./inte_qagil_gsl_ts$(EXEEXT) > inte_qagil_gsl.scr
inte_qagiu_gsl.scr: inte_qagiu_gsl_ts$(EXEEXT) 
	./inte_qagiu_gsl_ts$(EXEEXT) > inte_qagiu_gsl.scr
inte_qawc_gsl.scr: inte_qawc_gsl_ts$(EXEEXT) 
	./inte_qawc_gsl_ts$(EXEEXT) > inte_qawc_gsl.scr
inte_qawo_gsl.scr: inte_qawo_gsl_ts$(EXEEXT) 
	./inte_qawo_gsl_ts$(EXEEXT) > inte_qawo_gsl.scr
inte_qawf_gsl.scr: inte_qawf_gsl_ts$(EXEEXT) 
	./inte_qawf_gsl_ts$(EXEEXT) > inte_qawf_gsl.scr
inte_qaws_gsl.scr: inte_qaws_gsl_ts$(EXEEXT) 
	./inte_qaws_gsl_ts$(EXEEXT) > inte_qaws_gsl.scr
inte_gauss56_cern.scr: inte_gauss56_cern_ts$(EXEEXT) 
	./inte_gauss56_cern_ts$(EXEEXT) > inte_gauss56_cern.scr
inte_kronrod_boost.scr: inte_kronrod_boost_ts$(EXEEXT) 
	./inte_kronrod_boost_ts$(EXEEXT) | tee inte_kronrod_boost.scr
inte_double_exp_boost.scr: inte_double_exp_boost_ts$(EXEEXT) 
	./inte_double_exp_boost_ts$(EXEEXT) > inte_double_exp_boost.scr

inte_adapt_cern_ts_SOURCES = inte_adapt_cern_ts.cpp
inte_cauchy_cern_ts_SOURCES = inte_cauchy_cern_ts.cpp
inte_gauss_cern_ts_SOURCES = inte_gauss_cern_ts.cpp
inte_qng_gsl_ts_SOURCES = inte_qng_gsl_ts.cpp
inte_qag_gsl_ts_SOURCES = inte_qag_gsl_ts.cpp
inte_qags_gsl_ts_SOURCES = inte_qags_gsl_ts.cpp
inte_qagi_gsl_ts_SOURCES = inte_qagi_gsl_ts.cpp
inte_qagil_gsl_ts_SOURCES = inte_qagil_gsl_ts.cpp
inte_qagiu_gsl_ts_SOURCES = inte_qagiu_gsl_ts.cpp
inte_qawc_gsl_ts_SOURCES = inte_qawc_gsl_ts.cpp
inte_qawo_gsl_ts_SOURCES = inte_qawo_gsl_ts.cpp
inte_qawf_gsl_ts_SOURCES = inte_qawf_gsl_ts.cpp
inte_qaws_gsl_ts_SOURCES = inte_qaws_gsl_ts.cpp
inte_gauss56_cern_ts_SOURCES = inte_gauss56_cern_ts.cpp
inte_kronrod_boost_ts_SOURCES = inte_kronrod_boost_ts.cpp
inte_double_exp_boost_ts_SOURCES = inte_double_exp_boost_ts.cpp

# ------------------------------------------------------------
# Library o2scl_inte
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_inte.la

libo2scl_inte_la_SOURCES = $(INTE_SRCS) 

am__libo2scl_inte_la_SOURCES_DIST = $(INTE_SRCS) 

# ------------------------------------------------------------
# Misc
# ------------------------------------------------------------

EXTRA_DIST = *_ts.cpp

emacs-clean: 
	-rm *~

qbkup:
	cp *.cpp *.h Makefile.am bkup

# ------------------------------------------------------------
# Testing
# ------------------------------------------------------------

o2scl-test: $(TEST_VAR)
	echo $(TEST_VAR) > testlist
	../internal/tsumm$(EXEEXT)

o2scl-test-rec: $(TEST_VAR)
	echo "src/inte" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

