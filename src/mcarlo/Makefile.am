# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

TEST_VAR = mcarlo_plain.scr mcarlo_miser.scr mcarlo_vegas.scr \
	expval.scr rng_gsl.scr mcmc_para.scr rng.scr
#mcmc_para_new.scr

MCARLO_SRCS = expval.cpp rng_gsl.cpp

HEADER_VARS = mcarlo_miser.h mcarlo_plain.h mcarlo_vegas.h mcarlo.h \
	expval.h rng_gsl.h mcmc_para.h rng.h mcmc_para_new.h \
	emulator.h 

# fit_bayes.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -DO2SCL_OPENMP -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

if O2SCL_EIGEN
O2SCL_EIGEN_MVAR = -DO2SCL_EIGEN
else
O2SCL_EIGEN_MVAR =
endif

if O2SCL_ARMA
O2SCL_ARMA_MVAR = -DO2SCL_ARMA
else
O2SCL_ARMA_MVAR =
endif

if O2SCL_LINUX
O2SCL_LINUX_MVAR = -DO2SCL_LINUX
else 
O2SCL_LINUX_MVAR = 
endif

if O2SCL_OSX
O2SCL_OSX_MVAR = -DO2SCL_OSX
else 
O2SCL_OSX_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VARS)

pkginclude_HEADERS = $(HEADER_VARS)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

if O2SCL_OPENMP
check_PROGRAMS = mcarlo_miser_ts mcarlo_vegas_ts mcarlo_plain_ts expval_ts \
	rng_gsl_ts mcmc_para_ts mcmc_para_new_ts rng_ts
else
check_PROGRAMS = mcarlo_miser_ts mcarlo_vegas_ts mcarlo_plain_ts expval_ts \
	rng_gsl_ts mcmc_para_new_ts rng_ts
endif

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

mcarlo_miser_ts_LDADD = $(ADDL_TEST_LIBS)
rng_ts_LDADD = $(ADDL_TEST_LIBS)
rng_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
mcarlo_plain_ts_LDADD = $(ADDL_TEST_LIBS)
mcarlo_vegas_ts_LDADD = $(ADDL_TEST_LIBS)
expval_ts_LDADD = $(ADDL_TEST_LIBS)
mcmc_para_ts_LDADD = $(ADDL_TEST_LIBS)
mcmc_para_new_ts_LDADD = $(ADDL_TEST_LIBS)

mcarlo_miser_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
rng_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
rng_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mcarlo_plain_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mcarlo_vegas_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
expval_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mcmc_para_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mcmc_para_new_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

mcarlo_miser.scr: mcarlo_miser_ts$(EXEEXT)
	./mcarlo_miser_ts$(EXEEXT) > mcarlo_miser.scr
rng.scr: rng_ts$(EXEEXT)
	./rng_ts$(EXEEXT) > rng.scr
rng_gsl.scr: rng_gsl_ts$(EXEEXT)
	./rng_gsl_ts$(EXEEXT) > rng_gsl.scr
mcarlo_plain.scr: mcarlo_plain_ts$(EXEEXT)
	./mcarlo_plain_ts$(EXEEXT) > mcarlo_plain.scr
mcarlo_vegas.scr: mcarlo_vegas_ts$(EXEEXT)
	./mcarlo_vegas_ts$(EXEEXT) > mcarlo_vegas.scr
expval.scr: expval_ts$(EXEEXT)
	./expval_ts$(EXEEXT) > expval.scr
mcmc_para.scr: mcmc_para_ts$(EXEEXT)
	./mcmc_para_ts$(EXEEXT) -exit > mcmc_para.scr
mcmc_para_new.scr: mcmc_para_new_ts$(EXEEXT)
	./mcmc_para_new_ts$(EXEEXT) -exit > mcmc_para_new.scr

mcarlo_miser_ts_SOURCES = mcarlo_miser_ts.cpp
rng_ts_SOURCES = rng_ts.cpp
rng_gsl_ts_SOURCES = rng_gsl_ts.cpp
mcarlo_plain_ts_SOURCES = mcarlo_plain_ts.cpp
mcarlo_vegas_ts_SOURCES = mcarlo_vegas_ts.cpp
expval_ts_SOURCES = expval_ts.cpp
mcmc_para_ts_SOURCES = mcmc_para_ts.cpp
mcmc_para_new_ts_SOURCES = mcmc_para_new_ts.cpp

# ------------------------------------------------------------
# Library o2scl_mcarlo
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_mcarlo.la

libo2scl_mcarlo_la_SOURCES = $(MCARLO_SRCS)

am__libo2scl_mcarlo_la_SOURCES_DIST = $(MCARLO_SRCS)

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
	echo "src/mcarlo" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

