# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

TEST_VAR = mcarlo_plain.scr mcarlo_miser.scr mcarlo_vegas.scr \
	expval.scr rng_gsl.scr mcmc_para.scr

MCARLO_SRCS = expval.cpp rng_gsl.cpp

HEADER_VARS = mcarlo_miser.h mcarlo_plain.h mcarlo_vegas.h mcarlo.h \
	expval.h rng_gsl.h mcmc_para.h

# fit_bayes.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_PARTLIB
O2SCL_PART_MVAR = -DO2SCL_PART
else
O2SCL_PART_MVAR =
endif

if O2SCL_EOSLIB
O2SCL_EOS_MVAR = -DO2SCL_EOS
else
O2SCL_EOS_MVAR =
endif

if O2SCL_HDF
O2SCL_HDF_MVAR = -DO2SCL_HDF
else
O2SCL_HDF_MVAR =
endif

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

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG

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
	rng_gsl_ts mcmc_para_ts
else
check_PROGRAMS = mcarlo_miser_ts mcarlo_vegas_ts mcarlo_plain_ts expval_ts \
	rng_gsl_ts 
endif

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

mcarlo_miser_ts_LDADD = $(VCHECK_LIBS)
rng_gsl_ts_LDADD = $(VCHECK_LIBS)
mcarlo_plain_ts_LDADD = $(VCHECK_LIBS)
mcarlo_vegas_ts_LDADD = $(VCHECK_LIBS)
expval_ts_LDADD = $(VCHECK_LIBS)
mcmc_para_ts_LDFLAGS = -fopenmp
mcmc_para_ts_LDADD = $(VCHECK_LIBS)

mcarlo_miser.scr: mcarlo_miser_ts$(EXEEXT)
	./mcarlo_miser_ts$(EXEEXT) > mcarlo_miser.scr
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

mcarlo_miser_ts_SOURCES = mcarlo_miser_ts.cpp
rng_gsl_ts_SOURCES = rng_gsl_ts.cpp
mcarlo_plain_ts_SOURCES = mcarlo_plain_ts.cpp
mcarlo_vegas_ts_SOURCES = mcarlo_vegas_ts.cpp
expval_ts_SOURCES = expval_ts.cpp
mcmc_para_ts_SOURCES = mcmc_para_ts.cpp

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

