# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

MCARLO_SRCS = expval.cpp rng_gsl.cpp

HEADER_VARS = mcarlo_miser.h mcarlo_plain.h mcarlo_vegas.h mcarlo.h \
	expval.h rng_gsl.h

TEST_VAR = mcarlo_plain.scr mcarlo_miser.scr mcarlo_vegas.scr \
	expval.scr rng_gsl.scr

# fit_bayes.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_PARTLIB
O2SCL_PART_MVAR = 1
else
O2SCL_PART_MVAR = 0
endif

if O2SCL_EOSLIB
O2SCL_EOS_MVAR = 1
else
O2SCL_EOS_MVAR = 0
endif

if O2SCL_HDF
O2SCL_HDF_MVAR = 1
else
O2SCL_HDF_MVAR = 0
endif

if O2SCL_OPENMP
if O2SCL_EIGEN
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_ARMA -DO2SCL_EIGEN \
	-fopenmp -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_EIGEN -fopenmp \
	-DO2SCL_COND_FLAG
endif # ARMA
else
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_ARMA -fopenmp \
	-DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -fopenmp \
	-DO2SCL_COND_FLAG
endif # ARMA
endif # EIGEN
else
if O2SCL_EIGEN
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_ARMA -DO2SCL_EIGEN -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_EIGEN -DO2SCL_COND_FLAG
endif # ARMA
else
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_ARMA -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_COND_FLAG
endif # ARMA
endif # EIGEN
endif # OPENMP

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VARS)

pkginclude_HEADERS = $(HEADER_VARS)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = mcarlo_miser_ts mcarlo_vegas_ts mcarlo_plain_ts expval_ts \
	rng_gsl_ts

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl.la $(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl.la
endif

mcarlo_miser_ts_LDADD = $(VCHECK_LIBS)
rng_gsl_ts_LDADD = $(VCHECK_LIBS)
mcarlo_plain_ts_LDADD = $(VCHECK_LIBS)
mcarlo_vegas_ts_LDADD = $(VCHECK_LIBS)
expval_ts_LDADD = $(VCHECK_LIBS)

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

mcarlo_miser_ts_SOURCES = mcarlo_miser_ts.cpp
rng_gsl_ts_SOURCES = rng_gsl_ts.cpp
mcarlo_plain_ts_SOURCES = mcarlo_plain_ts.cpp
mcarlo_vegas_ts_SOURCES = mcarlo_vegas_ts.cpp
expval_ts_SOURCES = expval_ts.cpp

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

