# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

FIT_SRCS = 

HEADER_VAR = fit_base.h fit_nonlin.h fit_min.h fit_fix.h fit_linear.h \
	fit_bayes.h

TEST_VAR = fit_base.scr fit_min.scr fit_fix.scr fit_linear.scr \
	fit_nonlin.scr

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

if O2SCL_GSL2
O2SCL_GSL2_MVAR = -DO2SCL_GSL2
else
O2SCL_GSL2_MVAR =
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = fit_base_ts fit_min_ts fit_fix_ts fit_nonlin_ts \
	fit_linear_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

fit_base_ts_LDADD = $(VCHECK_LIBS)
fit_nonlin_ts_LDADD = $(VCHECK_LIBS)
fit_min_ts_LDADD = $(VCHECK_LIBS)
fit_fix_ts_LDADD = $(VCHECK_LIBS)
fit_linear_ts_LDADD = $(VCHECK_LIBS)
#fit_bayes_ts_LDADD = $(VCHECK_LIBS)

#fit_bayes.scr: fit_bayes_ts$(EXEEXT)
#	./fit_bayes_ts$(EXEEXT) > fit_bayes.scr
fit_base.scr: fit_base_ts$(EXEEXT)
	./fit_base_ts$(EXEEXT) > fit_base.scr
fit_nonlin.scr: fit_nonlin_ts$(EXEEXT)
	./fit_nonlin_ts$(EXEEXT) > fit_nonlin.scr
fit_min.scr: fit_min_ts$(EXEEXT)
	./fit_min_ts$(EXEEXT) > fit_min.scr
fit_fix.scr: fit_fix_ts$(EXEEXT)
	./fit_fix_ts$(EXEEXT) > fit_fix.scr
fit_linear.scr: fit_linear_ts$(EXEEXT)
	./fit_linear_ts$(EXEEXT) > fit_linear.scr

#fit_bayes_ts_SOURCES = fit_bayes_ts.cpp
fit_base_ts_SOURCES = fit_base_ts.cpp
fit_nonlin_ts_SOURCES = fit_nonlin_ts.cpp
fit_min_ts_SOURCES = fit_min_ts.cpp
fit_fix_ts_SOURCES = fit_fix_ts.cpp
fit_linear_ts_SOURCES = fit_linear_ts.cpp

# ------------------------------------------------------------
# No library o2scl_fit
# ------------------------------------------------------------

# lib_LTLIBRARIES = libo2scl_fit.la

# libo2scl_fit_la_SOURCES = $(FIT_SRCS)
# libo2scl_fit_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

# am__libo2scl_fit_la_SOURCES_DIST = $(FIT_SRCS)

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
	echo "src/fit" $(TEST_VAR) >> ../../testlist

test-clean: 
	rm -f *_ts.o *_ts *.scr

