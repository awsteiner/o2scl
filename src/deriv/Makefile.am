# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

DERIV_SRCS = 

HEADERS_VAR = deriv_cern.h deriv.h deriv_eqi.h deriv_gsl.h \
	vector_derint.h

TEST_VAR = deriv_cern.scr deriv_gsl.scr deriv_eqi.scr vector_derint.scr

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

AM_CPPFLAGS = -I@includedir@ \
	-I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADERS_VAR)

pkginclude_HEADERS = $(HEADERS_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = deriv_cern_ts deriv_eqi_ts deriv_gsl_ts vector_derint_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl.la $(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl.la
endif

deriv_cern_ts_LDADD = $(VCHECK_LIBS)
deriv_gsl_ts_LDADD = $(VCHECK_LIBS)
deriv_eqi_ts_LDADD = $(VCHECK_LIBS)
vector_derint_ts_LDADD = $(VCHECK_LIBS)

deriv_cern.scr: deriv_cern_ts$(EXEEXT)
	./deriv_cern_ts$(EXEEXT) > deriv_cern.scr
deriv_gsl.scr: deriv_gsl_ts$(EXEEXT)
	./deriv_gsl_ts$(EXEEXT) > deriv_gsl.scr
deriv_eqi.scr: deriv_eqi_ts$(EXEEXT)
	./deriv_eqi_ts$(EXEEXT) > deriv_eqi.scr
vector_derint.scr: vector_derint_ts$(EXEEXT)
	./vector_derint_ts$(EXEEXT) > vector_derint.scr

deriv_cern_ts_SOURCES = deriv_cern_ts.cpp
deriv_gsl_ts_SOURCES = deriv_gsl_ts.cpp
deriv_eqi_ts_SOURCES = deriv_eqi_ts.cpp
vector_derint_ts_SOURCES = vector_derint_ts.cpp

# ------------------------------------------------------------
# No library o2scl_deriv
# ------------------------------------------------------------

#am__libo2scl_deriv_la_SOURCES_DIST = $(DERIV_SRCS)

#noinst_LTLIBRARIES = libo2scl_deriv.la

#libo2scl_deriv_la_SOURCES = $(DERIV_SRCS)

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
	echo "src/deriv" $(TEST_VAR) >> ../../testlist

test-clean:
	rm -f *_ts.o *_ts *.scr

