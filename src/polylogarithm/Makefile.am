# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

POLYLOGARITHM_SRCS = Cl2.cpp Cl3.cpp Cl4.cpp Cl5.cpp \
	Cl6.cpp Li.cpp Li2.cpp Li3.cpp Li4.cpp Li5.cpp Li6.cpp

HEADERS_VAR = Cl2.h Cl3.h Cl4.h Cl5.h Cl6.h \
	Li.h Li2.h Li3.h Li4.h Li5.h Li6.h complex.h

TEST_VAR = polylogarithm.scr polylogarithm_orig.scr

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

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADERS_VAR)

pkginclude_HEADERS = $(HEADERS_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = polylogarithm_ts polylogarithm_orig_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

polylogarithm_ts_LDADD = $(VCHECK_LIBS)

polylogarithm.scr: polylogarithm_ts$(EXEEXT)
	./polylogarithm_ts$(EXEEXT) > polylogarithm.scr

polylogarithm_ts_SOURCES = polylogarithm_ts.cpp

polylogarithm_orig_ts_LDADD = $(VCHECK_LIBS)

polylogarithm_orig.scr: polylogarithm_orig_ts$(EXEEXT)
	./polylogarithm_orig_ts$(EXEEXT) > polylogarithm_orig.scr

polylogarithm_orig_ts_SOURCES = polylogarithm_orig_ts.cpp

if O2SCL_OPENMP
polylogarithm_ts_LDFLAGS = -fopenmp
polylogarithm_orig_ts_LDFLAGS = -fopenmp
endif

# ------------------------------------------------------------
# No library o2scl_polylogarithm
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_polylogarithm.la

libo2scl_polylogarithm_la_SOURCES = $(POLYLOGARITHM_SRCS)

am__libo2scl_polylogarithm_la_SOURCES_DIST = $(POLYLOGARITHM_SRCS)

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
	echo "src/polylogarithm" $(TEST_VAR) >> ../../testlist

test-clean:
	rm -f *_ts.o *_ts *.scr

