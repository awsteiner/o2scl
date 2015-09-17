# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

CUBATURE_SRCS = 

HEADERS_VAR = clencurt.h cubature.h cubature_new.h

TEST_VAR = cubature.scr 
#cubature_new.scr

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

am__pkginclude_HEADERS_DIST = $(HEADERS_VAR)

pkginclude_HEADERS = $(HEADERS_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = cubature_ts 
#cubature_new_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

cubature_ts_LDADD = $(VCHECK_LIBS)
#cubature_new_ts_LDADD = $(VCHECK_LIBS)

cubature.scr: cubature_ts$(EXEEXT)
	./cubature_ts$(EXEEXT) > cubature.scr
#cubature_new.scr: cubature_new_ts$(EXEEXT)
#	./cubature_new_ts$(EXEEXT) > cubature_new.scr

cubature_ts_SOURCES = cubature_ts.cpp
#cubature_new_ts_SOURCES = cubature_new_ts.cpp

# ------------------------------------------------------------
# No library o2scl_cubature
# ------------------------------------------------------------

# noinst_LTLIBRARIES = libo2scl_cubature.la

# libo2scl_cubature_la_SOURCES = $(CUBATURE_SRCS)

# am__libo2scl_cubature_la_SOURCES_DIST = $(CUBATURE_SRCS)

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
	echo "src/cubature" $(TEST_VAR) >> ../../testlist

test-clean:
	rm -f *_ts.o *_ts *.scr

