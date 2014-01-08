# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

HEADER_VAR = hdf_file.h hdf_io.h

HDF_SRCS = hdf_file.cpp hdf_io.cpp

TEST_VAR = hdf_file.scr hdf_io.scr

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

if O2SCL_CPP11
O2SCL_CPP11_MVAR = -std=c++0x
else
O2SCL_CPP11_MVAR = -DO2SCL_NO_CPP11
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_CPP11_MVAR) -DO2SCL_COND_FLAG

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

pkginclude_HEADERS = $(HEADER_VAR)

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

CPVAR = 

if O2SCL_HDF
check_PROGRAMS = hdf_file_ts hdf_io_ts $(CPVAR)
else
check_PROGRAMS = $(CPVAR)
endif

VCHECK_LIBS = $(libdir)/libo2scl_hdf.la $(libdir)/libo2scl.la 

hdf_file_ts_LDADD = $(VCHECK_LIBS)
hdf_io_ts_LDADD = $(VCHECK_LIBS)

hdf_file.scr: hdf_file_ts$(EXEEXT) 
	./hdf_file_ts$(EXEEXT) > hdf_file.scr
hdf_io.scr: hdf_io_ts$(EXEEXT) 
	./hdf_io_ts$(EXEEXT) > hdf_io.scr

hdf_file_ts_SOURCES = hdf_file_ts.cpp
hdf_io_ts_SOURCES = hdf_io_ts.cpp

# ------------------------------------------------------------
# Library o2scl_hdf
# ------------------------------------------------------------

lib_LTLIBRARIES = libo2scl_hdf.la

libo2scl_hdf_la_SOURCES = $(HDF_SRCS)
libo2scl_hdf_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

am__libo2scl_hdf_la_SOURCES_DIST = $(HDF_SRCS)

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
	echo "src/hdf" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

