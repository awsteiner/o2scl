# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

HEADER_VAR = hdf_file.h hdf_io.h cloud_file.h acolm.h hdf_python.h

HDF_SRCS = hdf_file.cpp hdf_io.cpp cloud_file.cpp acolm.cpp acolm_ac.cpp \
	acolm_o2graph.cpp acolm_df.cpp acolm_gi.cpp acolm_jo.cpp \
	acolm_ps.cpp acolm_tz.cpp hdf_python.cpp

TEST_VAR = hdf_file.scr hdf_io.scr

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

if O2SCL_READLINE
ADDINC = -DO2SCL_READLINE
else 
ADDINC = 
endif

if O2SCL_NCURSES
ADDINC = -DO2SCL_NCURSES
else 
ADDINC = 
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

if O2SCL_FFTW
O2SCL_FFTW_MVAR = -DO2SCL_FFTW
else 
O2SCL_FFTW_MVAR = 
endif

if O2SCL_PUGIXML
O2SCL_PUGIXML_MVAR = -DO2SCL_PUGIXML
else 
O2SCL_PUGIXML_MVAR = 
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(ADDINC) $(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR) \
	$(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(ADDINC) $(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR) \
	$(O2SCL_PYTHON_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

pkginclude_HEADERS = $(HEADER_VAR)

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

CPVAR = 

check_PROGRAMS = hdf_file_ts hdf_io_ts $(CPVAR)

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

hdf_file_ts_LDADD = $(ADDL_TEST_LIBS)
hdf_io_ts_LDADD = $(ADDL_TEST_LIBS)

hdf_file_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
hdf_io_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

hdf_file.scr: hdf_file_ts$(EXEEXT) 
	./hdf_file_ts$(EXEEXT) | tee hdf_file.scr
hdf_io.scr: hdf_io_ts$(EXEEXT) 
	./hdf_io_ts$(EXEEXT) | tee hdf_io.scr

hdf_file_ts_SOURCES = hdf_file_ts.cpp
hdf_io_ts_SOURCES = hdf_io_ts.cpp

# ------------------------------------------------------------
# Library o2scl_hdf
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_hdf.la

libo2scl_hdf_la_SOURCES = $(HDF_SRCS)

am__libo2scl_hdf_la_SOURCES_DIST = $(HDF_SRCS)

# ------------------------------------------------------------

interface: empty
	../../bin/yanic hdf.i hdf_python ~/wcs/o2sclpy/o2sclpy/hdf \
		~/wcs/o2sclpy/doc/hdf  ../../bin/yanic_header.txt \
		| tee hdf_python.yanic.out

empty:

# ------------------------------------------------------------
# Misc
# ------------------------------------------------------------

# These are the files included in the distribution
EXTRA_DIST = *_ts.cpp hdf_io_ts_table.txt

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

