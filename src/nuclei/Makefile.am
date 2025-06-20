# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

NUCLEI_SRCS = nucmass.cpp nucdist.cpp nucleus.cpp nucmass_fit.cpp \
	reaction_lib.cpp hdf_nucmass_io.cpp nucmass_frdm.cpp \
	nucmass_dz.cpp nucmass_ame.cpp nucmass_hfb.cpp nucmass_ktuy.cpp \
	nucmass_dglg.cpp nucmass_wlw.cpp nucmass_sdnp.cpp \
	nucmass_densmat.cpp nucmass_gen.cpp nuclei_python.cpp \
	part_funcs.cpp nucmass_two_interp.cpp

HEADER_VAR = nucmass.h nucdist.h nucleus.h nucmass_fit.h \
	reaction_lib.h hdf_nucmass_io.h nucmass_frdm.h \
	nucmass_dz.h nucmass_ame.h nucmass_hfb.h nucmass_ktuy.h \
	nucmass_dglg.h nucmass_wlw.h nucmass_sdnp.h \
	nucmass_densmat.h nucmass_gen.h nuclei_python.h \
	part_funcs.h nucmass_two_interp.h

TEST_VAR = nucmass.scr nucdist.scr nucleus.scr nucmass_fit.scr \
	reaction_lib.scr hdf_nucmass_io.scr nucmass_frdm.scr \
	nucmass_dz.scr nucmass_densmat.scr nucmass_gen.scr nucmass_ame.scr \
	part_funcs.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -fopenmp
else
O2SCL_OPENMP_MVAR =
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

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = nucmass_ts nucdist_ts nucleus_ts \
	nucmass_fit_ts reaction_lib_ts hdf_nucmass_io_ts nucmass_frdm_ts \
	nucmass_dz_ts nucmass_densmat_ts nucmass_gen_ts nucmass_ame_ts \
	part_funcs_ts

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

nucmass_ts_LDADD = $(ADDL_TEST_LIBS)
reaction_lib_ts_LDADD = $(ADDL_TEST_LIBS)
hdf_nucmass_io_ts_LDADD = $(ADDL_TEST_LIBS)
nucdist_ts_LDADD = $(ADDL_TEST_LIBS)
nucleus_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_fit_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_frdm_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_gen_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_ame_ts_LDADD = $(ADDL_TEST_LIBS)
part_funcs_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_dz_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_densmat_ts_LDADD = $(ADDL_TEST_LIBS)

nucmass_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
reaction_lib_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
hdf_nucmass_io_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucdist_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucleus_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_fit_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_frdm_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_gen_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_ame_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
part_funcs_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_dz_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_densmat_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

nucmass.scr: nucmass_ts$(EXEEXT)
	./nucmass_ts$(EXEEXT) > nucmass.scr
reaction_lib.scr: reaction_lib_ts$(EXEEXT)
	./reaction_lib_ts$(EXEEXT) > reaction_lib.scr
hdf_nucmass_io.scr: hdf_nucmass_io_ts$(EXEEXT)
	./hdf_nucmass_io_ts$(EXEEXT) > hdf_nucmass_io.scr
nucdist.scr: nucdist_ts$(EXEEXT)
	./nucdist_ts$(EXEEXT) > nucdist.scr
nucleus.scr: nucleus_ts$(EXEEXT)
	./nucleus_ts$(EXEEXT) > nucleus.scr
nucmass_fit.scr: nucmass_fit_ts$(EXEEXT)
	./nucmass_fit_ts$(EXEEXT) > nucmass_fit.scr
nucmass_frdm.scr: nucmass_frdm_ts$(EXEEXT)
	./nucmass_frdm_ts$(EXEEXT) > nucmass_frdm.scr
nucmass_gen.scr: nucmass_gen_ts$(EXEEXT)
	./nucmass_gen_ts$(EXEEXT) > nucmass_gen.scr
nucmass_ame.scr: nucmass_ame_ts$(EXEEXT)
	./nucmass_ame_ts$(EXEEXT) > nucmass_ame.scr
part_funcs.scr: part_funcs_ts$(EXEEXT)
	./part_funcs_ts$(EXEEXT) > part_funcs.scr
nucmass_dz.scr: nucmass_dz_ts$(EXEEXT)
	./nucmass_dz_ts$(EXEEXT) > nucmass_dz.scr
nucmass_densmat.scr: nucmass_densmat_ts$(EXEEXT)
	./nucmass_densmat_ts$(EXEEXT) > nucmass_densmat.scr

nucmass_ts_SOURCES = nucmass_ts.cpp
reaction_lib_ts_SOURCES = reaction_lib_ts.cpp
hdf_nucmass_io_ts_SOURCES = hdf_nucmass_io_ts.cpp
nucdist_ts_SOURCES = nucdist_ts.cpp
nucleus_ts_SOURCES = nucleus_ts.cpp
nucmass_fit_ts_SOURCES = nucmass_fit_ts.cpp
nucmass_frdm_ts_SOURCES = nucmass_frdm_ts.cpp
nucmass_gen_ts_SOURCES = nucmass_gen_ts.cpp
nucmass_ame_ts_SOURCES = nucmass_ame_ts.cpp
part_funcs_ts_SOURCES = part_funcs_ts.cpp
nucmass_dz_ts_SOURCES = nucmass_dz_ts.cpp
nucmass_densmat_ts_SOURCES = nucmass_densmat_ts.cpp

# ------------------------------------------------------------
# Library o2scl_nuclei
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_nuclei.la

libo2scl_nuclei_la_SOURCES = $(NUCLEI_SRCS)

am__libo2scl_nuclei_la_SOURCES_DIST = $(NUCLEI_SRCS)

# ------------------------------------------------------------

interface: empty
	../../bin/yanic nuclei.i nuclei_python ~/wcs/o2sclpy/o2sclpy/nuclei \
		~/wcs/o2sclpy/doc/nuclei  ../../bin/yanic_header.txt \
		| tee nuclei_python.yanic.out

empty:

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
	echo "src/nuclei" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

