# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

# The files nucmass, nucdist, nucmass_fit, reaction_lib, and
# hdf_nucmass_io can only be tested if HDF is enabled

if O2SCL_HDF

NUCLEI_SRCS = nucmass.cpp nucdist.cpp nucleus.cpp nucmass_fit.cpp \
	reaction_lib.cpp hdf_nucmass_io.cpp nucmass_frdm.cpp \
	nucmass_dz.cpp nucmass_ame.cpp nucmass_hfb.cpp nucmass_ktuy.cpp \
	nucmass_dglg.cpp nucmass_wlw.cpp nucmass_sdnp.cpp \
	nucmass_densmat.cpp

HEADER_VAR = nucmass.h nucdist.h nucleus.h nucmass_fit.h \
	reaction_lib.h hdf_nucmass_io.h nucmass_frdm.h \
	nucmass_dz.h nucmass_ame.h nucmass_hfb.h nucmass_ktuy.h \
	nucmass_dglg.h nucmass_wlw.h nucmass_sdnp.h \
	nucmass_densmat.h

TEST_VAR = nucmass.scr nucdist.scr nucleus.scr nucmass_fit.scr \
	reaction_lib.scr hdf_nucmass_io.scr nucmass_frdm.scr \
	nucmass_dz.scr nucmass_densmat.scr

else

NUCLEI_SRCS = nucmass.cpp nucdist.cpp nucleus.cpp nucmass_fit.cpp \
	reaction_lib.cpp nucmass_frdm.cpp nucmass_dz.cpp \
	nucmass_ame.cpp nucmass_hfb.cpp nucmass_ktuy.cpp \
	nucmass_dglg.cpp nucmass_wlw.cpp nucmass_sdnp.cpp \
	nucmass_densmat.cpp

HEADER_VAR = nucmass.h nucdist.h nucleus.h nucmass_fit.h reaction_lib.h \
	nucmass_frdm.h nucmass_dz.h nucmass_ame.h \
	nucmass_hfb.h nucmass_ktuy.h nucmass_dglg.h \
	nucmass_wlw.h nucmass_sdnp.h nucmass_densmat.h

TEST_VAR = nucleus.scr

endif

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
AM_CPPFLAGS = -I@top_srcdir@/include/
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

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = nucmass_ts nucdist_ts nucleus_ts \
	nucmass_fit_ts reaction_lib_ts hdf_nucmass_io_ts nucmass_frdm_ts \
	nucmass_dz_ts nucmass_densmat_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../libo2scl_part.la ../libo2scl.la \
	../hdf/libo2scl_hdf.la
else
VCHECK_LIBS = ../libo2scl_part.la ../libo2scl.la 
endif

nucmass_ts_LDADD = $(VCHECK_LIBS)
reaction_lib_ts_LDADD = $(VCHECK_LIBS)
hdf_nucmass_io_ts_LDADD = $(VCHECK_LIBS)
nucdist_ts_LDADD = $(VCHECK_LIBS)
nucleus_ts_LDADD = $(VCHECK_LIBS)
nucmass_fit_ts_LDADD = $(VCHECK_LIBS)
nucmass_frdm_ts_LDADD = $(VCHECK_LIBS)
nucmass_dz_ts_LDADD = $(VCHECK_LIBS)
nucmass_densmat_ts_LDADD = $(VCHECK_LIBS)

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
nucmass_dz_ts_SOURCES = nucmass_dz_ts.cpp
nucmass_densmat_ts_SOURCES = nucmass_densmat_ts.cpp

# ------------------------------------------------------------
# Library o2scl_nuclei
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_nuclei.la

libo2scl_nuclei_la_SOURCES = $(NUCLEI_SRCS)

am__libo2scl_nuclei_la_SOURCES_DIST = $(NUCLEI_SRCS)

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

