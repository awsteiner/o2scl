# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

# The files nuclear_mass, nuclear_dist, mass_fit, reaction_lib, and
# hdf_nucmass_io can only be tested if HDF is enabled

if O2SCL_HDF

NUCLEI_SRCS = nuclear_mass.cpp nuclear_dist.cpp nucleus.cpp mass_fit.cpp \
	reaction_lib.cpp hdf_nucmass_io.cpp frdm_mass.cpp dz_mass.cpp \
	ame_mass.cpp hfb_mass.cpp ktuy_mass.cpp

HEADER_VAR = nuclear_mass.h nuclear_dist.h nucleus.h mass_fit.h \
	reaction_lib.h hdf_nucmass_io.h frdm_mass.h arb_dist.h dz_mass.h \
	ame_mass.h hfb_mass.h ktuy_mass.h

TEST_VAR = nuclear_mass.scr nuclear_dist.scr nucleus.scr mass_fit.scr \
	reaction_lib.scr hdf_nucmass_io.scr frdm_mass.scr arb_dist.scr \
	dz_mass.scr

else

NUCLEI_SRCS = nuclear_mass.cpp nuclear_dist.cpp nucleus.cpp mass_fit.cpp \
	reaction_lib.cpp frdm_mass.cpp dz_mass.cpp ame_mass.cpp \
	hfb_mass.cpp ktuy_mass.cpp

HEADER_VAR = nuclear_mass.h nuclear_dist.h nucleus.h mass_fit.h \
	reaction_lib.h frdm_mass.h arb_dist.h dz_mass.h ame_mass.h \
	hfb_mass.h ktuy_mass.h

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

if O2SCL_CPP11
O2SCL_CPP11_MVAR = 
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

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = nuclear_mass_ts nuclear_dist_ts nucleus_ts \
	mass_fit_ts reaction_lib_ts hdf_nucmass_io_ts frdm_mass_ts \
	arb_dist_ts dz_mass_ts

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl_part.la $(libdir)/libo2scl.la \
	$(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl_part.la $(libdir)/libo2scl.la 
endif

nuclear_mass_ts_LDADD = $(VCHECK_LIBS)
reaction_lib_ts_LDADD = $(VCHECK_LIBS)
hdf_nucmass_io_ts_LDADD = $(VCHECK_LIBS)
nuclear_dist_ts_LDADD = $(VCHECK_LIBS)
nucleus_ts_LDADD = $(VCHECK_LIBS)
mass_fit_ts_LDADD = $(VCHECK_LIBS)
frdm_mass_ts_LDADD = $(VCHECK_LIBS)
arb_dist_ts_LDADD = $(VCHECK_LIBS)
dz_mass_ts_LDADD = $(VCHECK_LIBS)

nuclear_mass.scr: nuclear_mass_ts$(EXEEXT)
	./nuclear_mass_ts$(EXEEXT) > nuclear_mass.scr
reaction_lib.scr: reaction_lib_ts$(EXEEXT)
	./reaction_lib_ts$(EXEEXT) > reaction_lib.scr
hdf_nucmass_io.scr: hdf_nucmass_io_ts$(EXEEXT)
	./hdf_nucmass_io_ts$(EXEEXT) > hdf_nucmass_io.scr
nuclear_dist.scr: nuclear_dist_ts$(EXEEXT)
	./nuclear_dist_ts$(EXEEXT) > nuclear_dist.scr
nucleus.scr: nucleus_ts$(EXEEXT)
	./nucleus_ts$(EXEEXT) > nucleus.scr
mass_fit.scr: mass_fit_ts$(EXEEXT)
	./mass_fit_ts$(EXEEXT) > mass_fit.scr
frdm_mass.scr: frdm_mass_ts$(EXEEXT)
	./frdm_mass_ts$(EXEEXT) > frdm_mass.scr
arb_dist.scr: arb_dist_ts$(EXEEXT)
	./arb_dist_ts$(EXEEXT) > arb_dist.scr
dz_mass.scr: dz_mass_ts$(EXEEXT)
	./dz_mass_ts$(EXEEXT) > dz_mass.scr

nuclear_mass_ts_SOURCES = nuclear_mass_ts.cpp
reaction_lib_ts_SOURCES = reaction_lib_ts.cpp
hdf_nucmass_io_ts_SOURCES = hdf_nucmass_io_ts.cpp
nuclear_dist_ts_SOURCES = nuclear_dist_ts.cpp
nucleus_ts_SOURCES = nucleus_ts.cpp
mass_fit_ts_SOURCES = mass_fit_ts.cpp
frdm_mass_ts_SOURCES = frdm_mass_ts.cpp
arb_dist_ts_SOURCES = arb_dist_ts.cpp
dz_mass_ts_SOURCES = dz_mass_ts.cpp

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

