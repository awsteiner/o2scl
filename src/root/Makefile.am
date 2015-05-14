# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

ROOT_SRCS = 

HEADER_VAR = root_bkt_cern.h root.h root_cern.h mroot.h mroot_hybrids.h \
	root_stef.h root_brent_gsl.h mroot_cern.h \
	jacobian.h mroot_broyden.h root_toms748.h

TEST_VAR = root_bkt_cern.scr mroot_cern.scr mroot_hybrids.scr \
	root_stef.scr root_cern.scr root_brent_gsl.scr \
	jacobian.scr mroot_broyden.scr root_toms748.scr

SUBDIRS = arma eigen neither both

testlist:
	echo $(TEST_VAR) > testlist

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

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = root_bkt_cern_ts mroot_cern_ts mroot_hybrids_ts \
	root_stef_ts root_cern_ts root_brent_gsl_ts jacobian_ts \
	mroot_broyden_ts root_toms748_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl.la $(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl.la
endif

root_bkt_cern_ts_LDADD = $(VCHECK_LIBS)
mroot_cern_ts_LDADD = $(VCHECK_LIBS)
mroot_hybrids_ts_LDADD = $(VCHECK_LIBS)
mroot_broyden_ts_LDADD = $(VCHECK_LIBS)
root_stef_ts_LDADD = $(VCHECK_LIBS)
root_cern_ts_LDADD = $(VCHECK_LIBS)
root_brent_gsl_ts_LDADD = $(VCHECK_LIBS)
root_toms748_ts_LDADD = $(VCHECK_LIBS)
jacobian_ts_LDADD = $(VCHECK_LIBS)

root_bkt_cern.scr: root_bkt_cern_ts$(EXEEXT) 
	./root_bkt_cern_ts$(EXEEXT) > root_bkt_cern.scr
mroot_cern.scr: mroot_cern_ts$(EXEEXT) 
	./mroot_cern_ts$(EXEEXT) > mroot_cern.scr
mroot_hybrids.scr: mroot_hybrids_ts$(EXEEXT) 
	./mroot_hybrids_ts$(EXEEXT) > mroot_hybrids.scr
mroot_broyden.scr: mroot_broyden_ts$(EXEEXT) 
	./mroot_broyden_ts$(EXEEXT) > mroot_broyden.scr
root_stef.scr: root_stef_ts$(EXEEXT) 
	./root_stef_ts$(EXEEXT) > root_stef.scr
root_cern.scr: root_cern_ts$(EXEEXT) 
	./root_cern_ts$(EXEEXT) > root_cern.scr
root_brent_gsl.scr: root_brent_gsl_ts$(EXEEXT) 
	./root_brent_gsl_ts$(EXEEXT) > root_brent_gsl.scr
root_toms748.scr: root_toms748_ts$(EXEEXT) 
	./root_toms748_ts$(EXEEXT) > root_toms748.scr
jacobian.scr: jacobian_ts$(EXEEXT) 
	./jacobian_ts$(EXEEXT) > jacobian.scr

root_bkt_cern_ts_SOURCES = root_bkt_cern_ts.cpp
mroot_cern_ts_SOURCES = mroot_cern_ts.cpp
mroot_hybrids_ts_SOURCES = mroot_hybrids_ts.cpp
mroot_broyden_ts_SOURCES = mroot_broyden_ts.cpp
root_stef_ts_SOURCES = root_stef_ts.cpp
root_cern_ts_SOURCES = root_cern_ts.cpp
root_brent_gsl_ts_SOURCES = root_brent_gsl_ts.cpp
root_toms748_ts_SOURCES = root_toms748_ts.cpp
jacobian_ts_SOURCES = jacobian_ts.cpp

# ------------------------------------------------------------
# No library o2scl_root
# ------------------------------------------------------------

# noinst_LTLIBRARIES = libo2scl_root.la

# libo2scl_root_la_SOURCES = $(ROOT_SRCS)

# am__libo2scl_root_la_SOURCES_DIST = $(ROOT_SRCS)

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
	echo "src/root" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

