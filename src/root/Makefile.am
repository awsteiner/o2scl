# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

ROOT_SRCS = 

HEADER_VAR = root_bkt_cern.h root.h root_cern.h mroot.h mroot_hybrids.h \
	root_stef.h root_brent_gsl.h mroot_cern.h \
	jacobian.h mroot_broyden.h root_toms748.h root_robbins_monro.h

TEST_VAR = root_bkt_cern.scr mroot_cern.scr mroot_hybrids.scr \
	root_stef.scr root_cern.scr root_brent_gsl.scr \
	jacobian.scr mroot_broyden.scr root_toms748.scr

SUBDIRS = arma eigen neither both

testlist:
	echo $(TEST_VAR) > testlist

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR)
endif

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

if O2SCL_PYTHON

ADDL_TEST_LIBS = ../libo2scl.la -l$(PYTHON_LIB)
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
endif

else

ADDL_TEST_LIBS = ../libo2scl.la
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = 
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = 
endif

endif

root_bkt_cern_ts_LDADD = $(ADDL_TEST_LIBS)
mroot_cern_ts_LDADD = $(ADDL_TEST_LIBS)
mroot_hybrids_ts_LDADD = $(ADDL_TEST_LIBS)
mroot_broyden_ts_LDADD = $(ADDL_TEST_LIBS)
root_stef_ts_LDADD = $(ADDL_TEST_LIBS)
root_cern_ts_LDADD = $(ADDL_TEST_LIBS)
root_brent_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
root_toms748_ts_LDADD = $(ADDL_TEST_LIBS)
jacobian_ts_LDADD = $(ADDL_TEST_LIBS)

root_bkt_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mroot_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mroot_hybrids_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mroot_broyden_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
root_stef_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
root_cern_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
root_brent_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
root_toms748_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
jacobian_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

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

