# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

ANNEAL_SRCS = 

HEADER_VAR = anneal.h anneal_gsl.h anneal_para.h

if O2SCL_OPENMP
TEST_VAR = anneal_gsl.scr anneal_para.scr
else
TEST_VAR = anneal_gsl.scr
endif

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

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_ARMA
O2SCL_ARMA_MVAR = -DO2SCL_ARMA
else
O2SCL_ARMA_MVAR =
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

if O2SCL_OPENMP
check_PROGRAMS = anneal_gsl_ts anneal_para_ts
else
check_PROGRAMS = anneal_gsl_ts
endif

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


anneal_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
anneal_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
anneal_gsl_ts_CFLAGS = $(ADDL_TEST_FLGS)
anneal_gsl_ts_SOURCES = anneal_gsl_ts.cpp

anneal_gsl.scr: anneal_gsl_ts$(EXEEXT)
	./anneal_gsl_ts$(EXEEXT) > anneal_gsl.scr

anneal_para_ts_LDADD = $(ADDL_TEST_LIBS)
anneal_para_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
anneal_para_ts_CFLAGS = $(ADDL_TEST_FLGS)
anneal_para_ts_SOURCES = anneal_para_ts.cpp

anneal_para.scr: anneal_para_ts$(EXEEXT)
	./anneal_para_ts$(EXEEXT) > anneal_para.scr

# ------------------------------------------------------------
# No library o2scl_anneal
# ------------------------------------------------------------

# noinst_LTLIBRARIES = libo2scl_anneal.la

# libo2scl_anneal_la_SOURCES = $(ANNEAL_SRCS)

# am__libo2scl_anneal_la_SOURCES_DIST = $(ANNEAL_SRCS)

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
	echo "src/anneal" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

