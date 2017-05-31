# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

ANNEAL_SRCS = 

HEADER_VAR = anneal.h anneal_gsl.h anneal_para.h

TEST_VAR = anneal_gsl.scr anneal_para.scr

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

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

anneal_gsl_ts_LDADD = $(VCHECK_LIBS)

anneal_gsl.scr: anneal_gsl_ts$(EXEEXT)
	./anneal_gsl_ts$(EXEEXT) > anneal_gsl.scr

anneal_gsl_ts_SOURCES = anneal_gsl_ts.cpp

if O2SCL_OPENMP

anneal_para_ts_LDADD = $(VCHECK_LIBS)
anneal_para_ts_LDFLAGS = -fopenmp
anneal_para_ts_CFLAGS = -DO2SCL_OPENMP

anneal_para.scr: anneal_para_ts$(EXEEXT)
	./anneal_para_ts$(EXEEXT) > anneal_para.scr

anneal_para_ts_SOURCES = anneal_para_ts.cpp

endif

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

