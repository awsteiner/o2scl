# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

FIT_SRCS = 

HEADER_VAR = fit_base.h fit_nonlin.h fit_min.h fit_fix.h fit_linear.h \
	fit_bayes.h

TEST_VAR = fit_base.scr fit_min.scr fit_fix.scr fit_linear.scr \
	fit_nonlin.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -DO2SCL_OPENMP -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ -I$(PYTHON_INCLUDE_DIR) \
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

check_PROGRAMS = fit_base_ts fit_min_ts fit_fix_ts fit_nonlin_ts \
	fit_linear_ts

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

fit_base_ts_LDADD = $(ADDL_TEST_LIBS)
fit_nonlin_ts_LDADD = $(ADDL_TEST_LIBS)
fit_min_ts_LDADD = $(ADDL_TEST_LIBS)
fit_fix_ts_LDADD = $(ADDL_TEST_LIBS)
fit_linear_ts_LDADD = $(ADDL_TEST_LIBS)

fit_base_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fit_nonlin_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fit_min_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fit_fix_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fit_linear_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

fit_base.scr: fit_base_ts$(EXEEXT)
	./fit_base_ts$(EXEEXT) > fit_base.scr
fit_nonlin.scr: fit_nonlin_ts$(EXEEXT)
	./fit_nonlin_ts$(EXEEXT) > fit_nonlin.scr
fit_min.scr: fit_min_ts$(EXEEXT)
	./fit_min_ts$(EXEEXT) > fit_min.scr
fit_fix.scr: fit_fix_ts$(EXEEXT)
	./fit_fix_ts$(EXEEXT) > fit_fix.scr
fit_linear.scr: fit_linear_ts$(EXEEXT)
	./fit_linear_ts$(EXEEXT) > fit_linear.scr

fit_base_ts_SOURCES = fit_base_ts.cpp
fit_nonlin_ts_SOURCES = fit_nonlin_ts.cpp
fit_min_ts_SOURCES = fit_min_ts.cpp
fit_fix_ts_SOURCES = fit_fix_ts.cpp
fit_linear_ts_SOURCES = fit_linear_ts.cpp

# ------------------------------------------------------------
# No library o2scl_fit
# ------------------------------------------------------------

# lib_LTLIBRARIES = libo2scl_fit.la

# libo2scl_fit_la_SOURCES = $(FIT_SRCS)

# am__libo2scl_fit_la_SOURCES_DIST = $(FIT_SRCS)

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
	echo "src/fit" $(TEST_VAR) >> ../../testlist

test-clean: 
	rm -f *_ts.o *_ts *.scr

