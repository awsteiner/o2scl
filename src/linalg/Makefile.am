# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

LINALG_SRCS = permutation.cpp givens.cpp qr.cpp cholesky.cpp

HEADER_VAR = cblas.h lu.h lanczos.h tridiag.h permutation.h qr.h \
	householder.h givens.h hh.h cblas_base.h householder_base.h \
	qr_base.h hh_base.h tridiag_base.h givens_base.h lanczos_base.h \
	lu_base.h linear_solver.h svdstep_base.h svdstep.h cholesky.h \
	cholesky_base.h qrpt.h qrpt_base.h bidiag.h bidiag_base.h \
	svd.h svd_base.h

TEST_VAR = cblas.scr permutation.scr lanczos.scr tridiag.scr lu.scr \
	qr.scr qrpt.scr householder.scr hh.scr linear_solver.scr \
	cholesky.scr bidiag.scr svdstep.scr svd.scr

SUBDIRS = arma eigen both neither

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
AM_CPPFLAGS = -I@top_srcdir@/include/ \
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

check_PROGRAMS = cblas_ts permutation_ts lanczos_ts tridiag_ts lu_ts \
	qr_ts householder_ts hh_ts linear_solver_ts svdstep_ts \
	cholesky_ts qrpt_ts bidiag_ts svd_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

cblas_ts_LDADD = $(VCHECK_LIBS)
permutation_ts_LDADD = $(VCHECK_LIBS)
lanczos_ts_LDADD = $(VCHECK_LIBS)
tridiag_ts_LDADD = $(VCHECK_LIBS)
lu_ts_LDADD = $(VCHECK_LIBS)
qr_ts_LDADD = $(VCHECK_LIBS)
qrpt_ts_LDADD = $(VCHECK_LIBS)
householder_ts_LDADD = $(VCHECK_LIBS)
hh_ts_LDADD = $(VCHECK_LIBS)
linear_solver_ts_LDADD = $(VCHECK_LIBS)
svdstep_ts_LDADD = $(VCHECK_LIBS)
cholesky_ts_LDADD = $(VCHECK_LIBS)
svd_ts_LDADD = $(VCHECK_LIBS)
bidiag_ts_LDADD = $(VCHECK_LIBS)

if O2SCL_OPENMP
cblas_ts_LDFLAGS = -fopenmp
permutation_ts_LDFLAGS = -fopenmp
lanczos_ts_LDFLAGS = -fopenmp
tridiag_ts_LDFLAGS = -fopenmp
lu_ts_LDFLAGS = -fopenmp
qr_ts_LDFLAGS = -fopenmp
qrpt_ts_LDFLAGS = -fopenmp
householder_ts_LDFLAGS = -fopenmp
hh_ts_LDFLAGS = -fopenmp
linear_solver_ts_LDFLAGS = -fopenmp
svdstep_ts_LDFLAGS = -fopenmp
cholesky_ts_LDFLAGS = -fopenmp
svd_ts_LDFLAGS = -fopenmp
bidiag_ts_LDFLAGS = -fopenmp
else
cblas_ts_LDFLAGS =  
permutation_ts_LDFLAGS =  
lanczos_ts_LDFLAGS =  
tridiag_ts_LDFLAGS =  
lu_ts_LDFLAGS =  
qr_ts_LDFLAGS =  
qrpt_ts_LDFLAGS =  
householder_ts_LDFLAGS =  
hh_ts_LDFLAGS =  
linear_solver_ts_LDFLAGS =  
svdstep_ts_LDFLAGS =  
cholesky_ts_LDFLAGS =  
svd_ts_LDFLAGS =  
bidiag_ts_LDFLAGS =  
endif

cblas.scr: cblas_ts$(EXEEXT) 
	./cblas_ts$(EXEEXT) > cblas.scr
permutation.scr: permutation_ts$(EXEEXT) 
	./permutation_ts$(EXEEXT) > permutation.scr
lanczos.scr: lanczos_ts$(EXEEXT) 
	./lanczos_ts$(EXEEXT) > lanczos.scr
tridiag.scr: tridiag_ts$(EXEEXT) 
	./tridiag_ts$(EXEEXT) > tridiag.scr
lu.scr: lu_ts$(EXEEXT) 
	./lu_ts$(EXEEXT) > lu.scr
qr.scr: qr_ts$(EXEEXT) 
	./qr_ts$(EXEEXT) > qr.scr
qrpt.scr: qrpt_ts$(EXEEXT) 
	./qrpt_ts$(EXEEXT) > qrpt.scr
householder.scr: householder_ts$(EXEEXT) 
	./householder_ts$(EXEEXT) > householder.scr
hh.scr: hh_ts$(EXEEXT) 
	./hh_ts$(EXEEXT) > hh.scr
linear_solver.scr: linear_solver_ts$(EXEEXT) 
	./linear_solver_ts$(EXEEXT) > linear_solver.scr
svdstep.scr: svdstep_ts$(EXEEXT) 
	./svdstep_ts$(EXEEXT) > svdstep.scr
cholesky.scr: cholesky_ts$(EXEEXT) 
	./cholesky_ts$(EXEEXT) > cholesky.scr
svd.scr: svd_ts$(EXEEXT) 
	./svd_ts$(EXEEXT) > svd.scr
bidiag.scr: bidiag_ts$(EXEEXT) 
	./bidiag_ts$(EXEEXT) > bidiag.scr

cblas_ts_SOURCES = cblas_ts.cpp
permutation_ts_SOURCES = permutation_ts.cpp
lanczos_ts_SOURCES = lanczos_ts.cpp
tridiag_ts_SOURCES = tridiag_ts.cpp
lu_ts_SOURCES = lu_ts.cpp
qr_ts_SOURCES = qr_ts.cpp
qrpt_ts_SOURCES = qrpt_ts.cpp
householder_ts_SOURCES = householder_ts.cpp
hh_ts_SOURCES = hh_ts.cpp
linear_solver_ts_SOURCES = linear_solver_ts.cpp
svdstep_ts_SOURCES = svdstep_ts.cpp
cholesky_ts_SOURCES = cholesky_ts.cpp
svd_ts_SOURCES = svd_ts.cpp
bidiag_ts_SOURCES = bidiag_ts.cpp

# ------------------------------------------------------------
# Library o2scl_linalg
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_linalg.la

libo2scl_linalg_la_SOURCES = $(LINALG_SRCS)

am__libo2scl_linalg_la_SOURCES_DIST = $(LINALG_SRCS)

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

../internal/tsumm$(EXEEXT):
	cd ../internal; $(MAKE) tsumm

o2scl-test: $(TEST_VAR) ../internal/tsumm$(EXEEXT)
	echo $(TEST_VAR) > testlist
	../internal/tsumm$(EXEEXT)

o2scl-test-rec: $(TEST_VAR)
	echo "src/linalg" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

svdd:
	make install; make test-clean; make svd_ts; \
		svd_ts gsl > old.txt; svd_ts o2 > new.txt; \
		diff old.txt new.txt | more

