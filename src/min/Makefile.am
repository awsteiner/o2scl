# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

MIN_SRCS = mmin_bfgs2.cpp

HEADER_VAR = min.h min_cern.h min_brent_gsl.h min_brent_boost.h mmin_fix.h \
	mmin_simp.h mmin.h mmin_conf.h mmin_simp2.h \
	mmin_conp.h mmin_bfgs2.h mmin_constr.h mmin_constr_pgrad.h \
	mmin_constr_spg.h mmin_constr_gencan.h min_quad_golden.h diff_evo.h \
	diff_evo_adapt.h

TEST_VAR = min_cern.scr min_brent_gsl.scr min_brent_boost.scr mmin_simp.scr \
	mmin_conf.scr mmin_conp.scr mmin_bfgs2.scr \
	mmin_fix.scr mmin_constr_pgrad.scr mmin_constr_spg.scr \
	min.scr mmin_simp2.scr min_quad_golden.scr diff_evo.scr \
	diff_evo_adapt.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_PARTLIB
O2SCL_PART_MVAR = 1
else
O2SCL_PART_MVAR = 0
endif

if O2SCL_EOSLIB
O2SCL_EOS_MVAR = 1
else
O2SCL_EOS_MVAR = 0
endif

if O2SCL_HDF
O2SCL_HDF_MVAR = 1
else
O2SCL_HDF_MVAR = 0
endif

if O2SCL_OPENMP
if O2SCL_EIGEN
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_ARMA -DO2SCL_EIGEN \
	-fopenmp -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_EIGEN -fopenmp \
	-DO2SCL_COND_FLAG
endif # ARMA
else
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -DO2SCL_ARMA -fopenmp \
	-DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=1 -fopenmp \
	-DO2SCL_COND_FLAG
endif # ARMA
endif # EIGEN
else
if O2SCL_EIGEN
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_ARMA -DO2SCL_EIGEN -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_EIGEN -DO2SCL_COND_FLAG
endif # ARMA
else
if O2SCL_ARMA
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_ARMA -DO2SCL_COND_FLAG
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_PART_SVAR=$(O2SCL_PART_MVAR) \
	-DO2SCL_EOS_SVAR=$(O2SCL_EOS_MVAR) \
	-DO2SCL_HDF_SVAR=$(O2SCL_HDF_MVAR) \
	-DO2SCL_OPENMP_SVAR=0 -DO2SCL_COND_FLAG
endif # ARMA
endif # EIGEN
endif # OPENMP

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = min_cern_ts min_brent_gsl_ts mmin_simp_ts \
	mmin_conf_ts mmin_conp_ts mmin_bfgs2_ts \
	mmin_fix_ts mmin_constr_pgrad_ts mmin_constr_spg_ts \
	min_ts mmin_simp2_ts min_quad_golden_ts diff_evo_ts \
	diff_evo_adapt_ts min_brent_boost_ts

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl.la $(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl.la
endif

min_cern_ts_LDADD = $(VCHECK_LIBS)
min_brent_gsl_ts_LDADD = $(VCHECK_LIBS)
min_brent_boost_ts_LDADD = $(VCHECK_LIBS)
min_quad_golden_ts_LDADD = $(VCHECK_LIBS)
mmin_simp_ts_LDADD = $(VCHECK_LIBS)
mmin_simp2_ts_LDADD = $(VCHECK_LIBS)
mmin_conf_ts_LDADD = $(VCHECK_LIBS)
mmin_conp_ts_LDADD = $(VCHECK_LIBS)
mmin_bfgs2_ts_LDADD = $(VCHECK_LIBS)
mmin_fix_ts_LDADD = $(VCHECK_LIBS)
mmin_constr_pgrad_ts_LDADD = $(VCHECK_LIBS)
mmin_constr_spg_ts_LDADD = $(VCHECK_LIBS)
diff_evo_ts_LDADD = $(VCHECK_LIBS)
diff_evo_adapt_ts_LDADD = $(VCHECK_LIBS)
min_ts_LDADD = $(VCHECK_LIBS)

min_cern.scr: min_cern_ts$(EXEEXT) 
	./min_cern_ts$(EXEEXT) > min_cern.scr
min_brent_gsl.scr: min_brent_gsl_ts$(EXEEXT) 
	./min_brent_gsl_ts$(EXEEXT) > min_brent_gsl.scr
min_brent_boost.scr: min_brent_boost_ts$(EXEEXT) 
	./min_brent_boost_ts$(EXEEXT) > min_brent_boost.scr
min_quad_golden.scr: min_quad_golden_ts$(EXEEXT) 
	./min_quad_golden_ts$(EXEEXT) > min_quad_golden.scr
mmin_simp.scr: mmin_simp_ts$(EXEEXT) 
	./mmin_simp_ts$(EXEEXT) > mmin_simp.scr
mmin_simp2.scr: mmin_simp2_ts$(EXEEXT) 
	./mmin_simp2_ts$(EXEEXT) > mmin_simp2.scr
mmin_conf.scr: mmin_conf_ts$(EXEEXT) 
	./mmin_conf_ts$(EXEEXT) > mmin_conf.scr
mmin_conp.scr: mmin_conp_ts$(EXEEXT) 
	./mmin_conp_ts$(EXEEXT) > mmin_conp.scr
mmin_bfgs2.scr: mmin_bfgs2_ts$(EXEEXT) 
	./mmin_bfgs2_ts$(EXEEXT) > mmin_bfgs2.scr
mmin_fix.scr: mmin_fix_ts$(EXEEXT) 
	./mmin_fix_ts$(EXEEXT) > mmin_fix.scr
mmin_constr_pgrad.scr: mmin_constr_pgrad_ts$(EXEEXT) 
	./mmin_constr_pgrad_ts$(EXEEXT) > mmin_constr_pgrad.scr
mmin_constr_spg.scr: mmin_constr_spg_ts$(EXEEXT) 
	./mmin_constr_spg_ts$(EXEEXT) > mmin_constr_spg.scr
diff_evo.scr: diff_evo_ts$(EXEEXT) 
	./diff_evo_ts$(EXEEXT) > diff_evo.scr
diff_evo_adapt.scr: diff_evo_adapt_ts$(EXEEXT) 
	./diff_evo_adapt_ts$(EXEEXT) > diff_evo_adapt.scr
min.scr: min_ts$(EXEEXT) 
	./min_ts$(EXEEXT) > min.scr

min_cern_ts_SOURCES = min_cern_ts.cpp
min_brent_gsl_ts_SOURCES = min_brent_gsl_ts.cpp
min_brent_boost_ts_SOURCES = min_brent_boost_ts.cpp
min_quad_golden_ts_SOURCES = min_quad_golden_ts.cpp
mmin_simp_ts_SOURCES = mmin_simp_ts.cpp
mmin_simp2_ts_SOURCES = mmin_simp2_ts.cpp
mmin_conf_ts_SOURCES = mmin_conf_ts.cpp
mmin_conp_ts_SOURCES = mmin_conp_ts.cpp
mmin_bfgs2_ts_SOURCES = mmin_bfgs2_ts.cpp
mmin_fix_ts_SOURCES = mmin_fix_ts.cpp
mmin_constr_pgrad_ts_SOURCES = mmin_constr_pgrad_ts.cpp
mmin_constr_spg_ts_SOURCES = mmin_constr_spg_ts.cpp
diff_evo_ts_SOURCES = diff_evo_ts.cpp
diff_evo_adapt_ts_SOURCES = diff_evo_adapt_ts.cpp
min_ts_SOURCES = min_ts.cpp

# ------------------------------------------------------------
# Library o2scl_min
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_min.la

libo2scl_min_la_SOURCES = $(MIN_SRCS)

am__libo2scl_min_la_SOURCES_DIST = $(MIN_SRCS)

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
	echo "src/min" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

