# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

ODE_SRCS = 

HEADER_VAR = astep.h astep_gsl.h astep_nonadapt.h ode_rkck_gsl.h \
	ode_funct.h ode_step.h ode_bv_solve.h ode_iv_solve.h \
	ode_rk8pd_gsl.h ode_it_solve.h ode_bv_multishoot.h \
	ode_jac_funct.h ode_bsimp_gsl.h ode_rkf45_gsl.h \
	ode_iv_table.h ode_bv_mshoot.h

TEST_VAR = astep_gsl.scr ode_rkck_gsl.scr astep_nonadapt.scr \
	ode_rk8pd_gsl.scr ode_bsimp_gsl.scr \
	ode_rkf45_gsl.scr ode_iv_solve.scr ode_it_solve.scr

# ode_bv_mshoot.scr ode_iv_table.scr ode_bv_solve.scr 

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

check_PROGRAMS = astep_gsl_ts ode_rkck_gsl_ts astep_nonadapt_ts \
	ode_iv_solve_ts ode_bv_solve_ts ode_rk8pd_gsl_ts ode_it_solve_ts \
	ode_bsimp_gsl_ts ode_rkf45_gsl_ts ode_iv_table_ts ode_bv_mshoot_ts 

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl.la $(libdir)/libo2scl_hdf.la
else
VCHECK_LIBS = $(libdir)/libo2scl.la
endif

astep_gsl_ts_LDADD = $(VCHECK_LIBS)
ode_rkf45_gsl_ts_LDADD = $(VCHECK_LIBS)
ode_rkck_gsl_ts_LDADD = $(VCHECK_LIBS)
ode_bsimp_gsl_ts_LDADD = $(VCHECK_LIBS)
astep_nonadapt_ts_LDADD = $(VCHECK_LIBS)
ode_iv_solve_ts_LDADD = $(VCHECK_LIBS)
ode_iv_table_ts_LDADD = $(VCHECK_LIBS)
ode_bv_mshoot_ts_LDADD = $(VCHECK_LIBS)
ode_it_solve_ts_LDADD = $(VCHECK_LIBS)
ode_bv_solve_ts_LDADD = $(VCHECK_LIBS)
ode_rk8pd_gsl_ts_LDADD = $(VCHECK_LIBS)

astep_gsl.scr: astep_gsl_ts$(EXEEXT)
	./astep_gsl_ts$(EXEEXT) > astep_gsl.scr
ode_rkck_gsl.scr: ode_rkck_gsl_ts$(EXEEXT)
	./ode_rkck_gsl_ts$(EXEEXT) > ode_rkck_gsl.scr
ode_rkf45_gsl.scr: ode_rkf45_gsl_ts$(EXEEXT)
	./ode_rkf45_gsl_ts$(EXEEXT) > ode_rkf45_gsl.scr
ode_bsimp_gsl.scr: ode_bsimp_gsl_ts$(EXEEXT)
	./ode_bsimp_gsl_ts$(EXEEXT) > ode_bsimp_gsl.scr
ode_rk8pd_gsl.scr: ode_rk8pd_gsl_ts$(EXEEXT)
	./ode_rk8pd_gsl_ts$(EXEEXT) > ode_rk8pd_gsl.scr
astep_nonadapt.scr: astep_nonadapt_ts$(EXEEXT)
	./astep_nonadapt_ts$(EXEEXT) > astep_nonadapt.scr
ode_iv_solve.scr: ode_iv_solve_ts$(EXEEXT)
	./ode_iv_solve_ts$(EXEEXT) > ode_iv_solve.scr
ode_iv_table.scr: ode_iv_table_ts$(EXEEXT)
	./ode_iv_table_ts$(EXEEXT) > ode_iv_table.scr
ode_bv_mshoot.scr: ode_bv_mshoot_ts$(EXEEXT)
	./ode_bv_mshoot_ts$(EXEEXT) > ode_bv_mshoot.scr
ode_it_solve.scr: ode_it_solve_ts$(EXEEXT)
	./ode_it_solve_ts$(EXEEXT) > ode_it_solve.scr
ode_bv_solve.scr: ode_bv_solve_ts$(EXEEXT)
	./ode_bv_solve_ts$(EXEEXT) > ode_bv_solve.scr

astep_gsl_ts_SOURCES = astep_gsl_ts.cpp
astep_nonadapt_ts_SOURCES = astep_nonadapt_ts.cpp
ode_rkck_gsl_ts_SOURCES = ode_rkck_gsl_ts.cpp
ode_rkf45_gsl_ts_SOURCES = ode_rkf45_gsl_ts.cpp
ode_bsimp_gsl_ts_SOURCES = ode_bsimp_gsl_ts.cpp
ode_iv_solve_ts_SOURCES = ode_iv_solve_ts.cpp
ode_iv_table_ts_SOURCES = ode_iv_table_ts.cpp
ode_bv_mshoot_ts_SOURCES = ode_bv_mshoot_ts.cpp
ode_it_solve_ts_SOURCES = ode_it_solve_ts.cpp
ode_bv_solve_ts_SOURCES = ode_bv_solve_ts.cpp
ode_rk8pd_gsl_ts_SOURCES = ode_rk8pd_gsl_ts.cpp

# ------------------------------------------------------------
# No library o2scl_ode
# ------------------------------------------------------------

#noinst_LTLIBRARIES = libo2scl_ode.la

#libo2scl_ode_la_SOURCES = $(ODE_SRCS)

#am__libo2scl_ode_la_SOURCES_DIST = $(ODE_SRCS)

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
	echo "src/ode" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

