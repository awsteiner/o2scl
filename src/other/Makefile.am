# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

OTHER_SRCS = series_acc.cpp poly.cpp polylog.cpp \
	interp2_eqi.cpp pinside.cpp bin_size.cpp \
	contour.cpp smooth_gsl.cpp hist.cpp \
	hist_2d.cpp prob_dens_func.cpp

HEADER_VAR = contour.h cheb_approx.h \
	series_acc.h interp2_planar.h poly.h polylog.h \
	interp2_direct.h interp2_eqi.h pinside.h bin_size.h \
	vec_stats.h smooth_gsl.h hist.h \
	hist_2d.h prob_dens_func.h interp2_seq.h interp2_neigh.h \
	interpm_idw.h interp2.h

TEST_VAR = series_acc.scr interp2_planar.scr contour.scr \
	poly.scr polylog.scr cheb_approx.scr vec_stats.scr smooth_gsl.scr \
	hist.scr hist_2d.scr prob_dens_func.scr interp2_direct.scr \
	pinside.scr interp2_seq.scr interp2_neigh.scr \
	interpm_idw.scr

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

check_PROGRAMS = series_acc_ts interp2_planar_ts \
	contour_ts poly_ts polylog_ts pinside_ts \
	interp2_direct_ts vec_stats_ts cheb_approx_ts \
	smooth_gsl_ts hist_ts hist_2d_ts interp2_seq_ts \
	prob_dens_func_ts interp2_neigh_ts \
	interpm_idw_ts

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

smooth_gsl_ts_LDADD = $(VCHECK_LIBS)
series_acc_ts_LDADD = $(VCHECK_LIBS)
interp2_planar_ts_LDADD = $(VCHECK_LIBS)
interp2_neigh_ts_LDADD = $(VCHECK_LIBS)
interpm_idw_ts_LDADD = $(VCHECK_LIBS)
contour_ts_LDADD = $(VCHECK_LIBS)
poly_ts_LDADD = $(VCHECK_LIBS)
polylog_ts_LDADD = $(VCHECK_LIBS)
interp2_seq_ts_LDADD = $(VCHECK_LIBS)
interp2_direct_ts_LDADD = $(VCHECK_LIBS)
cheb_approx_ts_LDADD = $(VCHECK_LIBS)
pinside_ts_LDADD = $(VCHECK_LIBS)
hist_ts_LDADD = $(VCHECK_LIBS)
hist_2d_ts_LDADD = $(VCHECK_LIBS)
prob_dens_func_ts_LDADD = $(VCHECK_LIBS)
vec_stats_ts_LDADD = $(VCHECK_LIBS)

smooth_gsl.scr: smooth_gsl_ts$(EXEEXT) 
	./smooth_gsl_ts$(EXEEXT) > smooth_gsl.scr
series_acc.scr: series_acc_ts$(EXEEXT) 
	./series_acc_ts$(EXEEXT) > series_acc.scr
interp2_planar.scr: interp2_planar_ts$(EXEEXT) 
	./interp2_planar_ts$(EXEEXT) > interp2_planar.scr
interp2_neigh.scr: interp2_neigh_ts$(EXEEXT) 
	./interp2_neigh_ts$(EXEEXT) > interp2_neigh.scr
interpm_idw.scr: interpm_idw_ts$(EXEEXT) 
	./interpm_idw_ts$(EXEEXT) > interpm_idw.scr
contour.scr: contour_ts$(EXEEXT) 
	./contour_ts$(EXEEXT) > contour.scr
poly.scr: poly_ts$(EXEEXT) 
	./poly_ts$(EXEEXT) > poly.scr
polylog.scr: polylog_ts$(EXEEXT) 
	./polylog_ts$(EXEEXT) > polylog.scr
interp2_seq.scr: interp2_seq_ts$(EXEEXT) 
	./interp2_seq_ts$(EXEEXT) > interp2_seq.scr
interp2_direct.scr: interp2_direct_ts$(EXEEXT) 
	./interp2_direct_ts$(EXEEXT) > interp2_direct.scr
cheb_approx.scr: cheb_approx_ts$(EXEEXT) 
	./cheb_approx_ts$(EXEEXT) > cheb_approx.scr
pinside.scr: pinside_ts$(EXEEXT) 
	./pinside_ts$(EXEEXT) > pinside.scr
hist.scr: hist_ts$(EXEEXT) 
	./hist_ts$(EXEEXT) > hist.scr
hist_2d.scr: hist_2d_ts$(EXEEXT) 
	./hist_2d_ts$(EXEEXT) > hist_2d.scr
prob_dens_func.scr: prob_dens_func_ts$(EXEEXT) 
	./prob_dens_func_ts$(EXEEXT) > prob_dens_func.scr
vec_stats.scr: vec_stats_ts$(EXEEXT) 
	./vec_stats_ts$(EXEEXT) > vec_stats.scr

cheb_approx_ts_SOURCES = cheb_approx_ts.cpp
contour_ts_SOURCES = contour_ts.cpp
series_acc_ts_SOURCES = series_acc_ts.cpp
hist_ts_SOURCES = hist_ts.cpp
hist_2d_ts_SOURCES = hist_2d_ts.cpp
interp2_neigh_ts_SOURCES = interp2_neigh_ts.cpp
interpm_idw_ts_SOURCES = interpm_idw_ts.cpp
interp2_planar_ts_SOURCES = interp2_planar_ts.cpp
interp2_seq_ts_SOURCES = interp2_seq_ts.cpp
interp2_direct_ts_SOURCES = interp2_direct_ts.cpp
pinside_ts_SOURCES = pinside_ts.cpp
poly_ts_SOURCES = poly_ts.cpp
polylog_ts_SOURCES = polylog_ts.cpp
prob_dens_func_ts_SOURCES = prob_dens_func_ts.cpp
smooth_gsl_ts_SOURCES = smooth_gsl_ts.cpp
vec_stats_ts_SOURCES = vec_stats_ts.cpp

# ------------------------------------------------------------
# Library o2scl_other
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_other.la

libo2scl_other_la_SOURCES = $(OTHER_SRCS)

am__libo2scl_other_la_SOURCES_DIST = $(OTHER_SRCS)

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
	echo "src/other" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

