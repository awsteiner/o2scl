# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

OTHER_SRCS = series_acc.cpp poly.cpp polylog.cpp \
	interp2_eqi.cpp pinside.cpp other_python.cpp \
	contour.cpp smooth_gsl.cpp hist.cpp xml.cpp \
	hist_2d.cpp prob_dens_func.cpp vec_stats.cpp \
	interpm_python.cpp gmm_python.cpp

HEADER_VAR = contour.h cheb_approx.h other_python.h \
	series_acc.h interp2_planar.h poly.h polylog.h \
	interp2_direct.h interp2_eqi.h pinside.h \
	vec_stats.h smooth_gsl.h hist.h smooth_func.h \
	hist_2d.h prob_dens_func.h interp2_seq.h interp2_neigh.h \
	interpm_idw.h interp2.h interpm_krige.h prob_dens_mdim_amr.h \
	slack_messenger.h xml.h exp_max.h interpm_python.h \
	gmm_python.h kde_python.h

TEST_VAR = series_acc.scr interp2_planar.scr contour.scr \
	poly.scr polylog.scr cheb_approx.scr vec_stats.scr smooth_gsl.scr \
	hist.scr hist_2d.scr prob_dens_func.scr interp2_direct.scr \
	pinside.scr interp2_seq.scr interp2_neigh.scr \
	interpm_idw.scr interpm_krige.scr smooth_func.scr \
	prob_dens_mdim_amr.scr xml.scr exp_max.scr interpm_python.scr \
	gmm_python.scr kde_python.scr

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

if O2SCL_ARMA
O2SCL_ARMA_MVAR = -DO2SCL_ARMA
else
O2SCL_ARMA_MVAR =
endif

if O2SCL_LINUX
O2SCL_LINUX_MVAR = -DO2SCL_LINUX
else 
O2SCL_LINUX_MVAR = 
endif

if O2SCL_OSX
O2SCL_OSX_MVAR = -DO2SCL_OSX
else 
O2SCL_OSX_MVAR = 
endif

if O2SCL_MPFR
O2SCL_MPFR_MVAR = -DO2SCL_MPFR
else 
O2SCL_MPFR_MVAR = 
endif

if O2SCL_FFTW
O2SCL_FFTW_MVAR = -DO2SCL_FFTW
else 
O2SCL_FFTW_MVAR = 
endif

if O2SCL_PUGIXML
O2SCL_PUGIXML_MVAR = -DO2SCL_PUGIXML
else 
O2SCL_PUGIXML_MVAR = 
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_CUBATURE
O2SCL_CUBATURE_MVAR = -DO2SCL_CUBATURE
else 
O2SCL_CUBATURE_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_MPFR_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR) \
	$(O2SCL_PYTHON_MVAR) \
	$(O2SCL_CUBATURE_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_MPFR_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR) \
	$(O2SCL_PYTHON_MVAR) \
	$(O2SCL_CUBATURE_MVAR)
endif

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
	interpm_idw_ts interpm_krige_ts smooth_func_ts \
	prob_dens_mdim_amr_ts xml_ts exp_max_ts interpm_python_ts \
	gmm_python_ts kde_python_ts

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

smooth_gsl_ts_LDADD = $(ADDL_TEST_LIBS)
exp_max_ts_LDADD = $(ADDL_TEST_LIBS)
series_acc_ts_LDADD = $(ADDL_TEST_LIBS)
interp2_planar_ts_LDADD = $(ADDL_TEST_LIBS)
interp2_neigh_ts_LDADD = $(ADDL_TEST_LIBS)
interpm_idw_ts_LDADD = $(ADDL_TEST_LIBS)
interpm_python_ts_LDADD = $(ADDL_TEST_LIBS)
kde_python_ts_LDADD = $(ADDL_TEST_LIBS)
gmm_python_ts_LDADD = $(ADDL_TEST_LIBS)
smooth_func_ts_LDADD = $(ADDL_TEST_LIBS)
interpm_krige_ts_LDADD = $(ADDL_TEST_LIBS)
contour_ts_LDADD = $(ADDL_TEST_LIBS)
poly_ts_LDADD = $(ADDL_TEST_LIBS)
polylog_ts_LDADD = $(ADDL_TEST_LIBS)
interp2_seq_ts_LDADD = $(ADDL_TEST_LIBS)
prob_dens_mdim_amr_ts_LDADD = $(ADDL_TEST_LIBS)
interp2_direct_ts_LDADD = $(ADDL_TEST_LIBS)
cheb_approx_ts_LDADD = $(ADDL_TEST_LIBS)
pinside_ts_LDADD = $(ADDL_TEST_LIBS)
hist_ts_LDADD = $(ADDL_TEST_LIBS)
hist_2d_ts_LDADD = $(ADDL_TEST_LIBS)
prob_dens_func_ts_LDADD = $(ADDL_TEST_LIBS)
vec_stats_ts_LDADD = $(ADDL_TEST_LIBS)
xml_ts_LDADD = $(ADDL_TEST_LIBS)

smooth_gsl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
exp_max_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
series_acc_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interp2_planar_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interp2_neigh_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interpm_idw_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interpm_python_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
gmm_python_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
kde_python_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
smooth_func_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interpm_krige_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
contour_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
poly_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
polylog_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interp2_seq_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
prob_dens_mdim_amr_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interp2_direct_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
cheb_approx_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
pinside_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
hist_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
hist_2d_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
prob_dens_func_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
vec_stats_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
xml_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

smooth_gsl.scr: smooth_gsl_ts$(EXEEXT) 
	./smooth_gsl_ts$(EXEEXT) > smooth_gsl.scr
exp_max.scr: exp_max_ts$(EXEEXT) 
	./exp_max_ts$(EXEEXT) > exp_max.scr
series_acc.scr: series_acc_ts$(EXEEXT) 
	./series_acc_ts$(EXEEXT) > series_acc.scr
interp2_planar.scr: interp2_planar_ts$(EXEEXT) 
	./interp2_planar_ts$(EXEEXT) > interp2_planar.scr
interp2_neigh.scr: interp2_neigh_ts$(EXEEXT) 
	./interp2_neigh_ts$(EXEEXT) > interp2_neigh.scr
interpm_idw.scr: interpm_idw_ts$(EXEEXT) 
	./interpm_idw_ts$(EXEEXT) > interpm_idw.scr
interpm_python.scr: interpm_python_ts$(EXEEXT) 
	./interpm_python_ts$(EXEEXT) > interpm_python.scr
gmm_python.scr: gmm_python_ts$(EXEEXT) 
	./gmm_python_ts$(EXEEXT) > gmm_python.scr
kde_python.scr: kde_python_ts$(EXEEXT) 
	./kde_python_ts$(EXEEXT) > kde_python.scr
smooth_func.scr: smooth_func_ts$(EXEEXT) 
	./smooth_func_ts$(EXEEXT) > smooth_func.scr
interpm_krige.scr: interpm_krige_ts$(EXEEXT) 
	./interpm_krige_ts$(EXEEXT) > interpm_krige.scr
contour.scr: contour_ts$(EXEEXT) 
	./contour_ts$(EXEEXT) > contour.scr
poly.scr: poly_ts$(EXEEXT) 
	./poly_ts$(EXEEXT) > poly.scr
polylog.scr: polylog_ts$(EXEEXT) 
	./polylog_ts$(EXEEXT) > polylog.scr
interp2_seq.scr: interp2_seq_ts$(EXEEXT) 
	./interp2_seq_ts$(EXEEXT) > interp2_seq.scr
prob_dens_mdim_amr.scr: prob_dens_mdim_amr_ts$(EXEEXT) 
	./prob_dens_mdim_amr_ts$(EXEEXT) > prob_dens_mdim_amr.scr
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
xml.scr: xml_ts$(EXEEXT) 
	./xml_ts$(EXEEXT) > xml.scr

cheb_approx_ts_SOURCES = cheb_approx_ts.cpp
contour_ts_SOURCES = contour_ts.cpp
series_acc_ts_SOURCES = series_acc_ts.cpp
hist_ts_SOURCES = hist_ts.cpp
hist_2d_ts_SOURCES = hist_2d_ts.cpp
interp2_neigh_ts_SOURCES = interp2_neigh_ts.cpp
interpm_idw_ts_SOURCES = interpm_idw_ts.cpp
interpm_python_ts_SOURCES = interpm_python_ts.cpp
gmm_python_ts_SOURCES = gmm_python_ts.cpp
kde_python_ts_SOURCES = kde_python_ts.cpp
smooth_func_ts_SOURCES = smooth_func_ts.cpp
interpm_krige_ts_SOURCES = interpm_krige_ts.cpp
interp2_planar_ts_SOURCES = interp2_planar_ts.cpp
interp2_seq_ts_SOURCES = interp2_seq_ts.cpp
prob_dens_mdim_amr_ts_SOURCES = prob_dens_mdim_amr_ts.cpp
interp2_direct_ts_SOURCES = interp2_direct_ts.cpp
pinside_ts_SOURCES = pinside_ts.cpp
poly_ts_SOURCES = poly_ts.cpp
polylog_ts_SOURCES = polylog_ts.cpp
prob_dens_func_ts_SOURCES = prob_dens_func_ts.cpp
smooth_gsl_ts_SOURCES = smooth_gsl_ts.cpp
exp_max_ts_SOURCES = exp_max_ts.cpp
vec_stats_ts_SOURCES = vec_stats_ts.cpp
xml_ts_SOURCES = xml_ts.cpp

# ------------------------------------------------------------
# Library o2scl_other
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_other.la

libo2scl_other_la_SOURCES = $(OTHER_SRCS)

am__libo2scl_other_la_SOURCES_DIST = $(OTHER_SRCS)

# ------------------------------------------------------------

interface: empty
	../../bin/yanic other.i other_python ~/wcs/o2sclpy/o2sclpy/other \
		~/wcs/o2sclpy/doc/other ../../bin/yanic_header.txt \
		| tee other_python.yanic.out

empty:

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

