# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

# 5/19/11, AWS: Even though currently O2scl_eos requires HDF support,
# we leave this HDF conditional here, in case future versions of
# the EOS library will allow one to disable HDF support
if O2SCL_HDF

EOS_SRCS = apr_eos.cpp bag_eos.cpp \
	bps_eos.cpp cold_nstar.cpp eos.cpp gen_potential_eos.cpp \
	nambujl_eos.cpp ddc_eos.cpp ldrop_mass.cpp \
	quark_eos.cpp rmf_eos.cpp schematic_eos.cpp \
	skyrme_eos.cpp tov_solve.cpp nse_eos.cpp hdf_eos_io.cpp \
	rmf_delta_eos.cpp tov_eos.cpp sym4_eos.cpp hadronic_eos.cpp \
	rmf_nucleus.cpp gen_sn_eos.cpp \
	ldrop_shell.cpp cfl_njl_eos.cpp cfl6_eos.cpp 

HEADER_VAR = apr_eos.h bag_eos.h bps_eos.h ddc_eos.h \
	cold_nstar.h eos.h gen_potential_eos.h ldrop_mass.h \
	nambujl_eos.h quark_eos.h rmf_eos.h cfl_njl_eos.h cfl6_eos.h \
	schematic_eos.h skyrme_eos.h tabulated_eos.h nse_eos.h \
	tov_solve.h rmf_delta_eos.h sym4_eos.h tov_eos.h hadronic_eos.h \
	hdf_eos_io.h rmf_nucleus.h gen_sn_eos.h \
	ldrop_shell.h gogny_eos.h virial_eos.h

TEST_VAR = apr_eos.scr bag_eos.scr \
	cold_nstar.scr eos.scr gen_potential_eos.scr sym4_eos.scr \
	nambujl_eos.scr quark_eos.scr ldrop_mass.scr \
	rmf_eos.scr schematic_eos.scr skyrme_eos.scr tov_eos.scr \
	tabulated_eos.scr nse_eos.scr \
	rmf_delta_eos.scr bps_eos.scr ddc_eos.scr \
	hadronic_eos.scr rmf_nucleus.scr \
	virial_eos.scr gogny_eos.scr tov_solve.scr cfl6_eos.scr \
	cfl_njl_eos.scr ldrop_shell.scr 
# gen_sn_eos.scr 

else

EOS_SRCS = apr_eos.cpp bag_eos.cpp cfl_njl_eos.cpp cfl6_eos.cpp \
	bps_eos.cpp cold_nstar.cpp eos.cpp gen_potential_eos.cpp \
	nambujl_eos.cpp ddc_eos.cpp ldrop_mass.cpp \
	quark_eos.cpp rmf_eos.cpp schematic_eos.cpp \
	skyrme_eos.cpp tov_solve.cpp nse_eos.cpp \
	rmf_delta_eos.cpp tov_eos.cpp sym4_eos.cpp hadronic_eos.cpp \
	gen_sn_eos.cpp ldrop_shell.cpp

HEADER_VAR = apr_eos.h bag_eos.h bps_eos.h ddc_eos.h \
	cold_nstar.h eos.h gen_potential_eos.h ldrop_mass.h \
	nambujl_eos.h quark_eos.h rmf_eos.h cfl_njl_eos.h cfl6_eos.h \
	schematic_eos.h skyrme_eos.h tabulated_eos.h nse_eos.h \
	tov_solve.h rmf_delta_eos.h sym4_eos.h tov_eos.h hadronic_eos.h \
	gen_sn_eos.h ldrop_shell.h gogny_eos.h virial_eos.h

TEST_VAR = apr_eos.scr bag_eos.scr virial_eos.scr \
	cold_nstar.scr eos.scr gen_potential_eos.scr sym4_eos.scr \
	nambujl_eos.scr quark_eos.scr gogny_eos.scr \
	rmf_eos.scr schematic_eos.scr skyrme_eos.scr tov_eos.scr \
	tov_solve.scr cfl6_eos.scr tabulated_eos.scr nse_eos.scr \
	rmf_delta_eos.scr bps_eos.scr ddc_eos.scr cfl_njl_eos.scr \
	hadronic_eos.scr gen_sn_eos.scr ldrop_shell.scr 

endif

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

check_PROGRAMS = apr_eos_ts bag_eos_ts ldrop_mass_ts \
	cold_nstar_ts eos_ts gen_potential_eos_ts sym4_eos_ts \
	hadronic_eos_ts nambujl_eos_ts quark_eos_ts \
	rmf_eos_ts schematic_eos_ts skyrme_eos_ts tov_eos_ts \
	tov_solve_ts cfl6_eos_ts tabulated_eos_ts nse_eos_ts \
	rmf_delta_eos_ts bps_eos_ts ddc_eos_ts cfl_njl_eos_ts \
	rmf_nucleus_ts \
	ldrop_shell_ts gogny_eos_ts virial_eos_ts
# gen_sn_eos_ts 

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl_eos.la $(libdir)/libo2scl_part.la \
	$(libdir)/libo2scl_hdf.la $(libdir)/libo2scl.la 
else
VCHECK_LIBS = $(libdir)/libo2scl_eos.la $(libdir)/libo2scl_part.la \
	$(libdir)/libo2scl.la 
endif

apr_eos_ts_LDADD = $(VCHECK_LIBS)
bag_eos_ts_LDADD = $(VCHECK_LIBS)
cold_nstar_ts_LDADD = $(VCHECK_LIBS)
ldrop_mass_ts_LDADD = $(VCHECK_LIBS)
ldrop_shell_ts_LDADD = $(VCHECK_LIBS)
eos_ts_LDADD = $(VCHECK_LIBS)
gen_potential_eos_ts_LDADD = $(VCHECK_LIBS)
sym4_eos_ts_LDADD = $(VCHECK_LIBS)
hadronic_eos_ts_LDADD = $(VCHECK_LIBS)
nambujl_eos_ts_LDADD = $(VCHECK_LIBS)
quark_eos_ts_LDADD = $(VCHECK_LIBS)
rmf_eos_ts_LDADD = $(VCHECK_LIBS)
schematic_eos_ts_LDADD = $(VCHECK_LIBS)
skyrme_eos_ts_LDADD = $(VCHECK_LIBS)
tov_eos_ts_LDADD = $(VCHECK_LIBS)
tov_solve_ts_LDADD = $(VCHECK_LIBS)
cfl6_eos_ts_LDADD = $(VCHECK_LIBS)
tabulated_eos_ts_LDADD = $(VCHECK_LIBS)
rmf_delta_eos_ts_LDADD = $(VCHECK_LIBS)
bps_eos_ts_LDADD = $(VCHECK_LIBS)
ddc_eos_ts_LDADD = $(VCHECK_LIBS)
cfl_njl_eos_ts_LDADD = $(VCHECK_LIBS)
nse_eos_ts_LDADD = $(VCHECK_LIBS)
# gen_sn_eos_ts_LDADD = $(VCHECK_LIBS)
rmf_nucleus_ts_LDADD = $(VCHECK_LIBS)
gogny_eos_ts_LDADD = $(VCHECK_LIBS)
virial_eos_ts_LDADD = $(VCHECK_LIBS)

apr_eos.scr: apr_eos_ts$(EXEEXT) 
	./apr_eos_ts$(EXEEXT) > apr_eos.scr
bag_eos.scr: bag_eos_ts$(EXEEXT) 
	./bag_eos_ts$(EXEEXT) > bag_eos.scr
cold_nstar.scr: cold_nstar_ts$(EXEEXT) 
	./cold_nstar_ts$(EXEEXT) > cold_nstar.scr
ldrop_mass.scr: ldrop_mass_ts$(EXEEXT) 
	./ldrop_mass_ts$(EXEEXT) > ldrop_mass.scr
ldrop_shell.scr: ldrop_shell_ts$(EXEEXT) 
	./ldrop_shell_ts$(EXEEXT) > ldrop_shell.scr
eos.scr: eos_ts$(EXEEXT) 
	./eos_ts$(EXEEXT) > eos.scr
gen_potential_eos.scr: gen_potential_eos_ts$(EXEEXT) 
	./gen_potential_eos_ts$(EXEEXT) > gen_potential_eos.scr
sym4_eos.scr: sym4_eos_ts$(EXEEXT) 
	./sym4_eos_ts$(EXEEXT) > sym4_eos.scr
hadronic_eos.scr: hadronic_eos_ts$(EXEEXT) 
	./hadronic_eos_ts$(EXEEXT) > hadronic_eos.scr
nambujl_eos.scr: nambujl_eos_ts$(EXEEXT) 
	./nambujl_eos_ts$(EXEEXT) > nambujl_eos.scr
quark_eos.scr: quark_eos_ts$(EXEEXT) 
	./quark_eos_ts$(EXEEXT) > quark_eos.scr
rmf_eos.scr: rmf_eos_ts$(EXEEXT) 
	./rmf_eos_ts$(EXEEXT) > rmf_eos.scr
schematic_eos.scr: schematic_eos_ts$(EXEEXT) 
	./schematic_eos_ts$(EXEEXT) > schematic_eos.scr
skyrme_eos.scr: skyrme_eos_ts$(EXEEXT) 
	./skyrme_eos_ts$(EXEEXT) > skyrme_eos.scr
tov_eos.scr: tov_eos_ts$(EXEEXT) 
	./tov_eos_ts$(EXEEXT) > tov_eos.scr
tov_solve.scr: tov_solve_ts$(EXEEXT) 
	./tov_solve_ts$(EXEEXT) > tov_solve.scr
cfl6_eos.scr: cfl6_eos_ts$(EXEEXT) 
	./cfl6_eos_ts$(EXEEXT) > cfl6_eos.scr
tabulated_eos.scr: tabulated_eos_ts$(EXEEXT) 
	./tabulated_eos_ts$(EXEEXT) > tabulated_eos.scr
rmf_delta_eos.scr: rmf_delta_eos_ts$(EXEEXT) 
	./rmf_delta_eos_ts$(EXEEXT) > rmf_delta_eos.scr
bps_eos.scr: bps_eos_ts$(EXEEXT) 
	./bps_eos_ts$(EXEEXT) > bps_eos.scr
ddc_eos.scr: ddc_eos_ts$(EXEEXT) 
	./ddc_eos_ts$(EXEEXT) > ddc_eos.scr
cfl_njl_eos.scr: cfl_njl_eos_ts$(EXEEXT) 
	./cfl_njl_eos_ts$(EXEEXT) > cfl_njl_eos.scr
nse_eos.scr: nse_eos_ts$(EXEEXT) 
	./nse_eos_ts$(EXEEXT) > nse_eos.scr
# gen_sn_eos.scr: gen_sn_eos_ts$(EXEEXT) 
# 	./gen_sn_eos_ts$(EXEEXT) > gen_sn_eos.scr
rmf_nucleus.scr: rmf_nucleus_ts$(EXEEXT) 
	./rmf_nucleus_ts$(EXEEXT) > rmf_nucleus.scr
gogny_eos.scr: gogny_eos_ts$(EXEEXT) 
	./gogny_eos_ts$(EXEEXT) > gogny_eos.scr
virial_eos.scr: virial_eos_ts$(EXEEXT) 
	./virial_eos_ts$(EXEEXT) > virial_eos.scr

apr_eos_ts_SOURCES = apr_eos_ts.cpp
bag_eos_ts_SOURCES = bag_eos_ts.cpp
cold_nstar_ts_SOURCES = cold_nstar_ts.cpp
ldrop_mass_ts_SOURCES = ldrop_mass_ts.cpp
ldrop_shell_ts_SOURCES = ldrop_shell_ts.cpp
eos_ts_SOURCES = eos_ts.cpp
gen_potential_eos_ts_SOURCES = gen_potential_eos_ts.cpp
sym4_eos_ts_SOURCES = sym4_eos_ts.cpp
hadronic_eos_ts_SOURCES = hadronic_eos_ts.cpp
nambujl_eos_ts_SOURCES = nambujl_eos_ts.cpp
quark_eos_ts_SOURCES = quark_eos_ts.cpp
rmf_eos_ts_SOURCES = rmf_eos_ts.cpp
schematic_eos_ts_SOURCES = schematic_eos_ts.cpp
skyrme_eos_ts_SOURCES = skyrme_eos_ts.cpp
tov_eos_ts_SOURCES = tov_eos_ts.cpp
tov_solve_ts_SOURCES = tov_solve_ts.cpp
cfl6_eos_ts_SOURCES = cfl6_eos_ts.cpp
tabulated_eos_ts_SOURCES = tabulated_eos_ts.cpp
rmf_delta_eos_ts_SOURCES = rmf_delta_eos_ts.cpp
bps_eos_ts_SOURCES = bps_eos_ts.cpp
ddc_eos_ts_SOURCES = ddc_eos_ts.cpp
cfl_njl_eos_ts_SOURCES = cfl_njl_eos_ts.cpp
nse_eos_ts_SOURCES = nse_eos_ts.cpp
# gen_sn_eos_ts_SOURCES = gen_sn_eos_ts.cpp
rmf_nucleus_ts_SOURCES = rmf_nucleus_ts.cpp
gogny_eos_ts_SOURCES = gogny_eos_ts.cpp
virial_eos_ts_SOURCES = virial_eos_ts.cpp

ls_test: ls_ls_test ls_skm_test ls_ska_test ls_sk1_test

ls_ls_test:
	gen_sn_eos_ts ls $(HOME)/pkgs/ls_eos2/ls.dat

ls_skm_test:
	gen_sn_eos_ts ls $(HOME)/pkgs/ls_eos2/skm.dat

ls_ska_test:
	gen_sn_eos_ts ls $(HOME)/pkgs/ls_eos2/ska.dat

ls_sk1_test:
	gen_sn_eos_ts ls $(HOME)/pkgs/ls_eos2/sk1.dat

stos_test: stos_eos1_test stos_eos2_test stos_eos3_test

stos_eos1_test:
	gen_sn_eos_ts stos $(HOME)/pkgs/stos_eos/eos1

stos_eos2_test:
	gen_sn_eos_ts stos $(HOME)/pkgs/stos_eos/eos2

stos_eos3_test:
	gen_sn_eos_ts stos $(HOME)/pkgs/stos_eos/eos3

oo_test: oo_ls180_test oo_ls220_test oo_ls375_test

oo_ls180_test:
	gen_sn_eos_ts oo \
	$(HOME)/pkgs/oo_eos/LS180_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_ls220_test:
	gen_sn_eos_ts oo \
	$(HOME)/pkgs/oo_eos/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_ls375_test:
	gen_sn_eos_ts oo \
	$(HOME)/pkgs/oo_eos/LS375_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_stos_test:
	gen_sn_eos_ts oo \
	$(HOME)/pkgs/oo_eos/HShenEOS_rho220_temp180_ye65_version1.1_20120817.h5 \
	stos

oo_stos_hyp_test:
	gen_sn_eos_ts oo \
	$(HOME)/pkgs/oo_eos/HShen_HyperonEOS_rho220_temp180_ye65_version1.1_20131007.h5 \
	stos


x:
	gen_sn_eos_ts hfsl \
		$(HOME)/pkgs/hfsl_eos/dd2_frdm_eos_shen98format_v1.02.tab \
		 >> gen_sn_test.scr
	gen_sn_eos_ts hfsl \
		$(HOME)/pkgs/hfsl_eos/fsg_roca_eos_shen98format_v1.0.tab \
		 >> gen_sn_test.scr
	gen_sn_eos_ts hfsl \
		$(HOME)/pkgs/hfsl_eos/nl3_lala_eos_shen98format_v1.0.tab \
		 >> gen_sn_test.scr

# ------------------------------------------------------------
# Library o2scl_eos
# ------------------------------------------------------------

lib_LTLIBRARIES = libo2scl_eos.la

libo2scl_eos_la_SOURCES = $(EOS_SRCS)
libo2scl_eos_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

am__libo2scl_eos_la_SOURCES_DIST = $(EOS_SRCS)

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
	echo "src/eos" $(TEST_VAR) >> ../../testlist

test-clean: 
	rm -f *_ts.o *_ts *.scr

