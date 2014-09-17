# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

# 5/19/11, AWS: Even though currently O2scl_eos requires HDF support,
# we leave this HDF conditional here, in case future versions of
# the EOS library will allow one to disable HDF support
if O2SCL_HDF

EOS_SRCS = eos_had_apr.cpp eos_quark_bag.cpp \
	eos_crust.cpp nstar_cold.cpp eos_base.cpp eos_had_potential.cpp \
	eos_quark_njl.cpp eos_had_ddc.cpp nucmass_ldrop.cpp \
	eos_quark.cpp eos_had_rmf.cpp eos_had_schematic.cpp \
	eos_had_skyrme.cpp tov_solve.cpp eos_nse.cpp hdf_eos_io.cpp \
	eos_had_rmf_delta.cpp eos_tov.cpp eos_had_sym4.cpp \
	eos_had_base.cpp nucleus_rmf.cpp eos_sn.cpp \
	nucmass_ldrop_shell.cpp eos_quark_cfl.cpp eos_quark_cfl6.cpp \
	eos_had_hlps.cpp eos_nse_full.cpp

HEADER_VAR = eos_had_apr.h eos_quark_bag.h eos_crust.h eos_had_ddc.h \
	nstar_cold.h eos_base.h eos_had_potential.h nucmass_ldrop.h \
	eos_quark_njl.h eos_quark.h eos_had_rmf.h eos_quark_cfl.h \
	eos_quark_cfl6.h eos_had_schematic.h eos_had_skyrme.h \
	eos_had_tabulated.h eos_nse.h tov_solve.h eos_had_rmf_delta.h \
	eos_had_sym4.h eos_tov.h eos_had_base.h \
	hdf_eos_io.h nucleus_rmf.h eos_sn.h nucmass_ldrop_shell.h \
	eos_had_gogny.h eos_crust_virial.h eos_had_hlps.h \
	eos_nse_full.h

TEST_VAR = eos_had_apr.scr eos_quark_bag.scr nstar_cold.scr \
	eos_base.scr eos_had_potential.scr eos_had_sym4.scr \
	eos_quark_njl.scr eos_quark.scr nucmass_ldrop.scr \
	eos_had_rmf.scr eos_had_schematic.scr eos_had_skyrme.scr eos_tov.scr \
	eos_had_tabulated.scr eos_nse.scr eos_had_rmf_delta.scr \
	eos_crust.scr eos_had_ddc.scr eos_had_base.scr nucleus_rmf.scr \
	eos_crust_virial.scr eos_had_gogny.scr tov_solve.scr \
	eos_quark_cfl6.scr eos_quark_cfl.scr nucmass_ldrop_shell.scr \
	eos_nse_full.scr

else

EOS_SRCS = eos_had_apr.cpp eos_quark_bag.cpp eos_quark_cfl.cpp \
	eos_quark_cfl6.cpp \
	eos_crust.cpp nstar_cold.cpp eos_base.cpp eos_had_potential.cpp \
	eos_quark_njl.cpp eos_had_ddc.cpp nucmass_ldrop.cpp \
	eos_quark.cpp eos_had_rmf.cpp eos_had_schematic.cpp \
	eos_had_skyrme.cpp tov_solve.cpp eos_nse.cpp \
	eos_had_rmf_delta.cpp eos_tov.cpp eos_had_sym4.cpp \
	eos_had_base.cpp \
	eos_sn.cpp nucmass_ldrop_shell.cpp eos_had_hlps.cpp \
	eos_nse_full.cpp

HEADER_VAR = eos_had_apr.h eos_quark_bag.h eos_crust.h eos_had_ddc.h \
	nstar_cold.h eos_base.h eos_had_potential.h nucmass_ldrop.h \
	eos_quark_njl.h eos_quark.h eos_had_rmf.h eos_quark_cfl.h \
	eos_quark_cfl6.h eos_had_schematic.h eos_had_skyrme.h \
	eos_had_tabulated.h eos_nse.h tov_solve.h eos_had_rmf_delta.h \
	eos_had_sym4.h eos_tov.h eos_had_base.h \
	eos_sn.h nucmass_ldrop_shell.h eos_had_gogny.h \
	eos_crust_virial.h eos_had_hlps.h eos_nse_full.h

TEST_VAR = eos_had_apr.scr eos_quark_bag.scr eos_crust_virial.scr \
	nstar_cold.scr eos_base.scr eos_had_potential.scr \
	eos_had_sym4.scr \
	eos_quark_njl.scr eos_quark.scr eos_had_gogny.scr \
	eos_had_rmf.scr eos_had_schematic.scr eos_had_skyrme.scr eos_tov.scr \
	tov_solve.scr eos_quark_cfl6.scr eos_had_tabulated.scr eos_nse.scr \
	eos_had_rmf_delta.scr eos_crust.scr eos_had_ddc.scr eos_quark_cfl.scr \
	eos_had_base.scr eos_sn.scr nucmass_ldrop_shell.scr eos_nse_full.scr

endif

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

if O2SCL_CPP11
O2SCL_CPP11_MVAR = 
else
O2SCL_CPP11_MVAR = -DO2SCL_NO_CPP11
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_CPP11_MVAR) -DO2SCL_COND_FLAG

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = eos_had_apr_ts eos_quark_bag_ts nucmass_ldrop_ts \
	nstar_cold_ts eos_base_ts eos_had_potential_ts eos_had_sym4_ts \
	eos_had_base_ts eos_quark_njl_ts eos_quark_ts \
	eos_had_rmf_ts eos_had_schematic_ts eos_had_skyrme_ts eos_tov_ts \
	tov_solve_ts eos_quark_cfl6_ts eos_had_tabulated_ts eos_nse_ts \
	eos_had_rmf_delta_ts eos_crust_ts eos_had_ddc_ts eos_quark_cfl_ts \
	nucleus_rmf_ts eos_nse_full_ts \
	nucmass_ldrop_shell_ts eos_had_gogny_ts eos_crust_virial_ts
# eos_sn_ts 

if O2SCL_HDF
VCHECK_LIBS = $(libdir)/libo2scl_eos.la $(libdir)/libo2scl_part.la \
	$(libdir)/libo2scl_hdf.la $(libdir)/libo2scl.la 
else
VCHECK_LIBS = $(libdir)/libo2scl_eos.la $(libdir)/libo2scl_part.la \
	$(libdir)/libo2scl.la 
endif

eos_had_apr_ts_LDADD = $(VCHECK_LIBS)
eos_quark_bag_ts_LDADD = $(VCHECK_LIBS)
nstar_cold_ts_LDADD = $(VCHECK_LIBS)
nucmass_ldrop_ts_LDADD = $(VCHECK_LIBS)
nucmass_ldrop_shell_ts_LDADD = $(VCHECK_LIBS)
eos_base_ts_LDADD = $(VCHECK_LIBS)
eos_had_potential_ts_LDADD = $(VCHECK_LIBS)
eos_had_sym4_ts_LDADD = $(VCHECK_LIBS)
eos_had_base_ts_LDADD = $(VCHECK_LIBS)
eos_quark_njl_ts_LDADD = $(VCHECK_LIBS)
eos_quark_ts_LDADD = $(VCHECK_LIBS)
eos_had_rmf_ts_LDADD = $(VCHECK_LIBS)
eos_had_schematic_ts_LDADD = $(VCHECK_LIBS)
eos_had_skyrme_ts_LDADD = $(VCHECK_LIBS)
eos_tov_ts_LDADD = $(VCHECK_LIBS)
tov_solve_ts_LDADD = $(VCHECK_LIBS)
eos_quark_cfl6_ts_LDADD = $(VCHECK_LIBS)
eos_had_tabulated_ts_LDADD = $(VCHECK_LIBS)
eos_had_rmf_delta_ts_LDADD = $(VCHECK_LIBS)
eos_crust_ts_LDADD = $(VCHECK_LIBS)
eos_had_ddc_ts_LDADD = $(VCHECK_LIBS)
eos_quark_cfl_ts_LDADD = $(VCHECK_LIBS)
eos_nse_ts_LDADD = $(VCHECK_LIBS)
eos_nse_full_ts_LDADD = $(VCHECK_LIBS)
# eos_sn_ts_LDADD = $(VCHECK_LIBS)
nucleus_rmf_ts_LDADD = $(VCHECK_LIBS)
eos_had_gogny_ts_LDADD = $(VCHECK_LIBS)
eos_crust_virial_ts_LDADD = $(VCHECK_LIBS)

eos_had_apr.scr: eos_had_apr_ts$(EXEEXT) 
	./eos_had_apr_ts$(EXEEXT) > eos_had_apr.scr
eos_quark_bag.scr: eos_quark_bag_ts$(EXEEXT) 
	./eos_quark_bag_ts$(EXEEXT) > eos_quark_bag.scr
nstar_cold.scr: nstar_cold_ts$(EXEEXT) 
	./nstar_cold_ts$(EXEEXT) > nstar_cold.scr
nucmass_ldrop.scr: nucmass_ldrop_ts$(EXEEXT) 
	./nucmass_ldrop_ts$(EXEEXT) > nucmass_ldrop.scr
nucmass_ldrop_shell.scr: nucmass_ldrop_shell_ts$(EXEEXT) 
	./nucmass_ldrop_shell_ts$(EXEEXT) > nucmass_ldrop_shell.scr
eos_base.scr: eos_base_ts$(EXEEXT) 
	./eos_base_ts$(EXEEXT) > eos_base.scr
eos_had_potential.scr: eos_had_potential_ts$(EXEEXT) 
	./eos_had_potential_ts$(EXEEXT) > eos_had_potential.scr
eos_had_sym4.scr: eos_had_sym4_ts$(EXEEXT) 
	./eos_had_sym4_ts$(EXEEXT) > eos_had_sym4.scr
eos_had_base.scr: eos_had_base_ts$(EXEEXT) 
	./eos_had_base_ts$(EXEEXT) > eos_had_base.scr
eos_quark_njl.scr: eos_quark_njl_ts$(EXEEXT) 
	./eos_quark_njl_ts$(EXEEXT) > eos_quark_njl.scr
eos_quark.scr: eos_quark_ts$(EXEEXT) 
	./eos_quark_ts$(EXEEXT) > eos_quark.scr
eos_had_rmf.scr: eos_had_rmf_ts$(EXEEXT) 
	./eos_had_rmf_ts$(EXEEXT) > eos_had_rmf.scr
eos_had_schematic.scr: eos_had_schematic_ts$(EXEEXT) 
	./eos_had_schematic_ts$(EXEEXT) > eos_had_schematic.scr
eos_had_skyrme.scr: eos_had_skyrme_ts$(EXEEXT) 
	./eos_had_skyrme_ts$(EXEEXT) > eos_had_skyrme.scr
eos_tov.scr: eos_tov_ts$(EXEEXT) 
	./eos_tov_ts$(EXEEXT) > eos_tov.scr
tov_solve.scr: tov_solve_ts$(EXEEXT) 
	./tov_solve_ts$(EXEEXT) > tov_solve.scr
eos_quark_cfl6.scr: eos_quark_cfl6_ts$(EXEEXT) 
	./eos_quark_cfl6_ts$(EXEEXT) > eos_quark_cfl6.scr
eos_had_tabulated.scr: eos_had_tabulated_ts$(EXEEXT) 
	./eos_had_tabulated_ts$(EXEEXT) > eos_had_tabulated.scr
eos_had_rmf_delta.scr: eos_had_rmf_delta_ts$(EXEEXT) 
	./eos_had_rmf_delta_ts$(EXEEXT) > eos_had_rmf_delta.scr
eos_crust.scr: eos_crust_ts$(EXEEXT) 
	./eos_crust_ts$(EXEEXT) > eos_crust.scr
eos_had_ddc.scr: eos_had_ddc_ts$(EXEEXT) 
	./eos_had_ddc_ts$(EXEEXT) > eos_had_ddc.scr
eos_quark_cfl.scr: eos_quark_cfl_ts$(EXEEXT) 
	./eos_quark_cfl_ts$(EXEEXT) > eos_quark_cfl.scr
eos_nse.scr: eos_nse_ts$(EXEEXT) 
	./eos_nse_ts$(EXEEXT) > eos_nse.scr
eos_nse_full.scr: eos_nse_full_ts$(EXEEXT) 
	./eos_nse_full_ts$(EXEEXT) > eos_nse_full.scr
# eos_sn.scr: eos_sn_ts$(EXEEXT) 
# 	./eos_sn_ts$(EXEEXT) > eos_sn.scr
nucleus_rmf.scr: nucleus_rmf_ts$(EXEEXT) 
	./nucleus_rmf_ts$(EXEEXT) > nucleus_rmf.scr
eos_had_gogny.scr: eos_had_gogny_ts$(EXEEXT) 
	./eos_had_gogny_ts$(EXEEXT) > eos_had_gogny.scr
eos_crust_virial.scr: eos_crust_virial_ts$(EXEEXT) 
	./eos_crust_virial_ts$(EXEEXT) > eos_crust_virial.scr

eos_had_apr_ts_SOURCES = eos_had_apr_ts.cpp
eos_quark_bag_ts_SOURCES = eos_quark_bag_ts.cpp
nstar_cold_ts_SOURCES = nstar_cold_ts.cpp
nucmass_ldrop_ts_SOURCES = nucmass_ldrop_ts.cpp
nucmass_ldrop_shell_ts_SOURCES = nucmass_ldrop_shell_ts.cpp
eos_base_ts_SOURCES = eos_base_ts.cpp
eos_had_potential_ts_SOURCES = eos_had_potential_ts.cpp
eos_had_sym4_ts_SOURCES = eos_had_sym4_ts.cpp
eos_had_base_ts_SOURCES = eos_had_base_ts.cpp
eos_quark_njl_ts_SOURCES = eos_quark_njl_ts.cpp
eos_quark_ts_SOURCES = eos_quark_ts.cpp
eos_had_rmf_ts_SOURCES = eos_had_rmf_ts.cpp
eos_had_schematic_ts_SOURCES = eos_had_schematic_ts.cpp
eos_had_skyrme_ts_SOURCES = eos_had_skyrme_ts.cpp
eos_tov_ts_SOURCES = eos_tov_ts.cpp
tov_solve_ts_SOURCES = tov_solve_ts.cpp
eos_quark_cfl6_ts_SOURCES = eos_quark_cfl6_ts.cpp
eos_had_tabulated_ts_SOURCES = eos_had_tabulated_ts.cpp
eos_had_rmf_delta_ts_SOURCES = eos_had_rmf_delta_ts.cpp
eos_crust_ts_SOURCES = eos_crust_ts.cpp
eos_had_ddc_ts_SOURCES = eos_had_ddc_ts.cpp
eos_quark_cfl_ts_SOURCES = eos_quark_cfl_ts.cpp
eos_nse_ts_SOURCES = eos_nse_ts.cpp
eos_nse_full_ts_SOURCES = eos_nse_full_ts.cpp
# eos_sn_ts_SOURCES = eos_sn_ts.cpp
nucleus_rmf_ts_SOURCES = nucleus_rmf_ts.cpp
eos_had_gogny_ts_SOURCES = eos_had_gogny_ts.cpp
eos_crust_virial_ts_SOURCES = eos_crust_virial_ts.cpp

ls_test: ls_ls_test ls_skm_test ls_ska_test ls_sk1_test

ls_ls_test:
	eos_sn_ts ls $(HOME)/pkgs/ls_eos2/ls.dat

ls_skm_test:
	eos_sn_ts ls $(HOME)/pkgs/ls_eos2/skm.dat

ls_ska_test:
	eos_sn_ts ls $(HOME)/pkgs/ls_eos2/ska.dat

ls_sk1_test:
	eos_sn_ts ls $(HOME)/pkgs/ls_eos2/sk1.dat

stos_test: stos_eos1_test stos_eos2_test stos_eos3_test

stos_eos1_test:
	eos_sn_ts stos $(HOME)/pkgs/stos_eos/eos1

stos_eos2_test:
	eos_sn_ts stos $(HOME)/pkgs/stos_eos/eos2

stos_eos3_test:
	eos_sn_ts stos $(HOME)/pkgs/stos_eos/eos3

oo_test: oo_ls180_test oo_ls220_test oo_ls375_test

oo_ls180_test:
	eos_sn_ts oo \
	$(HOME)/pkgs/oo_eos/LS180_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_ls220_test:
	eos_sn_ts oo \
	$(HOME)/pkgs/oo_eos/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_ls375_test:
	eos_sn_ts oo \
	$(HOME)/pkgs/oo_eos/LS375_234r_136t_50y_analmu_20091212_SVNr26.h5 ls

oo_stos_test:
	eos_sn_ts oo \
	$(HOME)/pkgs/oo_eos/HShenEOS_rho220_temp180_ye65_version1.1_20120817.h5 \
	stos

oo_stos_hyp_test:
	eos_sn_ts oo \
	$(HOME)/pkgs/oo_eos/HShen_HyperonEOS_rho220_temp180_ye65_version1.1_20131007.h5 \
	stos


x:
	eos_sn_ts hfsl \
		$(HOME)/pkgs/hfsl_eos/dd2_frdm_eos_shen98format_v1.02.tab \
		 >> gen_sn_test.scr
	eos_sn_ts hfsl \
		$(HOME)/pkgs/hfsl_eos/fsg_roca_eos_shen98format_v1.0.tab \
		 >> gen_sn_test.scr
	eos_sn_ts hfsl \
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

