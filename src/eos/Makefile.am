# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

EOS_SRCS = eos_had_apr.cpp eos_quark_bag.cpp \
	eos_crust.cpp nstar_cold.cpp eos_base.cpp eos_had_potential.cpp \
	eos_quark_njl.cpp eos_had_ddc.cpp nucmass_ldrop.cpp \
	eos_quark.cpp eos_had_rmf.cpp eos_had_schematic.cpp \
	eos_had_skyrme.cpp tov_solve.cpp eos_nse.cpp hdf_eos_io.cpp \
	eos_had_rmf_delta.cpp eos_tov.cpp eos_had_sym4.cpp \
	eos_had_base.cpp nucleus_rmf.cpp eos_sn.cpp \
	nucmass_ldrop_shell.cpp eos_quark_cfl.cpp eos_quark_cfl6.cpp \
	eos_nse_full.cpp eos_crust_virial.cpp nucleus_bin.cpp \
	nstar_rot.cpp tov_love.cpp eos_had_rmf_hyp.cpp \
	eos_python.cpp eos_leptons.cpp 

HEADER_VAR = eos_had_apr.h eos_quark_bag.h eos_crust.h eos_had_ddc.h \
	nstar_cold.h eos_base.h eos_had_potential.h nucmass_ldrop.h \
	eos_quark_njl.h eos_quark.h eos_had_rmf.h eos_quark_cfl.h \
	eos_quark_cfl6.h eos_had_schematic.h eos_had_skyrme.h \
	eos_had_tabulated.h eos_nse.h tov_solve.h eos_had_rmf_delta.h \
	eos_had_sym4.h eos_tov.h eos_had_base.h eos_cs2_poly.h \
	hdf_eos_io.h nucleus_rmf.h eos_sn.h nucmass_ldrop_shell.h \
	eos_had_gogny.h eos_crust_virial.h \
	eos_nse_full.h nstar_rot.h tov_love.h eos_had_rmf_hyp.h \
	eos_had_virial.h eos_python.h eos_leptons.h nucleus_bin.h

TEST_VAR = eos_had_apr.scr eos_quark_bag.scr nstar_cold.scr \
	eos_base.scr eos_had_potential.scr eos_had_sym4.scr \
	eos_quark_njl.scr eos_quark.scr nucmass_ldrop.scr eos_cs2_poly.scr \
	eos_had_rmf.scr eos_had_schematic.scr eos_had_skyrme.scr eos_tov.scr \
	eos_had_tabulated.scr eos_nse.scr eos_had_rmf_delta.scr \
	eos_crust.scr eos_had_ddc.scr eos_had_base.scr nucleus_rmf.scr \
	eos_crust_virial.scr eos_had_gogny.scr tov_solve.scr \
	eos_quark_cfl6.scr eos_quark_cfl.scr nucmass_ldrop_shell.scr \
	eos_nse_full.scr nstar_rot.scr tov_love.scr \
	eos_had_rmf_hyp.scr eos_leptons.scr

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -fopenmp
else
O2SCL_OPENMP_MVAR =
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

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) \
	-DO2SCL_COND_FLAG $(O2SCL_PYTHON_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
endif

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
	nucleus_rmf_ts eos_nse_full_ts eos_leptons_ts \
	nucmass_ldrop_shell_ts eos_had_gogny_ts eos_crust_virial_ts \
	nstar_rot_ts tov_love_ts eos_cs2_poly_ts eos_had_rmf_hyp_ts

check_SCRIPTS = o2scl-test

if O2SCL_PYTHON

ADDL_TEST_LIBS = ../libo2scl.la -l$(PYTHON_LIB)
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
endif

else

ADDL_TEST_LIBS = ../libo2scl.la
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = 
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = 
endif

endif

eos_had_apr_ts_LDADD = $(ADDL_TEST_LIBS)
eos_quark_bag_ts_LDADD = $(ADDL_TEST_LIBS)
nstar_cold_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_ldrop_ts_LDADD = $(ADDL_TEST_LIBS)
nucmass_ldrop_shell_ts_LDADD = $(ADDL_TEST_LIBS)
eos_leptons_ts_LDADD = $(ADDL_TEST_LIBS)
eos_base_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_potential_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_sym4_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_base_ts_LDADD = $(ADDL_TEST_LIBS)
eos_quark_njl_ts_LDADD = $(ADDL_TEST_LIBS)
eos_quark_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_rmf_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_rmf_hyp_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_schematic_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_skyrme_ts_LDADD = $(ADDL_TEST_LIBS)
eos_tov_ts_LDADD = $(ADDL_TEST_LIBS)
eos_cs2_poly_ts_LDADD = $(ADDL_TEST_LIBS)
tov_solve_ts_LDADD = $(ADDL_TEST_LIBS)
tov_love_ts_LDADD = $(ADDL_TEST_LIBS)
nstar_rot_ts_LDADD = $(ADDL_TEST_LIBS)
eos_quark_cfl6_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_tabulated_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_rmf_delta_ts_LDADD = $(ADDL_TEST_LIBS)
eos_crust_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_ddc_ts_LDADD = $(ADDL_TEST_LIBS)
eos_quark_cfl_ts_LDADD = $(ADDL_TEST_LIBS)
eos_nse_ts_LDADD = $(ADDL_TEST_LIBS)
eos_nse_full_ts_LDADD = $(ADDL_TEST_LIBS)
nucleus_rmf_ts_LDADD = $(ADDL_TEST_LIBS)
eos_had_gogny_ts_LDADD = $(ADDL_TEST_LIBS)
eos_crust_virial_ts_LDADD = $(ADDL_TEST_LIBS)

eos_had_apr_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_quark_bag_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nstar_cold_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_ldrop_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucmass_ldrop_shell_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_base_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_leptons_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_potential_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_sym4_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_base_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_quark_njl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_quark_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_rmf_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_rmf_hyp_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_schematic_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_skyrme_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_tov_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_cs2_poly_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
tov_solve_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
tov_love_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nstar_rot_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_quark_cfl6_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_tabulated_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_rmf_delta_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_crust_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_ddc_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_quark_cfl_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_nse_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_nse_full_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
nucleus_rmf_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_had_gogny_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
eos_crust_virial_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

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
eos_leptons.scr: eos_leptons_ts$(EXEEXT) 
	./eos_leptons_ts$(EXEEXT) > eos_leptons.scr
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
eos_had_rmf_hyp.scr: eos_had_rmf_hyp_ts$(EXEEXT) 
	./eos_had_rmf_hyp_ts$(EXEEXT) > eos_had_rmf_hyp.scr
eos_had_schematic.scr: eos_had_schematic_ts$(EXEEXT) 
	./eos_had_schematic_ts$(EXEEXT) > eos_had_schematic.scr
eos_had_skyrme.scr: eos_had_skyrme_ts$(EXEEXT) 
	./eos_had_skyrme_ts$(EXEEXT) > eos_had_skyrme.scr
eos_tov.scr: eos_tov_ts$(EXEEXT) 
	./eos_tov_ts$(EXEEXT) > eos_tov.scr
eos_cs2_poly.scr: eos_cs2_poly_ts$(EXEEXT) 
	./eos_cs2_poly_ts$(EXEEXT) > eos_cs2_poly.scr
tov_solve.scr: tov_solve_ts$(EXEEXT) 
	./tov_solve_ts$(EXEEXT) > tov_solve.scr
tov_love.scr: tov_love_ts$(EXEEXT) 
	./tov_love_ts$(EXEEXT) > tov_love.scr
nstar_rot.scr: nstar_rot_ts$(EXEEXT) 
	./nstar_rot_ts$(EXEEXT) > nstar_rot.scr
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
eos_leptons_ts_SOURCES = eos_leptons_ts.cpp
eos_had_potential_ts_SOURCES = eos_had_potential_ts.cpp
eos_had_sym4_ts_SOURCES = eos_had_sym4_ts.cpp
eos_had_base_ts_SOURCES = eos_had_base_ts.cpp
eos_quark_njl_ts_SOURCES = eos_quark_njl_ts.cpp
eos_quark_ts_SOURCES = eos_quark_ts.cpp
eos_had_rmf_ts_SOURCES = eos_had_rmf_ts.cpp
eos_had_rmf_hyp_ts_SOURCES = eos_had_rmf_hyp_ts.cpp
eos_had_schematic_ts_SOURCES = eos_had_schematic_ts.cpp
eos_had_skyrme_ts_SOURCES = eos_had_skyrme_ts.cpp
eos_tov_ts_SOURCES = eos_tov_ts.cpp
eos_cs2_poly_ts_SOURCES = eos_cs2_poly_ts.cpp
tov_solve_ts_SOURCES = tov_solve_ts.cpp
tov_love_ts_SOURCES = tov_love_ts.cpp
nstar_rot_ts_SOURCES = nstar_rot_ts.cpp
eos_quark_cfl6_ts_SOURCES = eos_quark_cfl6_ts.cpp
eos_had_tabulated_ts_SOURCES = eos_had_tabulated_ts.cpp
eos_had_rmf_delta_ts_SOURCES = eos_had_rmf_delta_ts.cpp
eos_crust_ts_SOURCES = eos_crust_ts.cpp
eos_had_ddc_ts_SOURCES = eos_had_ddc_ts.cpp
eos_quark_cfl_ts_SOURCES = eos_quark_cfl_ts.cpp
eos_nse_ts_SOURCES = eos_nse_ts.cpp
eos_nse_full_ts_SOURCES = eos_nse_full_ts.cpp
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

interface: empty
	../../bin/yanic eos.i eos_python ~/wcs/o2sclpy/o2sclpy/eos \
		~/wcs/o2sclpy/doc/eos  ../../bin/yanic_header.txt \
		| tee eos_python.yanic.out

empty:

# ------------------------------------------------------------
# Library o2scl_eos
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_eos.la

libo2scl_eos_la_SOURCES = $(EOS_SRCS)

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
