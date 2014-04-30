snt:
	ex_gen_sn_eos -set directory ../../../../../pkgs/oo_eos/ \
		-oo fsu17 -slices
	ex_gen_sn_eos -set directory ../../../../../pkgs/sht_eos/ \
		-sht fsu17 -slices

# -------------------------------------------------------------
# Help
# -------------------------------------------------------------
help:
	@echo 
	@echo "This the O2scl makefile for the examples subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply. In addition there"
	@echo "are a couple extra targets which may be useful listed below"
	@echo
	@echo "Additional package-specific targets:"
	@echo "----------------------------------------------------------"
	@echo "help:             Show this help file"
	@echo "o2scl-examples:   Make the examples"
	@echo "o2scl-benchmarks: Make the benchmarks (time consuming)"
	@echo "examples-clean:   Clean this directory (remove *.o and exe's)"
	@echo 

# -------------------------------------------------------------
# Variables
# -------------------------------------------------------------

SUBDIRS = plot

if O2SCL_EOSLIB
BENCHMARK_PRGS = bm_poly.scr bm_mmin.scr 
	bm_mroot.scr bm_rkck.scr bm_mroot2.scr bm_lu.scr \
	bm_part.scr bm_part2.scr bm_root.scr bm_min.scr
EXAMPLE_PRGS = ex_diff.scr ex_fptr.scr ex_interp2.scr ex_table.scr \
	ex_cold_nstar.scr ex_mroot.scr ex_mmin.scr ex_apr_eos.scr \
	ex_inte.scr ex_deriv.scr ex_part.scr ex_ode.scr ex_anneal.scr \
	ex_mcarlo.scr ex_string.scr ex_contour.scr ex_bps_eos.scr \
	ex_mmin_fix.scr ex_nucmass_fit.scr \
	ex_stiff.scr ex_markov.scr ex_chebapp.scr ex_poly.scr ex_conmin.scr \
	ex_rmf_nuc.scr ex_methast.scr ex_skyrme_eos.scr ex_nucmass.scr \
	ex_minte.scr

# ex_ode_it.scr ex_gen_sn_eos.scr 

else

BENCHMARK_PRGS = bm_poly.scr bm_mmin.scr \
	bm_mroot.scr bm_rkck.scr bm_mroot2.scr bm_lu.scr \
	bm_root.scr bm_min.scr
EXAMPLE_PRGS = ex_diff.scr ex_fptr.scr ex_interp2.scr ex_table.scr \
	ex_mroot.scr ex_mmin.scr ex_inte.scr ex_deriv.scr \
	ex_ode.scr ex_anneal.scr ex_mcarlo.scr ex_string.scr ex_mmin_fix.scr \
	ex_stiff.scr ex_markov.scr \
	ex_chebapp.scr ex_poly.scr ex_conmin.scr ex_methast.scr \
	ex_minte.scr

# ex_ode_it.scr 

endif

if O2SCL_EOSLIB
if O2SCL_HDF
OOLIBSTWO = -lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl \
		-L../src -L../src/eos -L../src/hdf
else
OOLIBSTWO = -lo2scl_eos -lo2scl_part -lo2scl \
		-L../src -L../src/eos
endif
else
if O2SCL_HDF
OOLIBSTWO = -lo2scl_hdf -lo2scl -L../src -L../src/hdf
else 
OOLIBSTWO = -lo2scl -L../src
endif
endif

if O2SCL_CPP11
O2SCL_CPP11_MVAR = -DO2SCL_CPP11 -std=c++0x
else
O2SCL_CPP11_MVAR =
endif

if O2SCL_HDF
O2SCL_HDF_MVAR = -DO2SCL_HDF
else
O2SCL_HDF_MVAR =
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_HDF_MVAR) $(O2SCL_CPP11_MVAR) -DO2SCL_COND_FLAG

# -------------------------------------------------------------
# Distribution files
# -------------------------------------------------------------

EXTRA_DIST = ex*.cpp bm*.cpp Makefile* ex_bps_input.txt

# -------------------------------------------------------------
# Main targets
# -------------------------------------------------------------

../src/internal/tsumm$(EXEEXT): 
	cd ../src/internal; make tsumm

o2scl-examples: $(EXAMPLE_PRGS) ../src/internal/tsumm$(EXEEXT)
	echo $(EXAMPLE_PRGS) > testlist
	../src/internal/tsumm$(EXEEXT)

o2scl-benchmarks: $(BENCHMARK_PRGS)

o2scl-examples-doc: 

# -------------------------------------------------------------
# EXTRA programs target. This lists extra programs which aren't
# installed or included in the distribution. Some of these,
# like ex_hydrogen, are included here but not part of the 
# o2scl_examples target.
# -------------------------------------------------------------

if O2SCL_EOSLIB

EXTRA_PROGRAMS = \
	ex_fptr \
	ex_diff \
	ex_interp2 \
	ex_table \
	ex_cold_nstar \
	ex_mroot \
	ex_mmin \
	ex_inte \
	ex_deriv \
	ex_part \
	ex_ode \
	ex_ode_it \
	ex_anneal \
	ex_mcarlo \
	ex_minte \
	ex_string \
	ex_gen_sn_eos \
	ex_contour \
	ex_bps_eos \
	ex_apr_eos \
	ex_skyrme_eos \
	ex_nucmass \
	ex_hydrogen \
	ex_mmin_fix \
	ex_stiff \
	ex_nucmass_fit \
	ex_const_vec \
	ex_methast \
	ex_markov \
	ex_chebapp \
	ex_poly \
	ex_conmin \
	ex_rmf_nuc \
	ex_lambda \
	bm_lu \
	bm_deriv \
	bm_mmin \
	bm_mroot2 \
	bm_mroot \
	bm_root \
	bm_min \
	bm_part \
	bm_poly \
	bm_rk8pd \
	bm_rkck \
	bm_part2 

else

EXTRA_PROGRAMS = ex_fptr \
	ex_diff \
	ex_interp2 \
	ex_table \
	ex_mroot \
	ex_mmin \
	ex_inte \
	ex_deriv \
	ex_ode \
	ex_anneal \
	ex_mcarlo \
	ex_minte \
	ex_string \
	ex_contour \
	ex_hydrogen \
	ex_mmin_fix \
	ex_const_vec \
	ex_ode_it \
	ex_methast \
	ex_stiff \
	ex_markov \
	ex_chebapp \
	ex_poly \
	ex_conmin \
	ex_lambda \
	bm_lu \
	bm_deriv \
	bm_mmin \
	bm_mroot2 \
	bm_root \
	bm_min \
	bm_mroot \
	bm_poly \
	bm_rk8pd \
	bm_rkck 

endif

ex_diff_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_diff_SOURCES = ex_diff.cpp
ex_diff.scr: ex_diff
	./ex_diff > ex_diff.scr

ex_lambda_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_lambda_SOURCES = ex_lambda.cpp
ex_lambda.scr: ex_lambda
	./ex_lambda > ex_lambda.scr

ex_gen_sn_eos_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_gen_sn_eos_SOURCES = ex_gen_sn_eos.cpp
ex_gen_sn_eos.scr: ex_gen_sn_eos
	./ex_gen_sn_eos -help > ex_gen_sn_eos.scr

ex_hydrogen_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_hydrogen_SOURCES = ex_hydrogen.cpp
ex_hydrogen.scr: ex_hydrogen
	./ex_hydrogen > ex_hydrogen.scr

ex_fptr_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_fptr_SOURCES = ex_fptr.cpp
ex_fptr.scr: ex_fptr
	./ex_fptr > ex_fptr.scr

ex_interp2_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_interp2_SOURCES = ex_interp2.cpp
ex_interp2.scr: ex_interp2
	./ex_interp2 > ex_interp2.scr

ex_table_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_table_SOURCES = ex_table.cpp
ex_table.scr: ex_table
	./ex_table > ex_table.scr

ex_cold_nstar_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_cold_nstar_SOURCES = ex_cold_nstar.cpp
ex_cold_nstar.scr: ex_cold_nstar
	./ex_cold_nstar > ex_cold_nstar.scr

ex_mroot_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_mroot_SOURCES = ex_mroot.cpp
ex_mroot.scr: ex_mroot
	./ex_mroot > ex_mroot.scr

ex_mmin_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_mmin_SOURCES = ex_mmin.cpp
ex_mmin.scr: ex_mmin
	./ex_mmin > ex_mmin.scr

ex_mmin_fix_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_mmin_fix_SOURCES = ex_mmin_fix.cpp
ex_mmin_fix.scr: ex_mmin_fix
	./ex_mmin_fix > ex_mmin_fix.scr

ex_inte_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_inte_SOURCES = ex_inte.cpp
ex_inte.scr: ex_inte
	./ex_inte > ex_inte.scr

ex_deriv_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_deriv_SOURCES = ex_deriv.cpp
ex_deriv.scr: ex_deriv
	./ex_deriv > ex_deriv.scr

ex_part_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_part_SOURCES = ex_part.cpp
ex_part.scr: ex_part
	./ex_part > ex_part.scr

ex_nucmass_fit_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_nucmass_fit_SOURCES = ex_nucmass_fit.cpp
ex_nucmass_fit.scr: ex_nucmass_fit
	./ex_nucmass_fit > ex_nucmass_fit.scr

ex_ode_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_ode_SOURCES = ex_ode.cpp
ex_ode.scr: ex_ode
	./ex_ode > ex_ode.scr

ex_anneal_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_anneal_SOURCES = ex_anneal.cpp
ex_anneal.scr: ex_anneal
	./ex_anneal > ex_anneal.scr

ex_mcarlo_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_mcarlo_SOURCES = ex_mcarlo.cpp
ex_mcarlo.scr: ex_mcarlo
	./ex_mcarlo > ex_mcarlo.scr

ex_minte_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_minte_SOURCES = ex_minte.cpp
ex_minte.scr: ex_minte
	./ex_minte > ex_minte.scr

ex_string_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_string_SOURCES = ex_string.cpp
ex_string.scr: ex_string
	./ex_string > ex_string.scr

ex_contour_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_contour_SOURCES = ex_contour.cpp
ex_contour.scr: ex_contour
	./ex_contour > ex_contour.scr

ex_bps_eos_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_bps_eos_SOURCES = ex_bps_eos.cpp
ex_bps_eos.scr: ex_bps_eos
	./ex_bps_eos > ex_bps_eos.scr

ex_apr_eos_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_apr_eos_SOURCES = ex_apr_eos.cpp
ex_apr_eos.scr: ex_apr_eos
	./ex_apr_eos > ex_apr_eos.scr

ex_skyrme_eos_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_skyrme_eos_SOURCES = ex_skyrme_eos.cpp
ex_skyrme_eos.scr: ex_skyrme_eos
	./ex_skyrme_eos -test > ex_skyrme_eos.scr

ex_nucmass_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_nucmass_SOURCES = ex_nucmass.cpp
ex_nucmass.scr: ex_nucmass
	./ex_nucmass > ex_nucmass.scr

ex_const_vec_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_const_vec_SOURCES = ex_const_vec.cpp
ex_const_vec.scr: ex_const_vec
	./ex_const_vec > ex_const_vec.scr

ex_stiff_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_stiff_SOURCES = ex_stiff.cpp
ex_stiff.scr: ex_stiff
	./ex_stiff > ex_stiff.scr

ex_markov_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_markov_SOURCES = ex_markov.cpp
ex_markov.scr: ex_markov
	./ex_markov > ex_markov.scr

ex_methast_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_methast_SOURCES = ex_methast.cpp
ex_methast.scr: ex_methast
	./ex_methast > ex_methast.scr

ex_chebapp_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_chebapp_SOURCES = ex_chebapp.cpp
ex_chebapp.scr: ex_chebapp
	./ex_chebapp > ex_chebapp.scr

ex_poly_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_poly_SOURCES = ex_poly.cpp
ex_poly.scr: ex_poly
	./ex_poly > ex_poly.scr

ex_conmin_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_conmin_SOURCES = ex_conmin.cpp
ex_conmin.scr: ex_conmin
	./ex_conmin > ex_conmin.scr

ex_ode_it_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_ode_it_SOURCES = ex_ode_it.cpp
ex_ode_it.scr: ex_ode_it
	./ex_ode_it > ex_ode_it.scr

ex_rmf_nuc_LDADD = $(OOLIBS) $(OOLIBSTWO)
ex_rmf_nuc_SOURCES = ex_rmf_nuc.cpp
ex_rmf_nuc.scr: ex_rmf_nuc
	./ex_rmf_nuc > ex_rmf_nuc.scr

# -------------------------------------------------------------
# benchmark rule
# -------------------------------------------------------------

bm_deriv_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_deriv_SOURCES = bm_deriv.cpp
bm_deriv.scr: bm_deriv bm_deriv.cpp
	./bm_deriv > bm_deriv.scr

bm_mmin_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_mmin_SOURCES = bm_mmin.cpp
bm_mmin.scr: bm_mmin bm_mmin.cpp
	./bm_mmin > bm_mmin.scr

bm_lu_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_lu_SOURCES = bm_lu.cpp
bm_lu.scr: bm_lu bm_lu.cpp
	./bm_lu > bm_lu.scr

bm_mroot2_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_mroot2_SOURCES = bm_mroot2.cpp
bm_mroot2.scr: bm_mroot2 bm_mroot2.cpp
	./bm_mroot2 > bm_mroot2.scr

bm_root_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_root_SOURCES = bm_root.cpp
bm_root.scr: bm_root bm_root.cpp
	./bm_root > bm_root.scr

bm_min_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_min_SOURCES = bm_min.cpp
bm_min.scr: bm_min bm_min.cpp
	./bm_min > bm_min.scr

bm_mroot_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_mroot_SOURCES = bm_mroot.cpp
bm_mroot.scr: bm_mroot bm_mroot.cpp
	./bm_mroot > bm_mroot.scr

bm_part_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_part_SOURCES = bm_part.cpp
bm_part.scr: bm_part bm_part.cpp
	./bm_part > bm_part.scr

bm_part2_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_part2_SOURCES = bm_part2.cpp
bm_part2.scr: bm_part2 bm_part2.cpp
	./bm_part2 > bm_part2.scr

bm_poly_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_poly_SOURCES = bm_poly.cpp
bm_poly.scr: bm_poly bm_poly.cpp
	./bm_poly > bm_poly.scr

bm_rk8pd_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_rk8pd_SOURCES = bm_rk8pd.cpp
bm_rk8pd.scr: bm_rk8pd bm_rk8pd.cpp
	./bm_rk8pd > bm_rk8pd.scr

bm_rkck_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_rkck_SOURCES = bm_rkck.cpp
bm_rkck.scr: bm_rkck bm_rkck.cpp
	./bm_rkck > bm_rkck.scr

# clean

examples-clean: 
	-rm -f *.o *.out $(EXAMPLE_PRGS) $(BENCHMARK_PRGS) \
		$(EXTRA_PROGRAMS) examples-summary.txt

