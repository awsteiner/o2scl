emt:
	o2graph -read ex_mcmc_bimodal.o2 -histplot x bins=40 -create \
		table x "(-4)" 4.001 0.1 \
		-function "7800*exp(-x*x)*(sin(x-1.4)+1)" y \
		-plot x y -show

emt2:
	o2graph -read ex_mcmc_gaussian.o2 table2 -histplot x0 bins=40 -create \
		table x "(-1.5)" 1.5001 0.1 \
		-function "135*exp(-(x-0.2)*(x-0.2))" y \
		-plot x y -show

emt2b:
	o2graph -read ex_mcmc_gaussian.o2 table2 -histplot x0b bins=40 -create \
		table x "(-2)" 2.001 0.1 \
		-function "100*exp(-(x-0.2)*(x-0.2))" y \
		-plot x y -show

emt3:
	o2graph -read ex_mcmc_gaussian.o2 -histplot x1 bins=40 -create \
		table x "(-2)" 2.001 0.1 \
		-function "270*exp(-(x-0.5)*(x-0.5))" y \
		-plot x y -show

ein:
	o2graph -set xlo 0 -set xhi 4 -set ylo 0 -set yhi 4 \
		-read ex_interp2_nogrid_1000.o2 t3d -den-plot exact \
		-read ex_interp2_nogrid_1000.o2 tdata \
		-scatter x y None None color=black,marker=+ \
		-save exact.pdf 
	o2graph -set xlo 0 -set xhi 4 -set ylo 0 -set yhi 4 \
		-read ex_interp2_nogrid_1000.o2 t3d \
		-den-plot prob_dens_mdim_amr1 \
		-read ex_interp2_nogrid_1000.o2 tdata \
		-scatter x y None None color=black,marker=+ \
		-save pdma.pdf

nrapr:
	ex_eos_had_skyrme -set output_files 1 -summary NRAPR
	acol -read skyrme_NRAPR_eos.o2 -output skyrme_NRAPR_eos.txt
	acol -read skyrme_NRAPR_mvsr.o2 -output skyrme_NRAPR_mvsr.txt

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

BENCHMARK_PRGS = bm_poly.scr bm_root.scr bm_min.scr bm_polylog.scr
#	bm_mroot.scr bm_rkck.scr bm_mroot2.scr bm_lu.scr \
#	bm_part.scr bm_part2.scr 
# bm_mmin.scr
EXAMPLE_PRGS = ex_diff.scr ex_fptr.scr ex_interp2.scr ex_table.scr \
	ex_nstar_cold.scr ex_mroot.scr ex_mmin.scr ex_eos_had_apr.scr \
	ex_inte.scr ex_deriv.scr ex_part.scr ex_ode.scr ex_anneal.scr \
	ex_mcarlo.scr ex_string.scr ex_contour.scr ex_eos_crust.scr \
	ex_mmin_fix.scr ex_nucmass_fit.scr ex_eos_had_rmf.scr\
	ex_stiff.scr ex_chebapp.scr ex_poly.scr ex_conmin.scr \
	ex_nucleus_rmf.scr ex_methast.scr ex_eos_had_skyrme.scr \
	ex_nucmass.scr ex_tensor.scr ex_mcmc.scr ex_nstar_rot.scr \
	ex_interp.scr ex_eos_gibbs.scr \
	ex_fermion_summ.scr ex_lambda.scr

# ex_ode_it.scr ex_eos_sn.scr ex_interp2_nogrid.scr 

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

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -DO2SCL_OPENMP -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_MPFR_MVAR) \
	$(O2SCL_HDF_MVAR) -DO2SCL_COND_FLAG $(O2SCL_OPENMP_MVAR)

# -------------------------------------------------------------
# Distribution files
# -------------------------------------------------------------

EXTRA_DIST = ex*.cpp bm*.cpp Makefile* ex_bps_input.txt ex*.h

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

EXTRA_PROGRAMS = \
	ex_anneal \
	ex_chebapp \
	ex_conmin \
	ex_const_vec \
	ex_contour \
	ex_deriv \
	ex_diff \
	ex_eos_crust \
	ex_eos_gibbs \
	ex_eos_had_apr \
	ex_eos_had_skyrme \
	ex_eos_had_rmf \
	ex_eos_sn \
	ex_fermion_summ \
	ex_fptr \
	ex_hydrogen \
	ex_inte \
	ex_interp \
	ex_interp2 \
	ex_interp2_nogrid \
	ex_lambda \
	ex_mcarlo \
	ex_mcmc \
	ex_methast \
	ex_mmin \
	ex_mmin_fix \
	ex_mroot \
	ex_nucleus_rmf \
	ex_nucmass \
	ex_nucmass_fit \
	ex_nstar_rot \
	ex_nstar_cold \
	ex_ode \
	ex_ode_it \
	ex_part \
	ex_poly \
	ex_stiff \
	ex_string \
	ex_table \
	ex_tensor \
	bm_root \
	bm_min \
	bm_poly \
	bm_polylog

if O2SCL_PYTHON

ADDL_TEST_LIBS = -lo2scl -l$(PYTHON_LIB)
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp -L../src 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR) -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = -L../src 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
endif

else

ADDL_TEST_LIBS = -lo2scl
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp -L../src 
ADDL_TEST_FLGS = -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = -L../src 
ADDL_TEST_FLGS = 
endif

endif

ex_diff_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_diff_LDADD = $(ADDL_TEST_LIBS)
ex_diff_SOURCES = ex_diff.cpp
ex_diff.scr: ex_diff
	./ex_diff > ex_diff.scr

ex_lambda_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_lambda_LDADD = $(ADDL_TEST_LIBS)
ex_lambda_SOURCES = ex_lambda.cpp
ex_lambda.scr: ex_lambda
	./ex_lambda > ex_lambda.scr

ex_eos_sn_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_sn_LDADD = $(ADDL_TEST_LIBS)
ex_eos_sn_SOURCES = ex_eos_sn.cpp
ex_eos_sn.scr: ex_eos_sn
	./ex_eos_sn -help > ex_eos_sn.scr

ex_fermion_summ_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_fermion_summ_LDADD = $(ADDL_TEST_LIBS)
ex_fermion_summ_SOURCES = ex_fermion_summ.cpp
ex_fermion_summ.scr: ex_fermion_summ
	./ex_fermion_summ -help > ex_fermion_summ.scr

ex_eos_gibbs_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_gibbs_LDADD = $(ADDL_TEST_LIBS)
ex_eos_gibbs_SOURCES = ex_eos_gibbs.cpp
ex_eos_gibbs.scr: ex_eos_gibbs
	./ex_eos_gibbs -test > ex_eos_gibbs.scr

ex_hydrogen_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_hydrogen_LDADD = $(ADDL_TEST_LIBS)
ex_hydrogen_SOURCES = ex_hydrogen.cpp
ex_hydrogen.scr: ex_hydrogen
	./ex_hydrogen > ex_hydrogen.scr

ex_fptr_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_fptr_LDADD = $(ADDL_TEST_LIBS)
ex_fptr_SOURCES = ex_fptr.cpp
ex_fptr.scr: ex_fptr
	./ex_fptr > ex_fptr.scr

ex_mcmc_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_mcmc_LDADD = $(ADDL_TEST_LIBS)
ex_mcmc_SOURCES = ex_mcmc.cpp
ex_mcmc.scr: ex_mcmc
	./ex_mcmc > ex_mcmc.scr

ex_interp_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_interp_LDADD = $(ADDL_TEST_LIBS)
ex_interp_SOURCES = ex_interp.cpp
ex_interp.scr: ex_interp
	./ex_interp > ex_interp.scr

ex_interp2_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_interp2_LDADD = $(ADDL_TEST_LIBS)
ex_interp2_SOURCES = ex_interp2.cpp
ex_interp2.scr: ex_interp2
	./ex_interp2 > ex_interp2.scr

ex_interp2_nogrid_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_interp2_nogrid_LDADD = $(ADDL_TEST_LIBS)
ex_interp2_nogrid_SOURCES = ex_interp2_nogrid.cpp
ex_interp2_nogrid.scr: ex_interp2_nogrid
	./ex_interp2_nogrid > ex_interp2_nogrid.scr

ex_table_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_table_LDADD = $(ADDL_TEST_LIBS)
ex_table_SOURCES = ex_table.cpp
ex_table.scr: ex_table
	./ex_table > ex_table.scr

ex_nstar_cold_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_nstar_cold_LDADD = $(ADDL_TEST_LIBS)
ex_nstar_cold_SOURCES = ex_nstar_cold.cpp
ex_nstar_cold.scr: ex_nstar_cold
	./ex_nstar_cold > ex_nstar_cold.scr

ex_nstar_rot_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_nstar_rot_LDADD = $(ADDL_TEST_LIBS)
ex_nstar_rot_SOURCES = ex_nstar_rot.cpp
ex_nstar_rot.scr: ex_nstar_rot
	./ex_nstar_rot > ex_nstar_rot.scr

ex_mroot_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_mroot_LDADD = $(ADDL_TEST_LIBS)
ex_mroot_SOURCES = ex_mroot.cpp
ex_mroot.scr: ex_mroot
	./ex_mroot > ex_mroot.scr

ex_mmin_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_mmin_LDADD = $(ADDL_TEST_LIBS)
ex_mmin_SOURCES = ex_mmin.cpp
ex_mmin.scr: ex_mmin
	./ex_mmin > ex_mmin.scr

ex_mmin_fix_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_mmin_fix_LDADD = $(ADDL_TEST_LIBS)
ex_mmin_fix_SOURCES = ex_mmin_fix.cpp
ex_mmin_fix.scr: ex_mmin_fix
	./ex_mmin_fix > ex_mmin_fix.scr

ex_inte_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_inte_LDADD = $(ADDL_TEST_LIBS)
ex_inte_SOURCES = ex_inte.cpp
ex_inte.scr: ex_inte
	./ex_inte > ex_inte.scr

ex_deriv_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_deriv_LDADD = $(ADDL_TEST_LIBS)
ex_deriv_SOURCES = ex_deriv.cpp
ex_deriv.scr: ex_deriv
	./ex_deriv > ex_deriv.scr

ex_part_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_part_LDADD = $(ADDL_TEST_LIBS)
ex_part_SOURCES = ex_part.cpp
ex_part.scr: ex_part
	./ex_part > ex_part.scr

ex_nucmass_fit_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_nucmass_fit_LDADD = $(ADDL_TEST_LIBS)
ex_nucmass_fit_SOURCES = ex_nucmass_fit.cpp
ex_nucmass_fit.scr: ex_nucmass_fit
	./ex_nucmass_fit > ex_nucmass_fit.scr

ex_ode_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_ode_LDADD = $(ADDL_TEST_LIBS)
ex_ode_SOURCES = ex_ode.cpp
ex_ode.scr: ex_ode
	./ex_ode > ex_ode.scr

ex_anneal_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_anneal_LDADD = $(ADDL_TEST_LIBS)
ex_anneal_SOURCES = ex_anneal.cpp
ex_anneal.scr: ex_anneal
	./ex_anneal > ex_anneal.scr

ex_mcarlo_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_mcarlo_LDADD = $(ADDL_TEST_LIBS)
ex_mcarlo_SOURCES = ex_mcarlo.cpp
ex_mcarlo.scr: ex_mcarlo
	./ex_mcarlo > ex_mcarlo.scr

ex_string_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_string_LDADD = $(ADDL_TEST_LIBS)
ex_string_SOURCES = ex_string.cpp
ex_string.scr: ex_string
	./ex_string > ex_string.scr

ex_contour_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_contour_LDADD = $(ADDL_TEST_LIBS)
ex_contour_SOURCES = ex_contour.cpp
ex_contour.scr: ex_contour
	./ex_contour > ex_contour.scr

ex_tensor_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_tensor_LDADD = $(ADDL_TEST_LIBS)
ex_tensor_SOURCES = ex_tensor.cpp
ex_tensor.scr: ex_tensor
	./ex_tensor > ex_tensor.scr

ex_eos_crust_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_crust_LDADD = $(ADDL_TEST_LIBS)
ex_eos_crust_SOURCES = ex_eos_crust.cpp
ex_eos_crust.scr: ex_eos_crust
	./ex_eos_crust > ex_eos_crust.scr

ex_eos_had_apr_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_had_apr_LDADD = $(ADDL_TEST_LIBS)
ex_eos_had_apr_SOURCES = ex_eos_had_apr.cpp
ex_eos_had_apr.scr: ex_eos_had_apr
	./ex_eos_had_apr > ex_eos_had_apr.scr

ex_eos_had_skyrme_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_had_skyrme_LDADD = $(ADDL_TEST_LIBS)
ex_eos_had_skyrme_SOURCES = ex_eos_had_skyrme.cpp
ex_eos_had_skyrme.scr: ex_eos_had_skyrme
	./ex_eos_had_skyrme -test > ex_eos_had_skyrme.scr

ex_eos_had_rmf_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_eos_had_rmf_LDADD = $(ADDL_TEST_LIBS)
ex_eos_had_rmf_SOURCES = ex_eos_had_rmf.cpp
ex_eos_had_rmf.scr: ex_eos_had_rmf
	./ex_eos_had_rmf -test > ex_eos_had_rmf.scr

ex_nucmass_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_nucmass_LDADD = $(ADDL_TEST_LIBS)
ex_nucmass_SOURCES = ex_nucmass.cpp
ex_nucmass.scr: ex_nucmass
	./ex_nucmass > ex_nucmass.scr

ex_const_vec_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_const_vec_LDADD = $(ADDL_TEST_LIBS)
ex_const_vec_SOURCES = ex_const_vec.cpp
ex_const_vec.scr: ex_const_vec
	./ex_const_vec > ex_const_vec.scr

ex_stiff_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_stiff_LDADD = $(ADDL_TEST_LIBS)
ex_stiff_SOURCES = ex_stiff.cpp
ex_stiff.scr: ex_stiff
	./ex_stiff > ex_stiff.scr

ex_methast_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_methast_LDADD = $(ADDL_TEST_LIBS)
ex_methast_SOURCES = ex_methast.cpp
ex_methast.scr: ex_methast
	./ex_methast > ex_methast.scr

ex_chebapp_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_chebapp_LDADD = $(ADDL_TEST_LIBS)
ex_chebapp_SOURCES = ex_chebapp.cpp
ex_chebapp.scr: ex_chebapp
	./ex_chebapp > ex_chebapp.scr

ex_poly_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_poly_LDADD = $(ADDL_TEST_LIBS)
ex_poly_SOURCES = ex_poly.cpp
ex_poly.scr: ex_poly
	./ex_poly > ex_poly.scr

ex_conmin_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_conmin_LDADD = $(ADDL_TEST_LIBS)
ex_conmin_SOURCES = ex_conmin.cpp
ex_conmin.scr: ex_conmin
	./ex_conmin > ex_conmin.scr

ex_ode_it_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_ode_it_LDADD = $(ADDL_TEST_LIBS)
ex_ode_it_SOURCES = ex_ode_it.cpp
ex_ode_it.scr: ex_ode_it
	./ex_ode_it > ex_ode_it.scr

ex_nucleus_rmf_LDFLAGS = $(ADDL_TEST_LDFLGS)
ex_nucleus_rmf_LDADD = $(ADDL_TEST_LIBS)
ex_nucleus_rmf_SOURCES = ex_nucleus_rmf.cpp
ex_nucleus_rmf.scr: ex_nucleus_rmf
	./ex_nucleus_rmf > ex_nucleus_rmf.scr

# -------------------------------------------------------------
# benchmark rule
# -------------------------------------------------------------

# bm_deriv_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_deriv_SOURCES = bm_deriv.cpp
# bm_deriv.scr: bm_deriv bm_deriv.cpp
# 	./bm_deriv > bm_deriv.scr

# bm_mmin_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_mmin_SOURCES = bm_mmin.cpp
# bm_mmin.scr: bm_mmin bm_mmin.cpp
# 	./bm_mmin > bm_mmin.scr

# bm_lu_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_lu_SOURCES = bm_lu.cpp
# bm_lu.scr: bm_lu bm_lu.cpp
# 	./bm_lu > bm_lu.scr

# bm_mroot2_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_mroot2_SOURCES = bm_mroot2.cpp
# bm_mroot2.scr: bm_mroot2 bm_mroot2.cpp
# 	./bm_mroot2 > bm_mroot2.scr

bm_root_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_root_SOURCES = bm_root.cpp
bm_root.scr: bm_root bm_root.cpp
	./bm_root > bm_root.scr

bm_polylog_LDFLAGS = $(ADDL_TEST_LDFLGS)
bm_polylog_LDADD = $(ADDL_TEST_LIBS)
bm_polylog_SOURCES = bm_polylog.cpp
bm_polylog.scr: bm_polylog bm_polylog.cpp
	./bm_polylog > bm_polylog.scr

bm_min_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_min_SOURCES = bm_min.cpp
bm_min.scr: bm_min bm_min.cpp
	./bm_min > bm_min.scr

# bm_mroot_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_mroot_SOURCES = bm_mroot.cpp
# bm_mroot.scr: bm_mroot bm_mroot.cpp
# 	./bm_mroot > bm_mroot.scr

# bm_part_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_part_SOURCES = bm_part.cpp
# bm_part.scr: bm_part bm_part.cpp
# 	./bm_part > bm_part.scr

# bm_part2_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_part2_SOURCES = bm_part2.cpp
# bm_part2.scr: bm_part2 bm_part2.cpp
# 	./bm_part2 > bm_part2.scr

bm_poly_LDADD = $(OOLIBS) $(OOLIBSTWO)
bm_poly_SOURCES = bm_poly.cpp
bm_poly.scr: bm_poly bm_poly.cpp
	./bm_poly > bm_poly.scr

# bm_rk8pd_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_rk8pd_SOURCES = bm_rk8pd.cpp
# bm_rk8pd.scr: bm_rk8pd bm_rk8pd.cpp
# 	./bm_rk8pd > bm_rk8pd.scr

# bm_rkck_LDADD = $(OOLIBS) $(OOLIBSTWO)
# bm_rkck_SOURCES = bm_rkck.cpp
# bm_rkck.scr: bm_rkck bm_rkck.cpp
# 	./bm_rkck > bm_rkck.scr

# clean

examples-clean: 
	-rm -f *.o *.out $(EXAMPLE_PRGS) $(BENCHMARK_PRGS) \
		$(EXTRA_PROGRAMS) examples-summary.txt

slb00:
	ex_eos_gibbs -hadrons SLB00 -quarks SLB00_bag -mvsr slb00_bag.o2
	ex_eos_gibbs -hadrons SLB00 -quarks SLB00_njl -mvsr slb00_njl.o2
	ex_eos_gibbs -hadrons SLB00_hyp -quarks SLB00_bag \
		-mvsr slb00_hyp_bag.o2
	ex_eos_gibbs -hadrons SLB00_hyp -quarks SLB00_njl \
		-mvsr slb00_hyp_njl.o2
