# --------------------------------------------------------------------
# All figures should be generated with python or o2graph
#  
# Reworking June 2019
# 
# -------------------------------------------------------------------- 

help:
	@echo 
	@echo "This directory contains code which assists the developer"
	@echo "in generating the documentation plots and is not"
	@echo "for use by the casual end-user."
	@echo

# -------------------------------------------------------------------- 
# This directory generates the documentation images after all of
# the examples in the parent directory have been compiled and 
# executed. The main target, doc-extras, should create all of the files
# necessary to generate the documentation.
# -------------------------------------------------------------------- 

EXTRA_DIST = *.eps *.png

plot-all: ex_anneal_plot ex_chebapp_plot ex_contour_plots \
	ex_fptr_plot ex_nuc_prof ex_nucmass_plots \
	ex_ode_airy_plot ex_ode_bessel_plot ex_ode_bessel2_plot \
	ex_ode_bessel3_plot ex_stiff_plot 

if O2SCL_LINUX
IMAGE_VIEW = eog
IMAGE_VIEW2 = &
else
if O2SCL_OSX
IMAGE_VIEW = open
IMAGE_VIEW2 = 
else 
IMAGE_VIEW =
endif
endif

# UNAME_S := $(shell uname -s)

# ifeq ($(UNAME_S),Darwin)
#   IMAGE_VIEW = open
# else
#   IMAGE_VIEW = eog
# endif

# -------------------------------------------------------------------- 

ex_anneal_plot:
	o2graph -backend Agg -set colbar 1 \
		-read ../ex_anneal_plot.o2 -set xtitle x -set ytitle y \
		-den-plot f \
		-contours "(-0.9)" f ex_anneal_conts.o2 l1 \
		-contours "(-0.7)" f ex_anneal_conts.o2 l2 \
		-contours "(-0.5)" f ex_anneal_conts.o2 l3 \
		-contours "(-0.3)" f ex_anneal_conts.o2 l4 \
		-contours "(-0.2)" f ex_anneal_conts.o2 l5 \
		-contours "(-0.1)" f ex_anneal_conts.o2 l6 \
		-contours "(-0.05)" f ex_anneal_conts.o2 l7 \
		-contours "0.05" f ex_anneal_conts.o2 l8 \
		-contours "0.1" f ex_anneal_conts.o2 l9 \
		-read ex_anneal_conts.o2 l1 -plot color=black \
		-read ex_anneal_conts.o2 l2 -plot color=black \
		-read ex_anneal_conts.o2 l3 -plot color=black \
		-read ex_anneal_conts.o2 l4 -plot color=black \
		-read ex_anneal_conts.o2 l5 -plot color=black \
		-read ex_anneal_conts.o2 l6 -plot color=black \
		-read ex_anneal_conts.o2 l7 -plot color=black \
		-read ex_anneal_conts.o2 l8 -plot color=black \
		-read ex_anneal_conts.o2 l9 -plot color=black \
		-point 6 7 marker=s,color=black \
		-text 6 6.5 "Initial point" color=black \
		-point 2 "(-3)" marker=s,color=white \
		-text 2 "(-2.5)" "Global minimum" color=white \
		-save ex_anneal_plot_new.png
	rm -f ex_anneal_conts.o2
	convert +append ex_anneal_plot.png ex_anneal_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_contour_plot1:
	o2graph -backend Agg -read ../ex_contour.o2 c1_cl \
		-plot -save ex_contour_plot1_new.png
	convert +append ex_contour_plot.png ex_contour_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_contour_plot2:
	o2graph -backend Agg -read ../ex_contour.o2 c2_cl \
		-plot -save ex_contour_plot2_new.png
	convert +append ex_contour_plot.png ex_contour_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_chebapp_plot: 
	o2graph -backend Agg -generic table ../ex_chebapp.out \
		-plot c1 c2 color=black -plot c1 c4 color=red \
		-plot c1 c5 color=blue,ls=: \
		-line 0.8 "(-0.5)" 0.95 "(-0.5)" color=black \
		-text 0.78 "(-0.5)" Exact ha=right \
		-line 0.8 "(-0.65)" 0.95 "(-0.65)" color=red \
		-text 0.78 "(-0.65)" "Approx. (n=50)" ha=right,color=red \
		-line 0.8 "(-0.8)" 0.95 "(-0.8)" color=blue \
		-text 0.78 "(-0.8)" "Approx. (n=25)" ha=right,color=blue \
		-ttext 0.5 "(-0.1)" x \
		-save ex_chebapp_plot.png
	convert +append ex_chebapp_plot.png ex_chebapp_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_fptr_plot:
	cat ../ex_fptr.scr | grep -v t > exftemp1
	acol -generic table ../ex_fptr.out -internal exftemp2.o2
	o2graph -backend Agg -set xlo "(-1)" -set xhi 2 \
		-read exftemp2.o2 -plot x y color=red \
		-ttext 0.5 "(-0.1)" x -ttext "(-0.1)" 0.5 y \
		-generic table exftemp1 -plot c1 c2 marker=s \
		-save ex_fptr_plot.png
	rm -f exftemp1 exftemp2.o2
	convert +append ex_fptr_plot.png ex_fptr_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_mmin:
	o2graph -set xlo "(-2)" -set xhi 3 -set ylo "(-1)" -set yhi 15 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle z \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=green \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=red,ls=: \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=blue \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=magenta,ls=-. \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=brown \
		-generic table ../ex_mmin4.dat -internal ex_mmin.o2 \
		-scatter c1 c3 color=black \
		-show
	convert +append ex_mmin_plot.png ex_mmin_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_mmin2:
	o2graph -set xlo "(-3)" -set xhi 3 -set ylo "(-3)" -set yhi 3 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle y \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=green \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=red \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=blue \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=magenta \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=brown \
		-generic table ../ex_mmin4.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=black \
		-show
	convert +append ex_mmin_plot.png ex_mmin_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_nuc_prof:
	o2graph -read ../ex_nucleus_rmf_prof.o2 -set xlo 0 -set xhi 10 \
		-set ylo 0 -set yhi 0.1 \
		-set xtitle "$$ r~(\mathrm{fm})$$" -set ytitle \
		"$$ \mathrm{density}~(\mathrm{fm}^{-3}) $$" \
		-plot r rhon color=red -plot r rhop color=blue \
		-plot r chden color=black -plot r chden_exp \
		color=black,ls=-. \
		-text 4.5 0.09 "neutrons" color=red \
		-text 0.9 0.074 "protons" color=blue \
		-text 4.7 0.04 "charge" color=black \
		-line 5.4 0.04 6.0 0.0488 color=black \
		-text 3.5 0.066 "charge (expt.)" color=black \
		-save ex_nuc_prof.png
	convert +append ex_nuc_prof_plot.png ex_nuc_prof_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_nucmass_plots:
	o2graph -backend Agg \
		-set xtitle N -set ytitle Z -read ../ex_nucmass_table.o2 \
		-to-table3d N Z -set colbar 1 \
		-den-plot se cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Semi-empirical} $$" \
		-save ex_nucmass_se.png \
		-clf -den-plot mnmsk cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{M\ddot{o}ller~et~al.~(1995)} $$" \
		-save ex_nucmass_mnmsk.png \
		-clf -den-plot hfb14 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~14} $$" \
		-save ex_nucmass_hfb14.png \
		-clf -den-plot hfb21 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~21} $$" \
		-save ex_nucmass_hfb21.png \
		-clf -den-plot hfb27 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~27} $$" \
		-save ex_nucmass_hfb27.png \
		-clf -den-plot ame03 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{AME~(2003)} $$" \
		-save ex_nucmass_ame03.png \
		-clf -den-plot dz96 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Duflo~and~Zuker~(1996)} $$" \
		-save ex_nucmass_dz96.png \
		-clf -den-plot ktuy05 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Koura~et~al.~(1996)} $$" \
		-save ex_nucmass_ktuy05.png \
		-clf -den-plot dvi cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Dieperink~et~al.~(2009)} $$" \
		-save ex_nucmass_dvi.png \
		-clf -den-plot ws32 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Wang~et~al.~(2010)} $$" \
		-save ex_nucmass_ws32.png \
		-clf -den-plot ws36 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Liu~et~al.~(2011)} $$" \
		-save ex_nucmass_ws36.png 
	convert +append ex_anneal_plot.png ex_anneal_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_ode_airy: 
	o2graph -read ../ex_ode.o2 table_2 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err)" abserr \
		-plot x abserr -plot x diff \
		-read ../ex_ode.o2 table_3 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err)" abserr \
		-plot x abserr -plot x diff \
		-ttext 0.57 0.14 "Prince-Dormand act. error" \
		-ttext 0.25 0.96 "Cash-Karp act. error" \
		-ttext 0.60 0.86 "Cash-Karp est. error" \
		-ttext 0.5 0.21 "Prince-Dormand est. error" \
		-ttext 0.5 1.03 "Airy function w/non-adaptive steppers" \
		-save ex_ode_airy.png
	convert +append ex_ode_airy_plot.png ex_ode_airy_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_ode_bessel: 
	o2graph -backend Agg -read ../ex_ode.o2 table_0 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-read ../ex_ode.o2 table_1 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-ttext 0.46 0.87 "Prince-Dormand act. error" \
		-ttext 0.70 0.67 "Cash-Karp act. error" \
		-ttext 0.60 0.48 "Cash-Karp est. error" \
		-ttext 0.48 0.09 "Prince-Dormand est. error" \
		-ttext 0.5 1.03 "Bessel function w/non-adaptive steppers" \
		-save ex_ode_bessel.png
	convert +append ex_ode_bessel_plot.png ex_ode_bessel_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_ode_bessel2: 
	o2graph -backend Agg -read ../ex_ode.o2 table_4 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-11 -set yhi 2.5e-6 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err0)" abserr \
		-plot x abserr color=black,ls=-. -plot x diff color=black \
		-read ../ex_ode.o2 table_5 \
		-function "abs(calc-exact)" diff \
		-function "abs(err0)" abserr \
		-plot x abserr color=red,ls=-. -plot x diff color=red \
		-read ../ex_ode.o2 table_6 \
		-function "abs(calc-exact)" diff \
		-function "abs(err0)" abserr \
		-plot x abserr color=blue,ls=-. -plot x diff color=blue \
		-ttext 0.45 0.97 "Cash-Karp act. error" color=black \
		-ttext 0.22 0.92 "Cash-Karp est. error" color=black \
		-ttext 0.55 0.65 "Cash-Karp(2) act. error" color=red \
		-ttext 0.30 0.60 "Prince-Dormand act. error" color=blue \
		-ttext 0.70 0.10 "Cash-Karp(2) est. error" color=red \
		-ttext 0.35 0.05 "Prince-Dormand est. error" color=blue \
		-ttext 0.5 1.03 "Bessel function w/adaptive steppers" \
		-save ex_ode_bessel2.png
	convert +append ex_ode_bessel2.png ex_ode_bessel2_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_ode_bessel3: 
	o2graph -backend Agg -read ../ex_ode.o2 table_7 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-11 -set yhi 5.0e-8 \
		-set xlo 0 -set xhi 10 -set logy 1 \
		-function "abs(err0)" abserr \
		-plot x abserr color=red,ls=-. -plot x diff color=black \
		-ttext 0.45 0.95 "Act. error" color=black \
		-ttext 0.50 0.50 "Est. error" color=red \
		-ttext 0.5 1.03 \
		"Bessel function, high-level adaptive stepper" \
		-save ex_ode_bessel3.png
	convert +append ex_ode_bessel3.png ex_ode_bessel3_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

ex_stiff_plot:
	o2graph -backend Agg -read ../ex_stiff.o2 table_0 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-16 -set yhi 1.0e8 -set logy 1 \
		-set xlo 0 -set xhi 4 -plot x rel_err color=black \
		-plot x rel_diff color=red \
		-read ../ex_stiff.o2 table_1 -plot x rel_err color=blue \
		-plot x rel_diff color=green \
		-ttext 0.35 0.83 \
		"Explicit actual err." color=green \
		-ttext 0.5 0.53 \
		"Explicit estimated err." color=blue \
		-ttext 0.5 0.15 \
		"BS implicit actual err." color=red \
		-ttext 0.5 0.05 \
		"BS implicit estimated err." color=black \
		-ttext 0.5 1.03 \
		"Stiff exponential function, implicit vs. explicit" \
		-save ex_stiff_plot.png -show
	convert +append ex_stiff_plot.png ex_stiff_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png

# Old targets

# ex_mmin_plot.o: ex_mmin_plot.cpp
# 	$(CXX) $(CPPFLAGS) $(CPPFLAGS) -o ex_mmin_plot.o -c ex_mmin_plot.cpp

# ex_mmin_plot: ex_mmin_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_mmin_plot ex_mmin_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

# ex_stiff_plot.o: ex_stiff_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_stiff_plot.o -c ex_stiff_plot.cpp

# ex_stiff_plot: ex_stiff_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_stiff_plot ex_stiff_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl \
# 		-lo2scl_hdf

# ex_anneal_plot.o: ex_anneal_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_anneal_plot.o -c ex_anneal_plot.cpp

# ex_anneal_plot: ex_anneal_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_anneal_plot ex_anneal_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

#ex_chebapp_plot.eps: ex_chebapp_plot.py
#	python ex_chebapp_plot.py

#ex_fptr_plot.eps: ex_fptr_plot.py
#	python ex_fptr_plot.py

# ex_mandel_plot.o: ex_mandel_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_mandel_plot.o -c ex_mandel_plot.cpp

# ex_mandel_plot: ex_mandel_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_mandel_plot ex_mandel_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

# ex_mandel_plot_all: ex_mandel_plot
# 	./ex_mandel_plot
# 	convert ex_mandel_plot.png ex_mandel_plot.pdf

#ex_ode_bessel.eps: ex_ode_plot.py
#	python ex_ode_plot.py

# ex_planar_plot.o: ex_planar_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_planar_plot.o -c ex_planar_plot.cpp

# ex_planar_plot: ex_planar_plot.o ../../src/other/graph.o
# 	$(CXX) $(CPPFLAGS) -o ex_planar_plot ex_planar_plot.o \
# 		-lo2scl_hdf -lo2scl $(ROOT_LIB) $(LDFLAGS) 

# ex_planar_plot.pdf: ex_planar_plot
# 	./ex_planar_plot
# 	convert ex_planar_plot.png ex_planar_plot.pdf

# ex_ublas.o: ex_ublas.cpp
# 	$(CXX) -I/home/asteiner/install/include $(INCLUDES) \
# 		-o ex_ublas.o -c ex_ublas.cpp

# ex_ublas: ex_ublas.o  ../../src/other/graph.o
# 	$(CXX) -I/home/asteiner/install/include -L/usr/lib64/atlas \
# 		-L/home/asteiner/install/lib -o ex_ublas ex_ublas.o -lgsl \
# 		-lo2scl

