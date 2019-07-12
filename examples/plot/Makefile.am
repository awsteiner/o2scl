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

plot-all: ex_anneal_plot.png ex_chebapp_plot.png ex_contour_plot1.png \
	ex_contour_plot2.png ex_fptr_plot.png \
	ex_ode_airy.png 

# -------------------------------------------------------------------- 

ex_anneal_plot.png:
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
		-save ex_anneal_plot.png
	rm -f ex_anneal_conts.o2

ex_contour_plot1.png:
	o2graph -backend Agg -read ../ex_contour.o2 c1_cl \
		-plot -save ex_contour_plot1.png

ex_contour_plot2.png:
	o2graph -backend Agg -read ../ex_contour.o2 c2_cl \
		-plot -save ex_contour_plot2.png

ex_chebapp_plot.png: 
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

ex_fptr_plot.png:
	cat ../ex_fptr.scr | grep -v t > exftemp1
	acol -generic table ../ex_fptr.out -internal exftemp2.o2
	o2graph -backend Agg -set xlo "(-1)" -set xhi 2 \
		-read exftemp2.o2 -plot x y color=red \
		-ttext 0.5 "(-0.1)" x -ttext "(-0.1)" 0.5 y \
		-generic table exftemp1 -plot c1 c2 marker=s \
		-save ex_fptr_plot.png
	rm -f exftemp1 exftemp2.o2

ex_mmin_plot.o: ex_mmin_plot.cpp
	$(CXX) $(CPPFLAGS) $(CPPFLAGS) -o ex_mmin_plot.o -c ex_mmin_plot.cpp

ex_mmin_plot: ex_mmin_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_mmin_plot ex_mmin_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

ex_mmin.png:
	o2graph -set xlo "(-2)" -set xhi 3 -set ylo "(-1)" -set yhi 15 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle z \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=green \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=red \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=blue \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=magenta \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=brown \
		-generic table ../ex_mmin4.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=black \
		-show

ex_mmin.png:
	o2graph -set xlo "(-2)" -set xhi 3 -set ylo "(-1)" -set yhi 15 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle z \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=green \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=red \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=blue \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=magenta \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=brown \
		-generic table ../ex_mmin4.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=black \
		-show

ex_mmin2.png:
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

ex_nucmass_plot.o: ex_nucmass_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_nucmass_plot.o -c ex_nucmass_plot.cpp

ex_nucmass_plot: ex_nucmass_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_nucmass_plot ex_nucmass_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

ex_nuc_plot.o: ex_nuc_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_nuc_plot.o -c ex_nuc_plot.cpp

ex_nuc_plot: ex_nuc_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_nuc_plot ex_nuc_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

ex_ode_airy.png: 
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

ex_ode_bessel.png: 
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

ex_ode_bessel2.png: 
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

ex_ode_bessel3.png: 
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

ex_stiff.png:
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
		-save ex_stiff.png -show

ex_stiff_plot.o: ex_stiff_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_stiff_plot.o -c ex_stiff_plot.cpp

ex_stiff_plot: ex_stiff_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_stiff_plot ex_stiff_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl \
		-lo2scl_hdf

# Old targets

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

