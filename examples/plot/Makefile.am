# --------------------------------------------------------------------
# Over the long term, I'd like to make sure that all figures can be
# automatically generated with python
#  
# Working as of 12/28/13:
# ex_fptr*, ex_chebapp*, ex_ode_bessel, ex_ode_airy, ex_ode_bessel2
#
#   ex_anneal_plot.png
#   ex_anneal_plot.eps
#   ex_planar_intp.png
#   ex_planar_intp.pdf
#   ex_contour_plot1.png
#   ex_contour_plot1.eps
#   ex_contour_plot2.png
#   ex_contour_plot2.eps
#   ex_contour_plot3.png
#   ex_contour_plot3.eps
#   ex_ode_bessel3.png
#   ex_ode_bessel3.eps
#   ex_ode_stiff.png
#   ex_ode_stiff.eps
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

../../src/other/graph.o:
	cd ../../src/other; $(CXX) $(CPPFLAGS) -o graph.o -c graph.cpp

plot-all: ex_fptr_plot.eps ex_chebapp_plot.eps ex_ode_bessel.eps

# -------------------------------------------------------------------- 

# In approximate alphabetical order

ex_anneal_plot.png:
	o2graph -backend Agg -set colbar 1 \
		-read ../ex_anneal_plot.o2 -set xtitle x -set ytitle y \
		-den-plot f \
		-contours "(-0.9)" f conts.o2 l1 \
		-contours "(-0.7)" f conts.o2 l2 \
		-contours "(-0.5)" f conts.o2 l3 \
		-contours "(-0.3)" f conts.o2 l4 \
		-contours "(-0.2)" f conts.o2 l5 \
		-contours "(-0.1)" f conts.o2 l6 \
		-contours "(-0.05)" f conts.o2 l7 \
		-contours "0.05" f conts.o2 l8 \
		-contours "0.1" f conts.o2 l9 \
		-read conts.o2 l1 -plot color=black \
		-read conts.o2 l2 -plot color=black \
		-read conts.o2 l3 -plot color=black \
		-read conts.o2 l4 -plot color=black \
		-read conts.o2 l5 -plot color=black \
		-read conts.o2 l6 -plot color=black \
		-read conts.o2 l7 -plot color=black \
		-read conts.o2 l8 -plot color=black \
		-read conts.o2 l9 -plot color=black \
		-point 6 7 marker=s,color=black \
		-text 6 6.5 "Initial point" color=black \
		-point 2 "(-3)" marker=s,color=white \
		-text 2 "(-2.5)" "Global minimum" color=white \
		-save ex_anneal_plot.png
	rm -f conts.o2

# ex_anneal_plot.o: ex_anneal_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_anneal_plot.o -c ex_anneal_plot.cpp

# ex_anneal_plot: ex_anneal_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_anneal_plot ex_anneal_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

#ex_chebapp_plot.eps: ex_chebapp_plot.py
#	python ex_chebapp_plot.py

ex_contour_plot.o: ex_contour_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_contour_plot.o -c ex_contour_plot.cpp

ex_contour_plot: ex_contour_plot.o ../../src/other/graph.o
	$(CXX) $(CPPFLAGS) -o ex_contour_plot ex_contour_plot.o \
		-lo2scl_hdf -lo2scl $(ROOT_LIB) $(LDFLAGS) 

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

ex_mmin_plot.o: ex_mmin_plot.cpp
	$(CXX) $(CPPFLAGS) $(CPPFLAGS) -o ex_mmin_plot.o -c ex_mmin_plot.cpp

ex_mmin_plot: ex_mmin_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_mmin_plot ex_mmin_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

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
		-set xtitle x -set ytitle y -set top-margin 0.05 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-read ../ex_ode.o2 table_3 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-ttext 0.46 0.87 "Prince-Dormand act. error" \
		-ttext 0.70 0.67 "Cash-Karp act. error" \
		-ttext 0.60 0.48 "Cash-Karp est. error" \
		-ttext 0.48 0.09 "Prince-Dormand est. error" \
		-ttext 0.5 1.03 "Bessel function w/non-adaptive steppers" \
		-show

ex_ode_bessel.png: 
	o2graph -backend Agg -read ../ex_ode.o2 table_0 \
		-set xtitle x -set ytitle y -set top-margin 0.05 \
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

#ex_ode_bessel.eps: ex_ode_plot.py
#	python ex_ode_plot.py

ex_planar_plot.o: ex_planar_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_planar_plot.o -c ex_planar_plot.cpp

ex_planar_plot: ex_planar_plot.o ../../src/other/graph.o
	$(CXX) $(CPPFLAGS) -o ex_planar_plot ex_planar_plot.o \
		-lo2scl_hdf -lo2scl $(ROOT_LIB) $(LDFLAGS) 

ex_planar_plot.pdf: ex_planar_plot
	./ex_planar_plot
	convert ex_planar_plot.png ex_planar_plot.pdf

ex_stiff_plot.o: ex_stiff_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_stiff_plot.o -c ex_stiff_plot.cpp

ex_stiff_plot: ex_stiff_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_stiff_plot ex_stiff_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl \
		-lo2scl_hdf

ex_ublas.o: ex_ublas.cpp
	$(CXX) -I/home/asteiner/install/include $(INCLUDES) \
		-o ex_ublas.o -c ex_ublas.cpp

ex_ublas: ex_ublas.o  ../../src/other/graph.o
	$(CXX) -I/home/asteiner/install/include -L/usr/lib64/atlas \
		-L/home/asteiner/install/lib -o ex_ublas ex_ublas.o -lgsl \
		-lo2scl

