# --------------------------------------------------------------------
# Over the long term, I'd like to make sure that all figures can be
# automatically generated (using Imagemagick). As of 4/11/12, the
# following figures were successfully automatically generated:
#  
# complex.dox: 
#   ex_mandel_plot.png
#   ex_mandel_plot.pdf
#
# Others to fix
# funct.dox:
#   ex_fptr_plot.png
#   ex_fptr_plot.eps
# anneal.dox:
#   ex_anneal_plot.png
#   ex_anneal_plot.eps
# cheb.dox:
#   ex_chebapp_plot.png
#   ex_chebapp_plot.eps
# tintp.dox:
#   ex_planar_intp.png
#   ex_planar_intp.pdf
#   ex_contour_plot1.png
#   ex_contour_plot1.eps
#   ex_contour_plot2.png
#   ex_contour_plot2.eps
#   ex_contour_plot3.png
#   ex_contour_plot3.eps
# ode.dox:
#   ex_ode_bessel.png
#   ex_ode_bessel.eps
#   ex_ode_bessel2.png
#   ex_ode_bessel2.eps
#   ex_ode_bessel3.png
#   ex_ode_bessel3.eps
#   ex_ode_airy.png
#   ex_ode_airy.eps
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

test:
	echo $(CXX) 
	echo $(CPPFLAGS)

EXTRA_DIST = 

ROOT_LIB = ../../src/other/graph.o \
	-lgsl -lgslcblas -lreadline -lncurses \
	-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d \
	-lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
	-lMathCore -lThread -lz -pthread -lm -ldl

../../src/other/graph.o:
	cd ../../src/other; $(CXX) $(CPPFLAGS) -o graph.o -c graph.cpp

doc-extras: ex_fptr_plot.png ex_fptr_plot.eps \
	ex_mandel_plot_all ex_mandel_plot.eps 

# -------------------------------------------------------------------- 

# In approximate alphabetical order

ex_anneal_plot.o: ex_anneal_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_anneal_plot.o -c ex_anneal_plot.cpp

ex_anneal_plot: ex_anneal_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_anneal_plot ex_anneal_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

ex_chebapp_plot.o: ex_chebapp_plot.cpp
	$(CXX) -o ex_chebapp_plot.o -c ex_chebapp_plot.cpp

ex_chebapp_plot: ex_chebapp_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_chebapp_plot ex_chebapp_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl

ex_chebapp_plot.eps: ex_chebapp_plot
	./ex_chebapp_plot

ex_chebapp_plot.png: ex_chebapp_plot
	./ex_chebapp_plot

ex_contour_plot.o: ex_contour_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_contour_plot.o -c ex_contour_plot.cpp

ex_contour_plot: ex_contour_plot.o ../../src/other/graph.o
	$(CXX) $(CPPFLAGS) -o ex_contour_plot ex_contour_plot.o \
		-lo2scl_hdf -lo2scl $(ROOT_LIB) $(LDFLAGS) 

ex_fptr_plot.o: ex_fptr_plot.cpp
	$(CXX) -o ex_fptr_plot.o -c ex_fptr_plot.cpp

ex_fptr_plot: ex_fptr_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_fptr_plot ex_fptr_plot.o \
		$(ROOT_LIB) $(LDFLAGS) -lo2scl

ex_fptr_plot.eps: ex_fptr_plot
	./ex_fptr_plot

ex_fptr_plot.png: ex_fptr_plot
	./ex_fptr_plot

ex_mandel_plot.o: ex_mandel_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_mandel_plot.o -c ex_mandel_plot.cpp

ex_mandel_plot: ex_mandel_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_mandel_plot ex_mandel_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

ex_mandel_plot_all: ex_mandel_plot
	./ex_mandel_plot
	convert ex_mandel_plot.png ex_mandel_plot.pdf

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

ex_ode_plot.o: ex_ode_plot.cpp
	$(CXX) $(CPPFLAGS) -o ex_ode_plot.o -c ex_ode_plot.cpp

ex_ode_plot: ex_ode_plot.o  ../../src/other/graph.o
	$(CXX) -o ex_ode_plot ex_ode_plot.o \
		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl \
		-lo2scl_hdf

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

