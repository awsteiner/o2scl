
help:
	@echo "This directory contains code which assists the developer"
	@echo "and is not intended for use by the end-user."
	@echo
	@echo "File descriptions:"
	@echo "--------------------------------------------------------"

agraph.o: agraph.cpp
	$(CXX) $(CPPFLAGS) \
		-DO2SCL_COND_FLAG \
		-I../../include -o agraph.o -c agraph.cpp

agraphm.o: agraphm.cpp
	$(CXX) $(CPPFLAGS) \
		-DO2SCL_COND_FLAG \
	-I../../include -o agraphm.o -c agraphm.cpp

# Most dependent libraries first: ROOT, O2scl, then readline
AGRAPH_LIB = -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d \
        -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
        -lMathCore -lThread -lz -pthread -ldl -lo2scl_hdf \
	-lo2scl -lreadline -lncurses

agraph: agraph.o agraphm.o
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -L@libdir@ -I../../include \
		-DO2SCL_COND_FLAG \
		-o agraph agraph.o agraphm.o ../acolm.o \
		../other/graph.o $(LIBS) $(AGRAPH_LIB) 

update_lic: update_lic.cpp
	g++ -o update_lic update_lic.cpp

EXTRA_DIST = binsiz.f agraph.cpp agraphm.cpp agraphm.h \
	gesvd.cpp zheev.cpp ame_parse.cpp hfb_parse.cpp mnk_parse.cpp \
	o2comm.cpp osni88.cpp reparse.cpp make_units.cpp

# contour3.pdf contour4.pdf contour_gr.cpp
# update_lic.cpp verup.cpp

gesvd.o: gesvd.cpp
	$(CXX) -I$(HOME)/install/include $(INCLUDES) \
		-o gesvd.o -c gesvd.cpp

gesvd: gesvd.o 
	$(CXX) -I$(HOME)/install/include -L/usr/lib64/atlas \
		-L$(HOME)/install/lib -o gesvd gesvd.o -lgsl \
		-lgslcblas -llapack -lcblas -lblas -latlas -lf77blas 

zheev.o: zheev.cpp
	$(CXX) -I$(HOME)/install/include $(INCLUDES) \
		-o zheev.o -c zheev.cpp

zheev: zheev.o 
	$(CXX) -I$(HOME)/install/include -L/usr/lib64/atlas \
		-L$(HOME)/install/lib -o zheev zheev.o -lgsl \
		-llapack -lblas -lgslcblas -lgfortran -lo2scl

ex_ublas.o: ex_ublas.cpp
	$(CXX) -I$(HOME)/install/include $(INCLUDES) \
		-o ex_ublas.o -c ex_ublas.cpp

ex_ublas: ex_ublas.o 
	$(CXX) -I$(HOME)/install/include -L/usr/lib64/atlas \
		-L$(HOME)/install/lib -o ex_ublas ex_ublas.o -lgsl \
		-lo2scl

bm_lublas.o: bm_lublas.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include $(INCLUDES) \
		-o bm_lublas.o -c bm_lublas.cpp

bm_lublas: bm_lublas.o 
	$(CXX) -O3 -static -I$(HOME)/install/include \
		-L$(HOME)/install/lib -o bm_lublas bm_lublas.o \
		-lo2scl -lgsl -lgslcblas 

ame_parse.o: ame_parse.cpp
	$(CXX) $(CPPFLAGS) -o ame_parse.o -c ame_parse.cpp

ame_parse: ame_parse.o 
	$(CXX) $(CPPFLAGS) -o ame_parse ame_parse.o $(LDFLAGS) \
		-lo2scl_part -lo2scl -lgsl -lgslcblas -lo2scl_hdf -lhdf5 \
		-lhdf5_hl

hfb_parse.o: hfb_parse.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include/hdf5 \
		-I$(HOME)/install/include $(INCLUDES) \
		-o hfb_parse.o -c hfb_parse.cpp

hfb_parse: hfb_parse.o 
	$(CXX) -O3 -I$(HOME)/install/include \
		-I$(HOME)/install/include/hdf5 \
		-L$(HOME)/install/lib -o hfb_parse hfb_parse.o \
		-lo2scl_part -lo2scl -lgsl -lgslcblas -lo2scl_hdf -lhdf5 \
		-lhdf5_hl

mnk_parse.o: mnk_parse.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include/hdf5 \
		-I$(HOME)/install/include $(INCLUDES) \
		-o mnk_parse.o -c mnk_parse.cpp

mnk_parse: mnk_parse.o 
	$(CXX) -O3 -I$(HOME)/install/include \
		-I$(HOME)/install/include/hdf5 \
		-L$(HOME)/install/lib -o mnk_parse mnk_parse.o \
		-lo2scl_part -lo2scl -lgsl -lgslcblas -lo2scl_hdf -lhdf5 \
		-lhdf5_hl

reparse.o: reparse.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include/hdf5 \
		-I$(HOME)/install/include $(INCLUDES) \
		-o reparse.o -c reparse.cpp

reparse: reparse.o 
	$(CXX) -O3 -I$(HOME)/install/include \
		-I$(HOME)/install/include/hdf5 \
		-L$(HOME)/install/lib -o reparse reparse.o \
		-lo2scl_part -lo2scl -lgsl -lgslcblas -lhdf5

o2comm.o: o2comm.cpp
	$(CXX) $(CPPFLAGS) -I../../include -DO2SCL_COND_FLAG \
		-o o2comm.o -c o2comm.cpp

o2comm: o2comm.o 
	$(CXX) $(LDFLAGS) $(CPPFLAGS) -L@libdir@ -o o2comm o2comm.o \
		-lo2scl_hdf -lo2scl -lgsl -lgslcblas -lhdf5

osni88.o: osni88.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include/hdf5 \
		-I$(HOME)/install/include $(INCLUDES) \
		-o osni88.o -c osni88.cpp

osni88: osni88.o 
	$(CXX) $(LDFLAGS) -O3 -I$(HOME)/install/include \
		-I$(HOME)/install/include/hdf5 \
		-L$(HOME)/install/lib -o osni88 osni88.o \
		-lo2scl_hdf -lo2scl -lgsl -lgslcblas -lhdf5

make_units.o: make_units.cpp
	$(CXX) -O3 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
		-I$(HOME)/install/include/hdf5 \
		-I$(HOME)/install/include $(INCLUDES) \
		-o make_units.o -c make_units.cpp

make_units: make_units.o 
	$(CXX) -O3 -I$(HOME)/install/include \
		-I$(HOME)/install/include/hdf5 \
		-L$(HOME)/install/lib -o make_units make_units.o \
		-lo2scl_part -lo2scl -lo2scl_hdf -lgsl -lgslcblas -lhdf5

# ------------------------------------------------

check_PROGRAMS = tsumm

tsumm_SOURCES = tsumm.cpp

o2scl-test-rec: tsumm

# ------------------------------------------------
# Old code no longer used
# ------------------------------------------------

# partmath: partmath.cpp
# 	g++ -o partmath -I$(HOME)/install/include \
# 		-L$(HOME)/install/lib \
# 		partmath.cpp -lo2scl

# contour_gr.o: contour_gr.cpp
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include \
# 		-o contour_gr.o -c contour_gr.cpp

# contour_gr: contour_gr.o $(HOME)/svn/o2int/trunk/src/other/graph.o
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include -L$(HOME)/root/lib \
# 		-L$(HOME)/install/lib \
# 		-o contour_gr contour_gr.o \
# 		$(HOME)/svn/o2int/trunk/src/other/graph.o \
# 		$(ROOT_LIB) -lo2scl_base  \
# 		-lo2scl_other -lo2scl_rnga 

# contour_gr2.o: contour_gr2.cpp
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include \
# 		-o contour_gr2.o -c contour_gr2.cpp

# contour_gr2: contour_gr2.o $(HOME)/svn/o2int/trunk/src/other/graph.o
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include -L$(HOME)/root/lib \
# 		-L$(HOME)/install/lib \
# 		-o contour_gr2 contour_gr2.o \
# 		$(HOME)/svn/o2int/trunk/src/other/graph.o \
# 		$(ROOT_LIB) -lo2scl_base  \
# 		-lo2scl_other -lo2scl_rnga 

# mplot.o: mplot.cpp
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include \
# 		-o mplot.o -c mplot.cpp

# mplot: mplot.o $(HOME)/svn/o2int/trunk/src/other/graph.o
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include -L$(HOME)/root/lib \
# 		-L$(HOME)/install/lib \
# 		-o mplot mplot.o \
# 		$(HOME)/svn/o2int/trunk/src/other/graph.o \
# 		$(ROOT_LIB) -lo2scl  

# in_planar_intp.o: in_planar_intp.cpp
# 	$(CXX) -I$(HOME)/install/include \
# 		-I$(HOME)/root/include \
# 		-o in_planar_intp.o -c in_planar_intp.cpp

# in_planar_intp: in_planar_intp.o 
# 	$(CXX) -I$(HOME)/install/include \
# 		-L$(HOME)/install/lib -o in_planar_intp \
# 		in_planar_intp.o -L$(HOME)/root/lib \
# 		$(HOME)/svn/o2int/trunk/src/other/graph.o \
# 		-lo2scl_base -lo2scl_other -lo2scl_rnga  \
# 		$(ROOT_LIB)

