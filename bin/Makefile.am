help:
	@echo "This the O2scl makefile for the bin subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply."
	@echo 

EXTRA_DIST = yanic_header.txt

# ------------------------------------------------------------      
# Include statement
# ------------------------------------------------------------      

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

if O2SCL_FFTW
O2SCL_FFTW_MVAR = -DO2SCL_FFTW
else 
O2SCL_FFTW_MVAR = 
endif

if O2SCL_PUGIXML
O2SCL_PUGIXML_MVAR = -DO2SCL_PUGIXML
else 
O2SCL_PUGIXML_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(O2SCL_PUGIXML_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR)
endif

# ------------------------------------------------------------      

temp:
	./make_rst ../doc/o2scl/xml/index.xml 

# ------------------------------------------------------------      
# Program acol 
# ------------------------------------------------------------      

if O2SCL_ENABLE_ACOL

bin_PROGRAMS = acol yanic ame_parse make_rst xml_acol

acol_SOURCES = acol.cpp

ame_parse_SOURCES = ame_parse.cpp

yanic_SOURCES = yanic.cpp

make_rst_SOURCES = make_rst.cpp

xml_acol_SOURCES = xml_acol.cpp

# We look first for the libraries in the distribution directory.
# However, homebrew's superenv appears to move the libraries around,
# so we look for them in pkglibdir. Removing -fopenmp here since
# it's already taken care of above

if O2SCL_SNAP

#if O2SCL_OPENMP
#acol_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
#	-Wl,/snap/o2scl/current/lib -fopenmp
#ame_parse_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
#	-Wl,/snap/o2scl/current/lib -fopenmp
#yanic_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
#	-Wl,/snap/o2scl/current/lib -fopenmp
#else

acol_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
	-Wl,/snap/o2scl/current/lib
ame_parse_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
	-Wl,/snap/o2scl/current/lib
yanic_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
	-Wl,/snap/o2scl/current/lib
make_rst_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
	-Wl,/snap/o2scl/current/lib
xml_acol_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -Wl,-rpath \
	-Wl,/snap/o2scl/current/lib

#endif

else

# Removed -L${pkglibdir} on 5/11/21 because it was causing the following
# error, 
# 
# ld: warning: directory not found for option '-L/usr/local/lib/o2scl'
#
# but I bet it was there to assist with creating an ubuntu
# package...

if O2SCL_OPENMP

#acol_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -fopenmp
#ame_parse_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -fopenmp
#yanic_LDFLAGS = -L../src -L../src/hdf -L${pkglibdir} -fopenmp

if O2SCL_PYTHON

acol_LDFLAGS = -L../src -L../src/hdf -fopenmp -l$(PYTHON_LIB)
ame_parse_LDFLAGS = -L../src -L../src/hdf -fopenmp -l$(PYTHON_LIB)
yanic_LDFLAGS = -L../src -L../src/hdf -fopenmp -l$(PYTHON_LIB)
make_rst_LDFLAGS = -L../src -L../src/hdf -fopenmp -l$(PYTHON_LIB)
xml_acol_LDFLAGS = -L../src -L../src/hdf -fopenmp -l$(PYTHON_LIB)

else

acol_LDFLAGS = -L../src -L../src/hdf -fopenmp
ame_parse_LDFLAGS = -L../src -L../src/hdf -fopenmp
yanic_LDFLAGS = -L../src -L../src/hdf -fopenmp
make_rst_LDFLAGS = -L../src -L../src/hdf -fopenmp
xml_acol_LDFLAGS = -L../src -L../src/hdf -fopenmp

endif

else

acol_LDFLAGS = -L../src -L../src/hdf -l$(PYTHON_LIB)
ame_parse_LDFLAGS = -L../src -L../src/hdf -l$(PYTHON_LIB)
yanic_LDFLAGS = -L../src -L../src/hdf -l$(PYTHON_LIB)
make_rst_LDFLAGS = -L../src -L../src/hdf -l$(PYTHON_LIB)
xml_acol_LDFLAGS = -L../src -L../src/hdf -l$(PYTHON_LIB)

endif

endif

if O2SCL_PYTHON

acol_LDADD = -lo2scl -lgsl -lgslcblas -l$(PYTHON_LIB)
ame_parse_LDADD = -lo2scl -lgsl -lgslcblas -l$(PYTHON_LIB)
yanic_LDADD = -lo2scl -lgsl -lgslcblas -l$(PYTHON_LIB)
make_rst_LDADD = -lo2scl -lgsl -lgslcblas -l$(PYTHON_LIB)
xml_acol_LDADD = -lo2scl -lgsl -lgslcblas -l$(PYTHON_LIB)

else

acol_LDADD = -lo2scl -lgsl -lgslcblas
ame_parse_LDADD = -lo2scl -lgsl -lgslcblas
yanic_LDADD = -lo2scl -lgsl -lgslcblas
make_rst_LDADD = -lo2scl -lgsl -lgslcblas
xml_acol_LDADD = -lo2scl -lgsl -lgslcblas

endif

endif

