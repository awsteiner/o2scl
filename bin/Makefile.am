help:
	@echo "This the O2scl makefile for the bin subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply."
	@echo 

# ------------------------------------------------------------      
# Include statement
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

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_FFTW_MVAR)
endif

# ------------------------------------------------------------      
# Program acol 
# ------------------------------------------------------------      

if O2SCL_ENABLE_ACOL

if O2SCL_HDF

bin_PROGRAMS = acol yanic ame_parse

acol_SOURCES = acol.cpp

ame_parse_SOURCES = ame_parse.cpp

yanic_SOURCES = yanic.cpp

if O2SCL_READLINE

# AWS: 5/13/21: commeting this out since it's already taken care
# of by configure.ac

#RLIB = -lreadline -lncurses
RLIB =

else

RLIB =

endif

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

# AWS: 5/13/21: -fopenmp is both here and in RLIB above, but I've
# had errors on ubuntu if it's not here as well

acol_LDFLAGS = -L../src -L../src/hdf -fopenmp
ame_parse_LDFLAGS = -L../src -L../src/hdf -fopenmp
yanic_LDFLAGS = -L../src -L../src/hdf -fopenmp

else

acol_LDFLAGS = -L../src -L../src/hdf 
ame_parse_LDFLAGS = -L../src -L../src/hdf 
yanic_LDFLAGS = -L../src -L../src/hdf 

endif

endif

acol_LDADD = -lo2scl_hdf -lo2scl -lgsl -lgslcblas $(RLIB)
ame_parse_LDADD = -lo2scl_part -lo2scl_hdf -lo2scl -lgsl -lgslcblas $(RLIB)
yanic_LDADD = -lo2scl_hdf -lo2scl -lgsl -lgslcblas $(RLIB)

endif

endif

