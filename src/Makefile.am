help:
	@echo "This the O2scl makefile for the src subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply."
	@echo 

# ------------------------------------------------------------      
# Subdirectory and sublibrary lists
# -------------------------------------------------------------

# The ordering is important for these directories, since part, nuclei,
# and eos, depend on the base library. We also need to make sure
# 'internal' appears before the other subdirectories because 'tsumm'
# needs to be created first for 'make check'. We don't need to list
# libo2scl_eos and libo2scl_hdf here, since they're taken care of in
# their respective directories.

if O2SCL_PARTLIB
if O2SCL_EOSLIB
if O2SCL_HDF

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg \
	part nuclei eos hdf cubature
LIBLIST = libo2scl.la libo2scl_part.la

else

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg \
	part nuclei eos cubature
LIBLIST = libo2scl.la libo2scl_part.la

endif
else

if O2SCL_HDF

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg \
	nuclei part hdf cubature
LIBLIST = libo2scl.la libo2scl_part.la

else

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg \
	nuclei part cubature
LIBLIST = libo2scl.la libo2scl_part.la

endif
endif
else 

if O2SCL_HDF

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg hdf cubature
LIBLIST = libo2scl.la

else

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg cubature
LIBLIST = libo2scl.la

endif
endif

SUBLIBS = mcarlo/libo2scl_mcarlo.la other/libo2scl_other.la \
	min/libo2scl_min.la \
	inte/libo2scl_inte.la linalg/libo2scl_linalg.la \
	base/libo2scl_base.la

#	fit/libo2scl_fit.la

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

AM_CPPFLAGS = -I../include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	-DO2SCL_COND_FLAG

# ------------------------------------------------------------      
# Aggregate int libraries libo2scl and libo2scl_part
# -------------------------------------------------------------

lib_LTLIBRARIES = $(LIBLIST)

libo2scl_la_SOURCES = 
libo2scl_la_LIBADD = $(SUBLIBS)
libo2scl_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

libo2scl_part_la_SOURCES = 
libo2scl_part_la_LIBADD = nuclei/libo2scl_nuclei.la part/libo2scl_basepart.la
libo2scl_part_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

# ------------------------------------------------------------      
# Program acol 
# ------------------------------------------------------------      

if O2SCL_ENABLE_ACOL

if O2SCL_HDF

bin_PROGRAMS = acol 

acol_SOURCES = acol.cpp 

if O2SCL_READLINE
RLIB = -lreadline -lncurses
else 
RLIB = 
endif

# We look first for the libraries in the distribution directory.
# However, homebrew's superenv appears to move the libraries around,
# so we look for them in pkglibdir
if O2SCL_SNAP
acol_LDFLAGS = -L. -Lhdf -L${pkglibdir} -Wl,-rpath -Wl,/snap/o2scl/current/lib
else
acol_LDFLAGS = -L. -Lhdf -L${pkglibdir}
endif

acol_LDADD = -lo2scl_hdf -lo2scl -lgsl -lgslcblas $(RLIB)

endif

endif

# ------------------------------------------------------------      
# Testing targets
# -------------------------------------------------------------

if O2SCL_PARTLIB
if O2SCL_EOSLIB
if O2SCL_HDF
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_eos test_fit \
	test_inte test_mcarlo test_min test_nuclei test_ode \
	test_other test_part test_root test_linalg \
	test_internal test_hdf test_cubature
else
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_eos test_fit \
	test_inte test_mcarlo test_min test_nuclei test_ode \
	test_other test_part test_root test_linalg \
	test_internal test_cubature
endif
else
if O2SCL_HDF
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_fit \
	test_inte test_mcarlo test_min test_nuclei test_ode \
	test_other test_part test_root test_linalg \
	test_internal test_hdf test_cubature
else
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_fit \
	test_inte test_mcarlo test_min test_nuclei test_ode \
	test_other test_part test_root test_linalg \
	test_internal test_cubature
endif
endif
else 
if O2SCL_HDF
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_fit \
	test_inte test_mcarlo test_min test_ode \
	test_other test_root test_linalg \
	test_internal test_hdf test_cubature
else
O2SCL_TEST_ADD = test_anneal test_base test_deriv test_fit \
	test_inte test_mcarlo test_min test_ode \
	test_other test_root test_linalg \
	test_internal test_cubature
endif
endif

test_anneal:
	cd anneal && $(MAKE) o2scl-test-rec

test_base:
	cd base && $(MAKE) o2scl-test-rec

test_cubature:
	cd cubature && $(MAKE) o2scl-test-rec

test_deriv:
	cd deriv && $(MAKE) o2scl-test-rec

test_eos:
	cd eos && $(MAKE) o2scl-test-rec

test_fit:
	cd fit && $(MAKE) o2scl-test-rec

test_hdf:
	cd hdf && $(MAKE) o2scl-test-rec

test_inte:
	cd inte && $(MAKE) o2scl-test-rec

test_internal:
	cd internal && $(MAKE) o2scl-test-rec

test_linalg:
	cd linalg && $(MAKE) o2scl-test-rec

test_mcarlo:
	cd mcarlo && $(MAKE) o2scl-test-rec

test_min:
	cd min && $(MAKE) o2scl-test-rec

test_nuclei:
	cd nuclei && $(MAKE) o2scl-test-rec

test_ode:
	cd ode && $(MAKE) o2scl-test-rec

test_other:
	cd other && $(MAKE) o2scl-test-rec

test_part:
	cd part && $(MAKE) o2scl-test-rec

test_root:
	cd root && $(MAKE) o2scl-test-rec

o2scl-test-rec: $(O2SCL_TEST_ADD)
