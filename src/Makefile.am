help:
	@echo "This the O2scl makefile for the src subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply."
	@echo 

# ------------------------------------------------------------      
# Subdirectory and sublibrary lists
# -------------------------------------------------------------

# We also need to make sure
# 'internal' appears before the other subdirectories because 'tsumm'
# needs to be created first for 'make check'. 

SUBDIRS = internal anneal base deriv fit inte mcarlo \
	min ode other root linalg part nuclei eos hdf

LIBLIST = libo2scl.la

SUBLIBS = mcarlo/libo2scl_mcarlo.la other/libo2scl_other.la \
	min/libo2scl_min.la base/libo2scl_base.la \
	inte/libo2scl_inte.la linalg/libo2scl_linalg.la \
	hdf/libo2scl_hdf.la eos/libo2scl_eos.la \
	part/libo2scl_part.la nuclei/libo2scl_nuclei.la

#	fit/libo2scl_fit.la

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

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
endif

# ------------------------------------------------------------      
# Aggregate int libraries libo2scl and libo2scl_part
# -------------------------------------------------------------

lib_LTLIBRARIES = $(LIBLIST)

libo2scl_la_SOURCES = 
libo2scl_la_LIBADD = $(SUBLIBS)
libo2scl_la_LDFLAGS = -version-info $(OLIB_LT_VERSION)

# ------------------------------------------------------------

interface: empty
	cd base; $(MAKE) interface
	cd other; $(MAKE) interface
	cd hdf; $(MAKE) interface
	cd nuclei; $(MAKE) interface
	cd part; $(MAKE) interface
	cd eos; $(MAKE) interface

empty:


# ------------------------------------------------------------      
# Testing targets
# -------------------------------------------------------------

O2SCL_TEST_ADD = test_anneal test_base test_deriv test_eos test_fit \
	test_inte test_mcarlo test_min test_nuclei test_ode \
	test_other test_part test_root test_linalg \
	test_internal test_hdf

test_anneal:
	cd anneal && $(MAKE) o2scl-test-rec

test_base:
	cd base && $(MAKE) o2scl-test-rec

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
