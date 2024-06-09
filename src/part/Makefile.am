# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

PART_SRCS = fermion_eff.cpp boson_eff.cpp \
	part.cpp quark.cpp fermion_deriv_nr.cpp \
	fermion_mag_zerot.cpp \
	part_deriv.cpp part_pdg.cpp part_python.cpp
HEADER_VAR = boson.h classical.h boson_eff.h fermion_eff.h \
	fermion.h fermion_nonrel.h part.h quark.h fermion_rel.h \
	part_deriv.h fermion_deriv_rel.h fermion_deriv_nr.h \
	boson_rel.h fermion_mag_zerot.h part_pdg.h \
	part_python.h classical_deriv.h fermion_rel2.h
TEST_VAR = classical.scr fermion_eff.scr fermion_rel.scr boson.scr \
	fermion.scr fermion_nonrel.scr part.scr quark.scr \
	boson_eff.scr fermion_deriv_rel.scr fermion_deriv_nr.scr \
	classical_deriv.scr boson_rel.scr fermion_mag_zerot.scr \
	fermion_rel2.scr

# ------------------------------------------------------------
# Includes
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

if O2SCL_MPFR
O2SCL_MPFR_MVAR = -DO2SCL_MPFR
else 
O2SCL_MPFR_MVAR = 
endif

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON -I$(PYTHON_INCLUDE_DIR)
else 
O2SCL_PYTHON_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(ADDINC) $(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_MPFR_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(ADDINC) $(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) \
	$(O2SCL_ARMA_MVAR) -DO2SCL_COND_FLAG $(O2SCL_MPFR_MVAR) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR) $(O2SCL_PYTHON_MVAR)
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = classical_ts fermion_eff_ts fermion_rel_ts \
	boson_ts fermion_ts fermion_nonrel_ts fermion_rel2_ts \
	part_ts quark_ts fermion_deriv_rel_ts fermion_mag_zerot_ts \
	classical_deriv_ts fermion_deriv_nr_ts boson_rel_ts boson_eff_ts

check_SCRIPTS = o2scl-test

if O2SCL_PYTHON

ADDL_TEST_LIBS = ../libo2scl.la -l$(PYTHON_LIB)
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR) -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = $(PYTHON_INCLUDE_DIR)
endif

else

ADDL_TEST_LIBS = ../libo2scl.la
if O2SCL_OPENMP
ADDL_TEST_LDFLGS = -fopenmp 
ADDL_TEST_FLGS = -DO2SCL_OPENMP
else
ADDL_TEST_LDFLGS = 
ADDL_TEST_FLGS = 
endif

endif

classical_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_eff_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_rel_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_rel2_ts_LDADD = $(ADDL_TEST_LIBS)
boson_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_ts_LDADD = $(ADDL_TEST_LIBS)
part_ts_LDADD = $(ADDL_TEST_LIBS)
quark_ts_LDADD = $(ADDL_TEST_LIBS)
classical_deriv_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_nonrel_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_deriv_rel_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_deriv_nr_ts_LDADD = $(ADDL_TEST_LIBS)
boson_rel_ts_LDADD = $(ADDL_TEST_LIBS)
boson_eff_ts_LDADD = $(ADDL_TEST_LIBS)
fermion_mag_zerot_ts_LDADD = $(ADDL_TEST_LIBS)

classical_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_eff_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_rel_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_rel2_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
boson_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
part_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
quark_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
classical_deriv_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_nonrel_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_deriv_rel_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_deriv_nr_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
boson_rel_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
boson_eff_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
fermion_mag_zerot_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

classical.scr: classical_ts$(EXEEXT) 
	./classical_ts$(EXEEXT) > classical.scr
fermion_eff.scr: fermion_eff_ts$(EXEEXT) 
	./fermion_eff_ts$(EXEEXT) > fermion_eff.scr
fermion_rel.scr: fermion_rel_ts$(EXEEXT) 
	./fermion_rel_ts$(EXEEXT) > fermion_rel.scr
fermion_rel2.scr: fermion_rel2_ts$(EXEEXT) 
	./fermion_rel2_ts$(EXEEXT) > fermion_rel2.scr
boson.scr: boson_ts$(EXEEXT) 
	./boson_ts$(EXEEXT) > boson.scr
fermion.scr: fermion_ts$(EXEEXT) 
	./fermion_ts$(EXEEXT) > fermion.scr
part.scr: part_ts$(EXEEXT) 
	./part_ts$(EXEEXT) > part.scr
quark.scr: quark_ts$(EXEEXT) 
	./quark_ts$(EXEEXT) > quark.scr
classical_deriv.scr: classical_deriv_ts$(EXEEXT) 
	./classical_deriv_ts$(EXEEXT) > classical_deriv.scr
fermion_nonrel.scr: fermion_nonrel_ts$(EXEEXT) 
	./fermion_nonrel_ts$(EXEEXT) > fermion_nonrel.scr
fermion_deriv_rel.scr: fermion_deriv_rel_ts$(EXEEXT) 
	./fermion_deriv_rel_ts$(EXEEXT) > fermion_deriv_rel.scr
fermion_deriv_nr.scr: fermion_deriv_nr_ts$(EXEEXT) 
	./fermion_deriv_nr_ts$(EXEEXT) > fermion_deriv_nr.scr
boson_rel.scr: boson_rel_ts$(EXEEXT) 
	./boson_rel_ts$(EXEEXT) > boson_rel.scr
boson_eff.scr: boson_eff_ts$(EXEEXT) 
	./boson_eff_ts$(EXEEXT) > boson_eff.scr
fermion_mag_zerot.scr: fermion_mag_zerot_ts$(EXEEXT) 
	./fermion_mag_zerot_ts$(EXEEXT) > fermion_mag_zerot.scr

classical_ts_SOURCES = classical_ts.cpp
fermion_eff_ts_SOURCES = fermion_eff_ts.cpp
fermion_rel_ts_SOURCES = fermion_rel_ts.cpp
fermion_rel2_ts_SOURCES = fermion_rel2_ts.cpp
boson_ts_SOURCES = boson_ts.cpp
fermion_ts_SOURCES = fermion_ts.cpp
part_ts_SOURCES = part_ts.cpp
quark_ts_SOURCES = quark_ts.cpp
classical_deriv_ts_SOURCES = classical_deriv_ts.cpp
fermion_nonrel_ts_SOURCES = fermion_nonrel_ts.cpp
fermion_deriv_rel_ts_SOURCES = fermion_deriv_rel_ts.cpp
fermion_deriv_nr_ts_SOURCES = fermion_deriv_nr_ts.cpp
boson_rel_ts_SOURCES = boson_rel_ts.cpp
boson_eff_ts_SOURCES = boson_eff_ts.cpp
fermion_mag_zerot_ts_SOURCES = fermion_mag_zerot_ts.cpp

# ------------------------------------------------------------

interface: empty
	../../bin/yanic part.i part_python ~/wcs/o2sclpy/o2sclpy/part \
		~/wcs/o2sclpy/doc/part  ../../bin/yanic_header.txt \
		| tee part_python.yanic.out

empty:

# ------------------------------------------------------------
# Library o2scl_part
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_part.la

libo2scl_part_la_SOURCES = $(PART_SRCS)

am__libo2scl_part_la_SOURCES_DIST = $(PART_SRCS)

# ------------------------------------------------------------
# Misc
# ------------------------------------------------------------

EXTRA_DIST = *_ts.cpp

emacs-clean: 
	-rm *~

qbkup:
	cp *.cpp *.h Makefile.am bkup

# ------------------------------------------------------------
# Testing
# ------------------------------------------------------------

o2scl-test: $(TEST_VAR)
	echo $(TEST_VAR) > testlist
	../internal/tsumm$(EXEEXT)

o2scl-test-rec: $(TEST_VAR)
	echo "src/part" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm *_ts.o *_ts *.scr

