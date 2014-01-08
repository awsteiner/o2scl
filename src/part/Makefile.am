# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

PART_SRCS = boson.cpp classical.cpp eff_boson.cpp eff_fermion.cpp \
	fermion.cpp nonrel_fermion.cpp part.cpp quark.cpp \
	rel_fermion.cpp sn_fermion.cpp sn_nr_fermion.cpp \
	sn_classical.cpp rel_boson.cpp mag_fermion_zerot.cpp
HEADER_VAR = boson.h classical.h eff_boson.h eff_fermion.h \
	fermion.h nonrel_fermion.h part.h quark.h rel_fermion.h \
	part_deriv.h sn_fermion.h sn_nr_fermion.h sn_classical.h \
	rel_boson.h mag_fermion_zerot.h
TEST_VAR = classical.scr eff_fermion.scr rel_fermion.scr boson.scr \
	fermion.scr nonrel_fermion.scr part.scr quark.scr \
	eff_boson.scr sn_fermion.scr sn_nr_fermion.scr \
	sn_classical.scr rel_boson.scr mag_fermion_zerot.scr

# ------------------------------------------------------------
# Includes
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

if O2SCL_CPP11
O2SCL_CPP11_MVAR = -std=c++0x
else
O2SCL_CPP11_MVAR = -DO2SCL_NO_CPP11
endif

AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_CPP11_MVAR) -DO2SCL_COND_FLAG

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

pkginclude_HEADERS = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

check_PROGRAMS = classical_ts eff_fermion_ts rel_fermion_ts \
	boson_ts fermion_ts nonrel_fermion_ts \
	part_ts quark_ts sn_fermion_ts mag_fermion_zerot_ts \
	sn_classical_ts sn_nr_fermion_ts rel_boson_ts eff_boson_ts

VCHECK_LIBS = $(libdir)/libo2scl_part.la $(libdir)/libo2scl.la \
	$(libdir)/libo2scl_hdf.la 

classical_ts_LDADD = $(VCHECK_LIBS)
eff_fermion_ts_LDADD = $(VCHECK_LIBS)
rel_fermion_ts_LDADD = $(VCHECK_LIBS)
boson_ts_LDADD = $(VCHECK_LIBS)
fermion_ts_LDADD = $(VCHECK_LIBS)
part_ts_LDADD = $(VCHECK_LIBS)
quark_ts_LDADD = $(VCHECK_LIBS)
sn_classical_ts_LDADD = $(VCHECK_LIBS)
nonrel_fermion_ts_LDADD = $(VCHECK_LIBS)
sn_fermion_ts_LDADD = $(VCHECK_LIBS)
sn_nr_fermion_ts_LDADD = $(VCHECK_LIBS)
rel_boson_ts_LDADD = $(VCHECK_LIBS)
eff_boson_ts_LDADD = $(VCHECK_LIBS)
mag_fermion_zerot_ts_LDADD = $(VCHECK_LIBS)

classical.scr: classical_ts$(EXEEXT) 
	./classical_ts$(EXEEXT) > classical.scr
eff_fermion.scr: eff_fermion_ts$(EXEEXT) 
	./eff_fermion_ts$(EXEEXT) > eff_fermion.scr
rel_fermion.scr: rel_fermion_ts$(EXEEXT) 
	./rel_fermion_ts$(EXEEXT) > rel_fermion.scr
boson.scr: boson_ts$(EXEEXT) 
	./boson_ts$(EXEEXT) > boson.scr
fermion.scr: fermion_ts$(EXEEXT) 
	./fermion_ts$(EXEEXT) > fermion.scr
part.scr: part_ts$(EXEEXT) 
	./part_ts$(EXEEXT) > part.scr
quark.scr: quark_ts$(EXEEXT) 
	./quark_ts$(EXEEXT) > quark.scr
sn_classical.scr: sn_classical_ts$(EXEEXT) 
	./sn_classical_ts$(EXEEXT) > sn_classical.scr
nonrel_fermion.scr: nonrel_fermion_ts$(EXEEXT) 
	./nonrel_fermion_ts$(EXEEXT) > nonrel_fermion.scr
sn_fermion.scr: sn_fermion_ts$(EXEEXT) 
	./sn_fermion_ts$(EXEEXT) > sn_fermion.scr
sn_nr_fermion.scr: sn_nr_fermion_ts$(EXEEXT) 
	./sn_nr_fermion_ts$(EXEEXT) > sn_nr_fermion.scr
rel_boson.scr: rel_boson_ts$(EXEEXT) 
	./rel_boson_ts$(EXEEXT) > rel_boson.scr
eff_boson.scr: eff_boson_ts$(EXEEXT) 
	./eff_boson_ts$(EXEEXT) > eff_boson.scr
mag_fermion_zerot.scr: mag_fermion_zerot_ts$(EXEEXT) 
	./mag_fermion_zerot_ts$(EXEEXT) > mag_fermion_zerot.scr

classical_ts_SOURCES = classical_ts.cpp
eff_fermion_ts_SOURCES = eff_fermion_ts.cpp
rel_fermion_ts_SOURCES = rel_fermion_ts.cpp
boson_ts_SOURCES = boson_ts.cpp
fermion_ts_SOURCES = fermion_ts.cpp
part_ts_SOURCES = part_ts.cpp
quark_ts_SOURCES = quark_ts.cpp
sn_classical_ts_SOURCES = sn_classical_ts.cpp
nonrel_fermion_ts_SOURCES = nonrel_fermion_ts.cpp
sn_fermion_ts_SOURCES = sn_fermion_ts.cpp
sn_nr_fermion_ts_SOURCES = sn_nr_fermion_ts.cpp
rel_boson_ts_SOURCES = rel_boson_ts.cpp
eff_boson_ts_SOURCES = eff_boson_ts.cpp
mag_fermion_zerot_ts_SOURCES = mag_fermion_zerot_ts.cpp

# ------------------------------------------------------------
# Library o2scl_basepart
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_basepart.la

libo2scl_basepart_la_SOURCES = $(PART_SRCS)

am__libo2scl_basepart_la_SOURCES_DIST = $(PART_SRCS)

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

