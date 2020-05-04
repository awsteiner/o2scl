# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

BASE_HEADER_VAR =  err_hnd.h misc.h test_mgr.h table.h \
	constants.h funct.h multi_funct.h mm_funct.h \
	lib_settings.h search_vec.h interp.h \
	cli.h columnify.h convert_units.h string_conv.h \
	tensor.h vector.h table3d.h cli_readline.h tensor_grid.h \
	format_float.h table_units.h exception.h uniform_grid.h \
	shunting_yard.h interp_krige.h find_constants.h cursesw.h

HEADER_VAR = $(BASE_HEADER_VAR)

BASE_BASE_SRCS = err_hnd.cpp table3d.cpp \
	lib_settings.cpp misc.cpp cli.cpp \
	test_mgr.cpp vector.cpp \
	string_conv.cpp exception.cpp format_float.cpp \
	shunting_yard.cpp tensor.cpp find_constants.cpp cursesw.cpp

BASE_SRCS = $(BASE_BASE_SRCS)

BASE_TEST_VAR = err_hnd.scr funct.scr table3d.scr \
	lib_settings.scr misc.scr mm_funct.scr multi_funct.scr \
	search_vec.scr table.scr vector.scr interp_krige.scr \
	interp.scr columnify.scr convert_units.scr \
	string_conv.scr tensor.scr shunting_yard.scr \
	format_float.scr table_units.scr exception.scr uniform_grid.scr \
	tensor_grid.scr constants.scr cursesw.scr

TEST_VAR = $(BASE_TEST_VAR)

SUBDIRS = arma eigen both neither

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

if O2SCL_GSL2
O2SCL_GSL2_MVAR = -DO2SCL_GSL2
else
O2SCL_GSL2_MVAR =
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

if O2SCL_READLINE
ADDINC = -DO2SCL_READLINE
else 
ADDINC = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(ADDINC) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
else
AM_CPPFLAGS = -I@top_srcdir@/include/ -DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	$(O2SCL_PART_MVAR) $(O2SCL_EOS_MVAR) $(O2SCL_HDF_MVAR) \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(ADDINC) \
	$(O2SCL_OSX_MVAR) $(O2SCL_LINUX_MVAR)
endif


# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

pkginclude_HEADERS = $(HEADER_VAR)

am__pkginclude_HEADERS_DIST = $(HEADER_VAR)

# ------------------------------------------------------------
# libtool testing targets
# ------------------------------------------------------------

CPVAR = err_hnd_ts funct_ts convert_units_ts \
	lib_settings_ts misc_ts mm_funct_ts constants_ts \
	multi_funct_ts search_vec_ts table_ts interp_ts \
	columnify_ts shunting_yard_ts interp_krige_ts \
	string_conv_ts tensor_ts tensor_grid_ts vector_ts table3d_ts \
	format_float_ts table_units_ts exception_ts uniform_grid_ts \
	cursesw_ts

check_PROGRAMS = $(CPVAR)

check_SCRIPTS = o2scl-test

if O2SCL_HDF
VCHECK_LIBS = ../hdf/libo2scl_hdf.la ../libo2scl.la
else
VCHECK_LIBS = ../libo2scl.la
endif

if O2SCL_OPENMP
interp_krige_ts_LDFLAGS = -fopenmp
constants_ts_LDFLAGS = -fopenmp
err_hnd_ts_LDFLAGS = -fopenmp
convert_units_ts_LDFLAGS = -fopenmp
funct_ts_LDFLAGS = -fopenmp
format_float_ts_LDFLAGS = -fopenmp
lib_settings_ts_LDFLAGS = -fopenmp
misc_ts_LDFLAGS = -fopenmp
mm_funct_ts_LDFLAGS = -fopenmp
multi_funct_ts_LDFLAGS = -fopenmp
search_vec_ts_LDFLAGS = -fopenmp
table_ts_LDFLAGS = -fopenmp
table3d_ts_LDFLAGS = -fopenmp
vector_ts_LDFLAGS = -fopenmp
interp_ts_LDFLAGS = -fopenmp
columnify_ts_LDFLAGS = -fopenmp
cursesw_ts_LDFLAGS = -fopenmp
string_conv_ts_LDFLAGS = -fopenmp
tensor_ts_LDFLAGS = -fopenmp
tensor_grid_ts_LDFLAGS = -fopenmp
table_units_ts_LDFLAGS = -fopenmp
exception_ts_LDFLAGS = -fopenmp
uniform_grid_ts_LDFLAGS = -fopenmp
shunting_yard_ts_LDFLAGS = -fopenmp
endif

interp_krige_ts_LDADD = $(VCHECK_LIBS)
constants_ts_LDADD = $(VCHECK_LIBS)
err_hnd_ts_LDADD = $(VCHECK_LIBS)
convert_units_ts_LDADD = $(VCHECK_LIBS)
funct_ts_LDADD = $(VCHECK_LIBS)
format_float_ts_LDADD = $(VCHECK_LIBS)
lib_settings_ts_LDADD = $(VCHECK_LIBS)
misc_ts_LDADD = $(VCHECK_LIBS)
mm_funct_ts_LDADD = $(VCHECK_LIBS)
multi_funct_ts_LDADD = $(VCHECK_LIBS)
search_vec_ts_LDADD = $(VCHECK_LIBS)
table_ts_LDADD = $(VCHECK_LIBS)
table3d_ts_LDADD = $(VCHECK_LIBS)
vector_ts_LDADD = $(VCHECK_LIBS)
interp_ts_LDADD = $(VCHECK_LIBS)
columnify_ts_LDADD = $(VCHECK_LIBS)
cursesw_ts_LDADD = $(VCHECK_LIBS)
string_conv_ts_LDADD = $(VCHECK_LIBS)
tensor_ts_LDADD = $(VCHECK_LIBS)
tensor_grid_ts_LDADD = $(VCHECK_LIBS)
table_units_ts_LDADD = $(VCHECK_LIBS)
exception_ts_LDADD = $(VCHECK_LIBS)
uniform_grid_ts_LDADD = $(VCHECK_LIBS)
shunting_yard_ts_LDADD = $(VCHECK_LIBS)

interp_krige.scr: interp_krige_ts$(EXEEXT) 
	./interp_krige_ts$(EXEEXT) > interp_krige.scr
constants.scr: constants_ts$(EXEEXT) 
	./constants_ts$(EXEEXT) > constants.scr
err_hnd.scr: err_hnd_ts$(EXEEXT) 
	./err_hnd_ts$(EXEEXT) > err_hnd.scr
convert_units.scr: convert_units_ts$(EXEEXT) 
	./convert_units_ts$(EXEEXT) > convert_units.scr
funct.scr: funct_ts$(EXEEXT) 
	./funct_ts$(EXEEXT) > funct.scr
format_float.scr: format_float_ts$(EXEEXT) 
	./format_float_ts$(EXEEXT) > format_float.scr
lib_settings.scr: lib_settings_ts$(EXEEXT) 
	./lib_settings_ts$(EXEEXT) > lib_settings.scr
misc.scr: misc_ts$(EXEEXT) 
	./misc_ts$(EXEEXT) > misc.scr
mm_funct.scr: mm_funct_ts$(EXEEXT) 
	./mm_funct_ts$(EXEEXT) > mm_funct.scr
multi_funct.scr: multi_funct_ts$(EXEEXT) 
	./multi_funct_ts$(EXEEXT) > multi_funct.scr
search_vec.scr: search_vec_ts$(EXEEXT) 
	./search_vec_ts$(EXEEXT) > search_vec.scr
table.scr: table_ts$(EXEEXT) 
	./table_ts$(EXEEXT) > table.scr
table3d.scr: table3d_ts$(EXEEXT) 
	./table3d_ts$(EXEEXT) > table3d.scr
vector.scr: vector_ts$(EXEEXT) 
	./vector_ts$(EXEEXT) > vector.scr
interp.scr: interp_ts$(EXEEXT) 
	./interp_ts$(EXEEXT) > interp.scr
columnify.scr: columnify_ts$(EXEEXT) 
	./columnify_ts$(EXEEXT) > columnify.scr
cursesw.scr: cursesw_ts$(EXEEXT) 
	./cursesw_ts$(EXEEXT) > cursesw.scr
string_conv.scr: string_conv_ts$(EXEEXT) 
	./string_conv_ts$(EXEEXT) > string_conv.scr
tensor.scr: tensor_ts$(EXEEXT) 
	./tensor_ts$(EXEEXT) > tensor.scr
tensor_grid.scr: tensor_grid_ts$(EXEEXT) 
	./tensor_grid_ts$(EXEEXT) > tensor_grid.scr
table_units.scr: table_units_ts$(EXEEXT) 
	./table_units_ts$(EXEEXT) > table_units.scr
exception.scr: exception_ts$(EXEEXT) 
	./exception_ts$(EXEEXT) > exception.scr
uniform_grid.scr: uniform_grid_ts$(EXEEXT) 
	./uniform_grid_ts$(EXEEXT) > uniform_grid.scr
shunting_yard.scr: shunting_yard_ts$(EXEEXT) 
	./shunting_yard_ts$(EXEEXT) > shunting_yard.scr

interp_krige_ts_SOURCES = interp_krige_ts.cpp
constants_ts_SOURCES = constants_ts.cpp
err_hnd_ts_SOURCES = err_hnd_ts.cpp
convert_units_ts_SOURCES = convert_units_ts.cpp
funct_ts_SOURCES = funct_ts.cpp
format_float_ts_SOURCES = format_float_ts.cpp
lib_settings_ts_SOURCES = lib_settings_ts.cpp
misc_ts_SOURCES = misc_ts.cpp
mm_funct_ts_SOURCES = mm_funct_ts.cpp
multi_funct_ts_SOURCES = multi_funct_ts.cpp
search_vec_ts_SOURCES = search_vec_ts.cpp
table_ts_SOURCES = table_ts.cpp
table3d_ts_SOURCES = table3d_ts.cpp
vector_ts_SOURCES = vector_ts.cpp
interp_ts_SOURCES = interp_ts.cpp
columnify_ts_SOURCES = columnify_ts.cpp
cursesw_ts_SOURCES = cursesw_ts.cpp
string_conv_ts_SOURCES = string_conv_ts.cpp
tensor_ts_SOURCES = tensor_ts.cpp
tensor_grid_ts_SOURCES = tensor_grid_ts.cpp
table_units_ts_SOURCES = table_units_ts.cpp
exception_ts_SOURCES = exception_ts.cpp
uniform_grid_ts_SOURCES = uniform_grid_ts.cpp
shunting_yard_ts_SOURCES = shunting_yard_ts.cpp

# ------------------------------------------------------------
# Library o2scl_base
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_base.la

libo2scl_base_la_SOURCES = $(BASE_SRCS)

am__libo2scl_base_la_SOURCES_DIST = $(BASE_SRCS)

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
	echo "src/base" $(TEST_VAR) >> ../../testlist

test-clean:
	-rm -f *_ts.o *_ts *.scr

