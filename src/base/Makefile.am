# ------------------------------------------------------------
# Basic variables
# ------------------------------------------------------------

BASE_HEADER_VAR =  err_hnd.h misc.h test_mgr.h table.h \
	constants.h funct.h multi_funct.h mm_funct.h \
	lib_settings.h search_vec.h interp.h \
	cli.h columnify.h convert_units.h string_conv.h \
	tensor.h vector.h table3d.h cli_readline.h tensor_grid.h \
	format_float.h table_units.h exception.h uniform_grid.h \
	interp_krige.h find_constants.h cursesw.h \
	prev_commit.h auto_format.h base_python.h calc_utf8.h \
	funct_multip.h interp_vec.h funct_to_fp.h

HEADER_VAR = $(BASE_HEADER_VAR)

BASE_BASE_SRCS = err_hnd.cpp table3d.cpp \
	lib_settings.cpp misc.cpp cli.cpp \
	test_mgr.cpp vector.cpp auto_format.cpp \
	string_conv.cpp exception.cpp format_float.cpp \
	tensor.cpp cursesw.cpp \
	base_python.cpp funct.cpp funct_to_fp.cpp

BASE_SRCS = $(BASE_BASE_SRCS)

BASE_TEST_VAR = err_hnd.scr funct.scr table3d.scr \
	lib_settings.scr misc.scr mm_funct.scr multi_funct.scr \
	search_vec.scr table.scr vector.scr interp_krige.scr \
	interp.scr columnify.scr convert_units.scr \
	string_conv.scr tensor.scr \
	format_float.scr table_units.scr exception.scr uniform_grid.scr \
	tensor_grid.scr constants.scr cursesw.scr auto_format.scr \
	calc_utf8.scr

TEST_VAR = $(BASE_TEST_VAR)

SUBDIRS = arma eigen both neither python

# ------------------------------------------------------------
# Includes
# ------------------------------------------------------------

if O2SCL_OPENMP
O2SCL_OPENMP_MVAR = -DO2SCL_OPENMP -fopenmp
else
O2SCL_OPENMP_MVAR =
endif

if O2SCL_READLINE
O2SCL_READLINE_MVAR = -DO2SCL_READLINE
else
O2SCL_READLINE_MVAR =
endif

if O2SCL_NCURSES
O2SCL_NCURSES_MVAR = -DO2SCL_NCURSES
else
O2SCL_NCURSES_MVAR =
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

if O2SCL_CUBATURE
O2SCL_CUBATURE_MVAR = -DO2SCL_CUBATURE
else 
O2SCL_CUBATURE_MVAR = 
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

if O2SCL_READLINE
ADDINC = -DO2SCL_READLINE
else 
ADDINC = 
endif

if O2SCL_MODULES
O2SCL_MODULES_MVAR = -DO2SCL_MODULES -fmodules-ts
else 
O2SCL_MODULES_MVAR = 
endif

if O2SCL_PUGIXML
O2SCL_PUGIXML_MVAR = -DO2SCL_PUGIXML
else 
O2SCL_PUGIXML_MVAR = 
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	-DO2SCL_DOC_DIR=\"/snap/o2scl/current/share/doc/o2scl/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(ADDINC) \
	$(O2SCL_OSX_MVAR) $(O2SCL_FFTW_MVAR) $(O2SCL_LINUX_MVAR) \
	$(O2SCL_MPFR_MVAR) \
	$(O2SCL_PYTHON_MVAR) $(O2SCL_CUBATURE_MVAR) \
	$(O2SCL_PUGIXML_MVAR) 
else
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"${datadir}/o2scl/\" \
	-DO2SCL_DOC_DIR=\"${docdir}/\" \
	$(O2SCL_OPENMP_MVAR) $(O2SCL_EIGEN_MVAR) $(O2SCL_ARMA_MVAR) \
	$(O2SCL_GSL2_MVAR) -DO2SCL_COND_FLAG $(ADDINC) \
	$(O2SCL_OSX_MVAR) $(O2SCL_FFTW_MVAR) $(O2SCL_LINUX_MVAR) \
	$(O2SCL_MPFR_MVAR) \
	$(O2SCL_PYTHON_MVAR) $(O2SCL_CUBATURE_MVAR) \
	$(O2SCL_PUGIXML_MVAR) 
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
	columnify_ts interp_krige_ts \
	string_conv_ts tensor_ts tensor_grid_ts vector_ts table3d_ts \
	format_float_ts table_units_ts exception_ts uniform_grid_ts \
	cursesw_ts auto_format_ts calc_utf8_ts

check_PROGRAMS = $(CPVAR)

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

interp_krige_ts_LDADD = $(ADDL_TEST_LIBS)
constants_ts_LDADD = $(ADDL_TEST_LIBS)
err_hnd_ts_LDADD = $(ADDL_TEST_LIBS)
convert_units_ts_LDADD = $(ADDL_TEST_LIBS)
funct_ts_LDADD = $(ADDL_TEST_LIBS)
format_float_ts_LDADD = $(ADDL_TEST_LIBS)
lib_settings_ts_LDADD = $(ADDL_TEST_LIBS)
misc_ts_LDADD = $(ADDL_TEST_LIBS)
auto_format_ts_LDADD = $(ADDL_TEST_LIBS)
calc_utf8_ts_LDADD = $(ADDL_TEST_LIBS)
mm_funct_ts_LDADD = $(ADDL_TEST_LIBS)
multi_funct_ts_LDADD = $(ADDL_TEST_LIBS)
search_vec_ts_LDADD = $(ADDL_TEST_LIBS)
table_ts_LDADD = $(ADDL_TEST_LIBS)
table3d_ts_LDADD = $(ADDL_TEST_LIBS)
vector_ts_LDADD = $(ADDL_TEST_LIBS)
interp_ts_LDADD = $(ADDL_TEST_LIBS)
columnify_ts_LDADD = $(ADDL_TEST_LIBS)
cursesw_ts_LDADD = $(ADDL_TEST_LIBS)
string_conv_ts_LDADD = $(ADDL_TEST_LIBS)
tensor_ts_LDADD = $(ADDL_TEST_LIBS)
tensor_grid_ts_LDADD = $(ADDL_TEST_LIBS)
table_units_ts_LDADD = $(ADDL_TEST_LIBS)
exception_ts_LDADD = $(ADDL_TEST_LIBS)
uniform_grid_ts_LDADD = $(ADDL_TEST_LIBS)

interp_krige_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
constants_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
err_hnd_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
convert_units_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
funct_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
format_float_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
lib_settings_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
misc_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
auto_format_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
calc_utf8_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
mm_funct_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
multi_funct_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
search_vec_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
table_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
table3d_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
vector_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
interp_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
columnify_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
cursesw_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
string_conv_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
tensor_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
tensor_grid_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
table_units_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
exception_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)
uniform_grid_ts_LDFLAGS = $(ADDL_TEST_LDFLGS)

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
auto_format.scr: auto_format_ts$(EXEEXT) 
	./auto_format_ts$(EXEEXT) > auto_format.scr
calc_utf8.scr: calc_utf8_ts$(EXEEXT) 
	./calc_utf8_ts$(EXEEXT) > calc_utf8.scr
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

interp_krige_ts_SOURCES = interp_krige_ts.cpp
constants_ts_SOURCES = constants_ts.cpp
err_hnd_ts_SOURCES = err_hnd_ts.cpp
convert_units_ts_SOURCES = convert_units_ts.cpp
funct_ts_SOURCES = funct_ts.cpp
format_float_ts_SOURCES = format_float_ts.cpp
lib_settings_ts_SOURCES = lib_settings_ts.cpp
misc_ts_SOURCES = misc_ts.cpp
auto_format_ts_SOURCES = auto_format_ts.cpp
calc_utf8_ts_SOURCES = calc_utf8_ts.cpp
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

# ------------------------------------------------------------
# Library o2scl_base
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_base.la

libo2scl_base_la_SOURCES = $(BASE_SRCS)

am__libo2scl_base_la_SOURCES_DIST = $(BASE_SRCS)

# ------------------------------------------------------------

interface: empty
	../../bin/yanic base.i base_python ~/wcs/o2sclpy/o2sclpy/base \
		~/wcs/o2sclpy/doc/base ../../bin/yanic_header.txt \
		| tee base_python.yanic.out

empty:

# ------------------------------------------------------------
# Misc
# ------------------------------------------------------------

EXTRA_DIST = *_ts.cpp output_path.py

emacs-clean: 
	-rm *~

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

