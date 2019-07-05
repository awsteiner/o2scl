PYTHON_SRCS = python.cpp

if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON
else
O2SCL_PYTHON_MVAR =
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
	-DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
	$(O2SCL_PYTHON_MVAR) $(PYTHON_INCLUDE)
else
AM_CPPFLAGS = $(O2SCL_PYTHON_MVAR) $(PYTHON_INCLUDE) \
	-I../../include
endif

VCHECK_LIBS = ../libo2scl.la

python.cpp: python.i
	swig -c++ -python -o python.cpp python.i

#ifeq ($(MACHINE),hedgehog)
#FLAGS_SWIG = --library-dirs=/Users/awsteiner/install/o2scl-0.923/lib \
#        --include-dirs=/Users/awsteiner/install/o2scl-0.923/include
#else
#FLAGS_SWIG =
#endif

swig_build:
	python3 setup.py build_ext --inplace

#$(FLAGS_SWIG)

# ------------------------------------------------------------
# Library o2scl_python
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_python.la

libo2scl_python_la_SOURCES = $(PYTHON_SRCS)

am__libo2scl_python_la_SOURCES_DIST = $(PYTHON_SRCS)
