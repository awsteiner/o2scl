if O2SCL_PYTHON
O2SCL_PYTHON_MVAR = -DO2SCL_PYTHON $(PYTHON_INCLUDE)
else
O2SCL_PYTHON_MVAR =
endif

if O2SCL_SNAP
AM_CPPFLAGS = -I@top_srcdir@/include/ \
        -DO2SCL_DATA_DIR=\"/snap/o2scl/current/share/o2scl/\" \
        $(O2SCL_PYTHON_MVAR)
else
AM_CPPFLAGS = 
        $(O2SCL_PYTHON_MVAR)
endif

PYTHON_SRCS = python.cpp

python.cpp: python.i
	swig -c++ -python -o python.cpp python.i


# ------------------------------------------------------------
# Library o2scl_base
# ------------------------------------------------------------

noinst_LTLIBRARIES = libo2scl_python.la

libo2scl_python_la_SOURCES = $(PYTHON_SRCS)

am__libo2scl_python_la_SOURCES_DIST = $(PYTHON_SRCS)
