docextrasdir = $(prefix)/share/doc/o2scl/extras
dist_docextras_DATA = *.png *.cpp *.html *.txt *.odg

# *.eps 

help:
	@echo
	@echo "File descriptions:"
	@echo "--------------------------------------------------------"
	@echo "contourg2.eps: Picture for contour documentation "
	@echo "contourg2.png: Picture for contour documentation "
	@echo "contourg3.eps: Picture for contour documentation "
	@echo "contourg3.png: Picture for contour documentation "
	@echo "contourg.eps: Picture for contour documentation "
	@echo "contourg.png: Picture for contour documentation "
	@echo "desc.png: Description for HTML front page"
	@echo "dl.png: Download button for HTML front page"
	@echo "doc.png: Documentation button for HTML front page"
	@echo "fdl_license.txt: GNU Doc license"
	@echo "fptr_doc.cpp: Docs for function pointers used in main.h"
	@echo "gpl_license.txt: GNU GPL"
	@echo "gsl_vector_doc.cpp: Docs for gsl_vector used in main.h"
	@echo "index.html: HTML front page"
	@echo "io_doc1.cpp: Docs for I/O used in main.h"
	@echo "io_doc2.cpp: Docs for I/O used in main.h"
	@echo "io_doc4.cpp: Docs for I/O used in main.h"
	@echo "io_doc5.cpp: Docs for I/O used in main.h"
	@echo "io_doc6.cpp: Docs for I/O used in main.h"
	@echo "otitle.odg: OpenOffice base file for HTML images"
	@echo "otitle.png: Title for HTML front page"
	@echo "ovector_doc1.cpp: Docs for ovector used in main.h"
	@echo "ovector_doc2.cpp: Docs for ovector used in main.h"
	@echo "ovector_doc3.cpp: Docs for ovector used in main.h"
	@echo "ovector_doc4.cpp: Docs for ovector used in main.h"
	@echo

#ex_fptr.png: 
#	cd ../../../examples; make o2scl-examples
#	cat ../../../examples/ex_fptr.scr | grep -v brent | grep -v value \
#		> temp1

O2SCL_LIB_DIR = /home/asteiner/install/lib
O2SCL_INC_DIR = /home/asteiner/install/include

READLINE_LIBS = -lreadline -lncurses

FLAGS = -ggdb -Wreturn-type -Wparentheses -Wconversion -Wall \
	-Wno-unused -Wshadow -Wimplicit -I$(O2SCL_INC_DIR) \
	-DO2SCL_ARRAY_ABORT

# This ordering must be maintained for the static compilation to work
LIB =  -lm -L$(O2SCL_LIB_DIR) -lo2int_base -lo2scl_ext -lo2scl_eos \
	-lo2scl_part -lo2scl_root -lo2scl_other -lo2scl_nuclei \
	-lo2scl_minimize -lo2scl_inte -lo2scl_base -lgsl -lgslcblas \
	$(READLINE_LIBS)

ROOT_INC = -I/home/asteiner/root/include
ROOT_LIB = /home/asteiner/svn/o2int/trunk/src/other/graph.o \
	-L/home/asteiner/root/lib \
	-lgsl -lgslcblas -lm -lCore -lCint -lHist -lGraf \
	-lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lGeom \
	-lPhysics -lGui -lm -ldl -lpthread -rdynamic

