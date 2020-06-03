help:
	@echo "This the O2scl makefile for the doc/o2scl/sphinx subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply. In addition there"
	@echo "are a couple extra targets which may be useful listed below."
	@echo 
	@echo "Additional package-specific targets"
	@echo "-------------------------------------------------------------"
	@echo "help:                 Show this help file"
	@echo 
	@echo "Developer targets"
	@echo "-------------------------------------------------------------"
	@echo 

docdir = $(prefix)/share/doc/o2scl/part
docclassdir = $(prefix)/share/doc/o2scl/part/class
docfunctiondir = $(prefix)/share/doc/o2scl/part/function
dochtmldir = $(prefix)/share/doc/o2scl/part/html
dochtmlclassdir = $(prefix)/share/doc/o2scl/part/html/class
dochtmlfunctiondir = $(prefix)/share/doc/o2scl/part/html/function

dochtmlstaticdir = $(prefix)/share/doc/o2scl/part/html/_static
dochtmlsourcesdir = $(prefix)/share/doc/o2scl/part/html/_sources
dochtmlsourcesclassdir = $(prefix)/share/doc/o2scl/part/html/_sources/class
dochtmlsourcesfunctiondir = $(prefix)/share/doc/o2scl/part/html/_sources/function
dochtmlimagesdir = $(prefix)/share/doc/o2scl/part/html/_images

# These are the files installed by 'make install'
dist_doc_DATA = *.rst
dist_docclass_DATA = *.rst
dist_docfunction_DATA = *.rst
dist_dochtml_DATA = build/html/*.html build/html/*.js build/html/*.inv
dist_dochtmlclass_DATA = build/html/class/*.html
dist_dochtmlfunction_DATA = build/html/function/*.html

dist_dochtmlstatic_DATA = build/html/_static/*
dist_dochtmlimages_DATA = build/html/_images/*
dist_dochtmlsources_DATA = build/html/_sources/*.txt
dist_dochtmlsourcesclass_DATA = build/html/_sources/class/*.txt
dist_dochtmlsourcesfunction_DATA = build/html/_sources/function/*.txt

# These are the additional files included in the distribution
# 
# AWS 5/22/2020: I think all of the files listed for the make install
# target above are already included so all we need are the 'extra'
# files
EXTRA_DIST = conf.py class_list function_list

#SUBDOC_DIRS = eos part 

SUBDIRS = 

emacs-clean:
	rm -f *~
	rm -f */*~
	rm -f */*/*~
	rm -f */*/*/*~
	rm -f */*/*/*/*~

SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = o2scl_part
SOURCEDIR     = .
BUILDDIR      = build

sphinx-help:
	sphinx-build --help

html: Makefile
	sphinx-build -v -b html "." build/html

html-nit: Makefile
	sphinx-build -v -n -b html "." build/html

html-full:
	sphinx-build -v -a -E -b html "." build/html

html-full-nit:
	sphinx-build -v -a -E -n -b html "." build/html

rst-files: empty
	grep "kind=\"function\"" ../xml/index.xml | \
		grep "refid=\"namespaceo2scl__part" > function_list
	grep "kind=\"class\"" ../xml/index.xml | grep -v "::entry" \
		> class_list
	cp ref_header ref.rst
	btmanip -parse ../../../../data/o2scl/o2scl_part.bib -rst temp.rst
	cat temp.rst >> ref.rst
	rm -f temp.rst
	g++ -o ../../make_rst ../../make_rst.cpp
	rm -f class/*.rst
	rm -f function/*.rst
	../../make_rst part > make_rst.out

o2scl-doc: rst-files html

empty:
