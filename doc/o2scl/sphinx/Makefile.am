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

docdir = $(prefix)/share/doc/o2scl/
docclassdir = $(prefix)/share/doc/o2scl/class
docfunctiondir = $(prefix)/share/doc/o2scl/function
dochtmldir = $(prefix)/share/doc/o2scl/html
dochtmlclassdir = $(prefix)/share/doc/o2scl/html/class
dochtmlfunctiondir = $(prefix)/share/doc/o2scl/html/function

dochtmlstaticdir = $(prefix)/share/doc/o2scl/html/_static
dochtmlsourcesdir = $(prefix)/share/doc/o2scl/html/_sources
dochtmlsourcesclassdir = $(prefix)/share/doc/o2scl/html/_sources/class
dochtmlsourcesfunctiondir = $(prefix)/share/doc/o2scl/html/_sources/function
dochtmlimagesdir = $(prefix)/share/doc/o2scl/html/_images

# These are the files installed by 'make install'
dist_doc_DATA = *.rst
dist_docclass_DATA = *.rst
dist_docfunction_DATA = *.rst
dist_dochtml_DATA = ../html/*.html ../html/*.js ../html/*.inv
dist_dochtmlclass_DATA = ../html/class/*.html
dist_dochtmlfunction_DATA = ../html/function/*.html

dist_dochtmlstatic_DATA = ../html/_static/*
dist_dochtmlimages_DATA = ../html/_images/*
dist_dochtmlsources_DATA = ../html/_sources/*.txt
dist_dochtmlsourcesclass_DATA = ../html/_sources/class/*.txt
dist_dochtmlsourcesfunction_DATA = ../html/_sources/function/*.txt

# These are the additional files included in the distribution
# 
# AWS 5/22/2020: I think all of the files listed for the make install
# target above are already included so all we need are the 'extra'
# files
EXTRA_DIST = conf.py class_list function_list variable_list \
	enum_list 

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
SPHINXPROJ    = o2scl
SOURCEDIR     = .
BUILDDIR      = build

oinv:
	python3 -msphinx.ext.intersphinx ../html/objects.inv | more

sphinx-help:
	sphinx-build --help

html: Makefile
	sphinx-build -v -b html "." ../html 2>&1 | tee sphinx.log 

html-nit: Makefile
	sphinx-build -v -n -b html "." ../html

full-html:
	sphinx-build -v -a -E -b html "." ../html

rst-files: empty
	cp ref_header ref.rst
	btmanip -parse ../../../data/o2scl/o2scl.bib -rst temp.rst
	cat temp.rst >> ref.rst
	rm -f temp.rst
	rm -f class/*.rst
	rm -f function/*.rst
	../../../bin/make_rst ../xml | tee make_rst.out
	cd static; acol -run acol.scr > acol.out

o2scl-doc: rst-files html

empty:
