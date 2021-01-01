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
	python3 -msphinx.ext.intersphinx build/html/objects.inv | more

sphinx-help:
	sphinx-build --help

html: Makefile
	sphinx-build -v -b html "." build/html 2>&1 | tee sphinx.log 

html-nit: Makefile
	sphinx-build -v -n -b html "." build/html

full-html:
	sphinx-build -v -a -E -b html "." build/html

rst-files: empty
	grep "kind=\"function\"" ../xml/index.xml | \
		grep "refid=\"namespaceo2scl" | grep -v Eigen | \
		grep -v "refid=\"namespaceo2scl__const" | \
		grep -v arma > function_list
	grep "kind=\"class\"" ../xml/index.xml | grep -v "cli::parameter" \
		> class_list
	grep "kind=\"variable\"" ../xml/index.xml | \
		grep -v "<member refid=\"class" | \
		grep -v "<member refid=\"struct" > variable_list
	grep "kind=\"enum\"" ../xml/index.xml > enum_list
#	grep "kind=\"file\"" ../xml/index.xml > file_list
	cp ref_header ref.rst
	btmanip -parse ../../../data/o2scl/o2scl.bib -rst temp.rst
	cat temp.rst >> ref.rst
	rm -f temp.rst
	g++ -o ../make_rst ../make_rst.cpp
	rm -f class/*.rst
	rm -f function/*.rst
	../make_rst main | tee make_rst.out
	cd static; acol -run acol.scr > acol.out

o2scl-doc: rst-files html

empty:
