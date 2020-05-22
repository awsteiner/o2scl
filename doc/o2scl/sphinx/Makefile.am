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

# These are the files installed by 'make install'
dist_doc_DATA = *.rst

# These are the files included in the distribution
EXTRA_DIST = *.rst conf.py class_list function_list variable_list \
	enum_list 

#SUBDOC_DIRS = eos part 

SUBDIRS = class function

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

sphinx-help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

full-html:
	$(SPHINXBUILD) -M -a -E html "$(SOURCEDIR)" \
		"$(BUILDDIR)" $(SPHINXOPTS) 

function_list: empty
	grep "kind=\"function\"" ../xml/index.xml | \
		grep "refid=\"namespaceo2scl" | grep -v Eigen | \
		grep -v arma > function_list

class_list: empty
	grep "kind=\"class\"" ../xml/index.xml | grep -v "cli::parameter" \
		> class_list

variable_list: empty
	grep "kind=\"variable\"" ../xml/index.xml | \
		grep -v "<member refid=\"class" | \
		grep -v "<member refid=\"struct" > variable_list

enum_list: empty
	grep "kind=\"enum\"" ../xml/index.xml > enum_list

rst-files: empty
	cp ref_header ref.rst
	btmanip -parse ../../../data/o2scl/o2scl.bib -rst temp.rst
	cat temp.rst >> ref.rst
	rm -f temp.rst
	g++ -o ../make_rst ../make_rst.cpp
	../make_rst main > make_rst.out

empty:
