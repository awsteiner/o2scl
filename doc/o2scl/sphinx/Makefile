
# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line. 
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = o2scl
SOURCEDIR     = .
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

full_html:
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

make_rst: empty
	g++ -o ../make_rst ../make_rst.cpp

rst_files: empty
	cp ref_header ref.rst
	btmanip -parse ../../../data/o2scl/o2scl.bib -rst temp.rst
	cat temp.rst >> ref.rst
	rm -f temp.rst
	g++ -o ../make_rst ../make_rst.cpp
	../make_rst main > make_rst.out

clean: empty
	rm -rf rst

empty:
