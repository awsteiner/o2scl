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

function_list: empty
	grep "kind=\"function\"" ../xml/index.xml | \
		grep "refid=\"namespaceo2scl" > function_list

clean:
	rm -rf make_index

full_html:
	$(SPHINXBUILD) -M -a -E html "$(SOURCEDIR)" \
		"$(BUILDDIR)" $(SPHINXOPTS) 

class_list: empty
	grep "kind=\"class\"" ../xml/index.xml > class_list

make_index: empty
	g++ -o make_index make_index.cpp

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

empty:
