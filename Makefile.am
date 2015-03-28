#
# Makefile.am for O2scl.
#

ACLOCAL_AMFLAGS = -I m4

# Note that in SUBDIRS below, include must come before src
SUBDIRS = include src data doc examples

# Extra files to add to the distribution
EXTRA_DIST = autogen.sh TODOS README acinclude.m4 m4

# Output some assistance describing the usual make targets
help:
	@echo "This the O2scl root directory makefile."
	@echo
	@echo "The traditional GNU-style targets apply. In addition there"
	@echo "are a couple extra targets which may be useful listed below."
	@echo "The documentation is automatically installed with"
	@echo "'make install'."
	@echo 
	@echo "Additional package-specific targets"
	@echo "-------------------------------------------------------------"
	@echo "help:             Show this help file"
	@echo "install:          Install"
	@echo "sinstall:         Install headers and libraries only"
	@echo "o2scl-test:       Test installation"
	@echo "o2scl-examples:   Compile and run examples"
	@echo "o2scl-benchmarks: Compile and run benchmarks"
	@echo "test-clean:       Clean files created by tests"
	@echo 
	@echo "Developer targets"
	@echo "-------------------------------------------------------------"
	@echo "o2scl-doc:        Make HTML documentation"
#	@echo "o2scl-docp:       Make LaTeX documentation"
	@echo "clean-doc:        Remove generated doc files"
	@echo "dist-nodoc:       Make a distribution w/o the documentation"
	@echo "sf-sync-doc:      Upload current documentation to SourceForge"
	@echo "emacs-clean:      Remove files ending with a tilde"
	@echo 

#	@echo "pngfix:           Remove transparency in doc PNG images"
#	@echo "tharkad-web:      Copy current docs to tharkad"

#

# @echo "Please use 'make o2scl-test' to test the installation."
#	ls 

# Create a distribution without documentation (not recently tested)
dist-nodoc:
	cd doc/o2scl/html; mkdir ../htdir; mv -i *.html *.png ../htdir; 
	cd doc/o2scl/html; touch index.html
	cd doc/o2scl/latex; mkdir tmpdir; mv -i *.tex *.eps *.pdf tmpdir
	cd doc/o2scl/latex; touch refman.pdf
	cd doc/o2scl/eos/html; mkdir ../htdir; mv -i *.html *.png ../htdir; 
	cd doc/o2scl/eos/html; touch index.html
	cd doc/o2scl/eos/latex; mkdir tmpdir; mv -i *.tex *.eps *.pdf tmpdir
	cd doc/o2scl/eos/latex; touch refman.pdf
	cd doc/o2scl/part/html; mkdir ../htdir; mv -i *.html *.png ../htdir; 
	cd doc/o2scl/part/html; touch index.html
	cd doc/o2scl/part/latex; mkdir tmpdir; mv -i *.tex *.eps *.pdf tmpdir
	cd doc/o2scl/part/latex; touch refman.pdf
	make dist
	cd doc/o2scl/html; mv ../htdir/* .; rmdir ../htdir
	cd doc/o2scl/latex; mv tmpdir/* .; rmdir tmpdir
	cd doc/o2scl/eos/html; mv ../htdir/* .; rmdir ../htdir
	cd doc/o2scl/eos/latex; mv tmpdir/* .; rmdir tmpdir
	cd doc/o2scl/part/html; mv ../htdir/* .; rmdir ../htdir
	cd doc/o2scl/part/latex; mv tmpdir/* .; rmdir tmpdir

# Compile and run all examples
o2scl-examples:
	cd examples && $(MAKE) o2scl-examples

# Compile and run all benchmarks
o2scl-benchmarks:
	cd examples && $(MAKE) o2scl-benchmarks

# Run the tests and summarize using the tsumm program
#
# The command echo " " > ./testlist creates a blank testlist file
# which the various subdirectories can add to.
o2scl-test: 
	echo " " > ./testlist
	cd src && $(MAKE) o2scl-test-rec
	./src/internal/tsumm$(EXEEXT) summary

# Install source only
sinstall:
	cd include && $(MAKE) install
	cd src && $(MAKE) install

# Clean o2scl-test targets
test-clean:
	rm -f testlist
	rm -f src/*/*_ts
	rm -f src/*/*_ts.o
	rm -f src/*/*.scr

# Clean emacs checkpoint files
emacs-clean:
	rm -v -f *~
	rm -v -f */*~
	rm -v -f */*/*~
	rm -v -f */*/*/*~
	rm -v -f */*/*/*/*~

# Upload docs to sourceforge
sf-sync-doc:
	cd doc/o2scl && $(MAKE) sf-sync-doc

# Run doxygen to create documentation
o2scl-doc:
	cd doc/o2scl && $(MAKE) o2scl-doc-recursive

# Create .pdf file from doxygen latex files
o2scl-docp:
	echo "Disabled."
#	cd doc/o2scl && $(MAKE) o2scl-docp-recursive

# Create .pdf file from doxygen latex files verbosely
#o2scl-docp-verbose:
	echo "Disabled."
	cd doc/o2scl && $(MAKE) o2scl-docp-verbose-recursive

check_SCRIPTS = o2scl-test

# Clean all generated documentation
clean-doc:
	cd doc/o2scl && $(MAKE) clean-doc-recursive

pngfix:
	cd doc/o2scl; $(MAKE) pngfix-recursive
