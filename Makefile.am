# Makefile.am for O2scl.

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
	@echo "o2scl-test:       Test installation ('make check' works also)"
	@echo "o2scl-examples:   Compile and run examples"
	@echo "test-clean:       Clean files created by tests"
	@echo 
	@echo "Developer targets"
	@echo "-------------------------------------------------------------"
	@echo "sinstall:         Install headers and libraries only"
	@echo "o2scl-doc:        Make HTML documentation"
	@echo "clean-doc:        Remove generated documentation files"
	@echo "blank-doc:        Create blank documentation (used in "
	@echo "open-doc:         View HTML documentation"
	@echo "  the homebrew recipe and the travis script)"
	@echo "release-sync-doc:     Copy docs from /usr/local on isospin"
	@echo "  to the web page repository for uploading"
	@echo "prerelease-sync-doc: Copy docs from /usr/local on isospin"
	@echo "  to the the internal web directory."
	@echo "emacs-clean:      Remove files ending with a tilde"
	@echo "permfix:          Change permissions to current user with sudo"
	@echo

# -----------------------------------------------------------------
# GNU style variables

ACLOCAL_AMFLAGS = -I m4

# Note that in SUBDIRS below, include must come before src (so that
# the necessary include files can be found) and bin must come after
# src (so that the libraries can be found for the executables)
SUBDIRS = include src bin data doc examples 

# Extra files to add to the distribution
EXTRA_DIST = autogen.sh TODOS README acinclude.m4 m4 README.md 

check_SCRIPTS = o2scl-test

# -----------------------------------------------------------------
# Additional targets

# Compile and run all examples
o2scl-examples:
	cd examples && $(MAKE) o2scl-examples

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
	cd bin && $(MAKE) install

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

pack-sync:
	git archive master | tar -x -C /home/awsteiner/pack/o2scl

# Publish release documentation 
release-sync-doc:
	cd doc/o2scl && $(MAKE) release-sync-doc

# Publish prerelease documentation 
prerelease-sync-doc:
	cd doc/o2scl && $(MAKE) prerelease-sync-doc

# Run doxygen to create documentation
o2scl-doc:
	cd doc/o2scl && $(MAKE) o2scl-doc-recursive

# Clean all generated documentation
clean-doc:
	cd doc/o2scl && $(MAKE) clean-doc-recursive

#Change permissions to current user with sudo
permfix:
	sudo chown -R `whoami`:`whoami` * .git

if O2SCL_OSX
    BROWSER = open
else
    BROWSER = xdg-open
endif

open-doc:
	$(BROWSER) doc/o2scl/html/index.html

interface: empty
	cd bin; $(MAKE) yanic
	cd src; $(MAKE) interface

empty:

# Create empty documentation
blank-doc:
	mkdir -p doc/o2scl/html
	mkdir -p doc/o2scl/eos/html
	mkdir -p doc/o2scl/part/html
	mkdir -p doc/o2scl/html/class
	mkdir -p doc/o2scl/eos/html/class
	mkdir -p doc/o2scl/part/html/class
	mkdir -p doc/o2scl/html/function
	mkdir -p doc/o2scl/eos/html/function
	mkdir -p doc/o2scl/part/html/function
	mkdir -p doc/o2scl/html/_sources
	mkdir -p doc/o2scl/eos/html/_sources
	mkdir -p doc/o2scl/part/html/_sources
	mkdir -p doc/o2scl/html/_sources/class
	mkdir -p doc/o2scl/eos/html/_sources/class
	mkdir -p doc/o2scl/part/html/_sources/class
	mkdir -p doc/o2scl/html/_sources/function
	mkdir -p doc/o2scl/eos/html/_sources/function
	mkdir -p doc/o2scl/part/html/_sources/function
	mkdir -p doc/o2scl/html/_static
	mkdir -p doc/o2scl/eos/html/_static
	mkdir -p doc/o2scl/part/html/_static
	mkdir -p doc/o2scl/html/_images
	mkdir -p doc/o2scl/eos/html/_images
	mkdir -p doc/o2scl/part/html/_images
	@echo ""
	touch doc/o2scl/html/index.html
	touch doc/o2scl/part/html/index.html
	touch doc/o2scl/eos/html/index.html
	touch doc/o2scl/html/blank.js
	touch doc/o2scl/part/html/blank.js
	touch doc/o2scl/eos/html/blank.js
	touch doc/o2scl/html/blank.inv
	touch doc/o2scl/part/html/blank.inv
	touch doc/o2scl/eos/html/blank.inv
	touch doc/o2scl/html/class/index.html
	touch doc/o2scl/part/html/class/index.html
	touch doc/o2scl/eos/html/class/index.html
	touch doc/o2scl/html/function/index.html
	touch doc/o2scl/part/html/function/index.html
	touch doc/o2scl/eos/html/function/index.html
	touch doc/o2scl/html/_sources/blank.txt
	touch doc/o2scl/eos/html/_sources/blank.txt
	touch doc/o2scl/part/html/_sources/blank.txt
	touch doc/o2scl/html/_sources/class/blank.txt
	touch doc/o2scl/eos/html/_sources/class/blank.txt
	touch doc/o2scl/part/html/_sources/class/blank.txt
	touch doc/o2scl/html/_sources/function/blank.txt
	touch doc/o2scl/eos/html/_sources/function/blank.txt
	touch doc/o2scl/part/html/_sources/function/blank.txt
	touch doc/o2scl/html/_static/blank.png
	touch doc/o2scl/eos/html/_static/blank.png
	touch doc/o2scl/part/html/_static/blank.png
	touch doc/o2scl/html/_images/blank.png
	touch doc/o2scl/eos/html/_images/blank.png
	touch doc/o2scl/part/html/_images/blank.png

snapb:
	snapcraft build 2>&1 > ~/o2scl_snap.out &
