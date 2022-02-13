help:
	@echo "This the O2scl makefile for the doc/o2scl subdirectory."
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
	@echo "clean-doc:            Delete the base generated doc files."
	@echo "o2scl-doc:            Make the base documentation"
	@echo "o2scl-doc-recursive:  Make all of the documentation"
	@echo "emacs-clean:          Remove tilde files."
	@echo 

# -----------------------------------------------------------------
# GNU style variables

docdir = $(prefix)/share/doc/o2scl

# These are the files installed by 'make install'
dist_doc_DATA = 

# These are the files included in the distribution
EXTRA_DIST = o2scl.tag make_rst.cpp 

SUBDIRS = sphinx xml

emacs-clean:
	rm -f *~
	rm -f */*~
	rm -f */*/*~
	rm -f */*/*/*~
	rm -f */*/*/*/*~

# Publish release documentation 
release-sync-doc:
	rsync -Cavzu /usr/local/share/doc/o2scl/* \
		$(STATIC_DOC_DIR)/o2scl

prerelease-sync-doc:
	rsync -Cavzu -r /usr/local/share/doc/o2scl/* \
		$(STATIC_DOC_DIR)/o2scl-dev

o2scl-doc:
	doxygen doxyfile > doxygen.scr
	cd sphinx && $(MAKE) o2scl-doc

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
	rm -rf html
	rm -f o2scl.tag

clean-doc-noform:
	cd html; rm -f *.html *.png doxygen.css \
		installdox tab*.gif tabs.css
	rm -f o2scl.tag

refs:
	cat refs_head.txt > refs.dox
	btmanip -p ../../data/o2scl/o2scl.bib -dox /tmp/btmanip
	cat /tmp/btmanip >> refs.dox
	cat refs_foot.txt >> refs.dox
