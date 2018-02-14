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

docdir = $(prefix)/share/doc/o2scl

# These are the files installed by 'make install'
dist_doc_DATA = doxyfile Makefile.am Makefile.in *.html

# These are the files included in the distribution
EXTRA_DIST = *.dox o2scl.tag

SUBDOC_DIRS = eos part 

SUBDIRS = extras eos part html 

emacs-clean:
	rm -f *~
	rm -f */*~
	rm -f */*/*~
	rm -f */*/*/*~
	rm -f */*/*/*/*~

utk-sync-doc:
	sudo cp -r /usr/local/share/doc/o2scl/* $(STATIC_DOC_DIR)/o2scl
	sudo cp extras/index.html $(STATIC_DOC_DIR)/o2scl
	sudo cp extras/bootstrap.css $(STATIC_DOC_DIR)/o2scl
	sudo cp extras/bootstrap.js $(STATIC_DOC_DIR)/o2scl
	sudo cp extras/bootstrap-responsive.css $(STATIC_DOC_DIR)/o2scl
	sudo cp extras/jquery-2.0.3.min.js $(STATIC_DOC_DIR)/o2scl

prerelease-sync-doc:
	sudo cp -r /usr/local/share/doc/o2scl/* $(STATIC_DOC_DIR)/o2scl-dev
	sudo cp extras/index.html $(STATIC_DOC_DIR)/o2scl-dev
	sudo cp extras/bootstrap.css $(STATIC_DOC_DIR)/o2scl-dev
	sudo cp extras/bootstrap.js $(STATIC_DOC_DIR)/o2scl-dev
	sudo cp extras/bootstrap-responsive.css $(STATIC_DOC_DIR)/o2scl-dev
	sudo cp extras/jquery-2.0.3.min.js $(STATIC_DOC_DIR)/o2scl-dev

o2scl-doc-recursive: refs o2scl-doc
	cd eos && $(MAKE) refs o2scl-doc
	cd part && $(MAKE) refs o2scl-doc

o2scl-doc:
	cd extras; acol -run acol.scr > acol.out
	doxygen doxyfile > doxygen.scr
	cat doxygen.log

clean-doc-recursive: clean-doc
	cd eos && $(MAKE) clean-doc
	cd part && $(MAKE) clean-doc

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
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
