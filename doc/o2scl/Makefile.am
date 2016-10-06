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

docdir = $(prefix)/doc/o2scl

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
	cp -r $(HOME)/install/o2scl-0.920/doc/o2scl \
		$(HOME)/wcs/int3/web/utk
	cp extras/index.html $(HOME)/wcs/int3/web/utk/o2scl
	cp extras/bootstrap.css $(HOME)/wcs/int3/web/utk/o2scl
	cp extras/bootstrap.js $(HOME)/wcs/int3/web/utk/o2scl
	cp extras/bootstrap-responsive.css $(HOME)/wcs/int3/web/utk/o2scl
	cp extras/jquery-2.0.3.min.js $(HOME)/wcs/int3/web/utk/o2scl

isospin-sync-doc:
	cp -r /usr/local/doc/o2scl $(HOME)/data
	cp extras/index.html $(HOME)/data/o2scl
	cp extras/bootstrap.css $(HOME)/data/o2scl
	cp extras/bootstrap.js $(HOME)/data/o2scl
	cp extras/bootstrap-responsive.css $(HOME)/data/o2scl
	cp extras/jquery-2.0.3.min.js $(HOME)/data/o2scl

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
