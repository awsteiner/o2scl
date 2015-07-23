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
#	@echo "o2scl-docp:           Make the base PDF documentation"
	@echo "o2scl-doc-recursive:  Make all of the documentation"
#	@echo "o2scl-docp-recursive: Make all of the PDF documentation"
	@echo "pngfix:               Fix the base PNG images"
	@echo "pngfix-recursive:     Fix all of the PNG images"
	@echo "sf-sync-doc:          Upload documentation to SF."
	@echo "sf-sync-index:        Upload index HTML files to SF."
	@echo "emacs-clean:          Remove tilde files."
	@echo 

docdir = $(prefix)/doc/o2scl

# These are the files installed by 'make install'
dist_doc_DATA = doxyfile Makefile.am Makefile.in *.html

# These are the files included in the distribution
EXTRA_DIST = *.dox o2scl.tag

SUBDOC_DIRS = eos part 

SUBDIRS = extras eos part html 
#latex

emacs-clean:
	rm -f *~
	rm -f */*~
	rm -f */*/*~
	rm -f */*/*/*~
	rm -f */*/*/*/*~

sf-sync-index:
	rsync -Cavzu extras/index.html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/index.html
	rsync -Cavzu extras/bootstrap-responsive.css \
		awsteiner,o2scl@web.sourceforge.net:htdocs/boostrap-responsive.css
	rsync -Cavzu extras/bootstrap.css \
		awsteiner,o2scl@web.sourceforge.net:htdocs/bootstrap.css
	rsync -Cavzu extras/bootstrap.js \
		awsteiner,o2scl@web.sourceforge.net:htdocs/bootstrap.js
	rsync -Cavzu extras/jquery-2.0.3.min.js \
		awsteiner,o2scl@web.sourceforge.net:htdocs/jquery-2.0.3.min.js

sf-sync-doc:
	rsync -Cavzu --exclude=.svn html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl
#	rsync -Cavzu --exclude=.svn latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/latex/refman.pdf
	rsync -Cavzu --exclude=.svn eos/html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/eos
#	rsync -Cavzu --exclude=.svn eos/latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/eos/latex/refman.pdf
	rsync -Cavzu --exclude=.svn part/html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/part
#	rsync -Cavzu --exclude=.svn part/latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/part/latex/refman.pdf

utk-sync-doc:
	cp -r $(HOME)/install/o2scl-0.918/doc/o2scl \
		$(HOME)/svn/int3/web/utk
	cp extras/index.html $(HOME)/svn/int3/web/utk/o2scl
	cp extras/bootstrap.css $(HOME)/svn/int3/web/utk/o2scl
	cp extras/bootstrap.js $(HOME)/svn/int3/web/utk/o2scl
	cp extras/bootstrap-responsive.css $(HOME)/svn/int3/web/utk/o2scl
	cp extras/jquery-2.0.3.min.js $(HOME)/svn/int3/web/utk/o2scl

sf-test-doc:
	rsync -Cavzun --exclude=.svn html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl
#	rsync -Cavzun --exclude=.svn latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/latex/refman.pdf
	rsync -Cavzun --exclude=.svn eos/html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/eos
#	rsync -Cavzun --exclude=.svn eos/latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/eos/latex/refman.pdf
	rsync -Cavzun --exclude=.svn part/html \
		awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/part
#	rsync -Cavzun --exclude=.svn part/latex/refman.pdf \
#	awsteiner,o2scl@web.sourceforge.net:htdocs/o2scl/part/latex/refman.pdf

o2scl-doc-recursive: o2scl-doc
	cd eos && $(MAKE) o2scl-doc
	cd part && $(MAKE) o2scl-doc

# We have to move the Makefile out of the way because
# doxygen will overwrite it
o2scl-doc:
#	mv latex/Makefile latex/Makefile.tmp
	doxygen doxyfile > doxygen.scr
#	mv latex/Makefile latex/Makefile.dox
#	mv latex/Makefile.tmp latex/Makefile
#	perl rep.perl < html/search/search.css > temp.css
#	mv temp.css html/search/search.css
	cat doxygen.log

# o2scl-docp-recursive: o2scl-docp
# 	cd eos && $(MAKE) o2scl-docp
# 	cd part && $(MAKE) o2scl-docp

# o2scl-docp-verbose-recursive: o2scl-docp-verbose
# 	cd eos && $(MAKE) o2scl-docp-verbose
# 	cd part && $(MAKE) o2scl-docp-verbose

# o2scl-docp:
# 	cd latex && $(MAKE) -f Makefile.dox refman.pdf > makepdf.scr \
# 		2> makepdf.err

# o2scl-docp-verbose:
# 	cd latex && $(MAKE) -f Makefile.dox refman.pdf

clean-doc-recursive: clean-doc
	cd eos && $(MAKE) clean-doc
	cd part && $(MAKE) clean-doc

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
#	cd latex; rm -f *.tex *.eps *.pdf doxygen.sty *.ttf makepdf.* \
#		refman.*
	rm -f o2scl.tag

clean-doc-noform:
	cd html; rm -f *.html *.png doxygen.css \
		installdox tab*.gif tabs.css
#	cd latex; rm -f *.eps *.pdf doxygen.sty *.ttf makepdf.* \
#		refman.*
	rm -f o2scl.tag

pngfix:
	cd html; $(MAKE) pngfix

pngfix-recursive:
	cd html; $(MAKE) pngfix
	cd eos/html; $(MAKE) pngfix
	cd part/html; $(MAKE) pngfix


