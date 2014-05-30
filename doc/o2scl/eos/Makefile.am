docdir = $(prefix)/doc/o2scl/eos
dist_doc_DATA = doxyfile Makefile.am Makefile.in *.html

EXTRA_DIST = *.dox o2scl_eos.tag

SUBDIRS = html latex extras

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
#	cd latex; rm -f *.tex *.eps *.pdf doxygen.sty *.ttf makepdf.* \
#		refman.*
	rm -f o2scl_eos.tag

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

#o2scl-docp:
#	cd latex; make -f Makefile.dox refman.pdf > makepdf.scr 2> makepdf.err

#o2scl-docp-verbose:
#	cd latex; make -f Makefile.dox refman.pdf

pngfix:
	cd html; $(MAKE) pngfix
