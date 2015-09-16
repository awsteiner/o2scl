docdir = $(prefix)/doc/o2scl/eos
dist_doc_DATA = doxyfile Makefile.am Makefile.in *.html

EXTRA_DIST = *.dox o2scl_eos.tag

SUBDIRS = html extras

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
	rm -f o2scl_eos.tag

o2scl-doc:
	doxygen doxyfile > doxygen.scr
	cat doxygen.log
