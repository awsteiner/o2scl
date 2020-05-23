docdir = $(prefix)/share/doc/o2scl/part
dist_doc_DATA = doxyfile Makefile.am Makefile.in *.html

EXTRA_DIST = *.dox o2scl_part.tag refs_head.txt refs_foot.txt dlayout.xml

SUBDIRS = extras sphinx

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
	rm -f o2scl_part.tag

o2scl-doc:
	doxygen doxyfile > doxygen.scr
	cd sphinx && $(MAKE) o2scl-doc
#	doxygen doxyfile > doxygen.scr
#	cat doxygen.log

refs:
	cat refs_head.txt > refs.dox
	btmanip -p ../../../data/o2scl/o2scl_part.bib -dox /tmp/btmanip
	cat /tmp/btmanip >> refs.dox
	cat refs_foot.txt >> refs.dox
