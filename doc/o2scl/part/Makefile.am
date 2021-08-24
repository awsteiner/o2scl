docdir = $(prefix)/share/doc/o2scl/part
dist_doc_DATA = doxyfile Makefile.am Makefile.in

EXTRA_DIST = o2scl_part.tag

SUBDIRS = sphinx

clean-doc:
	cd html; rm -f *.html *.png formula.repository doxygen.css \
		_formulas* installdox tab*.gif tabs.css
	rm -rf html
	rm -f o2scl_part.tag

o2scl-doc:
	doxygen doxyfile > doxygen.scr
	cd sphinx && $(MAKE) o2scl-doc

refs:
	cat refs_head.txt > refs.dox
	btmanip -p ../../../data/o2scl/o2scl_part.bib -dox /tmp/btmanip
	cat /tmp/btmanip >> refs.dox
	cat refs_foot.txt >> refs.dox
