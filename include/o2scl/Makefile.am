help:
	@echo "This the O2scl makefile for the include/o2scl subdirectory."
	@echo
	@echo "The traditional GNU-style targets apply. Note that"
	@echo "this directory only contains links to the header"
	@echo "files which are found in the src subdirectories."
	@echo 

all: all-am header-links

clean: clean-am remove-links

distclean: distclean-am remove-links

install: install-am header-links

# Homebrew seems to delete the header links before running
# 'make check', so we add them back here
check_SCRIPTS = header-links

header-links: remove-links
# Links to normal header files
	HEADERLIST="$(top_srcdir)/src/*/*.h"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  test -r $$BASENAME || ln -s $$h $$BASENAME; \
	done
if O2SCL_OPENMP
	test -r set_openmp.h || ln -s $(top_srcdir)/src/base/openmp/set_openmp.h set_openmp.h
else
	test -r set_openmp.h || ln -s $(top_srcdir)/src/base/neither/set_openmp.h set_openmp.h
endif
if O2SCL_PYTHON
	test -r set_python.h || ln -s $(top_srcdir)/src/base/python/set_python.h set_python.h
else
	test -r set_python.h || ln -s $(top_srcdir)/src/base/neither/set_python.h set_python.h
endif
if O2SCL_PUGIXML
	test -r set_pugixml.h || ln -s $(top_srcdir)/src/base/pugixml/set_pugixml.h set_pugixml.h
else
	test -r set_pugixml.h || ln -s $(top_srcdir)/src/base/neither/set_pugixml.h set_pugixml.h
endif
if O2SCL_MPFR
	test -r set_mpfr.h || ln -s $(top_srcdir)/src/base/mpfr/set_mpfr.h set_mpfr.h
else
	test -r set_mpfr.h || ln -s $(top_srcdir)/src/base/neither/set_mpfr.h set_mpfr.h
endif
if O2SCL_FFTW
	test -r set_fftw.h || ln -s $(top_srcdir)/src/base/fftw/set_fftw.h set_fftw.h
else
	test -r set_fftw.h || ln -s $(top_srcdir)/src/base/neither/set_fftw.h set_fftw.h
endif
if O2SCL_READLINE
	test -r set_readline.h || ln -s $(top_srcdir)/src/base/readline/set_readline.h set_readline.h
else
	test -r set_readline.h || ln -s $(top_srcdir)/src/base/neither/set_readline.h set_readline.h
endif
if O2SCL_ARMA
if O2SCL_EIGEN
	test -r vector_special.h || ln -s $(top_srcdir)/src/base/both/vector_special.h vector_special.h
else
	test -r vector_special.h || ln -s $(top_srcdir)/src/base/arma/vector_special.h vector_special.h
endif
else
if O2SCL_EIGEN
	test -r vector_special.h || ln -s $(top_srcdir)/src/base/eigen/vector_special.h vector_special.h
else
	test -r vector_special.h || ln -s $(top_srcdir)/src/base/neither/vector_special.h vector_special.h
endif
endif

remove-links: 
	rm -f *.h

