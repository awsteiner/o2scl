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

remove-links: 
	rm -f *.h

