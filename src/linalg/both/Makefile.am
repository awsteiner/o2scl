
if O2SCL_EIGEN
if O2SCL_ARMA
HEADER_VAR = qr_special.h linear_special.h cholesky_special.h \
	invert_special.h
else
HEADER_VAR = 
endif
else
HEADER_VAR = 
endif

# ------------------------------------------------------------
# Install headers
# ------------------------------------------------------------

am__pkginclude_HEADERS_DIST = qr_special.h linear_special.h \
	cholesky_special.h invert_special.h

pkginclude_HEADERS = $(HEADER_VAR)

