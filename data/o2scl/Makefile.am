if O2SCL_PARTLIB
if O2SCL_EOSLIB

SUBDIRVAR = nucmass rmfdata skdata

DATAVAR = bbp71.o2 bps71.o2 fmt49.o2 hzd89.o2 dh01.o2 \
	litable nvbps.o2 apr98.o2 fp81.o2 hz90.o2 fermion_cal.o2 \
	odcgcp88.o2 newton_PNM.o2 newton_J35.o2 APR_cs01_feq.txt \
	SLy4_cs01_feq.txt Rs_cs01_feq.txt gogny.o2 newton_SL.o2 \
	eos19.o2 eos20.o2 eos21.o2 nv83.o2 o2scl.bib o2scl_eos.bib \
	o2scl_part.bib fermion_cal2.o2

else

SUBDIRVAR = nucmass

DATAVAR = litable

endif
else

SUBDIRVAR = nucmass

DATAVAR = litable

endif

SUBDIRS = $(SUBDIRVAR)

dist_pkgdata_DATA = $(DATAVAR)

