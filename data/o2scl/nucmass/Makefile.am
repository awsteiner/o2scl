ndatadir=$(datadir)/o2scl/nucmass

dist_ndata_DATA = ame03.o2 ame03round.o2 ame95exp.o2 hfb*.o2 \
			ame95rmd.o2 Makefile.am mnmsk.o2 \
			du_zu_95.o2 du_zu_96.o2 ktuy04.o2 ktuy05.o2 \
			ame12.o2 wl11.o2 wllw10.o2 wlw10.o2 \
			wlwm14.o2 lwdw11.o2 sd*.o2

WSall: WS3.2.o2 WS3.3.o2 WS3.6.o2 WS3_RBF.o2 WS4_RBF.o2

WS3.2.o2:
	echo "A Z Beta2 Beta4 Esh Eexp Eth Mexp Mth" > temp32
	tail -n +25 ~/svn/osf/data/nucmass/WS3.2.txt >> temp32
	acol -generic temp32 -internal WS3.2.o2

WS3.3.o2:
	echo "A Z Beta2 Beta4 Beta6 Esh Eexp Eth Mexp Mth" > temp33
	tail -n +27 ~/svn/osf/data/nucmass/WS3.3.txt >> temp33
	acol -generic temp33 -internal WS3.3.o2

WS3.6.o2:
	echo "A Z Beta2 Beta4 Beta6 Esh Dres Eexp Eth Mexp Mth" > temp36
	tail -n +33 ~/svn/osf/data/nucmass/WS3.6.txt >> temp36
	acol -generic temp36 -internal WS3.6.o2

WS3_RBF.o2:
	echo "A Z WS3 WS3_RBF" > temp3r
	tail -n +14 ~/svn/osf/data/nucmass/WS3_RBF.txt >> temp3r
	acol -generic temp3r -internal WS3_RBF.o2

WS4_RBF.o2:
	echo "A Z WS4 WS4_RBF" > temp4r
	tail -n +16 ~/svn/osf/data/nucmass/WS4_RBF.txt >> temp4r
	acol -generic temp4r -internal WS4_RBF.o2