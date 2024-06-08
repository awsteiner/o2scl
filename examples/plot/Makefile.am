# --------------------------------------------------------------------
# All figures should be generated with python or o2graph
#  
# Reworking June 2019
# 
# -------------------------------------------------------------------- 

help:
	@echo 
	@echo "This directory contains code which assists the developer"
	@echo "in generating the documentation plots and is not"
	@echo "for use by the casual end-user."
	@echo

# -------------------------------------------------------------------- 
# This directory generates the documentation images after all of
# the examples in the parent directory have been compiled and 
# executed. The main target, doc-extras, should create all of the files
# necessary to generate the documentation.
# -------------------------------------------------------------------- 

EXTRA_DIST = *.png

plot-all: ex_anneal_plot ex_chebapp_plot ex_contour_plots \
	ex_fptr_plot ex_mmin_plot ex_mmin_plot2 \
	ex_nuc_prof ex_nucmass_plots \
	ex_ode_airy_plot ex_ode_bessel_plot ex_ode_bessel2_plot \
	ex_ode_bessel3_plot ex_stiff_plot

# Use 'open' for OSX and 'eyes of gnome' for linux

if O2SCL_LINUX
IMAGE_VIEW = eog
IMAGE_VIEW2 = &
else
if O2SCL_OSX
IMAGE_VIEW = open
IMAGE_VIEW2 = 
else 
IMAGE_VIEW =
IMAGE_VIEW2 = 
endif
endif

# -------------------------------------------------------------------- 

ex_interp_plot:
	o2graph -read ../ex_interp.o2 tdata \
		-set ylo "(-2.0)" -set yhi 1.5 \
		-scatter x y \
		-read ../ex_interp.o2 tresult \
		-plot x ylin color=k -plot x ycsp color=r \
		-plot x ystef color=b \
		-plot x yiko color=r,ls=: \
		-plot x yiko_lml color=b,ls=: \
		-xtitle x -ytitle y \
		-line 1.0 0.6 1.15 0.6 color=k \
		-line 1.0 0.8 1.15 0.8 color=r \
		-line 1.0 1.0 1.15 1.0 color=b \
		-line 0.9 1.2 1.05 1.2 color=r,ls=: \
		-line 0.9 1.4 1.05 1.4 color=b,ls=: \
		-text 1.2 0.6 "Linear" ha=left \
		-text 1.2 0.8 "Cubic spline" ha=left \
		-text 1.2 1.0 "Steffen's" ha=left \
		-text 1.1 1.2 "Kriging (CV)" ha=left \
		-text 1.1 1.4 "Kriging (ML)" ha=left \
		-inset 0.35 0.20 0.5 0.4 \
		-set logy 1 -set ylo 1.0e-6 -set yhi 0.1 \
		-function "abs(ylin-y)" lin2 -plot x lin2 color=k \
		-function "abs(ycsp-y)" csp2 -plot x csp2 color=r \
		-function "abs(ystef-y)" stef2 -plot x stef2 color=b \
		-function "abs(yiko-y)" iko2 -plot x iko2 color=r,ls=: \
		-function "abs(yiko_lml-y)" iko_lml2 \
		-plot x iko_lml2 color=b,ls=: \
		-ytitle "$$ |\Delta y| $$" \
		-ttext 0.5 "(-0.15)" x \
		-save ex_interp.pdf -show

ex_anneal_plot:
	o2graph -backend Agg -set colbar 1 \
		-read ../ex_anneal_plot.o2 -set xtitle x -set ytitle y \
		-den-plot f \
		-contours "(-0.9)" f ex_anneal_conts.o2 l1 \
		-contours "(-0.7)" f ex_anneal_conts.o2 l2 \
		-contours "(-0.5)" f ex_anneal_conts.o2 l3 \
		-contours "(-0.3)" f ex_anneal_conts.o2 l4 \
		-contours "(-0.2)" f ex_anneal_conts.o2 l5 \
		-contours "(-0.1)" f ex_anneal_conts.o2 l6 \
		-contours "(-0.05)" f ex_anneal_conts.o2 l7 \
		-contours "0.05" f ex_anneal_conts.o2 l8 \
		-contours "0.1" f ex_anneal_conts.o2 l9 \
		-read ex_anneal_conts.o2 l1 -plot color=black \
		-read ex_anneal_conts.o2 l2 -plot color=black \
		-read ex_anneal_conts.o2 l3 -plot color=black \
		-read ex_anneal_conts.o2 l4 -plot color=black \
		-read ex_anneal_conts.o2 l5 -plot color=black \
		-read ex_anneal_conts.o2 l6 -plot color=black \
		-read ex_anneal_conts.o2 l7 -plot color=black \
		-read ex_anneal_conts.o2 l8 -plot color=black \
		-read ex_anneal_conts.o2 l9 -plot color=black \
		-point 6 7 marker=s,color=black \
		-text 6 6.5 "Initial point" color=black \
		-point 2 "(-3)" marker=s,color=white \
		-text 2 "(-2.5)" "Global minimum" color=white \
		-save ex_anneal_plot_new.png
	rm -f ex_anneal_conts.o2
	convert +append ex_anneal_plot.png ex_anneal_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_contour_plot1:
	o2graph -backend Agg -read ../ex_contour.o2 c1_cl \
		-set xtitle x -set ytitle y \
		-plot -save ex_contour_plot1_new.png
	convert +append ex_contour_plot1.png ex_contour_plot1_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_contour_plot2:
	o2graph -backend Agg -read ../ex_contour.o2 c2_cl \
		-set xtitle x -set ytitle y \
		-plot -save ex_contour_plot2_new.png
	convert +append ex_contour_plot2.png ex_contour_plot2_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_chebapp_plot: 
	o2graph -backend Agg -generic table ../ex_chebapp.out \
		-plot c1 c2 color=black -plot c1 c4 color=red \
		-plot c1 c5 color=blue,ls=: \
		-line 0.8 "(-0.5)" 0.95 "(-0.5)" color=black \
		-text 0.78 "(-0.5)" Exact ha=right \
		-line 0.8 "(-0.65)" 0.95 "(-0.65)" color=red \
		-text 0.78 "(-0.65)" "Approx. (n=50)" ha=right,color=red \
		-line 0.8 "(-0.8)" 0.95 "(-0.8)" color=blue \
		-text 0.78 "(-0.8)" "Approx. (n=25)" ha=right,color=blue \
		-ttext 0.5 "(-0.1)" x \
		-save ex_chebapp_plot_new.png
	convert +append ex_chebapp_plot.png ex_chebapp_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_eos_had_eos:
	o2graph -backend agg -read ../ex_eos_had_apr_nstar.o2 \
		-set xtitle \
		"$$ \varepsilon (\mathrm{MeV}/\mathrm{fm}^{3}) $$" \
		-set ytitle \
		"$$ P (\mathrm{MeV}/\mathrm{fm}^{3}) $$" \
		-function "ed*197.33" edhc -function "pr*197.33" prhc \
		-set xlo 0 -set xhi 2000 -set ylo 0.1 -set yhi 2000 \
		-set logy 1 -plot edhc prhc \
		-read ../ex_eos_had_skyrme_SLy4_eos.o2 \
		-function "ed*197.33" edhc -function "pr*197.33" prhc \
		-plot edhc prhc -save ex_eos_had_eos_new.png
	convert +append ex_eos_had_eos.png ex_eos_had_eos_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_eos_had_gibbs:
	o2graph -backend agg -read ../ex_eos_had_apr_nstar.o2 \
		-set xtitle \
		"$$ \varepsilon (\mathrm{MeV}/\mathrm{fm}^{3}) $$" \
		-set ytitle \
		"$$ P (\mathrm{MeV}/\mathrm{fm}^{3}) $$" \
		-function "ed*197.33" edhc -function "pr*197.33" prhc \
		-set xlo 0 -set xhi 2000 -set ylo 0.1 -set yhi 2000 \
		-set logy 1 -plot edhc prhc \
		-read ../ex_eos_had_skyrme_SLy4_eos.o2 \
		-function "ed*197.33" edhc -function "pr*197.33" prhc \
		-plot edhc prhc -save ex_eos_had_gibbs_new.png
	convert +append ex_eos_had_gibbs.png ex_eos_had_gibbs_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_eos_had_mvsr:
	o2graph -backend agg -read ../ex_eos_had_apr_mvsr.o2 \
		-set xtitle \
		"$$ R~(\mathrm{km}) $$" \
		-set ytitle \
		"$$ M~(\mathrm{M}_{\odot}) $$" \
		-set xlo 9 -set xhi 16 -set ylo 0 -set yhi 2.3 \
		-plot r gm \
		-read ../ex_eos_had_skyrme_SLy4_mvsr.o2 \
		-plot r gm -save ex_eos_had_mvsr_new.png
	convert +append ex_eos_had_mvsr.png ex_eos_had_mvsr_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_fptr_plot:
	cat ../ex_fptr.scr | grep -v t > exftemp1
	acol -generic table ../ex_fptr.out -internal exftemp2.o2
	o2graph -backend Agg -set xlo "(-1)" -set xhi 2 \
		-read exftemp2.o2 -plot x y color=red \
		-ttext 0.5 "(-0.1)" x -ttext "(-0.1)" 0.5 y \
		-generic table exftemp1 -plot c1 c2 marker=s \
		-save ex_fptr_plot_new.png
	rm -f exftemp1 exftemp2.o2
	convert +append ex_fptr_plot.png ex_fptr_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_mcmc_plot:
	o2graph -backend Agg \
		-set fig_dict "fig_size_x=6.0,fig_size_y=4.0" \
		-subplots 1 2 -subadj \
		"left=0.12,right=0.99,top=0.99,bottom=0.11,wspace=0.38" \
		-selax 0 -read ../ex_mcmc.o2 indep \
		-hist-plot x bins=20 -create table x "grid:-5,5,0.1" \
		-function "1000*exp(-x^2)*(sin(x-1.4)+1)" y \
		-plot x y "lw=2" -xtitle "x" -ytitle "y" \
		-selax 1 -read ../ex_mcmc_kde.o2 indep \
		-hist-plot x bins=20 -create table x "grid:-5,5,0.1" \
		-function "5500*exp(-x^2)*(sin(x-1.4)+1)" y \
		-plot x y "lw=2" -xtitle "x" -ytitle "y" \
		-save ex_mcmc_plot_new.png
	open ex_mcmc_plot_new.png

#convert +append ex_mcmc_plot.png ex_mcmc_plot_new.png \
#		temp.png
#	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

LEG_X1 = 2.0
LEG_X2 = 2.35
LEG_X3 = 2.45
LEG_Y1 = 7.5
LEG_Y2 = 8.25
LEG_Y3 = 9.0
LEG_Y4 = 9.75
LEG_Y5 = 10.5
LEG_Y6 = 11.25
LEG_Y6b = 6.75
LEG_Y6c = 6.0

ex_mmin_plot:
	o2graph -set xlo "(-2)" -set xhi 5.5 -set ylo "(-1)" -set yhi 15 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle z \
		-ttext 0.5 1.03 "Minimizer trajectories in x-z plane" \
		-generic table ../ex_mmin6.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=paleturquoise,ls=-- \
		-generic table ../ex_mmin5.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=navajowhite,ls=-- \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=green,ls=: \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=red \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=blue \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=magenta,ls=-. \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=brown,ls=-- \
		-generic table ../ex_mmin4g.dat -internal ex_mmin.o2 \
		-plot c1 c3 color=black,ls=-- \
		-rect 1.9 5.6 5.45 11.65 0 lw=0,fc=white,fill=True,zorder=20 \
		-line $(LEG_X1) $(LEG_Y1) $(LEG_X2) $(LEG_Y1) \
			color=black,ls=:,zorder=21 \
		-line $(LEG_X1) $(LEG_Y2) $(LEG_X2) $(LEG_Y2) \
			color=brown,ls=--,zorder=22 \
		-line $(LEG_X1) $(LEG_Y3) $(LEG_X2) $(LEG_Y3) \
			color=magenta,ls=-.,zorder=28 \
		-line $(LEG_X1) $(LEG_Y4) $(LEG_X2) $(LEG_Y4) \
			color=blue,zorder=27 \
		-line $(LEG_X1) $(LEG_Y5) $(LEG_X2) $(LEG_Y5) \
			color=red,zorder=26 \
		-line $(LEG_X1) $(LEG_Y6) $(LEG_X2) $(LEG_Y6) \
			color=green,ls=:,zorder=23 \
		-line $(LEG_X1) $(LEG_Y6b) $(LEG_X2) $(LEG_Y6b) \
			color=sandybrown,ls=:,zorder=24 \
		-line $(LEG_X1) $(LEG_Y6c) $(LEG_X2) $(LEG_Y6c) \
			color=mediumturquoise,ls=:,zorder=25 \
		-text $(LEG_X3) $(LEG_Y1) "BFGS with gradient" \
			color=black,ha=left,zorder=36 \
		-text $(LEG_X3) $(LEG_Y2) "PR with gradient" \
			color=brown,ha=left,zorder=35 \
		-text $(LEG_X3) $(LEG_Y3) "Polak-Ribere" \
			color=magenta,ha=left,zorder=34 \
		-text $(LEG_X3) $(LEG_Y4) "FR with gradient" \
			color=blue,ha=left,zorder=33 \
		-text $(LEG_X3) $(LEG_Y5) "Fletcher-Reeves" \
			color=red,ha=left,zorder=32 \
		-text $(LEG_X3) $(LEG_Y6) "Nelder-Mead" \
			color=green,ha=left,zorder=31 \
		-text $(LEG_X3) $(LEG_Y6b) "DE" \
			color=sandybrown,ha=left,zorder=30 \
		-text $(LEG_X3) $(LEG_Y6c) "Adaptive DE" \
			color=mediumturquoise,ha=left,zorder=29 \
		-point 1 0 marker=o,color=black \
		-text 0.8 0.0 "Actual minimum" ha=right \
		-save ex_mmin_plot_new.png
	convert +append ex_mmin_plot.png ex_mmin_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

LEG_X4 = 0.8
LEG_X5 = 1.15
LEG_X6 = 1.25

LEG_Y7 = 2.0
LEG_Y8 = 2.3
LEG_Y9 = 2.6
LEG_Y10 = 2.9
LEG_Y11 = 3.2
LEG_Y12 = 3.5

ex_mmin_plot2:
	o2graph -set xlo "(-2)" -set xhi 3.5 -set ylo "(-2)" -set yhi 3.8 \
		-set fig-dict "top_margin=0.05" -set xtitle x \
		-set ytitle y \
		-ttext 0.5 1.03 "Minimizer trajectories in x-y plane" \
		-generic table ../ex_mmin1.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=green,ls=: \
		-line $(LEG_X4) $(LEG_Y7) $(LEG_X5) $(LEG_Y7) color=green,ls=: \
		-text $(LEG_X6) $(LEG_Y7) "Nelder-Mead" color=green,ha=left \
		-generic table ../ex_mmin2.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=red \
		-line $(LEG_X4) $(LEG_Y8) $(LEG_X5) $(LEG_Y8) color=red \
		-text $(LEG_X6) $(LEG_Y8) "Fletcher-Reeves" color=red,ha=left \
		-generic table ../ex_mmin2g.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=blue \
		-line $(LEG_X4) $(LEG_Y9) $(LEG_X5) $(LEG_Y9) color=blue \
		-text $(LEG_X6) $(LEG_Y9) "FR with gradient" color=blue,ha=left \
		-generic table ../ex_mmin3.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=magenta,ls=-. \
		-line $(LEG_X4) $(LEG_Y10) $(LEG_X5) $(LEG_Y10) color=magenta,ls=-. \
		-text $(LEG_X6) $(LEG_Y10) "Polak-Ribere" color=magenta,ha=left \
		-generic table ../ex_mmin3g.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=brown,ls=-- \
		-line $(LEG_X4) $(LEG_Y11) $(LEG_X5) $(LEG_Y11) color=brown,ls=-- \
		-text $(LEG_X6) $(LEG_Y11) "PR with gradient" color=brown,ha=left \
		-generic table ../ex_mmin4g.dat -internal ex_mmin.o2 \
		-plot c1 c2 color=black,ls=-- \
		-line $(LEG_X4) $(LEG_Y12) $(LEG_X5) $(LEG_Y12) color=black,ls=-- \
		-text $(LEG_X6) $(LEG_Y12) "BFGS with gradient" color=black,ha=left \
		-point 1 0 marker=o,color=black \
		-point "(-1.6)" 3.2 marker=o,color=black \
		-text "(-1.4)" 3.2 "Actual minimum" ha=left \
		-save ex_mmin_plot2_new.png
	convert +append ex_mmin_plot2.png ex_mmin_plot2_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_nuc_prof:
	o2graph -read ../ex_nucleus_rmf_prof.o2 -set xlo 0 -set xhi 10 \
		-set ylo 0 -set yhi 0.1 \
		-set xtitle "$$ r~(\mathrm{fm})$$" -set ytitle \
		"$$ \mathrm{density}~(\mathrm{fm}^{-3}) $$" \
		-plot r rhon color=red -plot r rhop color=blue \
		-plot r chden color=black -plot r chden_exp \
		color=black,ls=-. \
		-text 4.5 0.09 "neutrons" color=red \
		-text 0.9 0.074 "protons" color=blue \
		-text 4.7 0.04 "charge" color=black \
		-line 5.4 0.04 6.0 0.0488 color=black \
		-text 3.5 0.066 "charge (expt.)" color=black \
		-save ex_nuc_prof_new.png
	convert +append ex_nuc_prof.png ex_nuc_prof_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_nucmass_plots:
	o2graph -backend Agg \
		-set xtitle N -set ytitle Z -read ../ex_nucmass_table.o2 \
		-to-table3d N Z -set colbar 1 \
		-den-plot se cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Semi-empirical} $$" \
		-save ex_nucmass_se.png \
		-clf -den-plot mnmsk cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{M\ddot{o}ller~et~al.~(1995)} $$" \
		-save ex_nucmass_mnmsk.png \
		-clf -den-plot hfb14 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~14} $$" \
		-save ex_nucmass_hfb14.png \
		-clf -den-plot hfb21 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~21} $$" \
		-save ex_nucmass_hfb21.png \
		-clf -den-plot hfb27 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{HFB~27} $$" \
		-save ex_nucmass_hfb27.png \
		-clf -den-plot ame03 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{AME~(2003)} $$" \
		-save ex_nucmass_ame03.png \
		-clf -den-plot dz96 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Duflo~and~Zuker~(1996)} $$" \
		-save ex_nucmass_dz96.png \
		-clf -den-plot ktuy05 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Koura~et~al.~(1996)} $$" \
		-save ex_nucmass_ktuy05.png \
		-clf -den-plot dvi cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Dieperink~et~al.~(2009)} $$" \
		-save ex_nucmass_dvi.png \
		-clf -den-plot ws32 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Wang~et~al.~(2010)} $$" \
		-save ex_nucmass_ws32.png \
		-clf -den-plot ws36 cmap=PuOr \
		-ttext 0.6 0.95 "$$ \mathrm{Liu~et~al.~(2011)} $$" \
		-save ex_nucmass_ws36.png 
	convert +append ex_anneal_plot.png ex_anneal_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_ode_airy: 
	o2graph -read ../ex_ode.o2 table_2 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err)" abserr \
		-plot x abserr -plot x diff \
		-read ../ex_ode.o2 table_3 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err)" abserr \
		-plot x abserr -plot x diff \
		-ttext 0.57 0.14 "Prince-Dormand act. error" \
		-ttext 0.25 0.96 "Cash-Karp act. error" \
		-ttext 0.60 0.86 "Cash-Karp est. error" \
		-ttext 0.5 0.21 "Prince-Dormand est. error" \
		-ttext 0.5 1.03 "Airy function w/non-adaptive steppers" \
		-save ex_ode_airy_plot_new.png
	convert +append ex_ode_airy.png ex_ode_airy_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_ode_bessel: 
	o2graph -backend Agg -read ../ex_ode.o2 table_0 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-read ../ex_ode.o2 table_1 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-plot x err -plot x diff \
		-ttext 0.46 0.87 "Prince-Dormand act. error" \
		-ttext 0.70 0.67 "Cash-Karp act. error" \
		-ttext 0.60 0.48 "Cash-Karp est. error" \
		-ttext 0.48 0.09 "Prince-Dormand est. error" \
		-ttext 0.5 1.03 "Bessel function w/non-adaptive steppers" \
		-save ex_ode_bessel_new.png
	convert +append ex_ode_bessel.png ex_ode_bessel_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_ode_bessel2: 
	o2graph -backend Agg -read ../ex_ode.o2 table_4 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-11 -set yhi 2.5e-6 \
		-function "abs(calc-exact)" diff -set logy 1 \
		-function "abs(err0)" abserr \
		-plot x abserr color=black,ls=-. -plot x diff color=black \
		-read ../ex_ode.o2 table_5 \
		-function "abs(calc-exact)" diff \
		-function "abs(err0)" abserr \
		-plot x abserr color=red,ls=-. -plot x diff color=red \
		-read ../ex_ode.o2 table_6 \
		-function "abs(calc-exact)" diff \
		-function "abs(err0)" abserr \
		-plot x abserr color=blue,ls=-. -plot x diff color=blue \
		-ttext 0.45 0.97 "Cash-Karp act. error" color=black \
		-ttext 0.22 0.92 "Cash-Karp est. error" color=black \
		-ttext 0.55 0.65 "Cash-Karp(2) act. error" color=red \
		-ttext 0.30 0.60 "Prince-Dormand act. error" color=blue \
		-ttext 0.70 0.10 "Cash-Karp(2) est. error" color=red \
		-ttext 0.35 0.05 "Prince-Dormand est. error" color=blue \
		-ttext 0.5 1.03 "Bessel function w/adaptive steppers" \
		-save ex_ode_bessel2_new.png
	convert +append ex_ode_bessel2.png ex_ode_bessel2_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_ode_bessel3: 
	o2graph -backend Agg -read ../ex_ode.o2 table_7 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-11 -set yhi 5.0e-8 \
		-set xlo 0 -set xhi 10 -set logy 1 \
		-function "abs(err0)" abserr \
		-plot x abserr color=red,ls=-. -plot x diff color=black \
		-ttext 0.45 0.95 "Act. error" color=black \
		-ttext 0.50 0.50 "Est. error" color=red \
		-ttext 0.5 1.03 \
		"Bessel function, high-level adaptive stepper" \
		-save ex_ode_bessel3_new.png
	convert +append ex_ode_bessel3.png ex_ode_bessel3_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

ex_stiff_plot:
	o2graph -backend Agg -read ../ex_stiff.o2 table_0 \
		-set xtitle x -set ytitle y -set fig-dict "top_margin=0.05" \
		-set ylo 1.0e-16 -set yhi 1.0e8 -set logy 1 \
		-set xlo 0 -set xhi 4 -plot x rel_err color=black \
		-plot x rel_diff color=red \
		-read ../ex_stiff.o2 table_1 -plot x rel_err color=blue \
		-plot x rel_diff color=green \
		-ttext 0.35 0.83 \
		"Explicit actual err." color=green \
		-ttext 0.5 0.53 \
		"Explicit estimated err." color=blue \
		-ttext 0.5 0.15 \
		"BS implicit actual err." color=red \
		-ttext 0.5 0.05 \
		"BS implicit estimated err." color=black \
		-ttext 0.5 1.03 \
		"Stiff exponential function, implicit vs. explicit" \
		-save ex_stiff_plot_new.png -show
	convert +append ex_stiff.png ex_stiff_plot_new.png \
		temp.png
	$(IMAGE_VIEW) temp.png $(IMAGE_VIEW2)

# Old targets

# ex_mmin_plot.o: ex_mmin_plot.cpp
# 	$(CXX) $(CPPFLAGS) $(CPPFLAGS) -o ex_mmin_plot.o -c ex_mmin_plot.cpp

# ex_mmin_plot: ex_mmin_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_mmin_plot ex_mmin_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

# ex_stiff_plot.o: ex_stiff_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_stiff_plot.o -c ex_stiff_plot.cpp

# ex_stiff_plot: ex_stiff_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_stiff_plot ex_stiff_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl \
# 		-lo2scl_hdf

# ex_anneal_plot.o: ex_anneal_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_anneal_plot.o -c ex_anneal_plot.cpp

# ex_anneal_plot: ex_anneal_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_anneal_plot ex_anneal_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

#ex_chebapp_plot.eps: ex_chebapp_plot.py
#	python ex_chebapp_plot.py

#ex_fptr_plot.eps: ex_fptr_plot.py
#	python ex_fptr_plot.py

# ex_mandel_plot.o: ex_mandel_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_mandel_plot.o -c ex_mandel_plot.cpp

# ex_mandel_plot: ex_mandel_plot.o  ../../src/other/graph.o
# 	$(CXX) -o ex_mandel_plot ex_mandel_plot.o \
# 		$(ROOT_LIB) $(CPPFLAGS) $(LDFLAGS) -lo2scl -lo2scl_hdf

# ex_mandel_plot_all: ex_mandel_plot
# 	./ex_mandel_plot
# 	convert ex_mandel_plot.png ex_mandel_plot.pdf

#ex_ode_bessel.eps: ex_ode_plot.py
#	python ex_ode_plot.py

# ex_planar_plot.o: ex_planar_plot.cpp
# 	$(CXX) $(CPPFLAGS) -o ex_planar_plot.o -c ex_planar_plot.cpp

# ex_planar_plot: ex_planar_plot.o ../../src/other/graph.o
# 	$(CXX) $(CPPFLAGS) -o ex_planar_plot ex_planar_plot.o \
# 		-lo2scl_hdf -lo2scl $(ROOT_LIB) $(LDFLAGS) 

# ex_planar_plot.pdf: ex_planar_plot
# 	./ex_planar_plot
# 	convert ex_planar_plot.png ex_planar_plot.pdf

# ex_ublas.o: ex_ublas.cpp
# 	$(CXX) -I/home/asteiner/install/include $(INCLUDES) \
# 		-o ex_ublas.o -c ex_ublas.cpp

# ex_ublas: ex_ublas.o  ../../src/other/graph.o
# 	$(CXX) -I/home/asteiner/install/include -L/usr/lib64/atlas \
# 		-L/home/asteiner/install/lib -o ex_ublas ex_ublas.o -lgsl \
# 		-lo2scl

