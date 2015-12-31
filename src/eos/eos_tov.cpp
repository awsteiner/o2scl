/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_tov.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

// ----------------------------------------------------------------
// eos_tov functions

eos_tov::eos_tov() {
  baryon_column=false;
  verbose=1;
}

void eos_tov::check_nb(double &avg_abs_dev, double &max_abs_dev) {
  if (!baryon_column) {
    O2SCL_ERR2("Variable 'baryon_column' false in",
	       "eos_tov::check_nb().",exc_einval);
  }
  std::vector<double> edv, prv, nbv, dedn;
  for (double pres=0.1;pres<3.0;pres*=1.001) {
    double eps, nb;
    ed_nb_from_pr(pres,eps,nb);
    edv.push_back(eps);
    prv.push_back(pres);
    nbv.push_back(nb);
  }
  dedn.resize(edv.size());
  vector_deriv_interp(edv.size(),nbv,edv,dedn,itp_linear);
  avg_abs_dev=0.0;
  max_abs_dev=0.0;
  for(size_t i=0;i<edv.size();i++) {
    double abs_dev=(fabs(edv[i]+prv[i]-dedn[i]*nbv[i])/
		    fabs(dedn[i]*nbv[i]));
    if (abs_dev>max_abs_dev) max_abs_dev=abs_dev;
    avg_abs_dev+=abs_dev;
  }
  avg_abs_dev/=((double)edv.size());
  return;
}

// ----------------------------------------------------------------
// eos_tov_buchdahl functions

eos_tov_buchdahl::eos_tov_buchdahl() {
  Pstar=3.2e-5;
}

void eos_tov_buchdahl::set_baryon_density(double nb, double ed) {
  if (nb<0.0 || ed<0.0) {
    O2SCL_ERR2("Negative densities not supported in ",
	       "eos_tov_buchdahl::set_coeff_index().",exc_einval);
  }
  baryon_column=true;
  nb1=nb;
  ed1=ed;
  if (36.0*Pstar*Pstar-5.0*Pstar*ed<0.0) {
    O2SCL_ERR2("Values of 'Pstar' and 'ed' incommensurate in ",
	       "eos_tov_buchdahl::set_baryon_density().",exc_einval);
  }
  pr1=0.04*(72.0*Pstar-5.0*ed+
	    12.0*sqrt(36.0*Pstar*Pstar-5.0*Pstar*ed));
  return;
}

void eos_tov_buchdahl::ed_nb_from_pr(double pr, double &ed, double &nb) {
  ed=12.0*sqrt(Pstar*pr)-5.0*pr;
  if (baryon_column) {
    double mu1=(pr1+ed1)/nb1;
    double t1=sqrt(pr/Pstar);
    double t2=sqrt(pr1/Pstar);
    double mu=mu1*pow((-pr1+9.0*Pstar)*(3.0+t1)*(3.0-t2)/
		      (-pr+9.0*Pstar)/(3.0-t1)/(3.0+t2),0.25);
    nb=(pr+ed)/mu;
  } else {
    nb=0.0;
  }
  return;
}

int eos_tov_buchdahl::solve_u_rp_fun
(size_t bv, const std::vector<double> &bx, std::vector<double> &by) {
  
  double u, rp;
  u=bx[0];
  rp=bx[1];
  //by[0]=rp*(1.0-beta+u)/(1.0-2.0*beta)-buchrad;
  //by[1]=beta/biga/rp*sin(biga*rp);
  return 0;
}

double eos_tov_buchdahl::ed_from_pr(double pr) {
  return 12.0*sqrt(Pstar*pr)-5.0*pr;
}

double eos_tov_buchdahl::pr_from_ed(double ed) {
  return 0.0;
}

double eos_tov_buchdahl::nb_from_ed(double ed) {
  return 0.0;
}

double eos_tov_buchdahl::nb_from_pr(double pr) {
  double ed=12.0*sqrt(Pstar*pr)-5.0*pr;
  double mu1=(pr1+ed1)/nb1;
  double t1=sqrt(pr/Pstar);
  double t2=sqrt(pr1/Pstar);
  double mu=mu1*pow((-pr1+9.0*Pstar)*(3.0+t1)*(3.0-t2)/
		    (-pr+9.0*Pstar)/(3.0-t1)/(3.0+t2),0.25);
  return (pr+ed)/mu;
}

double eos_tov_buchdahl::ed_from_nb(double nb) {
  return 0.0;
}

double eos_tov_buchdahl::pr_from_nb(double nb) {
  return 0.0;
}

// ----------------------------------------------------------------
// eos_tov_polytrope functions

eos_tov_polytrope::eos_tov_polytrope() {
  K=1.0;
  n=3.0;
}

void eos_tov_polytrope::set_coeff_index(double coeff, double index) {
  if (coeff<0.0 || index<0.0) {
    O2SCL_ERR2("Negative coefficients and indices not supported in ",
	       "eos_tov_polytrope::set_coeff_index().",exc_einval);
  }
  K=coeff;
  n=index;
  if (baryon_column) {
    pr1=K*pow(ed1,1.0+1.0/n);
  }
  return;
}

void eos_tov_polytrope::set_baryon_density(double nb, double ed) {
  if (nb<0.0 || ed<0.0) {
    O2SCL_ERR2("Negative densities not supported in ",
	       "eos_tov_polytrope::set_baryon_density().",exc_einval);
  }
  baryon_column=true;
  nb1=nb;
  ed1=ed;
  pr1=K*pow(ed1,1.0+1.0/n);
  return;
}

double eos_tov_polytrope::ed_from_pr(double pr) {
  return pow(pr/K,n/(1.0+n));
}

double eos_tov_polytrope::pr_from_ed(double ed) {
  return K*pow(ed,1.0+1.0/n);
}

double eos_tov_polytrope::nb_from_ed(double ed) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_polytrope::nb_from_ed().",exc_einval);
  }
#endif
  double pr=K*pow(ed,1.0+1.0/n);
  return nb1*pow(ed/ed1,1.0+n)/pow((ed+pr)/(ed1+pr1),n);
}

double eos_tov_polytrope::nb_from_pr(double pr) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_polytrope::nb_from_ed().",exc_einval);
  }
#endif
  double ed=pow(pr/K,n/(1.0+n));
  return nb1*pow(ed/ed1,1.0+n)/pow((ed+pr)/(ed1+pr1),n);
}

double eos_tov_polytrope::ed_from_nb(double nb) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_polytrope::nb_from_ed().",exc_einval);
  }
#endif
  return pow(pow(nb1/nb/ed1,1.0/n)*(1.0+pr1/ed1)-K,-n);
}

double eos_tov_polytrope::pr_from_nb(double nb) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_polytrope::nb_from_ed().",exc_einval);
  }
#endif
  return K*pow(pow(nb1/nb/ed1,1.0/n)*(1.0+pr1/ed1)-K,-(n+1.0));
}

void eos_tov_polytrope::ed_nb_from_pr(double pr, double &ed, double &nb) {
  ed=pow(pr/K,n/(1.0+n));
  if (baryon_column) {
    nb=nb1*pow(ed/ed1,1.0+n)/pow((ed+pr)/(ed1+pr1),n);
  } else {
    nb=0.0;
  }
  return;
}

// ----------------------------------------------------------------
// eos_tov_linear functions

eos_tov_linear::eos_tov_linear() {
  cs2=1.0;
  eps0=0.0;
  nb1=0.0;
}

void eos_tov_linear::set_cs2_eps0(double cs2_, double eps0_) {
  eps0=eps0_;
  cs2=cs2_;
  return;
}

void eos_tov_linear::set_baryon_density(double nb, double ed) {
  if (nb<=0.0 || ed<=0.0) {
    O2SCL_ERR2("Negative and zero densities not supported in ",
	       "eos_tov_linear::set_coeff_index().",exc_einval);
  }
  baryon_column=true;
  nb1=nb;
  ed1=ed;
  pr1=cs2*(ed1-eps0);
  return;
}

double eos_tov_linear::ed_from_pr(double pr) {
  return pr/cs2+eps0;
}

double eos_tov_linear::pr_from_ed(double ed) {
  return (ed-eps0)*cs2;
}

double eos_tov_linear::nb_from_ed(double ed) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_linear::nb_from_ed().",exc_einval);
  }
#endif
  return nb1*pow(ed+cs2*ed-cs2*eps0,1.0/(1.0+cs2))*
    pow(ed1+cs2*(-eps0+ed1),-1.0/(1.0+cs2));
}

double eos_tov_linear::nb_from_pr(double pr) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_linear::nb_from_ed().",exc_einval);
  }
#endif
  double ed=pr/cs2+eps0;
  return nb1*pow(ed+cs2*ed-cs2*eps0,1.0/(1.0+cs2))*
    pow(ed1+cs2*(-eps0+ed1),-1.0/(1.0+cs2));
}

double eos_tov_linear::ed_from_nb(double nb) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_linear::nb_from_ed().",exc_einval);
  }
#endif
  double ret=(pow(nb/pow(ed1+cs2*(-eps0+ed1),-1.0/(1.0+cs2))/nb1,
		  1.0+cs2)+cs2*eps0)/(1.0+cs2);
  return ret;
}

double eos_tov_linear::pr_from_nb(double nb) {
#if !O2SCL_NO_RANGE_CHECK
  if (nb1==0.0) {
    O2SCL_ERR2("Fiducial baryon density not specified in ",
	       "eos_tov_linear::nb_from_ed().",exc_einval);
  }
#endif
  double ed=(pow(nb/pow(ed1+cs2*(-eps0+ed1),-1.0/(1.0+cs2))/nb1,
		 1.0+cs2)+cs2*eps0)/(1.0+cs2);
  return (ed-eps0)*cs2;
}

void eos_tov_linear::ed_nb_from_pr(double pr, double &ed, double &nb) {
  ed=pr/cs2+eps0;
  if (baryon_column) {
    nb=nb1*pow(ed+cs2*ed-cs2*eps0,1.0/(1.0+cs2))*
      pow(ed1+cs2*(-eps0+ed1),-1.0/(1.0+cs2));
  } else {
    nb=0.0;
  }
  return;
}

// ----------------------------------------------------------------
// eos_tov_interp functions

eos_tov_interp::eos_tov_interp() {
  
  eos_read=false;
  use_crust=false;
  verbose=1;

  cole=0;
  colp=0;
  colnb=0;

  core_auxp=0;

  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;
  
  trans_width=1.0;
  core_table=0;
  core_set=false;

  transition_mode=smooth_trans;

  gen_int.set_type(itp_linear);
}

eos_tov_interp::~eos_tov_interp() {
}

void eos_tov_interp::get_names_units(size_t &np, 
				     std::vector<std::string> &pnames,
				     std::vector<std::string> &vs_units) {
  np=0;
  if (core_auxp>0) {
    for(int i=0;i<((int)core_table->get_ncolumns());i++) {
      if (i!=cole && i!=colp && i!=colnb) {
	np++;
	pnames.push_back(core_table->get_column_name(i));
	vs_units.push_back
	  (core_table->get_unit(core_table->get_column_name(i)));
      }
    }
  }
  return;
}

void eos_tov_interp::read_vectors(size_t n_core, std::vector<double> &core_ed, 
				  std::vector<double> &core_pr) {

  core_table=0;
  full_nlines=n_core;
  pfactor=1.0;
  efactor=1.0;
  std::swap(core_ed,full_vece);
  std::swap(core_pr,full_vecp);
  core_auxp=0;
  baryon_column=false;
  eos_read=true;

  pe_int.set(n_core,full_vecp,full_vece,itp_linear);

  return;
}

void eos_tov_interp::read_vectors(size_t n_core, std::vector<double> &core_ed, 
				  std::vector<double> &core_pr, 
				  std::vector<double> &core_nb) {

  read_vectors(n_core,core_ed,core_pr);
  efactor=1.0;
  nfactor=1.0;
  std::swap(core_nb,full_vecnb);
  baryon_column=true;

  pnb_int.set(n_core,full_vecp,full_vecnb,itp_linear);

  return;
}

void eos_tov_interp::read_table(table_units<> &eosat, string s_cole, 
				string s_colp, string s_colnb) {
  
  core_table=&eosat;
  size_t core_nlines=core_table->get_nlines();

  if (verbose>1) cout << "Lines read from EOS file: " << core_nlines << endl;
  
  // ---------------------------------------------------------------
  // Look for energy density, pressure, and baryon density columns

  cole=core_table->lookup_column(s_cole);
  colp=core_table->lookup_column(s_colp);
  if (s_colnb!="") {
    colnb=core_table->lookup_column(s_colnb);
    baryon_column=true;
  } else {
    baryon_column=false;
  }

  // ---------------------------------------------------------------
  // Set core_auxp

  if (baryon_column) core_auxp=core_table->get_ncolumns()-3;
  else core_auxp=core_table->get_ncolumns()-2;

  // ---------------------------------------------------------------
  // Determine conversion factors from units
 
  string eunits=core_table->get_unit(s_cole);
  string punits=core_table->get_unit(s_colp);
  if (eunits.length()==0 || eunits=="Msun/km^3" ||
      eunits=="solarmass/km^3") {
    efactor=1.0;
  } else {
    efactor=o2scl_settings.get_convert_units().convert
      (eunits,"Msun/km^3",1.0);
  }

  if (punits.length()==0 || punits=="Msun/km^3" ||
      punits=="solarmass/km^3") {
    pfactor=1.0;
  } else {
    pfactor=o2scl_settings.get_convert_units().convert
      (punits,"Msun/km^3",1.0);
  }

  nfactor=1.0;
  if (baryon_column) {
    string nunits=core_table->get_unit(s_colnb);
    if (nunits=="1/cm^3") {
      nfactor=1.0e-39;
    } else if (nunits=="1/m^3") {
      nfactor=1.0e-42;
    } else if (nunits=="1/fm^3") {
      nfactor=1.0;
    } else if (nunits=="") {
      nfactor=1.0;
    } else {
      nfactor=o2scl_settings.get_convert_units().convert
	(nunits,"1/fm^3",1.0);
    }
  }

  core_set=true;

  internal_read();

  return;
}

void eos_tov_interp::internal_read() {

  if (!core_set) {
    O2SCL_ERR("Core table not set in eos_tov_interp::internal_read().",
	      exc_esanity);
  }

  size_t core_nlines=core_table->get_nlines();

  // ---------------------------------------------------------------
  // Construct full vectors

  double pr_lo=trans_pres/trans_width;
  double pr_hi=trans_pres*trans_width;

  full_vece.clear();
  full_vecp.clear();
  if (baryon_column) {
    full_vecnb.clear();
  }

  // Add lines of crust before transition
  if (verbose>1) {
    cout << "Crust: " << endl;
  }
  if (use_crust) {
    for(size_t i=0;i<crust_nlines;i++) {
      if (crust_vecp[i]<pr_lo) {
	full_vece.push_back(crust_vece[i]);
	full_vecp.push_back(crust_vecp[i]);
	if (baryon_column) {
	  full_vecnb.push_back(crust_vecnb[i]);
	}
	if (verbose>1) {
	  cout << crust_vece[i] << " " << crust_vecp[i] << " ";
	  if (baryon_column) {
	    cout << crust_vecnb[i] << " ";
	  }
	  cout << crust_vece[i]/efactor << " " << crust_vecp[i]/pfactor << " ";
	  if (baryon_column) {
	    cout << crust_vecnb[i]/nfactor;
	  }
	  cout << endl;
	}
      }
    }
  }
  
  if (trans_width>1.0) {

    // Add transition lines
    if (verbose>1) {
      cout << "Transition: " << endl;
    }
    // Add transition lines
    double dpr=(pr_hi-pr_lo)/20.0;
    double ed_lo, ed_hi, nb_lo, nb_hi;

    if (transition_mode==match_line) {
      ed_lo=gen_int.eval(pr_lo,crust_nlines,crust_vecp,crust_vece);
      ed_hi=gen_int.eval(pr_hi/pfactor,core_nlines,(*core_table)[colp],
			 (*core_table)[cole])*efactor;
      nb_lo=gen_int.eval(pr_lo,crust_nlines,crust_vecp,crust_vecnb);
      nb_hi=gen_int.eval(pr_hi/pfactor,core_nlines,(*core_table)[colp],
			 (*core_table)[colnb])*nfactor;
    }

    double ed, nb=0.0;

    for(double pr=pr_lo;pr<pr_hi+dpr/100.0;pr+=dpr) {

      double chi=(pr-pr_lo)/(pr_hi-pr_lo);

      if (transition_mode==smooth_trans) {
	ed_lo=gen_int.eval(pr,crust_nlines,crust_vecp,crust_vece);
	ed_hi=gen_int.eval(pr/pfactor,core_nlines,(*core_table)[colp],
			   (*core_table)[cole])*efactor;
	ed=(1.0-chi)*ed_lo+chi*ed_hi;
      
	if (baryon_column) {
	  nb_lo=gen_int.eval(pr,crust_nlines,crust_vecp,crust_vecnb);
	  nb_hi=gen_int.eval(pr/pfactor,core_nlines,(*core_table)[colp],
			     (*core_table)[colnb])*nfactor;
	  nb=(1.0-chi)*nb_lo+chi*nb_hi;
	}
      
      } else {
      
	ed=(ed_hi-ed_lo)*chi+ed_lo;
	if (baryon_column) {
	  nb=(nb_hi-nb_lo)*chi+nb_lo;
	}
      
      }
    
      full_vecp.push_back(pr);
      full_vece.push_back(ed);
      if (baryon_column) {
	full_vecnb.push_back(nb);
      }
      if (verbose>1) {
	cout << ed << " " << pr << " ";
	if (baryon_column) {
	  cout << nb << " ";
	}
	cout << ed/efactor << " " << pr/pfactor << " ";
	if (baryon_column) {
	  cout << nb/nfactor;
	}
	cout << endl;
      }

    }    

  }

  // Add core lines
  if (verbose>1) {
    cout << "Core: " << endl;
  }
  for(size_t i=0;i<core_nlines;i++) {
    double pt=((*core_table)[colp])[i]*pfactor;
    double et=((*core_table)[cole])[i]*efactor;
    double nt=((*core_table)[colnb])[i]*nfactor;
    if (pt>pr_hi) {
      full_vece.push_back(et);
      full_vecp.push_back(pt);
      if (baryon_column) {
	full_vecnb.push_back(nt);
      }
      if (verbose>1) {
	cout << et << " " << pt << " ";
	if (baryon_column) {
	  cout << nt << " ";
	}
	cout << ((*core_table)[cole])[i] << " ";
	cout << ((*core_table)[colp])[i];
	if (baryon_column) {
	  cout << " " << ((*core_table)[colnb])[i];
	}
	cout << endl;
      }
    }
  }

  // ---------------------------------------------------------------
  // Set 'eos_read', interpolators, and 'full_nlines'

  full_nlines=full_vecp.size();
  pe_int.set(full_nlines,full_vecp,full_vece,itp_linear);
  if (baryon_column) {
    pnb_int.set(full_nlines,full_vecp,full_vecnb,itp_linear);
  }

  eos_read=true;

  return;
}

void eos_tov_interp::get_transition(double &p, double &w) {
  p=trans_pres/pfactor;
  w=trans_width;
  return;
}

void eos_tov_interp::set_transition(double p, double wid) {
  trans_pres=p*pfactor;
  trans_width=wid;
  if (trans_width<1.0) {
    O2SCL_ERR("Width less than 1 in eos_tov_interp::set_transition().",
	      exc_einval);
  }
  return;
}

void eos_tov_interp::default_low_dens_eos() {

  // Read default EOS. 
  static const size_t nlines=73;
  double ed_arr[nlines]=
    {3.89999984e-18,3.93000002e-18,3.95000001e-18,4.07499982e-18,
     5.80000020e-18,8.20000023e-18,2.25500006e-17,1.05999999e-16,
     5.75000005e-16,5.22000020e-15,1.31099999e-14,3.29349991e-14,
     8.27000008e-14,2.07800004e-13,5.21999978e-13,1.31100003e-12,
     3.29400006e-12,4.14649998e-12,8.27499961e-12,1.65100000e-11,
     3.29450009e-11,6.57500027e-11,1.31199995e-10,1.65200006e-10,
     2.61849986e-10,4.15049994e-10,5.22499988e-10,6.57999988e-10,
     8.28499991e-10,1.31300004e-09,2.08199991e-09,3.30050010e-09,
     4.15599999e-09,5.22999999e-09,6.59000010e-09,8.29500024e-09,
     1.04500000e-08,1.31549998e-08,1.65650000e-08,2.08599999e-08,
     2.62699995e-08,3.30850014e-08,4.16600017e-08,5.24500017e-08,
     6.61000001e-08,8.31999998e-08,9.21999970e-08,1.04800002e-07,
     1.31999997e-07,1.66250004e-07,2.09400000e-07,//2.14950006e-07,
     2.23000001e-07,2.61400004e-07,3.30500001e-07,3.98200001e-07,
     4.86399983e-07,5.97999986e-07,7.35500009e-07,//8.41528747e-07,
     //4.21516539e-06,
     8.43854053e-06,1.26672671e-05,1.69004320e-05,
     2.11374665e-05,2.53779855e-05,2.96217149e-05,3.38684539e-05,
     3.81180526e-05,4.23703981e-05,4.66254054e-05,5.08830106e-05,
     5.51431670e-05,5.94058410e-05,6.36710102e-05,6.79386612e-05};
  double pr_arr[nlines]=
    {5.64869998e-32,5.64869986e-31,5.64869986e-30,5.64870017e-29,
     6.76729990e-28,7.82989977e-27,9.50780029e-26,3.25500004e-24,
     1.06260006e-22,5.44959997e-21,2.77850000e-20,1.35959994e-19,
     6.43729996e-19,2.94519994e-18,1.29639997e-17,5.45580009e-17,
     2.18730006e-16,2.94129994e-16,8.02570017e-16,2.14369999e-15,
     5.62639983e-15,1.45639997e-14,3.73379984e-14,4.88699996e-14,
     9.11069993e-14,1.69410005e-13,2.30929994e-13,2.81650002e-13,
     3.83670011e-13,7.11400014e-13,1.31769996e-12,2.43959995e-12,
     3.16660001e-12,4.30759985e-12,5.86130016e-12,7.96970042e-12,
     1.08389998e-11,1.39989999e-11,1.90380003e-11,2.58830006e-11,
     3.32719997e-11,4.52399992e-11,6.15209966e-11,8.36119993e-11,
     1.13700001e-10,1.45250006e-10,1.61739996e-10,1.83999996e-10,
     2.50169996e-10,3.25279997e-10,4.38359987e-10,//4.36519987e-10,
     4.41269993e-10,4.67110017e-10,5.08829978e-10,5.49830015e-10,
     6.05699990e-10,6.81199985e-10,7.82430010e-10,//4.70257815e-10,
     //7.04778782e-09,
     1.45139718e-08,2.62697827e-08,4.05674724e-08,
     5.69532689e-08,7.52445638e-08,9.53657839e-08,1.17299621e-07,
     1.41064470e-07,1.66702091e-07,1.94270264e-07,2.23838165e-07,
     2.55483362e-07,2.89289800e-07,3.25346430e-07,3.63172009e-07};
  double nb_arr[nlines]=
    {4.00000001e-15,4.73000011e-15,4.75999990e-15,4.91000012e-15,
     6.99000006e-15,9.89999996e-15,2.71999999e-14,1.27000000e-13,
     6.93000019e-13,6.29500011e-12,1.58099991e-11,3.97200016e-11,
     9.97599989e-11,2.50600013e-10,6.29399977e-10,1.58100000e-09,
     3.97200006e-09,4.99999997e-09,9.97599958e-09,1.98999999e-08,
     3.97199997e-08,7.92400030e-08,1.58099994e-07,1.98999999e-07,
     3.15499989e-07,4.99999999e-07,6.29400006e-07,7.92399987e-07,
     9.97599955e-07,1.58099999e-06,2.50600010e-06,3.97199983e-06,
     4.99999987e-06,6.29399983e-06,7.92399987e-06,9.97600000e-06,
     1.25600000e-05,1.58099992e-05,1.99000006e-05,2.50600006e-05,
     3.15500001e-05,3.97199983e-05,4.99999987e-05,6.29400020e-05,
     7.92400024e-05,9.97600000e-05,1.10499997e-04,1.25599996e-04,
     1.58099996e-04,1.99000002e-04,2.50599987e-04,//2.57200008e-04,
     2.66999996e-04,3.12599994e-04,3.95100011e-04,4.75899986e-04,
     5.81200002e-04,7.14300026e-04,8.78599996e-04,//1.00000001e-03,
     //5.00000001e-03,
     1.00000000e-02,1.50000000e-02,2.00000000e-02,
     2.50000000e-02,3.00000000e-02,3.50000000e-02,4.00000000e-02,
     4.50000000e-02,5.00000000e-02,5.50000000e-02,6.00000000e-02,
     6.50000000e-02,7.00000000e-02,7.50000000e-02,8.00000000e-02};

  crust_nlines=nlines;
  crust_vece.resize(crust_nlines);
  crust_vecp.resize(crust_nlines);
  crust_vecnb.resize(crust_nlines);
  for(size_t i=0;i<nlines;i++) {
    crust_vece[i]=ed_arr[i];
    crust_vecp[i]=pr_arr[i];
    crust_vecnb[i]=nb_arr[i];
  }
  trans_pres=crust_vecp[crust_nlines-1];

  // The energy density and pressure are already in Msun/km^3 and the
  // baryon density is in fm^{-3}
  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;
    
  if (verbose>1) {
    cout << "Transition pressure: " << trans_pres << endl;
  }

  use_crust=true;
  
  if (core_set) internal_read();
    
  return;
}
    
void eos_tov_interp::sho11_low_dens_eos() {

  static const size_t nlines=98;

  double ed_arr[nlines]=
    {3.89999984e-18,3.93000002e-18,3.95000001e-18,4.07499982e-18,
     5.80000020e-18,8.20000023e-18,2.25500006e-17,1.05999999e-16,
     5.75000005e-16,5.22000020e-15,1.31099999e-14,3.29349991e-14,
     8.27000008e-14,2.07800004e-13,5.21999978e-13,1.31100003e-12,
     3.29400006e-12,4.14649998e-12,8.27499961e-12,1.65100000e-11,
     3.29450009e-11,6.57500027e-11,1.31199995e-10,1.65200006e-10,
     2.61849986e-10,4.15049994e-10,5.22499988e-10,6.57999988e-10,
     8.28499991e-10,1.31300004e-09,2.08199991e-09,3.30050010e-09,
     4.15599999e-09,5.22999999e-09,6.59000010e-09,8.29500024e-09,
     1.04500000e-08,1.31549998e-08,1.65650000e-08,2.08599999e-08,
     2.62699995e-08,3.30850014e-08,4.16600017e-08,5.24500017e-08,
     6.61000001e-08,8.31999998e-08,9.21999970e-08,1.04800002e-07,
     1.31999997e-07,
     1.778589e-07,1.995610e-07,2.239111e-07,2.512323e-07,
     2.818873e-07,3.162827e-07,3.548750e-07,5.624389e-07,6.310668e-07,
     7.080685e-07,7.944659e-07,8.914054e-07,1.000173e-06,
     1.122213e-06,1.259143e-06,1.412782e-06,1.585167e-06,
     1.778587e-06,1.995607e-06,2.239108e-06,2.512320e-06,
     2.818869e-06,3.162823e-06,3.548746e-06,3.981758e-06,
     4.467606e-06,5.012736e-06,5.624382e-06,6.310660e-06,
     7.080676e-06,7.944649e-06,8.914043e-06,1.000172e-05,
     1.122211e-05,1.259142e-05,1.412780e-05,1.585165e-05,
     1.778585e-05,1.995605e-05,2.239105e-05,2.512317e-05,
     2.818866e-05,3.162819e-05,3.548742e-05,3.981753e-05,
     4.467600e-05,5.012730e-05,5.624375e-05,6.310652e-05};
     
  double pr_arr[nlines]=
    {5.64869998e-32,5.64869986e-31,5.64869986e-30,5.64870017e-29,
     6.76729990e-28,7.82989977e-27,9.50780029e-26,3.25500004e-24,
     1.06260006e-22,5.44959997e-21,2.77850000e-20,1.35959994e-19,
     6.43729996e-19,2.94519994e-18,1.29639997e-17,5.45580009e-17,
     2.18730006e-16,2.94129994e-16,8.02570017e-16,2.14369999e-15,
     5.62639983e-15,1.45639997e-14,3.73379984e-14,4.88699996e-14,
     9.11069993e-14,1.69410005e-13,2.30929994e-13,2.81650002e-13,
     3.83670011e-13,7.11400014e-13,1.31769996e-12,2.43959995e-12,
     3.16660001e-12,4.30759985e-12,5.86130016e-12,7.96970042e-12,
     1.08389998e-11,1.39989999e-11,1.90380003e-11,2.58830006e-11,
     3.32719997e-11,4.52399992e-11,6.15209966e-11,8.36119993e-11,
     1.13700001e-10,1.45250006e-10,1.61739996e-10,1.83999996e-10,
     2.50169996e-10,
     3.205499e-10,3.459614e-10,3.645704e-10,3.776049e-10,
     3.860057e-10,3.949086e-10,4.207210e-10,4.652518e-10,5.200375e-10,
     5.845886e-10,6.601006e-10,7.587453e-10,9.152099e-10,
     1.076798e-09,1.228690e-09,1.415580e-09,1.640252e-09,
     1.908541e-09,2.225540e-09,2.602693e-09,3.047368e-09,
     3.578169e-09,4.205697e-09,4.948708e-09,5.821588e-09,
     6.891087e-09,8.186994e-09,9.481895e-09,1.116971e-08,
     1.209666e-08,1.418455e-08,1.678941e-08,1.995387e-08,
     2.365509e-08,2.801309e-08,3.321740e-08,3.951309e-08,
     4.721146e-08,5.671502e-08,6.851912e-08,8.332945e-08,
     1.021148e-07,1.261124e-07,1.570649e-07,1.969348e-07,
     2.473268e-07,3.106536e-07,3.919161e-07,4.982142e-07};
  
  double nb_arr[nlines]=
    {4.00000001e-15,4.73000011e-15,4.75999990e-15,4.91000012e-15,
     6.99000006e-15,9.89999996e-15,2.71999999e-14,1.27000000e-13,
     6.93000019e-13,6.29500011e-12,1.58099991e-11,3.97200016e-11,
     9.97599989e-11,2.50600013e-10,6.29399977e-10,1.58100000e-09,
     3.97200006e-09,4.99999997e-09,9.97599958e-09,1.98999999e-08,
     3.97199997e-08,7.92400030e-08,1.58099994e-07,1.98999999e-07,
     3.15499989e-07,4.99999999e-07,6.29400006e-07,7.92399987e-07,
     9.97599955e-07,1.58099999e-06,2.50600010e-06,3.97199983e-06,
     4.99999987e-06,6.29399983e-06,7.92399987e-06,9.97600000e-06,
     1.25600000e-05,1.58099992e-05,1.99000006e-05,2.50600006e-05,
     3.15500001e-05,3.97199983e-05,4.99999987e-05,6.29400020e-05,
     7.92400024e-05,9.97600000e-05,1.10499997e-04,1.25599996e-04,
     1.58099996e-04,
     2.113478e-04,2.371361e-04,2.660711e-04,2.985366e-04,
     3.349636e-04,3.758353e-04,4.216941e-04,6.683399e-04,7.498897e-04,
     8.413900e-04,9.440551e-04,1.059247e-03,1.188495e-03,
     1.333513e-03,1.496226e-03,1.678793e-03,1.883637e-03,
     2.113475e-03,2.371358e-03,2.660707e-03,2.985362e-03,
     3.349632e-03,3.758348e-03,4.216936e-03,4.731480e-03,
     5.308807e-03,5.956579e-03,6.683391e-03,7.498888e-03,
     8.413890e-03,9.440539e-03,1.059246e-02,1.188493e-02,
     1.333511e-02,1.496224e-02,1.678791e-02,1.883635e-02,
     2.113473e-02,2.371355e-02,2.660704e-02,2.985359e-02,
     3.349627e-02,3.758344e-02,4.216931e-02,4.731474e-02,
     5.308801e-02,5.956572e-02,6.683383e-02,7.498879e-02};

  crust_nlines=nlines;
  crust_vece.resize(crust_nlines);
  crust_vecp.resize(crust_nlines);
  crust_vecnb.resize(crust_nlines);
  for(size_t i=0;i<nlines;i++) {
    crust_vece[i]=ed_arr[i];
    crust_vecp[i]=pr_arr[i];
    crust_vecnb[i]=nb_arr[i];
  }
    
  // The energy density and pressure are already in Msun/km^3 and the
  // baryon density is in fm^{-3}
  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;

  trans_pres=crust_vecp[crust_nlines-1];
    
  if (verbose>1) {
    cout << "Transition pressure: " << endl;
    cout << trans_pres << endl;
  }

  use_crust=true;
    
  if (core_set) internal_read();
    
  return;
}

void eos_tov_interp::ngl13_low_dens_eos(double L, string model,
					bool external) {

  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    fname=dir+"newton_"+model+".o2";
  }
  
  if (L<25.0) L=25.0;
  if (L>115.0) L=115.0;
  
  // --------------------------------------------------------------
  // Load and process the data file containing the crusts

  table3d newton_eos;
  newton_eos.set_interp_type(itp_linear);

  hdf_file hf;
  string name;
  hf.open(fname);
  hdf_input(hf,newton_eos,name);
  hf.close();

  crust_nlines=newton_eos.get_ny();
  crust_vece.resize(crust_nlines);
  crust_vecp.resize(crust_nlines);
  crust_vecnb.resize(crust_nlines);
  for(size_t j=0;j<newton_eos.get_ny();j++) {
    double nbt=newton_eos.get_grid_y(j);
    crust_vece[j]=newton_eos.interp(L,nbt,"ed");
    crust_vecp[j]=newton_eos.interp(L,nbt,"pr");
    crust_vecnb[j]=nbt;
  }

  // --------------------------------------------------------------
  // Manually set the transition density by interpolating 

  std::vector<double> Lv, ntv;
  for(double Lt=25.0;Lt<116;Lt+=5.0) Lv.push_back(Lt);
  if (model=="PNM") {
    ntv.push_back(0.0898408);
    ntv.push_back(0.0862488);
    ntv.push_back(0.0831956);
    ntv.push_back(0.0805016);
    ntv.push_back(0.0781668);
    ntv.push_back(0.0760116);
    ntv.push_back(0.0743952);
    ntv.push_back(0.0727788);
    ntv.push_back(0.0713420);
    ntv.push_back(0.0700848);
    ntv.push_back(0.0688276);
    ntv.push_back(0.0673908);
    ntv.push_back(0.0666724);
    ntv.push_back(0.0663132);
    ntv.push_back(0.0654152);
    ntv.push_back(0.0641580);
    ntv.push_back(0.0645172);
    ntv.push_back(0.0641580);
    ntv.push_back(0.0636192);
  } else {
    ntv.push_back(0.113189);
    ntv.push_back(0.106646);
    ntv.push_back(0.0982820);
    ntv.push_back(0.0927144);
    ntv.push_back(0.0876856);
    ntv.push_back(0.0831956);
    ntv.push_back(0.0792444);
    ntv.push_back(0.0754728);
    ntv.push_back(0.0735992);
    ntv.push_back(0.0686480);
    ntv.push_back(0.0654152);
    ntv.push_back(0.0623620);
    ntv.push_back(0.0593088);
    ntv.push_back(0.0564352);
    ntv.push_back(0.0533820);
    ntv.push_back(0.0503288);
    ntv.push_back(0.0472756);
    ntv.push_back(0.0451204);
    ntv.push_back(0.0427856);
  }

  o2scl::interp<std::vector<double> > itp(itp_linear);
  double nt=itp.eval(L,19,Lv,ntv);
  trans_pres=itp.eval(nt,crust_nlines,crust_vecnb,crust_vecp);

  // The energy density and pressure are already in Msun/km^3 and the
  // baryon density is in fm^{-3}
  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;

  // --------------------------------------------------------------
  // Set columns and limiting values

  if (verbose>1) {
    cout << "Transition density and pressure: " << endl;
    cout << nt << " " << trans_pres << endl;
  }

  use_crust=true;
    
  if (core_set) internal_read();
    
  return;
}

/*
  S = 32 MeV, L = 105, 115 MeV
  S = 30 MeV, L = 95, 105, 115 MeV
  S = 28 MeV, L = 85, 95, 105, 115 MeV

  To be certain we aren't extrapolating too far, we need L<S*5-65.
  If S=28+e, L=75-e, then interpolate between (28,30) and (65,75)
  If S=29+e, L=80-e, then interpolate between (28,30) and (75,85)
  If S=30-e, L=85-e, then interpolate between (28,30) and (75,85)
*/
void eos_tov_interp::ngl13_low_dens_eos2(double S, double L, double nt,
					 string fname) {

  if (S<28.0 || S>38.0) {
    O2SCL_ERR("S out of range.",exc_efailed);
  }
  if (L<25.0 || L>115.0) {
    O2SCL_ERR("L out of range.",exc_efailed);
  }
  if (L>S*5.0-65.0) {
    O2SCL_ERR("S,L out of range.",exc_efailed);
  }
  if (nt<0.01 || nt>0.15) {
    O2SCL_ERR("nt out of range.",exc_efailed);
  }
  
  if (fname.length()==0) {
    std::string dir=o2scl::o2scl_settings.get_data_dir();
    fname=dir+"newton_SL.o2";
  }
 
  // --------------------------------------------------------------
  // Find limiting values of S
  
  size_t iSlow=((size_t)S);
  if (iSlow%2==1) iSlow--;
  size_t iShigh=iSlow+2;

  // Double check arithmetic
  if (iSlow<28 || iSlow>38) {
    O2SCL_ERR("iSlow out of range.",exc_efailed);
  }
  if (iShigh<28 || iShigh>38) {
    O2SCL_ERR("iSlow out of range.",exc_efailed);
  }
  
  // --------------------------------------------------------------
  // Weights for interpolation

  double weight_low=(2.0-(S-((double)iSlow)))/2.0;
  double weight_high=1.0-weight_low;
  if (weight_low<0.0 || weight_high<0.0) {
    O2SCL_ERR("Sanity check (1) failed in ngl13_low_dens_eos2().",
	      exc_esanity);
  }

  // --------------------------------------------------------------
  // Load and process the data file containing the crusts

  table3d nlow, nhigh;

  hdf_file hf;
  string name;
  hf.open(fname);
  hdf_input(hf,nlow,((string)"S")+itos(iSlow));
  hdf_input(hf,nhigh,((string)"S")+itos(iShigh));
  hf.close();

  nlow.set_interp_type(itp_linear);
  nhigh.set_interp_type(itp_linear);

  crust_nlines=nlow.get_ny();
  crust_vece.resize(crust_nlines);
  crust_vecp.resize(crust_nlines);
  crust_vecnb.resize(crust_nlines);

  for(size_t j=0;j<nlow.get_ny();j++) {

    double nbt=nlow.get_grid_y(j);
    double edval=nlow.interp(L,nbt,"ed")*weight_low+
      nhigh.interp(L,nbt,"ed")*weight_high;
    double prval=nlow.interp(L,nbt,"pr")*weight_low+
      nhigh.interp(L,nbt,"pr")*weight_high;
    
    if (edval<1.0e-100 || prval<1.0e-100 || nbt<1.0e-100) {
      O2SCL_ERR("Sanity check (2) failed in ngl13_low_dens_eos2().",
		exc_esanity);
    }
    
    crust_vece[j]=edval;
    crust_vecp[j]=prval;
    crust_vecnb[j]=nbt;
  }

  // --------------------------------------------------------------
  // Transition pressure

  o2scl::interp<std::vector<double> > itp(itp_linear);
  trans_pres=itp.eval(nt,crust_nlines,crust_vecnb,crust_vecp);
  
  // --------------------------------------------------------------
  // Unit system

  // The energy density and pressure are already in Msun/km^3 and the
  // baryon density is in fm^{-3}

  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;

  // --------------------------------------------------------------
  // Set columns and limiting values

  trans_pres=crust_vecp[crust_nlines-1];
    
  if (verbose>1) {
    cout << "Transition density and pressure: " << endl;
    cout << nt << " " << trans_pres << endl;
  }

  use_crust=true;
    
  if (core_set) internal_read();
    
  return;
}

void eos_tov_interp::s12_low_dens_eos(string model, bool external) {

  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    fname=dir+"/"+model+"_cs01_feq.txt";
  }
  
  // --------------------------------------------------------------
  // Load and process the data file containing the crusts

  ifstream fin;
  fin.open(fname.c_str());
  string stemp;
  fin >> stemp >> stemp >> stemp;
  double dtemp;

  double factor=o2scl_settings.get_convert_units().convert
    ("1/fm^4","Msun/km^3",1.0);
  
  crust_vece.clear();
  crust_vecp.clear();
  crust_vecnb.clear();

  while (fin >> dtemp) {
    double line[3];
    line[2]=dtemp;
    fin >> line[0] >> line[1];
    // Original text file is in 1/fm^4, so we convert to Msun/km^3
    line[0]*=factor;
    line[1]*=factor;
    crust_vece.push_back(line[0]);
    crust_vecp.push_back(line[1]);
    crust_vecnb.push_back(line[2]);
  }
  fin.close();

  // --------------------------------------------------------------
  // Manually set the transition density by interpolating 

  o2scl::interp<std::vector<double> > itp(itp_linear);
  trans_pres=itp.eval(0.08,crust_nlines,crust_vecnb,crust_vecp);

  // --------------------------------------------------------------
  // Set columns and limiting values

  if (verbose>1) {
    cout << "Transition pressure: " << endl;
    cout << trans_pres << endl;
  }

  use_crust=true;
    
  if (core_set) internal_read();
    
  return;
}

void eos_tov_interp::gcp10_low_dens_eos(string model, bool external) {

  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  cout << "dir: " << dir << endl;
  if (external) {
    fname=model;
  } else {
    if (model==((string)"BSk19")) {
      fname=dir+"/eos19.o2";
    } else if (model==((string)"BSk21")) {
      fname=dir+"/eos21.o2";
    } else {
      fname=dir+"/eos20.o2";
    }
  }
  
  // --------------------------------------------------------------
  // Load and process the data file

  table_units<> tab;
  hdf_file hf;
  hf.open(fname);
  std::string name;
  hdf_input(hf,tab,name);
  hf.close();

  /*
    We double check that nlines and maxlines are equal. This ensures
    that crust_vece is the same size as crust_nlines. This wouldn't be
    strictly necessary, but it might lead to some surprising behavior,
    so we enforce equality here.
  */
  if (tab.get_nlines()!=tab.get_maxlines()) {
    O2SCL_ERR2("Misalignment sanity check for GCP10 crust in ",
	       "eos_tov_interp::gcp10_low_dens_eos().",exc_esanity);
  }

  tab.convert_to_unit("rho","Msun/km^3");
  // I'm not sure why the units in the data file are erg/cm^2. Is that a typo?
  tab.set_unit("P","erg/cm^3");
  tab.convert_to_unit("P","Msun/km^3");

  crust_vece.resize(tab.get_maxlines());
  crust_vecp.resize(tab.get_maxlines());
  crust_vecnb.resize(tab.get_maxlines());

  tab.swap_column_data("rho",crust_vece);
  tab.swap_column_data("P",crust_vecp);
  tab.swap_column_data("nb",crust_vecnb);
  crust_nlines=tab.get_nlines();

  // --------------------------------------------------------------
  // Transition pressures from the table in Pearson12

  if (model==((string)"BSk19")) {
    trans_pres=o2scl_settings.get_convert_units().convert
      ("MeV/fm^3","Msun/km^3",0.428);
  } else if (model==((string)"BSk21")) {
    trans_pres=o2scl_settings.get_convert_units().convert
      ("MeV/fm^3","Msun/km^3",0.365);
  } else {
    trans_pres=o2scl_settings.get_convert_units().convert
      ("MeV/fm^3","Msun/km^3",0.268);
  }

  // --------------------------------------------------------------
  // Set columns and limiting values

  if (verbose>1) {
    cout << "Transition pressure: " << endl;
    cout << trans_pres << endl;
  }

  use_crust=true;
    
  if (core_set) internal_read();
    
  return;
}

void eos_tov_interp::ed_nb_from_pr(double pr, double &ed, double &nb) {

  if (!std::isfinite(pr)) {
    O2SCL_ERR("Pressure not finite in eos_tov_interp::ed_nb_from_pr().",
	      exc_efailed);
  }

  ed=pe_int.eval(pr);
  if (baryon_column) {
    nb=pnb_int.eval(pr);
  }

  if (!std::isfinite(ed) || (baryon_column && !std::isfinite(nb))) {
    string s="Energy density or baryon density not finite at pressure ";
    s+=dtos(pr)+" in eos_tov_interp::ed_nb_from_pr().";
    O2SCL_ERR(s.c_str(),exc_efailed);
  }
  
  return;
}

double eos_tov_interp::ed_from_pr(double pr) {
  return pe_int.eval(pr);
}

double eos_tov_interp::ed_from_nb(double nb) {
  return gen_int.eval(nb,full_nlines,full_vecnb,full_vece);
}

double eos_tov_interp::nb_from_pr(double pr) {
  return pnb_int.eval(pr);
}

double eos_tov_interp::nb_from_ed(double ed) {
  return gen_int.eval(ed,full_nlines,full_vece,full_vecnb);
}

double eos_tov_interp::pr_from_nb(double nb) {
  return gen_int.eval(nb,full_nlines,full_vecnb,full_vecp);
}

double eos_tov_interp::pr_from_ed(double ed) {
  return gen_int.eval(ed,full_nlines,full_vece,full_vecp);
}

void eos_tov_interp::get_eden_user(double pres, double &ed, double &nb) {
  ed=pe_int.eval(pres*pfactor)/efactor;
  if (baryon_column) {
    nb=pnb_int.eval(pres*pfactor)/nfactor;
  }
  return;
}

void eos_tov_interp::get_aux(double P, size_t &np, std::vector<double> &auxp) {
  np=0;
  if (core_auxp>0) {
    double first_pres=(*core_table)[colp][0]/pfactor;
    for(size_t i=0;i<((int)core_table->get_ncolumns());i++) {
      if (i!=cole && i!=colp && i!=colnb) {
	if (use_crust==true && P>=first_pres) {
	  auxp[np]=0.0;
	} else {
	  // Convert P to user-units by dividing by pfactor
	  auxp[np]=gen_int.eval(P/pfactor,core_table->get_nlines(),
				(*core_table)[colp],(*core_table)[i]);
	} 
	np++;
      }
    }
  }
  return;
}

