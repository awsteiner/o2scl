/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#include <o2scl/nstar_cold.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nstar_cold::nstar_cold() : eost(new table_units<>) {

  neut.init(o2scl_settings.get_convert_units().convert
            ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  prot.init(o2scl_settings.get_convert_units().convert
            ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  neut.non_interacting=false;
  prot.non_interacting=false;

  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  mu.init(o2scl_settings.get_convert_units().convert
	  ("kg","1/fm",o2scl_mks::mass_muon),2.0);
  
  def_tov.verbose=1;
  def_eos_tov.verbose=1;

  def_tov.set_units("1/fm^4","1/fm^4","1/fm^3");
  def_eos_tov.default_low_dens_eos();
  tp=&def_tov;

  // Initialize diagnostic quantities
  pressure_dec_nb=0.0;
  acausal_nb=0.0;
  acausal_pr=0.0;
  acausal_ed=0.0;
  eos_neg=false;
  deny_urca_nb=0.0;
  allow_urca_nb=0.0;
  max_row=0;

  nb_start=0.05;
  nb_end=2.0;
  dnb=0.01;

  include_muons=true;
  eos_set=false;

  verbose=0;

  err_nonconv=true;
  remove_rows=true;
  
  mh.err_nonconv=false;

  nb_last_min=0.48;
  
  // In case calc_nstar() is used without calc_eos()
  nb_last=2.0;
}

int nstar_cold::solve_fun(size_t nv, const ubvector &x, ubvector &y,
                          thermo &hb, double n_B) {
  
  prot.n=x[0];
  neut.n=n_B-prot.n;
  
  if (neut.n<0.0 || prot.n<0.0) return 1;
  int had_ret=hep->calc_e(neut,prot,hb);
  if (had_ret!=0) {
    // Most EOS failures result in exceptions, but in case there is
    // a failure which does not, then we return a non-zero value here
    return 2;
  }

  e.mu=neut.mu-prot.mu;
  fzt.calc_mu_zerot(e);
  y[0]=prot.n-e.n;
  
  if (include_muons) {
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);
    y[0]-=mu.n;
  }
  
  return 0;
}

int nstar_cold::calc_eos(double np_0) {

  if (verbose>0) {
    cout << "Starting calc_eos()." << endl;
  }

  double fac=(nb_end-nb_start)/dnb;
  if (fac<0.0 || fac>1.0e8) {
    std::cout << nb_end << " " << nb_start << " " << dnb << " " << fac << endl;
    O2SCL_ERR2("Invalid baryon density range in ",
               "nstar_cold::calc_eos().",o2scl::exc_einval);
  }
  
  if (eos_set==false) {
    O2SCL_ERR("EOS not set in calc_eos().",exc_efailed);
  }
  
  if (fabs(neut.g-2.0)>1.0e-10 || fabs(prot.g-2.0)>1.0e-10) {
    O2SCL_ERR((((std::string)"Neutron (")+std::to_string(neut.g)+
               ") or proton ("+std::to_string(prot.g)+") spin deg"+
               "eneracies wrong in "+
               "nstar_cold::calc_eos().").c_str(),
              exc_einval);
  }
  if (fabs(neut.m-4.5)>1.0 || fabs(prot.m-4.5)>1.0) {
    O2SCL_ERR((((std::string)"Neutron (")+std::to_string(neut.m)+
               ") or proton ("+std::to_string(prot.m)+") masses wrong "+
               "in nstar_cold::calc_eos().").c_str(),
              exc_einval);
  }
  
  eost->clear();
  eost->line_of_names(((string)"ed pr nb mun mup mue nn np ne kfn ")+
		      "kfp kfe");
  //eost->line_of_names(((string)"ed pr nb mun mup mue nn np ne kfn ")+
  //"kfp kfe fcs2 dednb_Ye dPdnb_Ye");
  eost->set_unit("ed","1/fm^4");
  eost->set_unit("pr","1/fm^4");
  eost->set_unit("nb","1/fm^3");
  eost->set_unit("mun","1/fm");
  eost->set_unit("mup","1/fm");
  eost->set_unit("mue","1/fm");
  eost->set_unit("nn","1/fm^3");
  eost->set_unit("np","1/fm^3");
  eost->set_unit("ne","1/fm^3");
  eost->set_unit("kfn","1/fm");
  eost->set_unit("kfp","1/fm");
  eost->set_unit("kfe","1/fm");
  //eost->set_unit("dednb_Ye","1/fm");
  //eost->set_unit("dPdnb_Ye","1/fm");
  if (include_muons) {
    eost->new_column("mumu");
    eost->new_column("nmu");
    eost->new_column("kfmu");
    eost->set_unit("mumu","1/fm");
    eost->set_unit("nmu","1/fm^3");
    eost->set_unit("kfmu","1/fm");
  }

  // Fix poor initial guesses for np_0
  if (fabs(np_0)<1.0e-12) {
    np_0=nb_start/3.0;
  }

  // Initialize diagnostic quantities
  pressure_dec_nb=0.0;
  acausal_nb=0.0;
  acausal_pr=0.0;
  acausal_ed=0.0;
  eos_neg=false;
  deny_urca_nb=0.0;
  allow_urca_nb=0.0;
  max_row=0;
  
  /// Thermodynamic quantities
  thermo h, hb;

  // Store the previous pressure so we can track if it decreases
  double oldpr=0.0;
  
  // The function object for the solver
  mm_funct sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                        ubvector &, thermo &,double)>
                        (&nstar_cold::solve_fun),
                        this,std::placeholders::_1,std::placeholders::_2,
                        std::placeholders::_3,std::ref(hb),nb_start);
  
  if (verbose>0) {
    cout << "n_B         n_n         n_p         n_e         "
         << "eps         P" << endl;
    cout << "[1/fm^3]    [1/fm^3]    [1/fm^3]    [1/fm^3]    "
	 << "[1/fm^4]    [1/fm^4]" << endl;
  }

  // Get a better initial guess for nb_start. This section of code was
  // put in to avoid negative densities for calc_eos from the solver
  // below. 
  double n_B=nb_start;
  
  ubvector ux(1), uy(1);
  ux[0]=np_0;
  int tret=mh.msolve(1,ux,sf);
  if (tret!=0) {
    O2SCL_CONV_RET("Initial solver failed.",
                   o2scl::exc_efailed,err_nonconv);
  }

  // Loop over the density range, and also determine pressure_dec_nb
  bool done=false;
  for(n_B=nb_start;(dnb>0.0 && n_B<=nb_end+dnb/10.0 && done==false) ||
        (dnb<0.0 && n_B>=nb_end-dnb/10.0);n_B+=dnb) {

    sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                 ubvector &, thermo &,double)>
                 (&nstar_cold::solve_fun),
                 this,std::placeholders::_1,std::placeholders::_2,
                 std::placeholders::_3,std::ref(hb),n_B);
    
    tret=mh.msolve(1,ux,sf);
    nb_last=n_B;
    if (tret!=0 && n_B<nb_last_min) {
      std::string str_err=((string)"Solver failed inside loop at ")+
        "density "+o2scl::dtos(n_B)+" with nb_last_min="+
        o2scl::dtos(nb_last_min)+" and nb_end="+o2scl::dtos(nb_end);
      O2SCL_CONV_RET(str_err.c_str(),o2scl::exc_efailed,err_nonconv);
    }
    
    // Compute the function at the final point
    sf(1,ux,uy);
      
    if (false) {
      // AWS: 3/9/21: this is nothing other than the speed of
      // sound, which can be easily approximated by dP/deps,
      // saving the trouble of computing these derivatives.
      // I'm commenting this out for now, but I may put it
      // back in later. Even if it should be put in later,
      // it's probably more efficient just to compute dP/depsilon
      // directly.
      
      // ------------------------------------------------------------
      // Compute dP/de at fixed Ye and Ymu. This code uses neut.n and
      // prot.n, so we'll have to recompute them below.
      
      // Compute the hadronic part
      double dednb_Yp, dPdnb_Yp;
      hep->const_pf_derivs(n_B,prot.n/n_B,dednb_Yp,dPdnb_Yp);
      // Compute the leptonic part
      double dne_dmue=sqrt(e.mu*e.mu-e.m*e.m)*e.mu/pi2;
      double dP_dne=e.n/dne_dmue;
      // Put them together
      double numer=dPdnb_Yp+dP_dne*e.n/n_B;
      double denom=dednb_Yp+e.mu*e.n/n_B;
      // Add the muon contribution
      if (include_muons && mu.n>0.0) {
        double dnmu_dmumu=sqrt(mu.mu*mu.mu-mu.m*mu.m)*mu.mu/pi2;
        double dP_dnmu=mu.n/dnmu_dmumu;
        numer+=dP_dnmu*mu.n/n_B;
        denom+=mu.mu*mu.n/n_B;
      }
      double fcs2=numer/denom;

      // ------------------------------------------------------------
      // Recompute neut.n and prot.n
      //y=solve_fun(x,hb);
      
    }

    if (neut.inc_rest_mass==false) {
      hb.ed+=neut.n*neut.m;
      neut.mu+=neut.m;
    }
    if (prot.inc_rest_mass==false) {
      hb.ed+=prot.n*prot.m;
      prot.mu+=prot.m;
    }

    if (include_muons) {

      h=hb+e+mu;
      
      //double line[18]={h.ed,h.pr,n_B,neut.mu,prot.mu,e.mu,neut.n,prot.n,e.n,
      //neut.kf,prot.kf,e.kf,fcs2,denom,numer,mu.mu,mu.n,mu.kf};
      //eost->line_of_data(18,line);

      vector<double> line={h.ed,h.pr,n_B,neut.mu,prot.mu,e.mu,neut.n,
        prot.n,e.n,neut.kf,prot.kf,e.kf,mu.mu,mu.n,mu.kf};
      eost->line_of_data(line.size(),line);

    } else {

      h=hb+e;

      vector<double> line={h.ed,h.pr,n_B,neut.mu,prot.mu,e.mu,neut.n,
        prot.n,e.n,neut.kf,prot.kf,e.kf};
      eost->line_of_data(line.size(),line);

    }
    
    if (verbose>0) {
      cout.precision(5);
      cout << n_B << " " << neut.n << " " << prot.n << " " << e.n << " " 
	   << h.ed << " " << h.pr << endl;
      cout.precision(6);
    }
    
    if (n_B>nb_start && pressure_dec_nb<=0.0 && h.pr<oldpr) {
      pressure_dec_nb=n_B;
    }
    oldpr=h.pr;

    // Proceed to next baryon density
  }

  // -----------------------------------------------------------------
  // Report if the pressure decreases

  if (pressure_dec_nb>nb_start) {
    if (verbose>0) {
      cout << "Pressure is flat near a baryon density of " << pressure_dec_nb
	   << " 1/fm^3." << endl;
    }
  }

  // -----------------------------------------------------------------
  // Compute the speed of sound
  
  eost->deriv("ed","pr","cs2");
  if (verbose>1) {
    cout << " ix ed           pr           cs2" << endl;
  }
  for(size_t i=0;i<eost->get_nlines()-1;i++) {
    double cs2_temp=eost->get("cs2",i);
    if (verbose>1) {
      cout.width(3);
      cout << i << " " << eost->get("ed",i) << " " << eost->get("pr",i) << " "
           << eost->get("cs2",i) << endl;
    }
    if (!std::isfinite(cs2_temp)) {
      O2SCL_CONV2_RET("Speed of sound not finite in ",
                      "nstar_cold::calc_eos().",
                      exc_efailed,err_nonconv);
    }
    if (cs2_temp>=1.0) {
      acausal_nb=eost->interp("cs2",1.0,"nb");
      acausal_pr=eost->interp("cs2",1.0,"pr");
      acausal_ed=eost->interp("cs2",1.0,"ed");
      i=eost->get_nlines();
    }
  }
  
  if (verbose>0) {
    if (acausal_nb>0.0) {
      cout << "The EOS is acausal at nb=" << acausal_nb << ", ed="
           << acausal_ed << ", pr="
           << acausal_pr << " ." << endl;
    } else {
      cout << "The EOS is causal." << endl;
    }
  }
    
  // -----------------------------------------------------------------
  // Calculate the adiabatic index and the Urca threshold

  eost->new_column("logp");
  eost->new_column("loge");
  eost->new_column("s");
  eost->new_column("urca");
  eost->set_unit("s","1/fm");
  eost->set_unit("urca","1/fm^4");

  if (verbose>0) {
    cout << "Going to Urca, pass 1: " << endl;
  }
  
  for(size_t i=0;i<eost->get_nlines();i++) {
    // Compute Urca threshold
    double stmp=(eost->get("kfn",i)+eost->get("kfp",i)+eost->get("kfe",i))/2.0;
    eost->set("s",i,stmp);
    double utmp=(stmp*(stmp-eost->get("kfn",i))*(stmp-eost->get("kfp",i))*
		 (stmp-eost->get("kfe",i)));
    eost->set("urca",i,utmp);
  }

  if (verbose>0) {
    cout << "Computing adiabatic index: " << endl;
  }
  
  for(size_t i=0;i<eost->get_nlines();i++) {
    if (eost->get("pr",i)<0.0) {
      eos_neg=true;
      eost->set("logp",i,0.0);
    } else {
      eost->set("logp",i,log(eost->get("pr",i)));
    }
    if (eost->get("ed",i)<0.0) {
      eos_neg=true;
      eost->set("loge",i,0.0);
    } else {
      eost->set("loge",i,log(eost->get("ed",i)));
    }
  }
  if (eos_neg==false) {
    eost->deriv("loge","logp","ad_index");
  } else {
    eost->new_column("ad_index");
    for(size_t i=0;i<eost->get_nlines();i++) {
      eost->set("ad_index",i,0.0);
    }
  }

  // -----------------------------------------------------------------
  // Determine the densities at which the Urca process turns on and
  // off
  
  if (verbose>0) {
    cout << "Going to Urca, pass 2: " << endl;
  }

  for(size_t i=0;i<eost->get_nlines()-1;i++) {
    if (deny_urca_nb<=0.0 && eost->get("urca",i)>0.0 && 
	eost->get("urca",i+1)<0.0) {
      deny_urca_nb=eost->get("nb",i);
    } else if (allow_urca_nb<=0.0 && eost->get("urca",i)>0.0) {
      /// Use linear interpolation to improve the urca estimate
      allow_urca_nb=eost->get("nb",i)+
	(-eost->get("urca",i))*(eost->get("nb",i+1)-eost->get("nb",i))/
	(eost->get("urca",i+1)-eost->get("urca",i));
      
    }
  }

  if (verbose>0) {
    if (allow_urca_nb>0.0) {
      if (deny_urca_nb>0.0) {
        cout << "Urca allowed at " << allow_urca_nb
             << " fm^{-3} and disallowed at " << deny_urca_nb
             << " fm^{-3}." << endl;
      } else {
        cout << "Urca allowed at " << allow_urca_nb
             << " fm^{-3}." << endl;
      }
    } else {
      cout << "Urca process is never allowed." << endl;
    }
    cout << "Done with calc_eos()." << endl;
  }
  
  return 0;
}

double nstar_cold::calc_urca(double np_0) {
  double old_urca=0.0, urca;
  
  // Fix poor initial guesses for np_0
  if (fabs(np_0)<1.0e-12) {
    np_0=nb_start/3.0;
  }

  thermo hb;

  bool success=true;

  ubvector ux(1), uy(1);
  ux[0]=np_0;
  
  for(double n_B=nb_start;n_B<=nb_end+dnb/10.0;n_B+=dnb) {
    
    mm_funct sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                          ubvector &,thermo &,double)>
                          (&nstar_cold::solve_fun),
                          this,std::placeholders::_1,
                          std::placeholders::_2,
                          std::placeholders::_3,
                          std::ref(hb),n_B);
    
    int ret=mh.msolve(1,ux,sf);
    sf(1,ux,uy);
    
    if (ret!=0) {
      success=false;
    }
    
    double s=(neut.kf+prot.kf+e.kf)/2.0;
    urca=s*(s-neut.kf)*(s-prot.kf)*(s-e.kf);
    if (n_B>nb_start && urca>0.0 && old_urca<0.0) {
      if (success==false) {
	O2SCL_ERR("Solution failed in calc_urca().",exc_efailed);
      }
      return n_B-dnb*urca/(urca-old_urca);
    }
    
    old_urca=urca;
  }

  if (success==false) {
    O2SCL_ERR("Solution failed in calc_urca().",exc_efailed);
  }
  
  return 0.0;
}

int nstar_cold::calc_nstar() {

  // Check to make sure required columns are present
  if (!eost->is_column("ed")) {
    O2SCL_ERR2("Column ed not found in table in ",
	       "eos_tov_interp::read_table().",o2scl::exc_einval);
  }
  if (!eost->is_column("pr")) {
    O2SCL_ERR2("Column pr not found in table in ",
	       "eos_tov_interp::read_table().",o2scl::exc_einval);
  }
  if (!eost->is_column("nb")) {
    O2SCL_ERR2("Column nb not found in table in ",
	       "eos_tov_interp::read_table().",o2scl::exc_einval);
  }

  // Read the table and set the EOS
  def_eos_tov.read_table(*eost,"ed","pr","nb");
  tp->set_units("1/fm^4","1/fm^4","1/fm^3");
  tp->set_eos(def_eos_tov);

  // Compute the M-R curve
  int mret=tp->mvsr();
  if (mret!=0) {
    O2SCL_CONV2_RET("TOV solver failed in ",
                    "nstar_cold::calc_nstar().",mret,err_nonconv);
  }

  // Find row containing maximum mass
  std::shared_ptr<table_units<> > mvsrt=tp->get_results();
  max_row=mvsrt->lookup("gm",mvsrt->max("gm"));

  // Compute central pressure and baryon density of maximum mass star
  double pr_max=mvsrt->get("pr",max_row);
  double nb_max=mvsrt->get("nb",max_row);

  if (nb_max>nb_last) {
    O2SCL_CONV2_RET("Central baryon density larger than last baryon ",
                    "density in nstar_cold::calc_nstar().",
                    o2scl::exc_einval,err_nonconv);
  }

  // Remove rows beyond maximum mass
  if (remove_rows) {
    mvsrt->set_nlines(max_row+1);
  }

  // Report problems
  if (acausal_nb>0.0 && nb_max>acausal_nb) {
    if (verbose>0) {
      cout << "Acausal (acausal_nb,nb_max): "
           << acausal_nb << " " << nb_max << endl;
    }
    O2SCL_CONV2_RET("EOS acausal at densities below central density ",
                    "of maximum mass star in nstar_cold::calc_nstar().",
                    o2scl::exc_einval,err_nonconv);
  }
  if (pressure_dec_nb>0.0 && nb_max>pressure_dec_nb) {
    if (verbose>0) {
      cout << "Pressure decreasing (pressure_dec_nb,nb_max): "
           << pressure_dec_nb << " " << nb_max << endl;
    }
    O2SCL_CONV2_RET("Pressure descreasing in ",
                    "maximum mass star in nstar_cold::calc_nstar().",
                    o2scl::exc_einval,err_nonconv);
  }
  
  return 0;
}

int nstar_cold::fixed(double target_mass) {
  def_eos_tov.read_table(*eost,"ed","pr","nb");
  
  tp->set_units("1/fm^4","1/fm^4","1/fm^3");
  tp->set_eos(def_eos_tov);
  
  int fret=tp->fixed(target_mass);
  std::shared_ptr<table_units<> > ft=tp->get_results();
  
  double nb_max=ft->get("nb",0);

  if (nb_max>nb_last) {
    O2SCL_CONV2_RET("Central baryon density larger than last baryon ",
                    "density in nstar_cold::fixed().",
                    o2scl::exc_einval,err_nonconv);
  }
  
  return fret;
}

int nstar_hot::solve_fun_T(size_t nv, const ubvector &x,
                           ubvector &y, thermo &hb, double T,
                           double n_B) {
  
  prot.n=x[0];
  neut.n=n_B-prot.n;

  if (neut.n<=0.0 || prot.n<=0.0) return 1;
  int had_ret=hepT->calc_temp_e(neut,prot,T,hb);
  if (had_ret!=0) {
    // Most EOS failures result in exceptions, but in case there is
    // a failure which does not, then we return a non-zero value here
    return 2;
  }
  
  e.mu=neut.mu-prot.mu;
  ft.calc_mu(e,T);
  y[0]=prot.n-e.n;
  
  if (include_muons) {
    mu.mu=e.mu;
    ft.calc_mu(mu,T);
    y[0]-=mu.n;
  }
  
  return 0;
}

int nstar_hot::solve_fun_s(size_t nv, const ubvector &x, ubvector &y,
			   thermo &hb, double s, double n_B) {
  neut.n=x[0];
  double T=x[1];
  
  prot.n=n_B-neut.n;
  int had_ret=hepT->calc_temp_e(neut,prot,T,hb);
  if (had_ret!=0) return had_ret;

  e.mu=neut.mu-prot.mu;
  ft.calc_mu(e,T);
  y[0]=prot.n-e.n;
  
  if (include_muons) {
    mu.mu=e.mu;
    ft.calc_mu(mu,T);
    y[0]-=mu.n;
  }

  y[1]=hb.en/n_B-s;
  
  return 0;
}

int nstar_hot::solve_fun_s_YLe(size_t nv, const ubvector &x, ubvector &y,
			       thermo &hb, double s, double YLe,
                               double n_B) {
  neut.n=x[0];
  double T=x[1];
  
  prot.n=n_B-neut.n;
  int had_ret=hepT->calc_temp_e(neut,prot,T,hb);
  if (had_ret!=0) return had_ret;

  e.mu=neut.mu-prot.mu;
  ft.calc_mu(e,T);
  y[0]=prot.n-e.n;
  
  if (include_muons) {
    mu.mu=e.mu;
    ft.calc_mu(mu,T);
    y[0]-=mu.n;
  }

  y[1]=hb.en/n_B-s;
  
  return 0;
}

int nstar_hot::calc_eos_point_beta_T(double nb, double T, double &np,
                                     thermo &hb) {

  if (eos_T_set==false) {
    O2SCL_ERR("EOS not set in calc_eos_T_point().",exc_efailed);
  }

  if (np<=0.0) np=nb/2.0;

  mm_funct sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                        ubvector &, thermo &, double,double)>
                        (&nstar_hot::solve_fun_T),
                        this,std::placeholders::_1,std::placeholders::_2,
                        std::placeholders::_3,std::ref(hb),T,nb);
  
  ubvector ux(1), uy(1);
  ux[0]=np;
  int tret=mh.msolve(1,ux,sf);
  np=ux[0];
  
  return tret;
}
  
int nstar_hot::calc_eos_T(double T, double np_0) {

  if (verbose>0) {
    cout << "Starting calc_eos_T()." << endl;
  }

  if (eos_T_set==false) {
    O2SCL_ERR("EOS not set in calc_eos_T().",exc_efailed);
  }

  double fac=(nb_end-nb_start)/dnb;
  if (fac<0.0 || fac>1.0e8) {
    O2SCL_ERR2("Invalid baryon density range in ",
               "nstar_hot::calc_eos_T().",o2scl::exc_einval);
  }
  
  eost->clear();
  eost->line_of_names(((string)"ed pr nb mun mup mue nn np ne kfn ")+
		      "kfp kfe");
  eost->set_unit("ed","1/fm^4");
  eost->set_unit("pr","1/fm^4");
  eost->set_unit("nb","1/fm^3");
  eost->set_unit("mun","1/fm");
  eost->set_unit("mup","1/fm");
  eost->set_unit("mue","1/fm");
  eost->set_unit("nn","1/fm^3");
  eost->set_unit("np","1/fm^3");
  eost->set_unit("ne","1/fm^3");
  eost->set_unit("kfn","1/fm");
  eost->set_unit("kfp","1/fm");
  eost->set_unit("kfe","1/fm");

  // Initial guess for the proton density
  double x;
  if (fabs(np_0)<1.0e-12) x=nb_start/3.0;
  else x=np_0;
  
  double oldpr=0.0;
  pressure_dec_nb=0.0;
  
  /// Thermodynamic quantities
  thermo h, hb;
  
  mm_funct sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                        ubvector &, thermo &, double,double)>
                        (&nstar_hot::solve_fun_T),
                        this,std::placeholders::_1,std::placeholders::_2,
                        std::placeholders::_3,std::ref(hb),T,nb_start);
  
  if (verbose>0) {
    cout << "n_B         n_n         n_p         n_e         "
         << "eps         P" << endl;
    cout << "[1/fm^3]    [1/fm^3]    [1/fm^3]    [1/fm^3]    "
	 << "[1/fm^4]    [1/fm^4]" << endl;
  }

  ubvector ux(1), uy(1);
  ux[0]=x;
  int tret=mh.msolve(1,ux,sf);
  if (tret!=0) {
    O2SCL_CONV2_RET("Solution at first baryon density failed ",
                    "in nstar_hot::calc_eos_T().",o2scl::exc_efailed,
                    err_nonconv);
  }
  x=ux[0];
    
  for(double n_B=nb_start;(dnb>0.0 && n_B<=nb_end+dnb/10.0) ||
        (dnb<0.0 && n_B>=nb_end-dnb/10.0);n_B+=dnb) {
    
    sf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                 ubvector &, thermo &,double,double)>
                 (&nstar_hot::solve_fun_T),
                 this,std::placeholders::_1,std::placeholders::_2,
                 std::placeholders::_3,std::ref(hb),T,n_B);

    double y;

    // In case the baryon density stepping led to a poor guess
    if (ux[0]>n_B) ux[0]=n_B*0.5;
    
    tret=mh.msolve(1,ux,sf);
    if (tret!=0) {
      std::cout << "Solver failed with error code " << tret << endl;
      cout << "x: " << ux[0] << " ret: " << sf(1,ux,uy)
           << " y: " << uy[0] << endl;
      mh.verbose=2;
      mh.msolve(1,ux,sf);
      O2SCL_CONV2_RET("Solver failed ",
                      "in nstar_hot::calc_eos_T().",o2scl::exc_efailed,
                      err_nonconv);
    }
    x=ux[0];
    sf(1,ux,uy);
    y=uy[0];
      
    // ------------------------------------------------------------
    
    if (neut.inc_rest_mass==false) {
      hb.ed+=neut.n*neut.m;
      neut.mu+=neut.m;
    }
    if (prot.inc_rest_mass==false) {
      hb.ed+=prot.n*prot.m;
      prot.mu+=prot.m;
    }

    if (include_muons) {

      h=hb+e+mu;
      
      vector<double> line={h.ed,h.pr,n_B,neut.mu,prot.mu,e.mu,neut.n,
                       prot.n,e.n,neut.kf,prot.kf,e.kf};
      eost->line_of_data(line.size(),line);

    } else {

      h=hb+e;

      vector<double> line={h.ed,h.pr,n_B,neut.mu,prot.mu,e.mu,neut.n,
                       prot.n,e.n,neut.kf,prot.kf,e.kf};
      eost->line_of_data(line.size(),line);
    }
    
    if (verbose>0) {
      cout.precision(5);
      cout << n_B << " " << neut.n << " " << prot.n << " " << e.n << " " 
	   << h.ed << " " << h.pr << endl;
      cout.precision(6);
    }
    
    if (n_B>nb_start && pressure_dec_nb<=0.0 && h.pr<oldpr) {
      pressure_dec_nb=n_B;
    }
    oldpr=h.pr;

    // Proceed to next baryon density
  }

  if (verbose>0) {
    cout << "Done with nstar_hot::calc_eos_T()." << endl;
  }
  
  return 0;
}

