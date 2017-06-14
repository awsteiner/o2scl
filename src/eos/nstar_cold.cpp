/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

  np.init(o2scl_settings.get_convert_units().convert
	     ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  pp.init(o2scl_settings.get_convert_units().convert
	     ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  np.non_interacting=false;
  pp.non_interacting=false;

  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  mu.init(o2scl_settings.get_convert_units().convert
	  ("kg","1/fm",o2scl_mks::mass_muon),2.0);
  
  def_tov.verbose=1;
  def_eos_tov.verbose=1;

  def_tov.set_units("1/fm^4","1/fm^4","1/fm^3");
  def_eos_tov.default_low_dens_eos();
  tp=&def_tov;

  acausal=0.0;
  pressure_dec=0.0;
  allow_urca=0.0;
  deny_urca=0.0;

  nb_start=0.05;
  nb_end=2.0;
  dnb=0.01;

  include_muons=false;
  eos_set=false;

  rp=&def_root;

  solver_tol=1.0e-4;

  verbose=0;

  err_nonconv=true;
  well_formed=true;
}

double nstar_cold::solve_fun(double x) {
  double y;
  
  np.n=x;
  
  pp.n=barn-np.n;
  int had_ret=hep->calc_e(np,pp,hb);
  if (had_ret!=0) return had_ret;

  e.mu=np.mu-pp.mu;
  fzt.calc_mu_zerot(e);
  y=pp.n-e.n;
  
  if (include_muons) {
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);
    y-=mu.n;
  }
  
  return y;
}

int nstar_cold::calc_eos(double np_0) {

  if (verbose>0) {
    cout << "Starting calc_eos()." << endl;
  }

  if (eos_set==false) {
    O2SCL_ERR("EOS not set in calc_eos().",exc_efailed);
  }
  
  eost->clear_table();
  eost->line_of_names(((string)"ed pr nb mun mup mue nn np ne kfn ")+
		      "kfp kfe fcs2 dednb_Ye dPdnb_Ye");
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
  eost->set_unit("dednb_Ye","1/fm");
  eost->set_unit("dPdnb_Ye","1/fm");
  if (include_muons) {
    eost->new_column("mumu");
    eost->new_column("nmu");
    eost->new_column("kfmu");
    eost->set_unit("mumu","1/fm");
    eost->set_unit("nmu","1/fm^3");
    eost->set_unit("kfmu","1/fm");
  }
  
  double x;
  if (fabs(np_0)<1.0e-12) x=nb_start/3.0;
  else x=np_0;
  double oldpr=0.0;
  pressure_dec=0.0;
  
  funct sf=std::bind(std::mem_fn<double(double)>(&nstar_cold::solve_fun),
		       this,std::placeholders::_1);
  
  if (verbose>0) {
    cout << "baryon dens neutrons    protons     electrons   " 
	 << "energy dens pressure" << endl;
    cout << "[1/fm^3]    [1/fm^3]    [1/fm^3]    [1/fm^3]    "
	 << "[1/fm^4]    [1/fm^4]" << endl;
  }

  // Get a better initial guess for nb_start. This section of code was
  // put in to avoid negative densities for calc_eos from the solver
  // below. 
  barn=nb_start;
  bool done=false;
  for(size_t i=0;done==false && i<20;i++) {
    if (solve_fun(x)>0.0) {
      x=sqrt(nb_start*x);
    } else {
      done=true;
    }
  }
  
  for(barn=nb_start;barn<=nb_end+dnb/10.0;barn+=dnb) {
    
    int ret=rp->solve(x,sf);
    double y=solve_fun(x);
    
    if (ret!=0 || fabs(y)>solver_tol) {
      // We don't use the CONV macro here because we
      // want to return only if err_nonconv is true
      if (err_nonconv) {
	O2SCL_ERR("Solver failed in nstar_cold::calc_eos().",
		      exc_efailed);
      }
      return o2scl::exc_efailed;
    }

    // ------------------------------------------------------------
    // Compute dP/de at fixed Ye and Ymu. This code uses np.n and
    // pp.n, so we'll have to recompute them below.

    // Compute the hadronic part
    double dednb_Yp, dPdnb_Yp;
    hep->const_pf_derivs(barn,pp.n/barn,dednb_Yp,dPdnb_Yp);
    // Compute the leptonic part
    double dne_dmue=sqrt(e.mu*e.mu-e.m*e.m)*e.mu/pi2;
    double dP_dne=e.n/dne_dmue;
    // Put them together
    double numer=dPdnb_Yp+dP_dne*e.n/barn;
    double denom=dednb_Yp+e.mu*e.n/barn;
    // Add the muon contribution
    if (include_muons && mu.n>0.0) {
      double dnmu_dmumu=sqrt(mu.mu*mu.mu-mu.m*mu.m)*mu.mu/pi2;
      double dP_dnmu=mu.n/dnmu_dmumu;
      numer+=dP_dnmu*mu.n/barn;
      denom+=mu.mu*mu.n/barn;
    }
    double fcs2=numer/denom;

    // ------------------------------------------------------------

    // Recompute np.n and pp.n
    y=solve_fun(x);

    if (include_muons) {

      h=hb+e+mu;
      
      double line[18]={h.ed,h.pr,barn,np.mu,pp.mu,e.mu,np.n,pp.n,e.n,
		       np.kf,pp.kf,e.kf,fcs2,denom,numer,mu.mu,mu.n,mu.kf};
      eost->line_of_data(18,line);

    } else {

      h=hb+e;

      double line[15]={h.ed,h.pr,barn,np.mu,pp.mu,e.mu,np.n,pp.n,e.n,
		       np.kf,pp.kf,e.kf,fcs2,denom,numer};
      eost->line_of_data(15,line);

    }
    
    if (verbose>0) {
      cout.precision(5);
      cout << barn << " " << np.n << " " << pp.n << " " << e.n << " " 
	   << h.ed << " " << h.pr << endl;
      cout.precision(6);
    }
    
    if (barn>nb_start && pressure_dec<=0.0 && h.pr<oldpr) {
      pressure_dec=barn;
    }
    oldpr=h.pr;

    // Proceed to next baryon density
  }

  // Calculate the squared speed of sound. If the EOS becomes acausal,
  // calculate the density and pressure at which this happens

  well_formed=true;
  if (eost->get("ed",0)<=0.0) {
    well_formed=false;
    if (verbose>0) {
      cout << "Initial energy density is negative." << endl;
    }
  }
  if (pressure_dec>nb_start) {
    if (verbose>0) {
      cout << "Pressure is flat near a baryon density of " << pressure_dec
	   << " 1/fm^3." << endl;
    }
    well_formed=false;
  }
  for(size_t i=0;i<eost->get_nlines()-1;i++) {
    if (eost->get("ed",i)>eost->get("ed",i+1) || eost->get("pr",i)<0.0) {
      if (verbose>0) {
	cout << "Pressure is negative or energy density is decreasing" << endl;
	cout << " near a baryon density of " << eost->get("nb",i) 
	     << " 1/fm^3." << endl;
      }
      well_formed=false;
    }
  }

  if (verbose>0) {
    cout << "well_formed: " << well_formed << endl;
  }

  if (well_formed==false) {
    // We don't use the CONV macro here because we
    // want to return only if err_nonconv is true
    if (err_nonconv) {
      O2SCL_ERR("EOS not well-formed in nstar_cold::calc_eos().",
		    exc_efailed);
    }
  }

  acausal=0.0;
  acausal_pr=0.0;
  acausal_ed=0.0;

  if (well_formed) {
    eost->deriv("ed","pr","cs2");
    for(size_t i=0;i<eost->get_nlines()-1;i++) {
      if (eost->get("cs2",i)>=0.99) {
	acausal=eost->interp("cs2",1.0,"nb");
	acausal_pr=eost->interp("cs2",1.0,"pr");
	acausal_ed=eost->interp("cs2",1.0,"ed");
	i=eost->get_nlines();
      }
    }
  } else {
    eost->new_column("cs2");
    for(size_t i=0;i<eost->get_nlines();i++) {
      eost->set("cs2",i,0.0);
    }
  }

  if (verbose>0) {
    cout << "The EOS is acausal at nb=" << acausal << ", ed="
	 << acausal_ed << ", pr="
	 << acausal_pr << " ." << endl;
  }
    
  // Naively fix the equation of state so that it is not acausal

  if (false) {
    for(size_t i=0;i<eost->get_nlines()-1;i++) {
      if (eost->get("ed",i)>acausal_ed) {
	eost->set("pr",i,eost->get("ed",i)+acausal_pr-acausal_ed);
      }
    }
  }

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
  
  if (well_formed) {
    for(size_t i=0;i<eost->get_nlines();i++) {
      eost->set("logp",i,log(eost->get("pr",i)));
      eost->set("loge",i,log(eost->get("ed",i))); 
    }
    eost->deriv("loge","logp","ad_index");
  } else {
    for(size_t i=0;i<eost->get_nlines();i++) {
      eost->set("logp",i,0.0);
      eost->set("loge",i,0.0);
    }
    eost->new_column("ad_index");
    for(size_t i=0;i<eost->get_nlines();i++) {
      eost->set("ad_index",i,0.0);
    }
  }

  if (verbose>0) {
    cout << "Going to Urca, pass 2: " << endl;
  }

  deny_urca=0.0;
  allow_urca=0.0;
  for(size_t i=0;i<eost->get_nlines()-1;i++) {
    if (deny_urca<0.0001 && eost->get("urca",i)>0.0 && 
	eost->get("urca",i+1)<0.0) {
      deny_urca=eost->get("nb",i);
    } else if (allow_urca<0.0001 && eost->get("urca",i)>0.0) {
      /// Use linear interpolation to improve the urca estimate
      allow_urca=eost->get("nb",i)+
	(-eost->get("urca",i))*(eost->get("nb",i+1)-eost->get("nb",i))/
	(eost->get("urca",i+1)-eost->get("urca",i));
      
    }
  }

  if (verbose>0) {
    cout << "Urca: " << allow_urca << " " << deny_urca << endl;
    cout << "Done with calc_eos()." << endl;
  }
  
  return 0;
}

double nstar_cold::calc_urca(double np_0) {
  int ret;
  double old_urca=0.0, urca;
  
  double x;
  if (fabs(np_0)<1.0e-12) x=nb_start/3.0;
  else x=np_0;

  funct sf=std::bind(std::mem_fn<double(double)>
		       (&nstar_cold::solve_fun),
		       this,std::placeholders::_1);
  
  bool success=true;
  for(barn=nb_start;barn<=nb_end+dnb/10.0;barn+=dnb) {
    
    ret=rp->solve(x,sf);
    double y=solve_fun(x);
    
    if (ret!=0 || fabs(y)>1.0e-4) {
      success=false;
    }
    
    double s=(np.kf+pp.kf+e.kf)/2.0;
    urca=s*(s-np.kf)*(s-pp.kf)*(s-e.kf);
    if (barn>nb_start && urca>0.0 && old_urca<0.0) {
      if (success==false) {
	O2SCL_ERR("Solution failed in calc_urca().",exc_efailed);
      }
      return barn-dnb*urca/(urca-old_urca);
    }
    
    old_urca=urca;
  }

  if (success==false) {
    O2SCL_ERR("Solution failed in calc_urca().",exc_efailed);
  }
  
  return 0.0;
}

int nstar_cold::calc_nstar() {
  def_eos_tov.read_table(*eost,"ed","pr","nb");
  
  tp->set_units("1/fm^4","1/fm^4","1/fm^3");
  tp->set_eos(def_eos_tov);
  
  return tp->mvsr();
}

int nstar_cold::fixed(double target_mass) {
  def_eos_tov.read_table(*eost,"ed","pr","nb");
  
  tp->set_units("1/fm^4","1/fm^4","1/fm^3");
  tp->set_eos(def_eos_tov);
  
  return tp->fixed(target_mass);
}

