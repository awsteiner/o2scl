/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

fermion_deriv_rel::fermion_deriv_rel() {
  
  deg_limit=2.0;
  upper_limit_fac=20.0;

  density_root=&def_density_root;
  nit=&def_nit;
  dit=&def_dit;
  
  method=automatic;
  intl_method=by_parts;

  exp_limit=200.0;

  err_nonconv=true;

  last_method=0;
}

fermion_deriv_rel::~fermion_deriv_rel() {
}

void fermion_deriv_rel::set_inte(inte<funct> &l_nit, inte<funct> &l_dit) {
  nit=&l_nit;
  dit=&l_dit;
  return;
}

int fermion_deriv_rel::calc_mu(fermion_deriv &f, double temper) {

  fr.calc_mu_tlate<fermion_deriv>(f,temper);
  last_method=fr.last_method*10;
  
  int iret;

  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in fermion_deriv_rel().",exc_eunimpl);
  }

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

  double prefac=f.g/2.0/pi2;

  // Compute the degeneracy parameter

  bool deg=false;
  double psi;
  if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
  else psi=(f.nu+f.m-f.ms)/temper;
  if (psi>deg_limit) deg=true;

  // Try the non-degenerate expansion if psi is small enough
  if (psi<-4.0) {
    bool acc=calc_mu_ndeg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      unc.dndT=f.dndT*1.0e-14;
      unc.dsdT=f.dsdT*1.0e-14;
      unc.dndmu=f.dndmu*1.0e-14;
      last_method+=1;
      return 0;
    }
  }

  // Try the degenerate expansion if psi is large enough
  if (psi>20.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      unc.dndT=f.dndT*1.0e-14;
      unc.dsdT=f.dsdT*1.0e-14;
      unc.dndmu=f.dndmu*1.0e-14;
      last_method+=2;
      return 0;
    }
  }

  if (deg==false) {
    
    // Set integration method
    if (method==automatic) {
      intl_method=by_parts;
      last_method+=3;
    } else {
      intl_method=method;
      last_method+=4;
    }

    // The non-degenerate case

    funct density_T_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::density_T_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(density_T_fun_f,0.0,0.0,f.dndT,unc.dndT);
    if (iret!=0) {
      O2SCL_ERR2("dndT integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",
		 exc_efailed);
    }
    f.dndT*=prefac;
    unc.dndT*=prefac;

    funct density_mu_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::density_mu_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(density_mu_fun_f,0.0,0.0,f.dndmu,unc.dndmu);
    if (iret!=0) {
      O2SCL_ERR2("dndmu integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",
		 exc_efailed);
    }
    f.dndmu*=prefac;
    unc.dndmu*=prefac;
    
    funct entropy_T_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::entropy_T_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(entropy_T_fun_f,0.0,0.0,f.dsdT,unc.dsdT);
    if (iret!=0) {
      O2SCL_ERR2("dsdT integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.dsdT*=prefac;
    unc.dsdT*=prefac;

  } else {

    // Compute the upper limit for degenerate integrals

    double arg;
    if (f.inc_rest_mass) {
      arg=pow(upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
    } else {
      arg=pow(upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
    }
    // Set to zero to avoid uninit'ed var. warnings
    double ul=0.0;
    if (arg>0.0) {
      ul=sqrt(arg);
    } else {
      O2SCL_ERR2("Zero density in degenerate limit in fermion_deriv_rel::",
		 "calc_mu(). Variable deg_limit set improperly?",
		 exc_efailed);
    }
    
    // Compute the lower limit for the entropy and derivative
    // integrations
    
    double ll;
    if (f.inc_rest_mass) {
      arg=pow(-upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
      if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
        ll=sqrt(arg);
      } else {
        ll=-1.0;
      }
    } else {
      arg=pow(-upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
      if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
        ll=sqrt(arg);
      } else {
        ll=-1.0;
      }
    }

    // Set integration method
    if (method==automatic) {
      if ((!f.inc_rest_mass && (f.nu+f.m-f.ms)/temper>1.0e3) ||
	  (f.inc_rest_mass && (f.nu-f.ms)/temper>1.0e3)) {
	intl_method=direct;
	last_method+=5;
      } else {
	intl_method=by_parts;
	last_method+=6;
      }
    } else {
      intl_method=method;
      last_method+=7;
    }

    funct deg_density_mu_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::deg_density_mu_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    if (intl_method==direct && ll>0.0) {
      iret=dit->integ_err(deg_density_mu_fun_f,ll,ul,
			  f.dndmu,unc.dndmu);
    } else {
      iret=dit->integ_err(deg_density_mu_fun_f,0.0,ul,
			  f.dndmu,unc.dndmu);
    }
    if (iret!=0) {
      O2SCL_ERR2("dndmu integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.dndmu*=prefac;
    unc.dndmu*=prefac;
    
    funct deg_density_T_fun_f=std::bind
      (std::mem_fn<double(double,fermion_deriv &,double)>
       (&fermion_deriv_rel::deg_density_T_fun),
       this,std::placeholders::_1,std::ref(f),temper);
    if (intl_method==direct && ll>0.0) {
      iret=dit->integ_err(deg_density_T_fun_f,ll,ul,f.dndT,unc.dndT);
    } else {
      iret=dit->integ_err(deg_density_T_fun_f,0.0,ul,f.dndT,unc.dndT);
    }
    if (iret!=0) {
      O2SCL_ERR2("dndT integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.dndT*=prefac;
    unc.dndT*=prefac;

    funct deg_entropy_T_fun_f=std::bind
      (std::mem_fn<double(double,fermion_deriv &,double)>
       (&fermion_deriv_rel::deg_entropy_T_fun),
       this,std::placeholders::_1,std::ref(f),temper);
    if (intl_method==direct && ll>0.0) {
      iret=dit->integ_err(deg_entropy_T_fun_f,ll,ul,f.dsdT,unc.dsdT);
    } else {
      iret=dit->integ_err(deg_entropy_T_fun_f,0.0,ul,f.dsdT,unc.dsdT);
    }
    if (iret!=0) {
      O2SCL_ERR2("dsdT integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.dsdT*=prefac;
    unc.dsdT*=prefac;

  }
  
  if (!std::isfinite(f.en)) {
    O2SCL_ERR2("Entropy not finite in ",
	       "fermion_deriv_rel::calc_mu().",exc_efailed);
  }
  f.pr=-f.ed+temper*f.en+f.nu*f.n;
  // Pressure uncertainties are not computed
  unc.pr=0.0;
  
  return 0;
}

int fermion_deriv_rel::nu_from_n(fermion_deriv &f, double temper) {
  int ret=fr.nu_from_n_tlate<fermion_deriv>(f,temper);
  last_method=fr.last_method*100;
  return ret;
}

int fermion_deriv_rel::calc_density(fermion_deriv &f, double temper) {
  
  if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
  nu_from_n(f,temper);
  int lm=last_method;
  
  if (f.non_interacting) { f.mu=f.nu; }

  calc_mu(f,temper);
  last_method+=lm;

  return 0;
}

double fermion_deriv_rel::deg_density_T_fun(double k, fermion_deriv &f, 
					    double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*(E-f.nu)/T/T*ff*(1.0-ff);
    } else {
      ret=(2.0*k*k/T+E*E/T-E*(f.nu)/T-k*k*(f.nu)/T/E)*
	fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*(E-f.nu)/T/T*ff*(1.0-ff);
    } else {
      ret=(2.0*k*k/T+E*E/T-E*(f.nu+f.m)/T-k*k*(f.nu+f.m)/T/E)*
	fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::deg_density_mu_fun(double k, fermion_deriv &f, 
					     double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k/T*ff*(1.0-ff);
    } else {
      ret=(E*E+k*k)/E*fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k/T*ff*(1.0-ff);
    } else {
      ret=(E*E+k*k)/E*fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::deg_entropy_T_fun(double k, fermion_deriv &f, 
					    double T) {
  double E, ret;
  E=gsl_hypot(k,f.ms);
  if (f.inc_rest_mass) {
    double ff=fermi_function(E,f.nu,T,exp_limit);
    if (intl_method==direct) {
      ret=k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.nu)/E/T/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*f.nu)*ff;
    }
  } else {
    double ff=fermi_function(E-f.m,f.nu,T,exp_limit);
    if (intl_method==direct) {
      ret=k*k*ff*(1.0-ff)*pow(E-f.nu-f.m,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.m-f.nu)/E/T/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu+f.m))*ff;
    }
  }
  return ret;
}

double fermion_deriv_rel::deg_density_ms_fun(double k, fermion_deriv &f, 
					     double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=-k*k*f.ms/E/T*ff*(1.0-ff);
    } else {
      ret=-f.ms*fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=-k*k*f.ms/(E+f.m)/T*ff*(1.0-ff);
    } else {
      ret=-f.ms*fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::density_T_fun(double u, fermion_deriv &f, 
					double T) {
  double k=u*(T), E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*(E-f.nu)/T*ff*(1.0-ff);
    } else {
      ret=(2.0*k*k/T+E*E/T-E*(f.nu)/T-k*k*(f.nu)/T/E)*
	T*fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*(E-f.nu)/T*ff*(1.0-ff);
    } else {
      ret=(2.0*k*k/T+E*E/T-E*(f.nu+f.m)/T-k*k*(f.nu+f.m)/T/E)*
	T*fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::density_mu_fun(double u, fermion_deriv &f, 
					 double T) {
  double k=u*(T), E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*ff*(1.0-ff);
    } else {
      ret=T*(E*E+k*k)/E*fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*ff*(1.0-ff);
    } else {
      ret=T*(E*E+k*k)/E*fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::entropy_T_fun(double u, fermion_deriv &f, 
					double T) {
  double k=u*T, E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=T*k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.nu)/E/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu))*
	fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=T*k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.m-f.nu)/E/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu+f.m))*
	fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  return ret;
}

double fermion_deriv_rel::density_ms_fun(double u, fermion_deriv &f, 
					 double T) {
  double k=u*T, E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=-k*k*f.ms/(E)/T*ff*(1.0-ff);
    } else {
      ret=-f.ms*fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=-k*k*f.ms/(E+f.m)/T*ff*(1.0-ff);
    } else {
      ret=-f.ms*fermi_function(E-f.m,f.nu,T,exp_limit);
    }
  }
  ret*=T;
  return ret;
}

int fermion_deriv_rel::pair_mu(fermion_deriv &f, double temper) {

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
  fermion_deriv antip(f.ms,f.g);
  f.anti(antip);

  calc_mu(f,temper);
  int lm=last_method*100;
  calc_mu(antip,temper);
  last_method+=lm;

  f.n-=antip.n;
  f.pr+=antip.pr;
  if (f.inc_rest_mass) {
    f.ed+=antip.ed;
  } else {
    f.ed=f.ed+antip.ed+2.0*antip.n*f.m;
  }
  f.en+=antip.en;
  f.dsdT+=antip.dsdT;
  f.dndT-=antip.dndT;
  f.dndmu+=antip.dndmu;
  
  return 0;
}

int fermion_deriv_rel::pair_density(fermion_deriv &f, double temper) {
  int ret=fr.pair_density_tlate<fermion_deriv>(f,temper);
  //cout << "ret: " << ret << endl;
  pair_mu(f,temper);
  return 0;
}
