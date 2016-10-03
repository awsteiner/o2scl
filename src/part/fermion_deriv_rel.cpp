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
}

fermion_deriv_rel::~fermion_deriv_rel() {
}

void fermion_deriv_rel::set_inte(inte<funct11> &l_nit, inte<funct11> &l_dit) {
  nit=&l_nit;
  dit=&l_dit;
  return;
}

int fermion_deriv_rel::calc_mu(fermion_deriv &f, double temper) {
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

  // Try the degenerate expansion if psi is large enough
  if (psi>20.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      return 0;
    }
  }

  if (deg==false) {
    
    // Set integration method
    if (method==automatic) {
      intl_method=by_parts;
    } else {
      intl_method=method;
    }

    // The non-degenerate case

    funct11 density_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::density_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(density_fun_f,0.0,0.0,f.n,unc.n);
    if (iret!=0) {
      O2SCL_ERR2("Density integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.n*=prefac;
    unc.n*=prefac;
    
    funct11 density_T_fun_f=
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

    funct11 density_mu_fun_f=
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
    
    funct11 energy_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::energy_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(energy_fun_f,0.0,0.0,f.ed,unc.ed);
    if (iret!=0) {
      O2SCL_ERR2("Energy integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.ed*=prefac;
    f.ed*=pow(temper,4.0);
    unc.ed*=prefac;
    
    funct11 entropy_fun_f=
      std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		(&fermion_deriv_rel::entropy_fun),
		this,std::placeholders::_1,std::ref(f),temper);
    iret=nit->integ_err(entropy_fun_f,0.0,0.0,f.en,unc.en);
    if (iret!=0) {
      O2SCL_ERR2("Entropy integration (ndeg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.en*=prefac;
    unc.en*=prefac;
    
    funct11 entropy_T_fun_f=
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
    double ul;
    if (arg>0.0) {
      ul=sqrt(arg);
    } else {
      f.n=0.0;
      f.ed=0.0;
      f.pr=0.0;
      f.en=0.0;
      unc.n=0.0;
      unc.ed=0.0;
      unc.pr=0.0;
      unc.en=0.0;
      O2SCL_ERR2("Zero density in degenerate limit in fermion_deriv_rel::",
		 "calc_mu(). Variable deg_limit set improperly?",
		 exc_efailed);
    }
    
    // Compute the lower limit
    
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
      } else {
	intl_method=by_parts;
      }
    } else {
      intl_method=method;
    }
    
    funct11 deg_density_fun_f=std::bind
      (std::mem_fn<double(double,fermion_deriv &,double)>
       (&fermion_deriv_rel::deg_density_fun),
       this,std::placeholders::_1,std::ref(f),temper);
    iret=dit->integ_err(deg_density_fun_f,0.0,ul,f.n,unc.n);
    if (iret!=0) {
      O2SCL_ERR2("Density integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.n*=prefac;
    unc.n*=prefac;
    
    funct11 deg_density_mu_fun_f=
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
    
    funct11 deg_density_T_fun_f=std::bind
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

    funct11 deg_energy_fun_f=std::bind
      (std::mem_fn<double(double,fermion_deriv &,double)>
       (&fermion_deriv_rel::deg_energy_fun),
       this,std::placeholders::_1,std::ref(f),temper);
    iret=dit->integ_err(deg_energy_fun_f,0.0,ul,f.ed,unc.ed);
    if (iret!=0) {
      O2SCL_ERR2("Energy integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.ed*=prefac;
    unc.ed*=prefac;

    funct11 deg_entropy_fun_f=std::bind
      (std::mem_fn<double(double,fermion_deriv &,double)>
       (&fermion_deriv_rel::deg_entropy_fun),
       this,std::placeholders::_1,std::ref(f),temper);
    if (ll>0.0) {
      iret=dit->integ_err(deg_entropy_fun_f,ll,ul,f.en,unc.en);
    } else {
      iret=dit->integ_err(deg_entropy_fun_f,0.0,ul,f.en,unc.en);
    }
    if (iret!=0) {
      O2SCL_ERR2("Entropy integration (deg) failed in ",
		 "fermion_deriv_rel::calc_mu().",exc_efailed);
    }
    f.en*=prefac;
    unc.en*=prefac;
    
    funct11 deg_entropy_T_fun_f=std::bind
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
  
  return 0;
}

int fermion_deriv_rel::nu_from_n(fermion_deriv &f, double temper) {
  double nex;

  // Check that initial guess is close enough

  nex=f.nu/temper;
  double y=solve_fun(nex,f,temper);

  // If that didn't work, try a different guess
  if (y==1.0) {
    if (f.inc_rest_mass) {
      nex=f.ms/temper;
    } else {
      nex=(f.ms-f.m)/temper;
    }
    y=solve_fun(nex,f,temper);
  }

  // If nothing worked, call the error handler
  if (y==1.0) {
    O2SCL_ERR2("Couldn't find reasonable initial guess in ",
	       "fermion_deriv_rel::nu_from_n().",exc_efailed);
  }

  // Perform full solution
  funct11 mf=std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		       (&fermion_deriv_rel::solve_fun),
		       this,std::placeholders::_1,std::ref(f),temper);
  int ret=density_root->solve(nex,mf);
  if (ret!=0) {
    O2SCL_ERR("Solver failed in fermion_deriv_rel::nu_from_n().",exc_efailed);
  }

  f.nu=nex*temper;
  
  return 0;
}

int fermion_deriv_rel::calc_density(fermion_deriv &f, double temper) {

  if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
  nu_from_n(f,temper);
  
  if (f.non_interacting) { f.mu=f.nu; }
  
  calc_mu(f,temper);

  return 0;
}

double fermion_deriv_rel::deg_density_fun(double k, fermion_deriv &f, 
					  double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }
  ret=k*k*fermi_function(E,f.nu,T,exp_limit);
  return ret;
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

double fermion_deriv_rel::deg_energy_fun(double k, fermion_deriv &f, 
					 double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }
  ret=k*k*E*fermi_function(E,f.nu,T,exp_limit);
  return ret;
}

double fermion_deriv_rel::deg_entropy_fun(double k, fermion_deriv &f, 
					  double T) {

  double E, ret, nx, nm1;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }
  
  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-f.nu)/(T))<-200.0) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-f.nu)/T)<-exp_limit) {
    // Should this be E/T-nu/T or (E-nu)/T ?
    nm1=-exp(E/T-f.nu/T);
    ret=k*k*nm1*log(-nm1);
  } else {
    nx=fermi_function(E,f.nu,T,exp_limit);
    if (nx==0.0) ret=0.0;
    else ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx));
  }
  if (!std::isfinite(ret)) {
    ret=0.0;
    //O2SCL_ERR("Entropy not finite in fermion_deriv_rel::deg_entropy_fun().",
    //exc_efailed);
  }

  return ret;
}
  
double fermion_deriv_rel::deg_entropy_T_fun(double k, fermion_deriv &f, 
					    double T) {
  double E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.nu)/E/T/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu))*
	fermi_function(E,f.nu,T,exp_limit);
    }
  } else {
    E=gsl_hypot(k,f.ms);
    if (intl_method==direct) {
      E-=f.m;
      double ff=fermi_function(E,f.nu,T,exp_limit);
      ret=k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
    } else {
      ret=(E-f.m-f.nu)/E/T/T*
	(pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu+f.m))*
	fermi_function(E-f.m,f.nu,T,exp_limit);
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
  return ret;
}

double fermion_deriv_rel::density_fun(double u, fermion_deriv &f, 
				      double T) {
  double k=u*(T), E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }
  ret=T*k*k*fermi_function(E,f.nu,T,exp_limit);
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

double fermion_deriv_rel::energy_fun(double u, fermion_deriv &f, double T) {
  double k=u*(T), E, ret;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }
  ret=u*u*E*fermi_function(E,f.nu,T,exp_limit)/T;
  return ret;
}

double fermion_deriv_rel::entropy_fun(double u, fermion_deriv &f, double T) {
  double k=u*(T), E, ret, nx, nm1;
  if (f.inc_rest_mass) {
    E=gsl_hypot(k,f.ms);
  } else {
    E=gsl_hypot(k,f.ms)-f.m;
  }

  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-f.nu)/(T))<-200.0) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-f.nu)/T)<-30.0) {
    nm1=-exp(E/(T)-f.nu/(T));
    ret=k*k*nm1*log(-nm1)*(T);
  } else {
    nx=fermi_function(E,f.nu,T,exp_limit);
    if (nx==0.0) ret=0.0;
    else ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx))*T;
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

double fermion_deriv_rel::solve_fun(double x, fermion_deriv &f, double T) {
  double nden, yy;
  
  f.nu=T*x;
  
  // 6/6/03 - Should this be included? I think yes!
  if (f.non_interacting) f.mu=f.nu;
  
  bool deg=true;
  if (f.inc_rest_mass) {
    if ((f.nu-f.ms)/T<deg_limit) deg=false;
  } else {
    if ((f.nu+f.m-f.ms)/T<deg_limit) deg=false;
  }
  
  funct11 density_fun_f=
    std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
	      (&fermion_deriv_rel::density_fun),
	      this,std::placeholders::_1,std::ref(f),T);
  funct11 deg_density_fun_f=
    std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
	      (&fermion_deriv_rel::deg_density_fun),
	      this,std::placeholders::_1,std::ref(f),T);
  
  // Set integration method
  if (method==automatic) {
    if ((!f.inc_rest_mass && (f.nu+f.m-f.ms)/T>1.0e3) ||
	(f.inc_rest_mass && (f.nu-f.ms)/T>1.0e3)) {
      intl_method=direct;
    } else {
      intl_method=by_parts;
    }
  } else {
    intl_method=method;
  }
  
  if (!deg) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=f.g/2.0/pi2;
    
    yy=(f.n-nden)/f.n;

  } else {
    
    double arg;
    if (f.inc_rest_mass) {
      arg=pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms;
    } else {
      arg=pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms;
    }

    double ul;

    if (arg>0.0) {
      
      ul=sqrt(arg);
      nden=dit->integ(deg_density_fun_f,0.0,ul);
      nden*=f.g/2.0/pi2;

    } else {
      nden=0.0;
    }
    
    yy=(f.n-nden)/f.n;
  }

  return yy;
}

int fermion_deriv_rel::pair_mu(fermion_deriv &f, double temper) {

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
  fermion_deriv antip(f.ms,f.g);
  f.anti(antip);
  
  calc_mu(f,temper);
  calc_mu(antip,temper);

  f.n-=antip.n;
  f.pr+=antip.pr;
  f.ed+=antip.ed;
  f.en+=antip.en;
  f.dsdT+=antip.dsdT;
  f.dndT+=antip.dndT;
  f.dndmu+=antip.dndmu;
  
  return 0;
}

int fermion_deriv_rel::pair_density(fermion_deriv &f, double temper) {
  double nex;
  int ret;
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in fermion_deriv_rel().",exc_eunimpl);
  }

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  nex=f.nu/temper;
  funct11 mf=std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		       (&fermion_deriv_rel::pair_fun),
		       this,std::placeholders::_1,std::ref(f),temper);
  ret=density_root->solve(nex,mf);
  f.nu=nex*temper;

  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  return 0;
}

double fermion_deriv_rel::pair_fun(double x, fermion_deriv &f, double T) {
  double nden, yy;
  
  f.nu=T*x;
  
  if (f.non_interacting) f.mu=f.nu;
  
  // Set integration method
  if (method==automatic) {
    if ((!f.inc_rest_mass && (f.nu+f.m-f.ms)/T>1.0e3) ||
	(f.inc_rest_mass && (f.nu-f.ms)/T>1.0e3)) {
      intl_method=direct;
    } else {
      intl_method=by_parts;
    }
  } else {
    intl_method=method;
  }

  funct11 density_fun_f=
    std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
	      (&fermion_deriv_rel::density_fun),
	      this,std::placeholders::_1,std::ref(f),T);
  funct11 deg_density_fun_f=
    std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
	      (&fermion_deriv_rel::deg_density_fun),
	      this,std::placeholders::_1,std::ref(f),T);
  
  if (f.nu/T<deg_limit) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=f.g/2.0/pi2;
    yy=nden;

  } else {
    
    double ulimit;
    if (f.inc_rest_mass) {
      ulimit=sqrt(pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms);
    } else {
      ulimit=sqrt(pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms);
    }
    nden=dit->integ(deg_density_fun_f,0.0,ulimit);
    nden*=f.g/2.0/pi2;
    yy=nden;
  }
  
  if (f.inc_rest_mass) f.nu=-T*x;
  else f.nu=-T*x-2.0*f.m;
  
  if (f.nu/T < deg_limit) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=f.g/2.0/pi2;
    yy-=nden;
    
  } else {
    
    double ulimit;
    if (f.inc_rest_mass) {
      ulimit=sqrt(pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms);
    } else {
      ulimit=sqrt(pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms);
    }
    nden=dit->integ(deg_density_fun_f,0.0,ulimit);
    nden*=f.g/2.0/pi2;
    yy-=nden;
  }
  
  yy=yy/f.n-1.0;

  return yy;
}

double fermion_deriv_rel::deriv_calibrate(fermion_deriv &f, int verbose, 
					  std::string fname) {

  double ret=0;
  
  // ----------------------------------------------------------------
  // Will return to these original values afterwards

  part orig;
  orig.mu=f.mu;
  orig.m=f.m;
  orig.ms=f.ms;
  orig.g=f.g;
  orig.non_interacting=f.non_interacting;
  orig.inc_rest_mass=f.inc_rest_mass;

  // ----------------------------------------------------------------
  // Read data file
  
  if (fname=="") {
    fname=o2scl_settings.get_data_dir()+"fermion_cal.o2";
  }

  if (verbose>1) {
    cout << "In fermion_deriv_rel::deriv_calibrate(), loading file named\n\t'" 
	 << fname << "'.\n" << endl;
  }
  table<> tab;
  hdf_file hf;
  hf.open(fname);
  string name;
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input(hf,tab,name);
#endif
  hf.close();
  
  if (tab.get_nlines()==0) {
    string str="Failed to load data from file '"+fname+
      "' in fermion_T::calibrate(). Bad filename?";
    O2SCL_ERR(str.c_str(),exc_efilenotfound);
  }
  
  // ----------------------------------------------------------------

  f.g=2.0;
  
  size_t cnt=0;
  fermion_deriv bad, dev, exact;
  double m_bad=0.0, mu_bad=0.0, T_bad=0.0, psi_bad=0.0, mot_bad=0.0;
  f.non_interacting=true;

  // ----------------------------------------------------------------
  // First pass, test calc_mu() 

  // k=0 is with rest mass, k=1 is without
  for(size_t k=1;k<2;k++) {

    // Initialize storage
    dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
    dev.dndT=0.0; dev.dndmu=0.0; dev.dsdT=0.0; 
    bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    bad.dndT=0.0; bad.dndmu=0.0; bad.dsdT=0.0; 
    
    // Temperature loop
    for(double T2=1.0e-2;T2<=1.001e2;T2*=1.0e2) {

      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("mot",i);
	double psi=tab.get("psi",i);
	exact.n=tab.get("n",i);
	exact.ed=tab.get("ed",i);
	exact.pr=tab.get("pr",i);
	exact.en=tab.get("en",i);
	exact.dndT=tab.get("dndT",i);
	exact.dndmu=tab.get("dndmu",i);
	exact.dsdT=tab.get("dsdT",i);
      
	if (k==0) {
	  
	  f.inc_rest_mass=true;
          
	  f.m=mot*T2;
	  f.mu=f.m+T2*psi;
	  
	} else {
	  
	  f.inc_rest_mass=false;
	  
	  f.m=mot*T2;
	  f.mu=T2*psi;
	  
	}
	
	calc_mu(f,T2);
	
	exact.n*=pow(T2,3.0);
	if (k==0) {
	  exact.ed*=pow(T2,4.0);
	} else {
	  exact.ed=exact.ed*pow(T2,4.0)-exact.n*f.m;
	}
	exact.pr*=pow(T2,4.0);
	exact.en*=pow(T2,3.0);
	exact.dndT*=pow(T2,2.0);
	exact.dndmu*=pow(T2,2.0);
	exact.dsdT*=pow(T2,2.0);
	
	if (verbose>1) {
	  cout << "T,m,mu: " << T2 << " " << f.m << " " << f.mu << endl;
	  cout << "n,ed,pr,en: " << endl;
	  cout << "approx: " << f.n << " " << f.ed << " " << f.pr << " " 
	       << f.en << endl;
	  cout << "\t" << f.dndT << " " << f.dndmu << " " << f.dsdT << endl;
	  cout << "exact : " << exact.n << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "\t" << exact.dndT << " " << exact.dndmu << " " 
	       << exact.dsdT << endl;
	  cout << "bad   : " << bad.n << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << "\t" << bad.dndT << " " << bad.dndmu << " " 
	       << bad.dsdT << endl;
	  cout << endl;
	  if (verbose>2) {
	    char ch;
	    cin >> ch;
	  }
	}

	dev.n+=fabs((f.n-exact.n)/exact.n);
	dev.ed+=fabs((f.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((f.pr-exact.pr)/exact.pr);
	dev.en+=fabs((f.en-exact.en)/exact.en);
	dev.dndT+=fabs((f.dndT-exact.dndT)/exact.dndT);
	dev.dndmu+=fabs((f.dndmu-exact.dndmu)/exact.dndmu);
	dev.dsdT+=fabs((f.dsdT-exact.dsdT)/exact.dsdT);
	
	cnt++;
	if (fabs((f.n-exact.n)/exact.n)>bad.n) {
	  bad.n=fabs((f.n-exact.n)/exact.n);
	  if (bad.n>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.n;
	  }
	}
	if (fabs((f.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((f.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((f.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((f.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((f.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((f.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}
	if (fabs((f.dndT-exact.dndT)/exact.dndT)>bad.dndT) {
	  bad.dndT=fabs((f.dndT-exact.dndT)/exact.dndT);
	  if (bad.dndT>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dndT;
	  }
	}
	if (fabs((f.dndmu-exact.dndmu)/exact.dndmu)>bad.dndmu) {
	  bad.dndmu=fabs((f.dndmu-exact.dndmu)/exact.dndmu);
	  if (bad.dndmu>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dndmu;
	  }
	}
	if (fabs((f.dsdT-exact.dsdT)/exact.dsdT)>bad.dsdT) {
	  bad.dsdT=fabs((f.dsdT-exact.dsdT)/exact.dsdT);
	  if (bad.dsdT>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dsdT;
	  }
	}

	// End of loop over points in data file
      }
      // End of temperature loop
    }

    dev.n/=cnt;
    dev.ed/=cnt;
    dev.pr/=cnt;
    dev.en/=cnt;
    dev.dndT/=cnt;
    dev.dndmu/=cnt;
    dev.dsdT/=cnt;

    if (verbose>0) {
      if (k==0) {
	cout << "Function calc_mu(), include rest mass:" << endl;
      } else {
	cout << "Function calc_mu(), without rest mass:" << endl;
      }

      cout << "Average performance: " << endl;
      cout << "n: " << dev.n << " ed: " << dev.ed << " pr: " 
	   << dev.pr << " en: " << dev.en << endl;
      cout << "dndT: " << dev.dndT << " dndmu: " << dev.dndmu 
	   << " dsdT: " << dev.dsdT << endl;
      cout << "Worst case: " << endl;
      cout << "mu: " << mu_bad << " m: " << m_bad << " T: " << T_bad 
	   << " mot: " << mot_bad << "\n\tpsi: " << psi_bad << endl;
      cout << "n: " << bad.n << " ed: " << bad.ed << " pr: " 
	   << bad.pr << " en: " << bad.en << endl;
      cout << "dndT: " << bad.dndT << " dndmu: " << bad.dndmu 
	   << " dsdT: " << bad.dsdT << endl;
      cout << endl;
      if (verbose>2) {
	char ch;
	cin >> ch;
      }
    }

    // End of k loop
  }

  // ----------------------------------------------------------------
  // Second pass, test calc_density()

  // k=0 is with rest mass, k=1 is without
  for(size_t k=1;k<2;k++) {

    // Initialize storage
    dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
    dev.dndT=0.0; dev.dndmu=0.0; dev.dsdT=0.0; 
    bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    bad.dndT=0.0; bad.dndmu=0.0; bad.dsdT=0.0; 
    
    // Temperature loop
    for(double T2=1.0e-2;T2<=1.001e2;T2*=1.0e2) {
      
      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("mot",i);
	double psi=tab.get("psi",i);
	f.n=tab.get("n",i);	
	exact.ed=tab.get("ed",i);
	exact.pr=tab.get("pr",i);
	exact.en=tab.get("en",i);
	exact.dndT*=pow(T2,2.0);
	exact.dndmu*=pow(T2,2.0);
	exact.dsdT*=pow(T2,2.0);

	f.m=mot*T2;
	if (k==0) {
	  f.inc_rest_mass=true;
	  exact.mu=f.m+T2*psi;
	} else {
	  f.inc_rest_mass=false;
	  exact.mu=T2*psi;
	}

	f.n*=pow(T2,3.0);
	if (k==0) {
	  exact.ed*=pow(T2,4.0);
	} else {
	  exact.ed=exact.ed*pow(T2,4.0)-f.n*f.m;
	}
	exact.pr*=pow(T2,4.0);
	exact.en*=pow(T2,3.0);

	// Give it a guess for the chemical potential
	f.mu=f.m;

	calc_density(f,T2);
	
	if (verbose>1) {
	  cout << "T2, m, n: " << T2 << " " << f.m << " " << f.n << endl;
	  cout << "mu,ed,pr,en: " << endl;
	  cout << "approx: " << f.mu << " " << f.ed << " " << f.pr << " " 
	       << f.en << endl;
	  cout << "\t" << f.dndT << " " << f.dndmu << " " << f.dsdT << endl;
	  cout << "exact : " << exact.mu << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "\t" << exact.dndT << " " << exact.dndmu << " " 
	       << exact.dsdT << endl;
	  cout << "bad   : " << bad.mu << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << "\t" << bad.dndT << " " << bad.dndmu << " " 
	       << bad.dsdT << endl;
	  cout << endl;
	  if (verbose>2) {
	    char ch;
	    cin >> ch;
	  }
	}

	dev.mu+=fabs((f.mu-exact.mu)/exact.mu);
	dev.ed+=fabs((f.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((f.pr-exact.pr)/exact.pr);
	dev.en+=fabs((f.en-exact.en)/exact.en);
	dev.dndT+=fabs((f.dndT-exact.dndT)/exact.dndT);
	dev.dndmu+=fabs((f.dndmu-exact.dndmu)/exact.dndmu);
	dev.dsdT+=fabs((f.dsdT-exact.dsdT)/exact.dsdT);
	
	cnt++;
	if (fabs((f.mu-exact.mu)/exact.mu)>bad.mu) {
	  bad.mu=fabs((f.mu-exact.mu)/exact.mu);
	  if (bad.n>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.n;
	  }
	}
	if (fabs((f.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((f.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((f.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((f.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((f.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((f.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}
	if (fabs((f.dndT-exact.dndT)/exact.dndT)>bad.dndT) {
	  bad.dndT=fabs((f.dndT-exact.dndT)/exact.dndT);
	  if (bad.dndT>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dndT;
	  }
	}
	if (fabs((f.dndmu-exact.dndmu)/exact.dndmu)>bad.dndmu) {
	  bad.dndmu=fabs((f.dndmu-exact.dndmu)/exact.dndmu);
	  if (bad.dndmu>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dndmu;
	  }
	}
	if (fabs((f.dsdT-exact.dsdT)/exact.dsdT)>bad.dsdT) {
	  bad.dsdT=fabs((f.dsdT-exact.dsdT)/exact.dsdT);
	  if (bad.dsdT>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T2;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.dsdT;
	  }
	}

	// End of loop over points in data file
      }
      // End of temperature loop
    }

    dev.mu/=cnt;
    dev.ed/=cnt;
    dev.pr/=cnt;
    dev.en/=cnt;
    dev.dndT/=cnt;
    dev.dndmu/=cnt;
    dev.dsdT/=cnt;

    if (verbose>0) {
      if (k==0) {
	cout << "Function calc_density(), include rest mass:" << endl;
      } else {
	cout << "Function calc_density(), without rest mass:" << endl;
      }

      cout << "Average performance: " << endl;
      cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
	   << dev.pr << " en: " << dev.en << endl;
      cout << "dndT: " << dev.dndT << " dndmu: " << dev.dndmu 
	   << " dsdT: " << dev.dsdT << endl;
      cout << "Worst case: " << endl;
      cout << "mu: " << mu_bad << " m: " << m_bad << " T: " << T_bad 
	   << " mot: " << mot_bad << "\n\tpsi: " << psi_bad << endl;
      cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
	   << bad.pr << " en: " << bad.en << endl;
      cout << "dndT: " << bad.dndT << " dndmu: " << bad.dndmu 
	   << " dsdT: " << bad.dsdT << endl;
      cout << endl;
      if (verbose>2) {
	char ch;
	cin >> ch;
      }
    }

    // End of k loop
  }

  // ----------------------------------------------------------------
  // Return to the original values 

  f.mu=orig.mu;
  f.m=orig.m;
  f.ms=orig.ms;
  f.g=orig.g;
  f.non_interacting=orig.non_interacting;
  f.inc_rest_mass=orig.inc_rest_mass;
  
  return ret;
}

