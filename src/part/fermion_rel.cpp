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

#include <o2scl/fermion_rel.h>
#include <o2scl/root_cern.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

fermion_rel::fermion_rel() : nit(new inte_qagiu_gsl<>), 
			     dit(new inte_qag_gsl<>), 
			     density_root(new root_cern<>) {
  deg_limit=2.0;
  
  exp_limit=200.0;
  upper_limit_fac=20.0;
  deg_entropy_fac=30.0;
  min_psi=-4.0;
  err_nonconv=true;
  use_expansions=true;
}

fermion_rel::~fermion_rel() {
}

void fermion_rel::calc_mu(fermion &f, double temper) {

  // -----------------------------------------------------------------
  // Handle T<=0

  if (temper<0.0) {
    O2SCL_ERR("Temperature less than zero in fermion_rel::calc_density().",
	      exc_einval);
  }
  if (temper==0.0) {
    calc_density_zerot(f);
    return;
  }

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

  // Compute the degeneracy parameter
  
  bool deg=true;
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/temper;
  } else {
    psi=(f.nu+(f.m-f.ms))/temper;
  }
  if (psi<deg_limit) deg=false;
  
  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    bool acc=calc_mu_ndeg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      return;
    }
  }

  // Try the degenerate expansion if psi is large enough
  if (use_expansions && psi>20.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      return;
    }
  }

  if (!deg) {
    
    // If the temperature is large enough, perform the full integral
    
    funct mfd=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::density_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::energy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::entropy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
      
    double prefac=f.g*pow(temper,3.0)/2.0/pi2;

    // Compute the number density
    
    f.n=nit->integ(mfd,0.0,0.0);
    f.n*=prefac;
    unc.n=nit->get_error()*prefac;

    // Compute the energy density

    f.ed=nit->integ(mfe,0.0,0.0);
    f.ed*=prefac*temper;
    if (!f.inc_rest_mass) f.ed-=f.n*f.m;
    unc.ed=nit->get_error()*prefac*temper;
    
    // Compute the entropy

    f.en=nit->integ(mfs,0.0,0.0);
    f.en*=prefac;
    unc.en=nit->get_error()*prefac;

  } else {
    
    // Otherwise, apply a degenerate approximation, by making the
    // upper integration limit finite
    
    funct mfd=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_density_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_energy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_entropy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);

    double prefac=f.g/2.0/pi2;
    
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
      O2SCL_ERR2("Zero density in degenerate limit in fermion_rel::",
		 "calc_mu(). Variable deg_limit set improperly?",
		 exc_efailed);
      return;
    }
    
    // Compute the number density

    f.n=dit->integ(mfd,0.0,ul);
    f.n*=prefac;
    unc.n=dit->get_error()*prefac;
    
    // Compute the energy density

    f.ed=dit->integ(mfe,0.0,ul);
    f.ed*=prefac;
    unc.ed=dit->get_error()*prefac;

    // Compute the lower limit for the entropy integration

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

    // Compute the entropy

    if (ll>0.0) {
      f.en=dit->integ(mfs,ll,ul);
    } else {
      f.en=dit->integ(mfs,0.0,ul);
    }
    f.en*=prefac;
    unc.en=dit->get_error()*prefac;
    
  }

  // Compute the pressure using the thermodynamic identity

  f.pr=-f.ed+temper*f.en+f.nu*f.n;
  unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
	      f.nu*unc.n*f.nu*unc.n);

  return;
}

int fermion_rel::nu_from_n(fermion &f, double temper) {
  double nex;

  // Try to ensure a good initial guess

  nex=f.nu/temper;
  double y=solve_fun(nex,f,temper);

  if (y>1.0-1.0e-6) {
    double scale=f.ms;
    if (temper>scale) scale=temper;
    for(size_t i=0;i<10;i++) {
      if (nex<0.0) nex+=scale*1.0e5;
      else nex*=10.0;
      y=solve_fun(nex,f,temper);
      if (y<1.0-1.0e-6) i=10;
    }
  }

  // If that didn't work, try a different guess
  if (y>1.0-1.0e-6) {
    if (f.inc_rest_mass) {
      nex=f.ms/temper;
    } else {
      nex=(f.ms-f.m)/temper;
    }
    y=solve_fun(nex,f,temper);
  }
  
  // If neither worked, call the error handler
  if (y==1.0 || !std::isfinite(y)) {
    O2SCL_CONV2_RET("Couldn't find reasonable initial guess in ",
		    "fermion_rel::nu_from_n().",exc_einval,this->err_nonconv);
  }

  // Perform full solution
  funct mf=std::bind(std::mem_fn<double(double,fermion &,double)>
		       (&fermion_rel::solve_fun),
		       this,std::placeholders::_1,std::ref(f),temper);

  bool drec=density_root->err_nonconv;
  density_root->err_nonconv=false;
  int ret=density_root->solve(nex,mf);

  if (ret!=0) {

    // If it fails, try to make the integrators more accurate
    double tol1=dit->tol_rel, tol2=dit->tol_abs;
    double tol3=nit->tol_rel, tol4=nit->tol_abs;
    dit->tol_rel/=1.0e2;
    dit->tol_abs/=1.0e2;
    nit->tol_rel/=1.0e2;
    nit->tol_abs/=1.0e2;
    ret=density_root->solve(nex,mf);

    // Return tolerances to their original values
    dit->tol_rel=tol1;
    dit->tol_abs=tol2;
    nit->tol_rel=tol3;
    nit->tol_abs=tol4;
  }

  density_root->err_nonconv=drec;

  if (ret!=0) {
    O2SCL_CONV2_RET("Density solver failed in ",
		    "fermion_rel::nu_from_n().",exc_efailed,this->err_nonconv);
  }
  f.nu=nex*temper;

  return success;
}

int fermion_rel::calc_density(fermion &f, double temper) {

  // The function pair_mu() can modify the density, which is
  // confusing to the user, so we return it to the user-specified
  // value.
  double density_temp=f.n;
  
  // -----------------------------------------------------------------
  // Handle T<=0

  if (temper<0.0) {
    O2SCL_ERR("Temperature less than zero in fermion_rel::calc_density().",
	      exc_einval);
  }
  if (temper==0.0) {
    calc_density_zerot(f);
    return 0;
  }

#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the density is not
  // finite, but I've found this extra checking of the inputs useful
  // for debugging.
  if (!std::isfinite(f.n)) {
    O2SCL_ERR2("Density not finite in ",
	       "fermion_rel::calc_density().",exc_einval);
  }
#endif

  // -----------------------------------------------------------------
  // First determine the chemical potential by solving for the density

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  int ret=nu_from_n(f,temper);
  if (ret!=0) {
    O2SCL_CONV2_RET("Function calc_density() failed in fermion_rel::",
		    "calc_density().",exc_efailed,this->err_nonconv);
  }

  if (f.non_interacting) { f.mu=f.nu; }

  // -----------------------------------------------------------------
  // Now use the chemical potential to compute the energy density,
  // pressure, and entropy

  bool deg=true;
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/temper;
  } else {
    psi=(f.nu+(f.m-f.ms))/temper;
  }
  if (psi<deg_limit) deg=false;

  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    bool acc=calc_mu_ndeg(f,temper,1.0e-14);
    if (acc) {
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      f.n=density_temp;
      return 0;
    }
  }
  
  // Try the degenerate expansion if psi is large enough
  if (use_expansions && psi>20.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      f.n=density_temp;
      return 0;
    }
  }

  if (!deg) {
    
    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::energy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::entropy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    
    f.ed=nit->integ(mfe,0.0,0.0);
    f.ed*=f.g*pow(temper,4.0)/2.0/pi2;
    if (!f.inc_rest_mass) f.ed-=f.n*f.m;
    unc.ed=nit->get_error()*f.g*pow(temper,4.0)/2.0/pi2;
    
    f.en=nit->integ(mfs,0.0,0.0);
    f.en*=f.g*pow(temper,3.0)/2.0/pi2;
    unc.en=nit->get_error()*f.g*pow(temper,3.0)/2.0/pi2;

  } else {

    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_energy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
    funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_entropy_fun),
			  this,std::placeholders::_1,std::ref(f),temper);
      
    double arg;
    if (f.inc_rest_mass) {
      arg=pow(upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
    } else {
      arg=pow(upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
    }
    double ul;
    if (arg>0.0) {
      
      ul=sqrt(arg);
      
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
      
      f.ed=dit->integ(mfe,0.0,ul);
      f.ed*=f.g/2.0/pi2;
      unc.ed=dit->get_error()*f.g/2.0/pi2;
      
      if (ll>0.0) {
	f.en=dit->integ(mfs,ll,ul);
      } else {
	f.en=dit->integ(mfs,0.0,ul);
      }
      f.en*=f.g/2.0/pi2;
      unc.en=dit->get_error()*f.g/2.0/pi2;
      
    } else {

      f.ed=0.0;
      f.en=0.0;
      unc.ed=0.0;
      unc.en=0.0;
      O2SCL_ERR2("Zero density in degenerate limit in fermion_rel::",
		 "calc_mu(). Variable deg_limit set improperly?",exc_efailed);
      
    }
  }

  f.n=density_temp;
  f.pr=-f.ed+temper*f.en+f.mu*f.n;
  unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
	      f.mu*unc.n*f.mu*unc.n);
  
  return 0;
}

double fermion_rel::deg_density_fun(double k, fermion &f, double T) {
  
  double E=gsl_hypot(k,f.ms), ret;
  if (!f.inc_rest_mass) E-=f.m;

  ret=k*k/(1.0+exp((E-f.nu)/T));

  if (!std::isfinite(ret)) {
    O2SCL_ERR2("Returned not finite result ",
	       "in fermion_rel::deg_density_fun().",exc_einval);
  }

  return ret;
}

double fermion_rel::deg_energy_fun(double k, fermion &f, double T) {

  double E=gsl_hypot(k,f.ms), ret;
  if (!f.inc_rest_mass) E-=f.m;

  ret=k*k*E/(1.0+exp((E-f.nu)/T));

  if (!std::isfinite(ret)) {
    O2SCL_ERR2("Returned not finite result ",
	       "in fermion_rel::deg_energy_fun().",exc_einval);
  }
  
  return ret;
}
  
double fermion_rel::deg_entropy_fun(double k, fermion &f, double T) {
  
  double E=gsl_hypot(k,f.ms), ret;
  if (!f.inc_rest_mass) E-=f.m;
  
  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-f.nu)/T)<-exp_limit) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-f.nu)/T)<-deg_entropy_fac) {
    double arg=E/T-f.nu/T;
    ret=-k*k*arg*exp(arg);
  } else {
    double nx=1.0/(1.0+exp(E/T-f.nu/T));
    ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx));
  }

  if (!std::isfinite(ret)) {
    O2SCL_ERR2("Returned not finite result ",
	       "in fermion_rel::deg_entropy_fun().",exc_einval);
  }

  return ret;
}
  
double fermion_rel::density_fun(double u, fermion &f, double T) {

  double ret, y, mx;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  mx=f.ms/T;
  
  if (y-mx-u>exp_limit) {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u);
  } else if (y>u+exp_limit && mx>u+exp_limit) {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u)/(exp(mx+u-y)+1.0);
  } else {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u)*exp(y)/(exp(mx+u)+exp(y));
  }

  if (!std::isfinite(ret)) {
    ret=0.0;
  }

  return ret;
}

double fermion_rel::energy_fun(double u, fermion &f, double T) {
  double ret, y, mx;

  mx=f.ms/T;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  if (y>u+exp_limit && mx>u+exp_limit) {
    ret=(mx+u)*(mx+u)*sqrt(u*u+2.0*mx*u)/(exp(mx+u-y)+1.0);
  } else {
    ret=(mx+u)*(mx+u)*sqrt(u*u+2.0*mx*u)*exp(y)/(exp(mx+u)+exp(y));
  }
 
  if (!std::isfinite(ret)) {
    return 0.0;
  }

  return ret;
}

double fermion_rel::entropy_fun(double u, fermion &f, double T) {
  double ret, y, mx, term1, term2;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  mx=f.ms/T;

  term1=log(1.0+exp(y-mx-u))/(1.0+exp(y-mx-u));
  term2=log(1.0+exp(mx+u-y))/(1.0+exp(mx+u-y));
  ret=(mx+u)*sqrt(u*u+2.0*mx*u)*(term1+term2);
  
  if (!std::isfinite(ret)) {
    return 0.0;
  }

  return ret;
}

double fermion_rel::solve_fun(double x, fermion &f, double T) {
  double nden, yy;
  
  f.nu=T*x;

  if (f.non_interacting) f.mu=f.nu;

  bool deg=true;
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/T;
  } else {
    psi=(f.nu+(f.m-f.ms))/T;
  }
  if (psi<deg_limit) deg=false;

  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    double ntemp=f.n;
    bool acc=calc_mu_ndeg(f,T,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      yy=(ntemp-f.n)/ntemp;
      f.n=ntemp;
      return yy;
    }
    f.n=ntemp;
  }

  // Try the degenerate expansion if psi is large enough
  if (use_expansions && psi>20.0) {
    double ntemp=f.n;
    bool acc=calc_mu_deg(f,T,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      yy=(ntemp-f.n)/ntemp;
      f.n=ntemp;
      return yy;
    }
    f.n=ntemp;
  }

  // Otherwise, directly perform the integration
  if (!deg) {

    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::density_fun),
			  this,std::placeholders::_1,std::ref(f),T);
    
    nden=nit->integ(mfe,0.0,0.0);
    nden*=f.g*pow(T,3.0)/2.0/pi2;
    unc.n=nit->get_error()*f.g*pow(T,3.0)/2.0/pi2;
    
    yy=(f.n-nden)/f.n;

  } else {
    
    funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			  (&fermion_rel::deg_density_fun),
			  this,std::placeholders::_1,std::ref(f),T);
    
    double arg;
    if (f.inc_rest_mass) {
      arg=pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms;
    } else {
      arg=pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms;
    }

    double ul;

    if (arg>0.0) {

      ul=sqrt(arg);
      
      nden=dit->integ(mfe,0.0,ul);
      nden*=f.g/2.0/pi2;
      unc.n=dit->get_error()*f.g/2.0/pi2;
      
    } else {

      nden=0.0;

    }

    yy=(f.n-nden)/f.n;
  }
  
  return yy;
}

void fermion_rel::pair_mu(fermion &f, double temper) {

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

  if (use_expansions) {
    if (calc_mu_ndeg(f,temper,1.0e-14,true)) {
      unc.n=1.0e-14*f.n;
      unc.ed=1.0e-14*f.ed;
      unc.en=1.0e-14*f.en;
      unc.pr=1.0e-14*f.pr;
      return;
    }
  }

  fermion antip(f.ms,f.g);
  f.anti(antip);

  // Particles
  calc_mu(f,temper);
  double unc_n=unc.n;
  double unc_pr=unc.pr;
  double unc_ed=unc.ed;
  double unc_en=unc.en;

  // Antiparticles
  calc_mu(antip,temper);

  // Add up thermodynamic quantities
  f.n-=antip.n;
  f.pr+=antip.pr;
  f.ed+=antip.ed;
  f.en+=antip.en;

  // Add up uncertainties
  unc.n=gsl_hypot(unc.n,unc_n);
  unc.ed=gsl_hypot(unc.ed,unc_ed);
  unc.pr=gsl_hypot(unc.pr,unc_pr);
  unc.en=gsl_hypot(unc.ed,unc_en);

  return;
}

int fermion_rel::pair_density(fermion &f, double temper) {

  if (f.n==0.0) {
    O2SCL_ERR("Zero density sent to fermion_rel::pair_density().",
	      exc_einval);
  }

  // -----------------------------------------------------------------
  // Handle T<=0

  if (temper<=0.0) {
    calc_density_zerot(f);
    return success;
  }

  double density_temp=f.n;
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  double nex=f.nu/temper;

  // -----------------------------------------------------------------
  // If the chemical potential is too small, then the solver
  // will fail
  
  // Find the larger of either the temperature or the mass
  double lg=temper;
  if (f.ms>lg) lg=temper;
  // Try increasing the chemical potential
  double y=pair_fun(nex,f,temper,false);
  for(size_t i=0;i<10 && fabs(y+1.0)<1.0e-6;i++) {
    nex+=lg/temper;
    y=pair_fun(nex,f,temper,false);
  }

  funct mf=std::bind(std::mem_fn<double(double,fermion &,double,bool)>
		       (&fermion_rel::pair_fun),
		       this,std::placeholders::_1,std::ref(f),temper,false);

  bool drec=density_root->err_nonconv;
  density_root->err_nonconv=false;
  int ret=density_root->solve(nex,mf);

  if (ret!=0) {

    // If it fails, try to make the integrators more accurate
    double tol1=dit->tol_rel, tol2=dit->tol_abs;
    double tol3=nit->tol_rel, tol4=nit->tol_abs;
    dit->tol_rel/=1.0e2;
    dit->tol_abs/=1.0e2;
    nit->tol_rel/=1.0e2;
    nit->tol_abs/=1.0e2;
    ret=density_root->solve(nex,mf);

    // Function in log units
    funct lmf=std::bind(std::mem_fn<double(double,fermion &,double,bool)>
			  (&fermion_rel::pair_fun),
			  this,std::placeholders::_1,std::ref(f),temper,true);
    
    if (ret!=0) {
      // If that failed, try working in log units
      nex=log(nex);
      ret=density_root->solve(nex,lmf);
      nex=exp(nex);
    }
    
    if (ret!=0) {
      // If that failed, try a different solver
      root_brent_gsl<> rbg;
      rbg.err_nonconv=false;
      nex=log(nex);
      ret=rbg.solve(nex,lmf);
      nex=exp(nex);
    }

    // Return tolerances to their original values
    dit->tol_rel=tol1;
    dit->tol_abs=tol2;
    nit->tol_rel=tol3;
    nit->tol_abs=tol4;
  }

  density_root->err_nonconv=drec;

  if (ret!=0) {
    O2SCL_CONV2_RET("Density solver failed in fermion_rel::",
		    "pair_density().",exc_efailed,this->err_nonconv);
  }

  f.nu=nex*temper;
  
  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  // The function pair_mu() can modify the density, which is
  // confusing to the user, so we return it to the user-specified
  // value.
  f.n=density_temp;

  return success;
}

double fermion_rel::pair_fun(double x, fermion &f, double T, bool log_mode) { 

  // Temporary storage for density to match
  double nn_match=f.n;
  // Temporary storage for integration results
  double nden;
  // The return value, n/(n') -1
  double yy;

  // -----------------------------------------------------------------

  f.nu=T*x;
  if (log_mode) f.nu=T*exp(x);

  //if (!std::isfinite(f.nu)) {
  //O2SCL_ERR("Chemical potential not finite in fermion_rel::pair_fun().",
  //exc_einval);
  //}

  if (f.non_interacting) f.mu=f.nu;

  // -----------------------------------------------------------------
  // First, try the non-degenerate expansion with both particles and
  // antiparticles together
  
  if (use_expansions) {
    if (calc_mu_ndeg(f,T,1.0e-8,true) && std::isfinite(f.n)) {
      yy=f.n;
      f.n=nn_match;
      yy=yy/nn_match-1.0;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (10) in fermion_rel::pair_fun().",
		  exc_einval);
      }
      return yy;
    }
  }

  // -----------------------------------------------------------------
  // If that doesn't work, evaluate particles and antiparticles 
  // separately. This is the contribution for particles

  bool deg=true;
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/T;
  } else {
    psi=(f.nu+(f.m-f.ms))/T;
  }
  if (psi<deg_limit) deg=false;

  bool particles_done=false;

  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    if (calc_mu_ndeg(f,T,1.0e-8) && std::isfinite(f.n)) {
      particles_done=true;
      yy=f.n;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (1) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }
  
  // Try the degenerate expansion if psi is large enough
  if (use_expansions && particles_done==false && psi>20.0) {
    if (calc_mu_deg(f,T,1.0e-8) && std::isfinite(f.n)) {
      particles_done=true;
      yy=f.n;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (2) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }

  // If neither expansion worked, use direct integration
  if (particles_done==false) {
    
    if (!deg) {
      
      // Nondegenerate case
      
      funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::density_fun),
			    this,std::placeholders::_1,std::ref(f),T);
      
      nden=nit->integ(mfe,0.0,0.0);
      nden*=f.g*pow(T,3.0)/2.0/pi2;
      yy=nden;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (3) in fermion_rel::pair_fun().",
		  exc_einval);
      }
      
    } else {
      
      // Degenerate case
      
      funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_density_fun),
			    this,std::placeholders::_1,std::ref(f),T);
      
      double arg;
      if (f.inc_rest_mass) {
	arg=pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms;
      } else {
	arg=pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms;
      }
      
      double ul;
      if (arg>0.0) {
	ul=sqrt(arg);
	nden=dit->integ(mfe,0.0,ul);
	nden*=f.g/2.0/pi2;
      } else {
	nden=0.0;
      }
      
      yy=nden;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (4) in fermion_rel::pair_fun().",
		  exc_einval);
      }

    }

    particles_done=true;

    // End of 'if (particles_done==false)'
  }

  // -----------------------------------------------------------------
  // Compute the contribution from the antiparticles

  f.nu=-T*x;
  if (log_mode) f.nu=-T*exp(x);
  if (f.non_interacting) f.mu=f.nu;

  bool antiparticles_done=false;

  // Evaluate the degeneracy parameter
  deg=true;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/T;
  } else {
    psi=(f.nu+(f.m-f.ms))/T;
  }
  if (psi<deg_limit) deg=false;
  
  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    if (calc_mu_ndeg(f,T,1.0e-8)) {
      antiparticles_done=true;
      yy-=f.n;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (5) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }

  // Try the degenerate expansion if psi is large enough
  if (use_expansions && antiparticles_done==false && psi>20.0) {
    if (calc_mu_deg(f,T,1.0e-8)) {
      antiparticles_done=true;
      yy-=f.n;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (6) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }

  // If neither expansion worked, use direct integration
  if (antiparticles_done==false) {
    
    if (!deg) {
      
      // Nondegenerate case
      
      funct mf=std::bind(std::mem_fn<double(double,fermion &,double)>
			   (&fermion_rel::density_fun),
			   this,std::placeholders::_1,std::ref(f),T);
      
      nden=nit->integ(mf,0.0,0.0);
      nden*=f.g*pow(T,3.0)/2.0/pi2;
      yy-=nden;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (7) in fermion_rel::pair_fun().",
		  exc_einval);
      }
      
    } else {
      
      // Degenerate case
      
      funct mf=std::bind(std::mem_fn<double(double,fermion &,double)>
			   (&fermion_rel::deg_density_fun),
			   this,std::placeholders::_1,std::ref(f),T);
      
      double arg;
      if (f.inc_rest_mass) {
	arg=pow(upper_limit_fac*T+f.nu,2.0)-f.ms*f.ms;
      } else {
	arg=pow(upper_limit_fac*T+f.nu+f.m,2.0)-f.ms*f.ms;
      }
      
      double ul;
      if (arg>0.0) {
	ul=sqrt(arg);
	nden=dit->integ(mf,0.0,ul);
	nden*=f.g/2.0/pi2;
      } else {
	nden=0.0;
      }
      yy-=nden;
      if (!std::isfinite(yy)) {
	O2SCL_ERR("Value 'yy' not finite (8) in fermion_rel::pair_fun().",
		  exc_einval);
      }

    }

    antiparticles_done=true;
  }

  // Finish computing the function value
  f.n=nn_match;
  yy=yy/nn_match-1.0;

  if (!std::isfinite(yy)) {
    O2SCL_ERR("Value 'yy' not finite (9) in fermion_rel::pair_fun().",
	      exc_einval);
  }

  return yy;
}

