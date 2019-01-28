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
  density_root->tol_rel=4.0e-7;
  verbose=0;
}

fermion_rel::~fermion_rel() {
}

void fermion_rel::calc_mu(fermion &f, double temper) {
  calc_mu_tlate<fermion>(f,temper);
  return;
}

int fermion_rel::nu_from_n(fermion &f, double temper) {
  return nu_from_n_tlate<fermion>(f,temper);
}

int fermion_rel::calc_density(fermion &f, double temper) {
  return calc_density_tlate<fermion>(f,temper);
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
    ret=-k*k*(-1.0+arg)*exp(arg);
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

  double ret, y, eta;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  eta=f.ms/T;
  
  if (y-eta-u>exp_limit) {
    ret=(eta+u)*sqrt(u*u+2.0*eta*u);
  } else if (y>u+exp_limit && eta>u+exp_limit) {
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)+1.0);
  } else {
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)+exp(y));
  }

  if (!std::isfinite(ret)) {
    ret=0.0;
  }

  return ret;
}

double fermion_rel::energy_fun(double u, fermion &f, double T) {
  double ret, y, eta;

  eta=f.ms/T;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  if (y>u+exp_limit && eta>u+exp_limit) {
    ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)+1.0);
  } else {
    ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)+exp(y));
  }
 
  if (!std::isfinite(ret)) {
    return 0.0;
  }

  return ret;
}

double fermion_rel::entropy_fun(double u, fermion &f, double T) {
  double ret, y, eta, term1, term2;

  if (f.inc_rest_mass) {
    y=f.nu/T;
  } else {
    y=(f.nu+f.m)/T;
  }
  eta=f.ms/T;

  term1=log(1.0+exp(y-eta-u))/(1.0+exp(y-eta-u));
  term2=log(1.0+exp(eta+u-y))/(1.0+exp(eta+u-y));
  ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2);
  
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
  pair_mu_tlate<fermion>(f,temper);
  return;
}

int fermion_rel::pair_density(fermion &f, double temper) {
  return pair_density_tlate<fermion>(f,temper);
}

double fermion_rel::pair_fun(double x, fermion &f, double T, bool log_mode) { 

  // Temporary storage for density to match
  double nn_match=f.n;
  // Number density of particles and antiparticles
  double nden_p, nden_ap;

  // -----------------------------------------------------------------

  f.nu=T*x;
  if (log_mode) {
    f.nu=T*exp(x);
  }

  // Sometimes the exp() call above causes an overflow, so
  // we avoid extreme values
  if (!std::isfinite(f.nu)) return 3;

  if (f.non_interacting) f.mu=f.nu;

  // -----------------------------------------------------------------
  // First, try the non-degenerate expansion with both particles and
  // antiparticles together

  // AWS: 7/25/18: I'm taking this section out because it doesn't seem
  // to make sense to me, it apparently uses calc_mu_ndeg() which is
  // for particles only, and I'm not sure that's sufficient here. This
  // section also caused problems for the n=0, T!=0 case.

  if (false && use_expansions) {
    if (calc_mu_ndeg(f,T,1.0e-8,true) && std::isfinite(f.n)) {
      double y1=f.n/nn_match-1.0;
      if (!std::isfinite(y1)) {
	O2SCL_ERR("Value 'y1' not finite (10) in fermion_rel::pair_fun().",
		  exc_einval);
      }
      // Make sure to restore the value of f.n to it's original value,
      // nn_match
      f.n=nn_match;
      return y1;
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
      nden_p=f.n;
      if (!std::isfinite(nden_p)) {
	O2SCL_ERR("Value 'nden_p' not finite (1) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }
  
  // Try the degenerate expansion if psi is large enough
  if (use_expansions && particles_done==false && psi>20.0) {
    if (calc_mu_deg(f,T,1.0e-8) && std::isfinite(f.n)) {
      particles_done=true;
      nden_p=f.n;
      if (!std::isfinite(nden_p)) {
	O2SCL_ERR("Value 'nden_p' not finite (2) in fermion_rel::pair_fun().",
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
      
      nden_p=nit->integ(mfe,0.0,0.0);
      nden_p*=f.g*pow(T,3.0)/2.0/pi2;
      if (!std::isfinite(nden_p)) {
	O2SCL_ERR("Value 'nden_p' not finite (3) in fermion_rel::pair_fun().",
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
	nden_p=dit->integ(mfe,0.0,ul);
	nden_p*=f.g/2.0/pi2;
      } else {
	nden_p=0.0;
      }
      
      if (!std::isfinite(nden_p)) {
	O2SCL_ERR("Value 'nden_p' not finite (4) in fermion_rel::pair_fun().",
		  exc_einval);
      }

    }

    particles_done=true;

    // End of 'if (particles_done==false)'
  }

  // -----------------------------------------------------------------
  // Compute the contribution from the antiparticles

  if (f.inc_rest_mass) {
    f.nu=-T*x;
    if (log_mode) f.nu=-T*exp(x);
  } else {
    f.nu=-T*x-2.0*f.m;
    if (log_mode) f.nu=-T*exp(x)-2.0*f.m;
  }
  if (f.non_interacting) f.mu=f.nu;

  bool antiparticles_done=false;

  // Evaluate the degeneracy parameter
  deg=true;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/T;
  } else {
    psi=(f.nu+f.m-f.ms)/T;
  }
  if (psi<deg_limit) deg=false;

  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions && psi<min_psi) {
    if (calc_mu_ndeg(f,T,1.0e-8)) {
      antiparticles_done=true;
      nden_ap=f.n;
      if (!std::isfinite(nden_ap)) {
	O2SCL_ERR("Value 'nden_ap' not finite (5) in fermion_rel::pair_fun().",
		  exc_einval);
      }
    }
  }

  // Try the degenerate expansion if psi is large enough
  if (use_expansions && antiparticles_done==false && psi>20.0) {
    if (calc_mu_deg(f,T,1.0e-8)) {
      antiparticles_done=true;
      nden_ap=f.n;
      if (!std::isfinite(nden_ap)) {
	O2SCL_ERR("Value 'nden_ap' not finite (6) in fermion_rel::pair_fun().",
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
      
      nden_ap=nit->integ(mf,0.0,0.0);
      nden_ap*=f.g*pow(T,3.0)/2.0/pi2;
      if (!std::isfinite(nden_ap)) {
	O2SCL_ERR("Value 'nden_ap' not finite (7) in fermion_rel::pair_fun().",
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
	nden_ap=dit->integ(mf,0.0,ul);
	nden_ap*=f.g/2.0/pi2;
      } else {
	nden_ap=0.0;
      }
      if (!std::isfinite(nden_ap)) {
	O2SCL_ERR("Value 'nden_ap' not finite (8) in fermion_rel::pair_fun().",
		  exc_einval);
      }

    }

    antiparticles_done=true;
  }

  double y2;
  // Finish computing the function value
  if (nn_match==0.0) {
    y2=fabs(nden_p-nden_ap)/fabs(nden_p);
  } else {
    y2=(nden_p-nden_ap)/nn_match-1.0;
  }

  if (!std::isfinite(y2)) {
    O2SCL_ERR("Value 'y2' not finite (9) in fermion_rel::pair_fun().",
	      exc_einval);
  }
  
  // Make sure to restore the value of f.n to it's original value,
  // nn_match
  f.n=nn_match;
  return y2;
}

