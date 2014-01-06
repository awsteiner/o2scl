/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#include <o2scl/rel_fermion.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

rel_fermion::rel_fermion() : nit(new inte_qagiu_gsl<>), 
			     dit(new inte_qag_gsl<>), 
			     density_root(new root_cern<>) {
  deg_limit=2.0;
  
  exp_limit=200.0;
  upper_limit_fac=20.0;
  deg_entropy_fac=30.0;
  min_psi=-4.0;
  err_nonconv=true;
}

rel_fermion::~rel_fermion() {
}

void rel_fermion::calc_mu(fermion &f, double temper) {

  if (temper<=0.0) {
    calc_mu_zerot(f);
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
  if (psi<min_psi) {
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
  if (psi>20.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-14);
    if (acc) {
      unc.n=f.n*1.0e-14;
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      return;
    }
  }

  T=temper;
  fp=&f;
  
  if (!deg) {
    
    // If the temperature is large enough, perform the full integral
    
    funct_mfptr<rel_fermion> mfd(this,&rel_fermion::density_fun);
    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::energy_fun);
    funct_mfptr<rel_fermion> mfs(this,&rel_fermion::entropy_fun);

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
    
    funct_mfptr<rel_fermion> mfd(this,&rel_fermion::deg_density_fun);
    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::deg_energy_fun);
    funct_mfptr<rel_fermion> mfs(this,&rel_fermion::deg_entropy_fun);

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
      O2SCL_ERR2("Zero density in degenerate limit in rel_fermion::",
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

int rel_fermion::nu_from_n(fermion &f, double temper) {
  double nex;

  T=temper;
  fp=&f;

  // Check that initial guess is close enough

  nex=f.nu/temper;
  double y=solve_fun(nex);

  // If that didn't work, try a different guess
  if (y==1.0) {
    if (f.inc_rest_mass) {
      nex=f.ms/temper;
    } else {
      nex=(f.ms-f.m)/temper;
    }
  }

  y=solve_fun(nex);

  // If neither worked, call the error handler
  if (y==1.0 || !o2scl::is_finite(y)) {
    O2SCL_CONV2_RET("Couldn't find reasonable initial guess in ",
		    "rel_fermion::nu_from_n().",exc_einval,this->err_nonconv);
  }
  
  // Perform full solution
  funct_mfptr<rel_fermion> mf(this,&rel_fermion::solve_fun);
  int ret=density_root->solve(nex,mf);
  if (ret!=0) {
    O2SCL_CONV2_RET("Density solver failed in ",
		    "rel_fermion::nu_from_n().",exc_efailed,this->err_nonconv);
  }
  f.nu=nex*temper;
  
  return success;
}

void rel_fermion::calc_density(fermion &f, double temper) {

#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the density is not
  // finite, but I've found this extra checking of the inputs useful
  // for debugging.
  if (!o2scl::is_finite(f.n)) {
    O2SCL_ERR2("Density not finite in ",
	       "rel_fermion::calc_density().",exc_einval);
  }
#endif

  fp=&f;
  T=temper;

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  nu_from_n(f,temper);

  if (f.non_interacting) { f.mu=f.nu; }

  bool deg=true;
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/temper;
  } else {
    psi=(f.nu+(f.m-f.ms))/temper;
  }
  if (psi<deg_limit) deg=false;

  if (psi<min_psi) {
    bool acc=calc_mu_ndeg(f,temper,1.0e-14);
    if (acc) {
      unc.ed=f.ed*1.0e-14;
      unc.pr=f.pr*1.0e-14;
      unc.en=f.en*1.0e-14;
      return;
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
      return;
    }
  }
  if (!deg) {
    
    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::energy_fun);
    funct_mfptr<rel_fermion> mfs(this,&rel_fermion::entropy_fun);
    
    f.ed=nit->integ(mfe,0.0,0.0);
    f.ed*=f.g*pow(temper,4.0)/2.0/pi2;
    if (!f.inc_rest_mass) f.ed-=f.n*f.m;
    unc.ed=nit->get_error()*f.g*pow(temper,4.0)/2.0/pi2;
    
    f.en=nit->integ(mfs,0.0,0.0);
    f.en*=f.g*pow(temper,3.0)/2.0/pi2;
    unc.en=nit->get_error()*f.g*pow(temper,3.0)/2.0/pi2;

  } else {

    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::deg_energy_fun);
    funct_mfptr<rel_fermion> mfs(this,&rel_fermion::deg_entropy_fun);
      
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
      O2SCL_ERR2("Zero density in degenerate limit in rel_fermion::",
		 "calc_mu(). Variable deg_limit set improperly?",exc_efailed);
      
    }
  }

  f.pr=-f.ed+temper*f.en+f.mu*f.n;
  unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
	      f.mu*unc.n*f.mu*unc.n);
  
  return;
}

double rel_fermion::deg_density_fun(double k) {

  double E=gsl_hypot(k,fp->ms), ret;
  if (!fp->inc_rest_mass) E-=fp->m;
  
  ret=k*k/(1.0+exp((E-fp->nu)/T));
  
  return ret;
}

double rel_fermion::deg_energy_fun(double k) {

  double E=gsl_hypot(k,fp->ms), ret;
  if (!fp->inc_rest_mass) E-=fp->m;

  ret=k*k*E/(1.0+exp((E-fp->nu)/T));
  
  return ret;
}
  
double rel_fermion::deg_entropy_fun(double k) {
  
  double E=gsl_hypot(k,fp->ms), ret;
  if (!fp->inc_rest_mass) E-=fp->m;
  
  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-fp->nu)/T)<-exp_limit) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-fp->nu)/T)<-deg_entropy_fac) {
    double arg=E/T-fp->nu/T;
    ret=-k*k*arg*exp(arg);
  } else {
    double nx=1.0/(1.0+exp(E/T-fp->nu/T));
    ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx));
  }

  return ret;
}
  
double rel_fermion::density_fun(double u) {
  double ret, y, mx;

  if (fp->inc_rest_mass) {
    y=fp->nu/T;
  } else {
    y=(fp->nu+fp->m)/T;
  }
  mx=fp->ms/T;
  
  if (y-mx-u>exp_limit) {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u);
  } else if (y>u+exp_limit && mx>u+exp_limit) {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u)/(exp(mx+u-y)+1.0);
  } else {
    ret=(mx+u)*sqrt(u*u+2.0*mx*u)*exp(y)/(exp(mx+u)+exp(y));
  }

  if (!o2scl::is_finite(ret)) {
    ret=0.0;
  }

  return ret;
}

double rel_fermion::energy_fun(double u) {
  double ret, y, mx;

  mx=fp->ms/T;

  if (fp->inc_rest_mass) {
    y=fp->nu/T;
  } else {
    y=(fp->nu+fp->m)/T;
  }
  if (y>u+exp_limit && mx>u+exp_limit) {
    ret=(mx+u)*(mx+u)*sqrt(u*u+2.0*mx*u)/(exp(mx+u-y)+1.0);
  } else {
    ret=(mx+u)*(mx+u)*sqrt(u*u+2.0*mx*u)*exp(y)/(exp(mx+u)+exp(y));
  }
 
  if (!o2scl::is_finite(ret)) {
    return 0.0;
  }

  return ret;
}

double rel_fermion::entropy_fun(double u) {
  double ret, y, mx, term1, term2;

  if (fp->inc_rest_mass) {
    y=fp->nu/T;
  } else {
    y=(fp->nu+fp->m)/T;
  }
  mx=fp->ms/T;

  term1=log(1.0+exp(y-mx-u))/(1.0+exp(y-mx-u));
  term2=log(1.0+exp(mx+u-y))/(1.0+exp(mx+u-y));
  ret=(mx+u)*sqrt(u*u+2.0*mx*u)*(term1+term2);
  
  if (!o2scl::is_finite(ret)) {
    return 0.0;
  }

  return ret;
}

double rel_fermion::solve_fun(double x) {
  double nden, yy;
  
  fp->nu=T*x;

  if (fp->non_interacting) fp->mu=fp->nu;

  bool deg=true;
  double psi;
  if (fp->inc_rest_mass) {
    psi=(fp->nu-fp->ms)/T;
  } else {
    psi=(fp->nu+(fp->m-fp->ms))/T;
  }
  if (psi<deg_limit) deg=false;
  
  if (psi<min_psi) {
    double ntemp=fp->n;
    bool acc=calc_mu_ndeg(*fp,T,1.0e-14);
    if (acc) {
      unc.n=fp->n*1.0e-14;
      yy=(ntemp-fp->n)/ntemp;
      fp->n=ntemp;
      return yy;
    }
    fp->n=ntemp;
  }

  // Try the degenerate expansion if psi is large enough
  if (psi>20.0) {
    double ntemp=fp->n;
    bool acc=calc_mu_deg(*fp,T,1.0e-14);
    if (acc) {
      unc.n=fp->n*1.0e-14;
      yy=(ntemp-fp->n)/ntemp;
      fp->n=ntemp;
      return yy;
    }
    fp->n=ntemp;
  }


  if (!deg) {

    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::density_fun);
    
    nden=nit->integ(mfe,0.0,0.0);
    nden*=fp->g*pow(T,3.0)/2.0/pi2;
    unc.n=nit->get_error()*fp->g*pow(T,3.0)/2.0/pi2;
    
    yy=(fp->n-nden)/fp->n;

  } else {
    
    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::deg_density_fun);
    
    double arg;
    if (fp->inc_rest_mass) {
      arg=pow(upper_limit_fac*T+fp->nu,2.0)-fp->ms*fp->ms;
    } else {
      arg=pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms;
    }

    double ul;

    if (arg>0.0) {

      ul=sqrt(arg);
      
      nden=dit->integ(mfe,0.0,ul);
      nden*=fp->g/2.0/pi2;
      unc.n=err_hnd->get_errno()*fp->g/2.0/pi2;
      
    } else {

      nden=0.0;

    }

    yy=(fp->n-nden)/fp->n;
  }
  
  return yy;
}

void rel_fermion::pair_mu(fermion &f, double temper) {

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

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
  unc.n=sqrt(unc.n*unc.n+unc_n*unc_n);
  unc.ed=sqrt(unc.ed*unc.ed+unc_ed*unc_ed);
  unc.pr=sqrt(unc.pr*unc.pr+unc_pr*unc_pr);
  unc.en=sqrt(unc.en*unc.en+unc_en*unc_en);

  return;
}

void rel_fermion::pair_density(fermion &f, double temper) {
  double nex;

  T=temper;
  fp=&f;

  if (temper<=0.0) {
    calc_density_zerot(f);
    return;
  }
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  nex=f.nu/temper;
  funct_mfptr<rel_fermion> mf(this,&rel_fermion::pair_fun);
  density_root->solve(nex,mf);
  f.nu=nex*temper;
  
  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  return;
}

double rel_fermion::pair_fun(double x) { 

  double nden, yy;

  // -------------------------------------------------------------
  // Compute the contribution from the particles

  fp->nu=T*x;

  if (fp->non_interacting) fp->mu=fp->nu;

  // Evaluate the degeneracy parameter
  bool deg=true;
  if (fp->inc_rest_mass) {
    if ((fp->nu-fp->ms)/T<deg_limit) deg=false;
  } else {
    if ((fp->nu+(fp->m-fp->ms))/T<deg_limit) deg=false;
  }

  if (!deg) {
    
    // Nondegenerate case

    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::density_fun);

    nden=nit->integ(mfe,0.0,0.0);
    nden*=fp->g*pow(T,3.0)/2.0/pi2;
    yy=nden;

  } else {

    // Degenerate case

    funct_mfptr<rel_fermion> mfe(this,&rel_fermion::deg_density_fun);

    double arg;
    if (fp->inc_rest_mass) {
      arg=pow(upper_limit_fac*T+fp->nu,2.0)-fp->ms*fp->ms;
    } else {
      arg=pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms;
    }

    double ul;
    if (arg>0.0) {
      ul=sqrt(arg);
      nden=dit->integ(mfe,0.0,ul);
      nden*=fp->g/2.0/pi2;
    } else {
      nden=0.0;
    }

    yy=nden;

  }

  // -------------------------------------------------------------
  // Compute the contribution from the antiparticles

  fp->nu=-T*x;

  // Evaluate the degeneracy parameter
  deg=true;
  if (fp->inc_rest_mass) {
    if ((fp->nu-fp->ms)/T<deg_limit) deg=false;
  } else {
    if ((fp->nu+(fp->m-fp->ms))/T<deg_limit) deg=false;
  }
  
  if (!deg) {
    
    // Nondegenerate case

    funct_mfptr<rel_fermion> mf(this,&rel_fermion::density_fun);
    
    nden=nit->integ(mf,0.0,0.0);
    nden*=fp->g*pow(T,3.0)/2.0/pi2;
    yy-=nden;

  } else {

    // Degenerate case

    funct_mfptr<rel_fermion> mf(this,&rel_fermion::deg_density_fun);
    
    double arg;
    if (fp->inc_rest_mass) {
      arg=pow(upper_limit_fac*T+fp->nu,2.0)-fp->ms*fp->ms;
    } else {
      arg=pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms;
    }

    double ul;
    if (arg>0.0) {
      ul=sqrt(arg);
      nden=dit->integ(mf,0.0,ul);
      nden*=fp->g/2.0/pi2;
    } else {
      nden=0.0;
    }
    yy-=nden;

  }

  // Construct the function value
  yy=yy/fp->n-1.0;

  return yy;
}

