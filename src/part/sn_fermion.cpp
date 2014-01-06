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

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/sn_fermion.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

sn_fermion::sn_fermion() {
  
  deg_limit=2.0;
  upper_limit_fac=20.0;

  density_root=&def_density_root;
  nit=&def_nit;
  dit=&def_dit;
  
  method=byparts;

  exp_limit=200.0;
}

sn_fermion::~sn_fermion() {
}

void sn_fermion::set_inte(inte<funct> &l_nit, inte<funct> &l_dit) {
  nit=&l_nit;
  dit=&l_dit;
  return;
}

void sn_fermion::calc_mu(fermion_deriv &f, double temper) {
  int ret=success, iret;
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in sn_fermion().",exc_eunimpl);
  }

  T=temper;
  fp=&f;

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

  double prefac=f.g/2.0/pi2;

  // Compute the degeneracy parameter

  bool deg=false;
  if (f.inc_rest_mass && (f.nu-f.ms)/temper>deg_limit) deg=true;
  else if (!f.inc_rest_mass && (f.nu+f.m-f.ms)/temper>deg_limit) deg=true;

  if (deg==false) {
    
    // The non-degenerate case

    funct_mfptr<sn_fermion> density_fun_f(this,&sn_fermion::density_fun);
    iret=nit->integ_err(density_fun_f,0.0,0.0,f.n,unc.n);
    if (iret!=0) {
      O2SCL_ERR2("Density integration (ndeg) failed in ",
		 "sn_fermion::calc_mu().",exc_efailed);
      ret=exc_efailed;
    }
    f.n*=prefac;
    unc.n*=prefac;
    
    funct_mfptr<sn_fermion> 
      density_T_fun_f(this,&sn_fermion::density_T_fun);
    iret=nit->integ_err(density_T_fun_f,0.0,0.0,f.dndT,unc.dndT);
    if (iret!=0) {
      O2SCL_ERR("dndT integration (ndeg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndT*=prefac;
    unc.dndT*=prefac;

    funct_mfptr<sn_fermion> 
      density_mu_fun_f(this,&sn_fermion::density_mu_fun);
    iret=nit->integ_err(density_mu_fun_f,0.0,0.0,f.dndmu,unc.dndmu);
    if (iret!=0) {
      O2SCL_ERR("dndmu integration (ndeg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndmu*=prefac;
    unc.dndmu*=prefac;
    
    funct_mfptr<sn_fermion> 
      energy_fun_f(this,&sn_fermion::energy_fun);
    iret=nit->integ_err(energy_fun_f,0.0,0.0,f.ed,unc.ed);
    if (iret!=0) {
      O2SCL_ERR2("Energy integration (ndeg) failed in ",
		 "sn_fermion::calc_mu().",exc_efailed);
      ret=exc_efailed;
    }
    f.ed*=prefac;
    f.ed*=pow(temper,4.0);
    unc.ed*=prefac;
    
    funct_mfptr<sn_fermion> 
      entropy_fun_f(this,&sn_fermion::entropy_fun);
    iret=nit->integ_err(entropy_fun_f,0.0,0.0,f.en,unc.en);
    if (iret!=0) {
      O2SCL_ERR2("Entropy integration (ndeg) failed in ",
		 "sn_fermion::calc_mu().",exc_efailed);
      ret=exc_efailed;
    }
    f.en*=prefac;
    unc.en*=prefac;
    
    funct_mfptr<sn_fermion> 
      entropy_T_fun_f(this,&sn_fermion::entropy_T_fun);
    iret=nit->integ_err(entropy_T_fun_f,0.0,0.0,f.dsdT,unc.dsdT);
    if (iret!=0) {
      O2SCL_ERR("dsdT integration (ndeg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dsdT*=prefac;
    unc.dsdT*=prefac;

    funct_mfptr<sn_fermion> 
      density_ms_fun_f(this,&sn_fermion::density_ms_fun);
    iret=nit->integ_err(density_ms_fun_f,0.0,0.0,f.dndm,unc.dndm);
    if (iret!=0) {
      O2SCL_ERR("dndm integration (ndeg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndm*=prefac;
    unc.dndm*=prefac;

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
      O2SCL_ERR2("Zero density in degenerate limit in sn_fermion::",
		 "calc_mu(). Variable deg_limit set improperly?",
		 exc_efailed);
      return;
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

    int old_meth=method;
    if ((f.nu+f.m-f.ms)/temper>1.0e3) {
      method=direct;
    }

    funct_mfptr<sn_fermion> 
      deg_density_fun_f(this,&sn_fermion::deg_density_fun);
    iret=dit->integ_err(deg_density_fun_f,0.0,ul,f.n,unc.n);
    if (iret!=0) {
      O2SCL_ERR2("Density integration (deg) failed in ",
		 "sn_fermion::calc_mu().",exc_efailed);
      ret=exc_efailed;
    }
    f.n*=prefac;
    unc.n*=prefac;
    
    funct_mfptr<sn_fermion> 
      deg_density_mu_fun_f(this,&sn_fermion::deg_density_mu_fun);
    if (method==direct && ll>0.0) {
      iret=dit->integ_err(deg_density_mu_fun_f,ll,ul,
			  f.dndmu,unc.dndmu);
    } else {
      iret=dit->integ_err(deg_density_mu_fun_f,0.0,ul,
			  f.dndmu,unc.dndmu);
    }
    if (iret!=0) {
      O2SCL_ERR("dndmu integration (deg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndmu*=prefac;
    unc.dndmu*=prefac;
    
    funct_mfptr<sn_fermion> 
      deg_density_T_fun_f(this,&sn_fermion::deg_density_T_fun);
    if (method==direct && ll>0.0) {
      iret=dit->integ_err(deg_density_T_fun_f,ll,ul,f.dndT,unc.dndT);
    } else {
      iret=dit->integ_err(deg_density_T_fun_f,0.0,ul,f.dndT,unc.dndT);
    }
    if (iret!=0) {
      O2SCL_ERR("dndT integration (deg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndT*=prefac;
    unc.dndT*=prefac;

    funct_mfptr<sn_fermion> 
      deg_energy_fun_f(this,&sn_fermion::deg_energy_fun);
    iret=dit->integ_err(deg_energy_fun_f,0.0,ul,f.ed,unc.ed);
    if (iret!=0) {
      O2SCL_ERR("Energy integration (deg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.ed*=prefac;
    unc.ed*=prefac;

    funct_mfptr<sn_fermion> 
      deg_entropy_fun_f(this,&sn_fermion::deg_entropy_fun);
    if (ll>0.0) {
      iret=dit->integ_err(deg_entropy_fun_f,ll,ul,f.en,unc.en);
    } else {
      iret=dit->integ_err(deg_entropy_fun_f,0.0,ul,f.en,unc.en);
    }
    if (iret!=0) {
      O2SCL_ERR2("Entropy integration (deg) failed in ",
		 "sn_fermion::calc_mu().",exc_efailed);
      ret=exc_efailed;
    }
    f.en*=prefac;
    unc.en*=prefac;
    
    funct_mfptr<sn_fermion> 
      deg_entropy_T_fun_f(this,&sn_fermion::deg_entropy_T_fun);
    if (method==direct && ll>0.0) {
      iret=dit->integ_err(deg_entropy_T_fun_f,ll,ul,f.dsdT,unc.dsdT);
    } else {
      iret=dit->integ_err(deg_entropy_T_fun_f,0.0,ul,f.dsdT,unc.dsdT);
    }
    if (iret!=0) {
      O2SCL_ERR("dsdT integration (deg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dsdT*=prefac;
    unc.dsdT*=prefac;

    funct_mfptr<sn_fermion> 
      deg_density_ms_fun_f(this,&sn_fermion::deg_density_ms_fun);
    if (method==direct && ll>0.0) {
      iret=dit->integ_err(deg_density_ms_fun_f,ll,ul,f.dndm,unc.dndm);
    } else {
      iret=dit->integ_err(deg_density_ms_fun_f,0.0,ul,f.dndm,unc.dndm);
    }
    if (iret!=0) {
      O2SCL_ERR("dndm integration (deg) failed in sn_fermion::calc_mu().",
		exc_efailed);
      ret=exc_efailed;
    }
    f.dndm*=prefac;
    unc.dndm*=prefac;

    if ((f.nu+f.m-f.ms)/temper>1.0e3) {
      method=old_meth;
    }
    
    
  }
  
  if (!o2scl::is_finite(f.en)) {
    O2SCL_ERR("Entropy not finite in sn_fermion::calc_mu().",exc_efailed);
    ret=exc_efailed;
  }
  f.pr=-f.ed+temper*f.en+f.nu*f.n;
  
  return;
}

void sn_fermion::nu_from_n(fermion_deriv &f, double temper) {
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
    y=solve_fun(nex);
  }

  if (y==1.0) {
    /*
      This section, which might be unnecessary, is to resolve a
      strange heisenbug. The 'cout' statement below seems to help. If
      someone gets an exception that throws here please let me know.
    */
    for(double nex2=1.0e-6;nex2<1.0e6;nex2*=100.0) {
      y=solve_fun(nex2);
      cout << f.ms << " " << f.m << " " << temper << endl;
      cout << f.nu << " " << f.mu << " " << f.n << endl;
      cout << nex2 << " " << y << endl;
      O2SCL_ERR2("At location 'X' in ",
		 "sn_fermion::nu_from_n().",exc_efailed);
      if (y!=1.0) {
	nex=nex2;
	nex2=1.0e10;
      }
    }
  }

  // If nothing worked, call the error handler
  if (y==1.0) {
    O2SCL_ERR2("Couldn't find reasonable initial guess in ",
	       "sn_fermion::nu_from_n().",exc_efailed);
  }

  // Perform full solution
  funct_mfptr<sn_fermion> mf(this,&sn_fermion::solve_fun);
  int ret=density_root->solve(nex,mf);
  f.nu=nex*temper;
  if (ret!=0) {
    O2SCL_ERR("Solver failed in sn_fermion::nu_from_n().",exc_efailed);
  }
  
  return;
}

void sn_fermion::calc_density(fermion_deriv &f, double temper) {

  T=temper;
  fp=&f;
  
  if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
  nu_from_n(f,temper);
  
  if (f.non_interacting) { f.mu=f.nu; }
  
  calc_mu(f,temper);

  return;
}

double sn_fermion::deg_density_fun(double k) {
  double E, ret;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }
  ret=k*k*fermi_function(E,fp->nu,T,exp_limit);
  return ret;
}

double sn_fermion::deg_density_T_fun(double k) {
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=k*k*(E-fp->nu)/T/T*
      fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=(2.0*k*k/T+E*E/T-E*(fp->nu+fp->m)/T-k*k*(fp->nu+fp->m)/T/E)*
      fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::deg_density_mu_fun(double k) {
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=k*k/T*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=(E*E+k*k)/E*fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::deg_energy_fun(double k) {
  double E, ret;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }
  ret=k*k*E*fermi_function(E,fp->nu,T,exp_limit);
  return ret;
}

double sn_fermion::deg_entropy_fun(double k) {

  double E, ret, nx, nm1;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }
  
  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-fp->nu)/(T))<-200.0) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-fp->nu)/T)<-exp_limit) {
    // Should this be E/T-nu/T or (E-nu)/T ?
    nm1=-exp(E/T-fp->nu/T);
    ret=k*k*nm1*log(-nm1);
  } else {
    nx=fermi_function(E,fp->nu,T,exp_limit);
    if (nx==0.0) ret=0.0;
    else ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx));
  }
  if (!o2scl::is_finite(ret)) {
    ret=0.0;
    //O2SCL_ERR("Entropy not finite in sn_fermion::deg_entropy_fun().",
    //exc_efailed);
  }

  return ret;
}
  
double sn_fermion::deg_entropy_T_fun(double k) {
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=k*k*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit))*
      pow(E-fp->nu,2.0)/pow(T,3.0);
  } else {
    ret=(E-fp->m-fp->nu)/E/T/T*
      (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(fp->nu+fp->m))*
      fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::deg_density_ms_fun(double k) {
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=-k*k*fp->ms/(E+fp->m)/T*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=-fp->ms*fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::density_fun(double u) {
  double k=u*(T), E, ret;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }
  ret=(T)*k*k*fermi_function(E,fp->nu,T,exp_limit);
  return ret;
}

double sn_fermion::density_T_fun(double u) {
  double k=u*(T);
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=k*k*(E-fp->nu)/T*
      fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=(2.0*k*k/T+E*E/T-E*(fp->nu+fp->m)/T-k*k*(fp->nu+fp->m)/T/E)*
      T*fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::density_mu_fun(double u) {
  double k=u*(T);
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=k*k*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=T*(E*E+k*k)/E*fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::energy_fun(double u) {
  double k=u*(T), E, ret;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }
  ret=u*u*E*fermi_function(E,fp->nu,T,exp_limit)/T;
  return ret;
}

double sn_fermion::entropy_fun(double u) {
  double k=u*(T), E, ret, nx, nm1;
  if (fp->inc_rest_mass) {
    E=gsl_hypot(k,fp->ms);
  } else {
    E=gsl_hypot(k,fp->ms)-fp->m;
  }

  // If the argument to the exponential is really small, then the
  // value of the integrand is just zero
  if (((E-fp->nu)/(T))<-200.0) {
    ret=0.0;
    // Otherwise, if the argument to the exponential is still small,
    // then addition of 1 makes us lose precision, so we use an
    // alternative:
  } else if (((E-fp->nu)/T)<-30.0) {
    nm1=-exp(E/(T)-fp->nu/(T));
    ret=k*k*nm1*log(-nm1)*(T);
  } else {
    nx=fermi_function(E,fp->nu,T,exp_limit);
    if (nx==0.0) ret=0.0;
    else ret=-k*k*(nx*log(nx)+(1.0-nx)*log(1.0-nx))*T;
  }

  return ret;
}

double sn_fermion::entropy_T_fun(double u) {
  double k=u*T;
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=T*k*k*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit))*
      pow(E-fp->nu,2.0)/pow(T,3.0);
  } else {
    ret=(E-fp->m-fp->nu)/E/T*
      (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(fp->nu+fp->m))*
      fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  return ret;
}

double sn_fermion::density_ms_fun(double u) {
  double k=u*T;
  double E=gsl_hypot(k,fp->ms), ret;
  if (method==direct) {
    E-=fp->m;
    ret=-k*k*fp->ms/(E+fp->m)/T*fermi_function(E,fp->nu,T,exp_limit)*
      (1.0-fermi_function(E,fp->nu,T,exp_limit));
  } else {
    ret=-fp->ms*fermi_function(E-fp->m,fp->nu,T,exp_limit);
  }
  ret*=T;
  return ret;
}

double sn_fermion::solve_fun(double x) {
  double nden, yy;
  
  fp->nu=T*x;
  
  // 6/6/03 - Should this be included? I think yes!
  if (fp->non_interacting) fp->mu=fp->nu;
  
  bool deg=true;
  if (fp->inc_rest_mass) {
    if ((fp->nu-fp->ms)/T<deg_limit) deg=false;
  } else {
    if ((fp->nu+(fp->m-fp->ms))/T<deg_limit) deg=false;
  }
  
  funct_mfptr<sn_fermion> 
    density_fun_f(this,&sn_fermion::density_fun);
  funct_mfptr<sn_fermion> 
    deg_density_fun_f(this,&sn_fermion::deg_density_fun);
  
  if (!deg) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=fp->g/2.0/pi2;
    
    yy=(fp->n-nden)/fp->n;

  } else {
    
    double arg;
    if (fp->inc_rest_mass) {
      arg=pow(upper_limit_fac*T+fp->nu,2.0)-fp->ms*fp->ms;
    } else {
      arg=pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms;
    }

    double ul;

    if (arg>0.0) {
      
      ul=sqrt(arg);
      nden=dit->integ(deg_density_fun_f,0.0,ul);
      nden*=fp->g/2.0/pi2;

    } else {
      nden=0.0;
    }
    
    yy=(fp->n-nden)/fp->n;
  }

  return yy;
}

void sn_fermion::pair_mu(fermion_deriv &f, double temper) {

  T=temper;
  fp=&f;

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
  
  return;
}

void sn_fermion::pair_density(fermion_deriv &f, double temper) {
  double nex;
  int ret;
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in sn_fermion().",exc_eunimpl);
  }

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  T=temper;
  fp=&f;

  nex=f.nu/temper;
  funct_mfptr<sn_fermion> mf(this,&sn_fermion::pair_fun);
  ret=density_root->solve(nex,mf);
  f.nu=nex*temper;

  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  return;
}

double sn_fermion::pair_fun(double x) {
  double nden, yy;
  
  fp->nu=T*x;
  
  if (fp->non_interacting) fp->mu=fp->nu;
  
    funct_mfptr<sn_fermion> density_fun_f(this,&sn_fermion::density_fun);
    funct_mfptr<sn_fermion> 
      deg_density_fun_f(this,&sn_fermion::deg_density_fun);

  if (fp->nu/T < deg_limit) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=fp->g/2.0/pi2;
    yy=nden;

  } else {
    
    double ulimit=sqrt(pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms);
    nden=dit->integ(deg_density_fun_f,0.0,ulimit);
    nden*=fp->g/2.0/pi2;
    yy=nden;
  }
  
  if (fp->inc_rest_mass) fp->nu=-T*x;
  else fp->nu=-T*x-2.0*fp->m;
  
  if (fp->nu/T < deg_limit) {
    
    nden=nit->integ(density_fun_f,0.0,0.0);
    nden*=fp->g/2.0/pi2;
    yy-=nden;
    
  } else {
    
    double ulimit=sqrt(pow(upper_limit_fac*T+fp->nu+fp->m,2.0)-fp->ms*fp->ms);
    nden=dit->integ(deg_density_fun_f,0.0,ulimit);
    nden*=fp->g/2.0/pi2;
    yy-=nden;
  }
  
  yy=yy/fp->n-1.0;

  return yy;
}

double sn_fermion::deriv_calibrate(fermion_deriv &f, int verbose, 
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
    cout << "In fermion_T::calibrate(), loading file named\n\t" 
	 << fname << "\n" << endl;
  }
  table<> tab;
  hdf_file hf;
  hf.open(fname);
  string name;
  hdf_input(hf,tab,name);
  hf.close();
  
  if (tab.get_nlines()==0) {
    string str="Failed to load data from file '"+fname+
      "' in fermion_T::calibrate(). Bad filename?";
    O2SCL_ERR_RET(str.c_str(),exc_efilenotfound);
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
    bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    bad.dndT=0.0; bad.dndmu=0.0; bad.dsdT=0.0; bad.dndm=0.0;
    
    // Temperature loop
    for(double T2=1.0e-2;T2<=1.001e2;T2*=1.0e2) {

      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("c1",i);
	double psi=tab.get("c2",i);
	exact.n=tab.get("c3",i);
	exact.ed=tab.get("c4",i);
	exact.pr=tab.get("c5",i);
	exact.en=tab.get("c6",i);
	exact.dndT=tab.get("c7",i);
	exact.dndmu=tab.get("c8",i);
	exact.dsdT=tab.get("c9",i);
	exact.dndm=0.0;
      
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
	  cout << "\t" << f.dndT << " " << f.dndmu << " " << f.dsdT << " "
	       << f.dndm << endl;
	  cout << "exact : " << exact.n << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "\t" << exact.dndT << " " << exact.dndmu << " " 
	       << exact.dsdT << " " << exact.dndm << endl;
	  cout << "bad   : " << bad.n << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << "\t" << bad.dndT << " " << bad.dndmu << " " 
	       << bad.dsdT << " " << bad.dndm << endl;
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
    dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
    bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
    // Temperature loop
    for(double T2=1.0e-2;T2<=1.001e2;T2*=1.0e2) {
      
      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("c1",i);
	double psi=tab.get("c2",i);
	f.n=tab.get("c3",i);	
	exact.ed=tab.get("c4",i);
	exact.pr=tab.get("c5",i);
	exact.en=tab.get("c6",i);
	exact.dndT*=pow(T2,2.0);
	exact.dndmu*=pow(T2,2.0);
	exact.dsdT*=pow(T2,2.0);
	exact.dndm=0.0;

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
	  cout << "\t" << f.dndT << " " << f.dndmu << " " << f.dsdT << " "
	       << f.dndm << endl;
	  cout << "exact : " << exact.mu << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "\t" << exact.dndT << " " << exact.dndmu << " " 
	       << exact.dsdT << " " << exact.dndm << endl;
	  cout << "bad   : " << bad.mu << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << "\t" << bad.dndT << " " << bad.dndmu << " " 
	       << bad.dsdT << " " << bad.dndm << endl;
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

