/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

#include <o2scl/fermion_deriv_nr.h>
// Used to get initial guesses for chemical potential if necessary
#include <o2scl/classical.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// fermion_deriv_nr class

fermion_deriv_nr::fermion_deriv_nr() {
  
  flimit=20.0;

  density_root=&def_density_root;
}

fermion_deriv_nr::~fermion_deriv_nr() {
}

int fermion_deriv_nr::calc_mu(fermion_deriv &f, double temper) {
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in fermion_deriv_nr().",exc_eunimpl);
  }
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  double pfac2=f.g*pow(f.ms*temper/2.0/pi,1.5)/temper, y;
  
  if (f.inc_rest_mass) {
    y=(f.nu-f.m)/temper;
  } else {
    y=f.nu/temper;
  }

  double half=gsl_sf_fermi_dirac_half(y);
  double mhalf=gsl_sf_fermi_dirac_mhalf(y);
  double thalf=gsl_sf_fermi_dirac_3half(y);
  
  // Number density:
  f.n=pfac2*half*temper;
  
  // Energy density:
  f.ed=pfac2*temper*temper*1.5*thalf;

  if (f.inc_rest_mass) {
    
    // Finish energy density
    f.ed+=f.n*f.m;
    
    // entropy density
    f.en=((f.ed-f.n*f.m)/0.6-(f.nu-f.m)*f.n)/temper;
    
    // pressure
    f.pr=(f.ed-f.n*f.m)/1.5;
    
  } else {
    
    // entropy density
    f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
    
    // pressure
    f.pr=2.0*f.ed/3.0;
    
  }
    
  f.dndT=pfac2*(1.5*half-y*mhalf);
  f.dndmu=pfac2*mhalf;
  f.dsdT=pfac2*(3.75*thalf-3.0*y*half+y*y*mhalf);
  
  return 0;
}

int fermion_deriv_nr::nu_from_n(fermion_deriv &f, double temper) {

  // Use initial value of nu for initial guess
  double nex;
  if (f.inc_rest_mass) {
    nex=-(f.nu-f.m)/temper;
  } else {
    nex=f.nu/temper;
  }

  // Make a correction if nex is too small and negative
  // (Note GSL_LOG_DBL_MIN is about -708)
  if (nex>-GSL_LOG_DBL_MIN*0.9) nex=-GSL_LOG_DBL_MIN/2.0;
  
  funct mf=std::bind(std::mem_fn<double(double,double,double)>
		     (&fermion_deriv_nr::solve_fun),
		     this,std::placeholders::_1,f.n/f.g,f.ms*temper);
    
  // Turn off convergence errors temporarily, since we'll
  // try again if it fails
  bool enc=density_root->err_nonconv;
  density_root->err_nonconv=false;
  int ret=density_root->solve(nex,mf);
  density_root->err_nonconv=enc;

  if (ret!=0) {

    // If it failed, try to get a guess from classical particle
    
    classical cl;
    cl.calc_density(f,temper);
    if (f.inc_rest_mass) {
      nex=-(f.nu-f.m)/temper;
    } else {
      nex=-f.nu/temper;
    } 
    ret=density_root->solve(nex,mf);
    
    // If it failed again, add error information
    if (ret!=0) {
      O2SCL_ERR("Solver failed in fermion_nonrel::nu_from_n().",ret);
    }
  }

  if (f.inc_rest_mass) {
    f.nu=-nex*temper+f.m;
  } else {
    f.nu=-nex*temper;
  }

  return 0;
}

int fermion_deriv_nr::calc_density(fermion_deriv &f, double temper) {

  if (f.m<0.0 || (f.non_interacting==false && f.ms<0.0)) {
    O2SCL_ERR2("Mass negative in ",
	       "fermion_deriv_nr::calc_density().",exc_einval);
  }
  if (temper<0.0) {
    O2SCL_ERR2("Temperature less than zero in ",
	       "fermion_deriv_nr::calc_density().",exc_einval);
  }
  if (f.n<=0.0) {
    O2SCL_ERR2("Density less than or equal to zero in ",
	       "fermion_deriv_nr::calc_density().",exc_einval);
  }
  
  if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }

  nu_from_n(f,temper);
  
  if (f.non_interacting) { f.mu=f.nu; }
  
  calc_mu(f,temper);
  
  return 0;
}

double fermion_deriv_nr::solve_fun(double x, double nog, double msT) {

  double nden;
  
  // If the argument to gsl_sf_fermi_dirac_half() is less
  // than GSL_LOG_DBL_MIN (which is about -708), then 
  // an underflow occurs. We just set nden to zero in this 
  // case, as this helps the solver find the right root.
  
  if (((-x)<GSL_LOG_DBL_MIN) || !std::isfinite(x)) nden=0.0;
  else nden=gsl_sf_fermi_dirac_half(-x)*sqrt(pi)/2.0;
  
  nden*=pow(2.0*msT,1.5)/4.0/pi2;
  return nden/nog-1.0;
}

int fermion_deriv_nr::pair_mu(fermion_deriv &f, double temper) {
  
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

int fermion_deriv_nr::pair_density(fermion_deriv &f, double temper) {
  double nex;
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in fermion_deriv_nr().",exc_eunimpl);
  }
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  nex=f.nu/temper;
  funct mf=std::bind(std::mem_fn<double(double,fermion_deriv &,double)>
		       (&fermion_deriv_nr::pair_fun),
		       this,std::placeholders::_1,std::ref(f),temper);

  density_root->solve(nex,mf);
  f.nu=nex*temper;

  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  return 0;
}

double fermion_deriv_nr::pair_fun(double x, fermion_deriv &f, double T) {
  double nden, y, yy;

  f.nu=T*x;

  // 6/6/03 - Should this be included? I think yes!
  if (f.non_interacting) f.mu=f.nu;

  if (f.inc_rest_mass) {
    y=(f.nu-f.m)/T;
  } else {
    y=f.nu/T;
  }

  nden=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  nden*=f.g*pow(2.0*f.ms*T,1.5)/4.0/pi2;
  
  yy=nden;

  if (f.inc_rest_mass) {
    y=-(f.nu-f.m)/T;
  } else {
    y=-f.nu/T;
  }
  
  nden=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  nden*=f.g*pow(2.0*f.ms*T,1.5)/4.0/pi2;
  
  yy-=nden;
  
  yy=yy/f.n-1.0;

  return yy;
}

