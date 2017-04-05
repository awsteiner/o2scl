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

#include <o2scl/fermion_nonrel.h>
#include <o2scl/classical.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

fermion_nonrel::fermion_nonrel() {
  density_root=&def_density_root;
}

fermion_nonrel::~fermion_nonrel() {
}

void fermion_nonrel::calc_mu_zerot(fermion &f) {
  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  if (f.inc_rest_mass) {
    f.kf=sqrt(2.0*f.ms*(f.nu-f.m));
  } else {
    f.kf=sqrt(2.0*f.ms*f.nu);
  }
  f.n=f.kf*f.kf*f.kf*f.g/6.0/pi2;
  f.ed=f.g*pow(f.kf,5.0)/20.0/pi2/f.ms;
  if (f.inc_rest_mass) f.ed+=f.n*f.m;
  f.pr=-f.ed+f.n*f.nu;
  f.en=0.0;
  return;
}

void fermion_nonrel::calc_density_zerot(fermion &f) {
  if (f.non_interacting) { f.ms=f.m; }
  kf_from_density(f);
  f.nu=f.kf*f.kf/2.0/f.ms;
  f.ed=f.g*pow(f.kf,5.0)/20.0/pi2/f.ms;
  if (f.inc_rest_mass) {
    f.ed+=f.n*f.m;
    f.nu+=f.m;
  }
  f.pr=-f.ed+f.n*f.nu;
  f.en=0.0;
  
  if (f.non_interacting) { f.mu=f.nu; }
  return;
}

void fermion_nonrel::calc_mu(fermion &f, double temper) {

  if (temper<0.0) {
    O2SCL_ERR("Temperature less than zero in fermion_nonrel::calc_mu().",
	      exc_einval);
  }
  if (temper==0.0) {
    calc_mu_zerot(f);
    return;
  }

  double y, sy, spi, ey, int1, int2;

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

  if (f.ms<0.0) {
    O2SCL_ERR2("Mass negative in ",
	       "fermion_nonrel::calc_mu().",exc_einval);
  }

  if (temper<=0.0) {
    calc_mu_zerot(f);
    return;
  }

  if (f.inc_rest_mass) {
    y=(f.nu-f.m)/temper;
  } else {
    y=f.nu/temper;
  }

  // Number density
  f.n=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  f.n*=f.g*pow(2.0*f.ms*temper,1.5)/4.0/pi2;

  // Energy density:
  f.ed=gsl_sf_fermi_dirac_3half(y)*0.75*sqrt(pi);

  if (f.inc_rest_mass) {
    
    // Finish energy density
    f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/pi2/f.ms;
    f.ed+=f.n*f.m;
    
    // entropy density
    f.en=(5.0*(f.ed-f.n*f.m)/3.0-(f.nu-f.m)*f.n)/temper;
    
    // pressure
    f.pr=2.0*(f.ed-f.n*f.m)/3.0;
    
  } else {

    // Finish energy density
    f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/pi2/f.ms;

    // entropy density
    f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
    
    // pressure
    f.pr=2.0*f.ed/3.0;
    
  }

  if (!std::isfinite(f.nu) || !std::isfinite(f.n)) {
    O2SCL_ERR2("Chemical potential or density in ",
		   "fermion_nonrel::calc_mu().",exc_efailed);
  }
  
  return;
}

void fermion_nonrel::nu_from_n(fermion &f, double temper) {

  // Use initial value of nu for initial guess
  double nex;
  if (f.inc_rest_mass) {
    nex=-(f.nu-f.m)/temper;
  } else {
    nex=-f.nu/temper;
  } 
  
  // Make a correction if nex is too small and negative
  // (Note GSL_LOG_DBL_MIN is about -708)
  if (nex>-GSL_LOG_DBL_MIN*0.9) nex=-GSL_LOG_DBL_MIN/2.0;
  
  funct mf=std::bind(std::mem_fn<double(double,double,double)>
		       (&fermion_nonrel::solve_fun),
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

  return;
}

int fermion_nonrel::calc_density(fermion &f, double temper) {

  if (f.ms<0.0) {
    O2SCL_ERR2("Mass negative in ",
	       "fermion_nonrel::calc_density().",exc_einval);
  }
  if (temper<0.0) {
    O2SCL_ERR2("Temperature less than zero in ",
	       "fermion_nonrel::calc_density().",exc_einval);
  }
  if (temper==0.0) {
    calc_density_zerot(f);
    return 0;
  }

  if (f.n<=0.0) {
    O2SCL_ERR2("Density less than or equal to zero in ",
	       "fermion_nonrel::calc_density().",exc_einval);
  }
  
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  double guess=f.nu;
  nu_from_n(f,temper);
  
  if (f.non_interacting) { f.mu=f.nu; }

  double y, spi, ey, sy;
  if (f.inc_rest_mass) {
    y=-(f.nu-f.m)/temper;
  } else {
    y=-f.nu/temper;
  }

  // energy density
  f.ed=gsl_sf_fermi_dirac_3half(-y)*sqrt(pi)*0.75;

  if (f.inc_rest_mass) {
    
    // Finish energy density
    f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/pi2/f.ms;
    f.ed+=f.n*f.m;
    
    // entropy density
    f.en=(5.0*(f.ed-f.n*f.m)/3.0-(f.nu-f.m)*f.n)/temper;
    
    // pressure
    f.pr=2.0*(f.ed-f.n*f.m)/3.0;
    
  } else {

    // Finish energy density
    f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/pi2/f.ms;

    // entropy density
    f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
    
    // pressure
    f.pr=2.0*f.ed/3.0;
    
  }
  
  if (f.non_interacting==true) { f.mu=f.nu; }
  
  return 0;
}

double fermion_nonrel::solve_fun(double x, double nog, double msT) {

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

