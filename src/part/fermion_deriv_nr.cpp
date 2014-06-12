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

#include <o2scl/fermion_deriv_nr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// fermion_deriv_nr class

fermion_deriv_nr::fermion_deriv_nr() {
  
  flimit=20.0;
  guess_from_nu=true;

  density_root=&def_density_root;
}

fermion_deriv_nr::~fermion_deriv_nr() {
}

void fermion_deriv_nr::calc_mu(fermion_deriv &f, double temper) {
  
  T=temper;
  fp=&f;

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
  f.dndm=1.5*pfac2*temper/f.ms*half;
  
  return;
}

void fermion_deriv_nr::nu_from_n(fermion_deriv &f, double temper) {
  double nex;
  
  T=temper;
  fp=&f;

  if (guess_from_nu) {
    nex=f.nu/temper;
  } else {
    O2SCL_ERR("guess_from_nu==false not implemented in fermion_deriv_nr.",
	      exc_eunimpl);
  }
  funct11 mf=std::bind(std::mem_fn<double(double)>
		       (&fermion_deriv_nr::solve_fun),
		       this,std::placeholders::_1);
    
  density_root->solve(nex,mf);
  f.nu=nex*temper;

  return;
}

void fermion_deriv_nr::calc_density(fermion_deriv &f, double temper) {

  T=temper;
  fp=&f;
  
  if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
  nu_from_n(f,temper);
  
  if (f.non_interacting) { f.mu=f.nu; }
  
  calc_mu(f,temper);
  
  return;
}

double fermion_deriv_nr::solve_fun(double x) {
  double nden, y, yy;
  
  fp->nu=T*x;
  
  // 6/6/03 - I think this should this be included.
  if (fp->non_interacting) fp->mu=fp->nu;

  if (fp->inc_rest_mass) {
    y=(fp->nu-fp->m)/T;
  } else {
    y=fp->nu/T;
  }

  nden=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  nden*=fp->g*pow(2.0*fp->ms*T,1.5)/4.0/pi2;

  yy=(fp->n-nden)/fp->n;
  
  return yy;
}

void fermion_deriv_nr::pair_mu(fermion_deriv &f, double temper) {
  
  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

  fermion_deriv antip(f.ms,f.g);
  f.anti(antip);

  calc_mu(f,temper);

  calc_mu(antip,temper);
  f.n-=antip.n;
  f.pr+=antip.pr;
  f.ed+=antip.ed;
  f.en+=antip.en;
  
  return;
}

void fermion_deriv_nr::pair_density(fermion_deriv &f, double temper) {
  double nex;
  
  if (temper<=0.0) {
    O2SCL_ERR("T=0 not implemented in fermion_deriv_nr().",exc_eunimpl);
  }
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  T=temper;
  fp=&f;

  nex=f.nu/temper;
  funct11 mf=std::bind(std::mem_fn<double(double)>
		       (&fermion_deriv_nr::pair_fun),
		       this,std::placeholders::_1);

  density_root->solve(nex,mf);
  f.nu=nex*temper;

  if (f.non_interacting==true) { f.mu=f.nu; }
  
  pair_mu(f,temper);

  return;
}

double fermion_deriv_nr::pair_fun(double x) {
  double nden, y, yy;

  fp->nu=T*x;

  // 6/6/03 - Should this be included? I think yes!
  if (fp->non_interacting) fp->mu=fp->nu;

  if (fp->inc_rest_mass) {
    y=(fp->nu-fp->m)/T;
  } else {
    y=fp->nu/T;
  }

  nden=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  nden*=fp->g*pow(2.0*fp->ms*T,1.5)/4.0/pi2;
  
  yy=nden;

  if (fp->inc_rest_mass) {
    y=-(fp->nu-fp->m)/T;
  } else {
    y=-fp->nu/T;
  }
  
  nden=gsl_sf_fermi_dirac_half(y)*sqrt(pi)/2.0;
  nden*=fp->g*pow(2.0*fp->ms*T,1.5)/4.0/pi2;
  
  yy-=nden;
  
  yy=yy/fp->n-1.0;

  return yy;
}

