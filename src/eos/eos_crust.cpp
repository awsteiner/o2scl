/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

#include <o2scl/eos_crust.h>
// For unit conversions
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_crust::eos_crust() {
  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  // Put an initial guess for the electron density - Is this necessary?
  // A new initial guess is given in calc_pressure()
  e.n=0.4;
  def_mass.B=-15.76;
  def_mass.Sv=94.8/4.0;
  def_mass.Ss=17.81;
  def_mass.Ec=0.71;
  def_mass.Epair=39.0;
  nmp=&def_mass;
}

double eos_crust::lattice_energy(int Z) {
  return -1.444*pow(Z,2.0/3.0)*o2scl_const::fine_structure_f<double>()*
    pow(e.n,4.0/3.0);
}

double eos_crust::mass_formula(int Z, int A) {
  return nmp->total_mass(Z,A-Z)/hc_mev_fm;
}

int eos_crust::eq274(size_t nv, const ubvector &nx, ubvector &ny, 
		     int Zt) {
  
  // Lattice energy density and pressure
  double PL, epsL;

  e.n=nx[0];

  // Equations 2.3.5-8:
  fzt.calc_density_zerot(e);
  
  epsL=lattice_energy(Zt);
  
  // Equation 2.7.5
  PL=epsL/3.0;

  // Pressure is electron Fermi gas pressure plus pressure due
  // to the lattice
  ny[0]=log(eos_thermo->pr)-log(e.pr+PL);

  if (!std::isfinite(ny[0])) {
    return 1;
  }
  
  return 0;
}

double eos_crust::gibbs(int Z, int A) {
  double g;
  ubvector nx(1);
  
  nx[0]=e.n;
  mm_funct mff=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &, int)>
     (&eos_crust::eq274),this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,Z);
  gs.msolve(1,nx,mff);
  e.n=nx[0];
  fzt.kf_from_density(e);
  e.mu=sqrt(e.kf*e.kf+e.m*e.m);
  fzt.calc_mu_zerot(e);

  // Equation 2.7.7:
  g=(mass_formula(Z,A)+Z*e.mu+4.0/3.0*Z*lattice_energy(Z)/e.n)/A;

  return g;
}

double eos_crust::energy(double barn, int Z, int A) {
  double ed;

  fzt.calc_density_zerot(e);
  
  // Eq. 2.7.1
  ed=barn/A*mass_formula(Z,A)+e.ed+lattice_energy(Z);

  return ed;
}

int eos_crust::calc_pressure(thermo &th, double &barn, int &minz, 
			     int &mina) {
  double gb, mingb=10.0, MAZ, epsL, ne=0.0;
  int ret;

  set_thermo(th);
  mingb=10.0;

  // Provide an initial guess for gibbs() which solves eq274().
  e.n=th.pr/10.0;
  
  for(int A=50;A<140;A++) {
    for(int Z=26;Z<50;Z++) {
      
      gb=gibbs(Z,A);
      
      if (gb<mingb) {
	mina=A;
	minz=Z;
	mingb=gb;
	ne=e.n;
      }
      
    }
  }
  
  // Get the right electron density for the favored nucleus
  gb=gibbs(minz,mina);
    
  // Right before Eq. 2.7.6 (this is just charge neutrality)
  barn=e.n*mina/minz;
  
  epsL=lattice_energy(minz);

  MAZ=mass_formula(minz,mina);
  
  // Eq. 2.7.1
  th.ed=barn/mina*MAZ+e.ed+epsL;
  th.en=0.0;

  return 0;
}

int eos_crust::calc_density(double barn, thermo &th, int &Z, int &A) {
  double ed, mined=10.0;
  int mina=0, minz=0;
  
  set_thermo(th);

  mined=10.0;

  for(A=50;A<140;A++) {
    for(Z=26;Z<50;Z++) {
      if (nmp->is_included(Z,A-Z)) {
	
	e.n=barn*Z/A;
	ed=energy(barn,Z,A);
	
	if (ed<mined) {
	  mina=A;
	  minz=Z;
	  mined=ed;
	}
      }
    }
  }
    
  A=mina;
  Z=minz;
  e.n=barn*Z/A;
  ed=energy(barn,Z,A);
  
  // Eq. 2.7.1
  th.ed=barn/A*mass_formula(Z,A)+e.ed+lattice_energy(Z);
  th.pr=e.pr+lattice_energy(Z)/3.0;
  th.en=0.0;

  return 0;
}

int eos_crust::calc_density_fixedA(double barn, thermo &th, int &Z, int A) {
  double ed, mined=10.0;
  int mina=0, minz=0;

  set_thermo(th);

  mined=10.0;

  for(Z=2;Z<50;Z++) {
    if (nmp->is_included(Z,A-Z)) {
      
      e.n=barn*Z/A;
      ed=energy(barn,Z,A);
      
      if (ed<mined) {
	mina=A;
	minz=Z;
	mined=ed;
      }
    }
  }
    
  A=mina;
  Z=minz;
  e.n=barn*Z/A;
  ed=energy(barn,Z,A);
  
  // Eq. 2.7.1
  th.ed=barn/A*mass_formula(Z,A)+e.ed+lattice_energy(Z);
  th.pr=e.pr+lattice_energy(Z)/3.0;
  th.en=0.0;

  return 0;
}

