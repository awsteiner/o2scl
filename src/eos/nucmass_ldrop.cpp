/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/nucmass_ldrop.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_ldrop::nucmass_ldrop() {

  def_neutron.init(o2scl_settings.get_convert_units().convert
		   ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  def_proton.init(o2scl_settings.get_convert_units().convert
		  ("kg","1/fm",o2scl_mks::mass_proton),2.0);

  def_neutron.non_interacting=false;
  def_proton.non_interacting=false;
  def_neutron.inc_rest_mass=true;
  def_proton.inc_rest_mass=true;
  n=&def_neutron;
  p=&def_proton;
  heos=&def_had_eos;

  // Load NL4 EOS
  def_had_eos.ms=508.194;
  def_had_eos.mw=782.501;
  def_had_eos.mr=763.0;
  def_had_eos.mnuc=939.0;
  def_had_eos.ms/=o2scl_const::hc_mev_fm; 
  def_had_eos.mw/=o2scl_const::hc_mev_fm; 
  def_had_eos.mr/=o2scl_const::hc_mev_fm; 
  def_had_eos.mnuc/=o2scl_const::hc_mev_fm;
	
  double gs, gw, gr;
  gs=10.217;
  gw=12.868;
  gr=4.474;
  def_had_eos.b=-10.431;
  def_had_eos.c=-28.885;
  def_had_eos.b/=-def_had_eos.mnuc*pow(fabs(gs),3.0);
  def_had_eos.c/=pow(gs,4.0);
  gr*=2.0;
  def_had_eos.cs=gs/def_had_eos.ms;
  def_had_eos.cw=gw/def_had_eos.mw;
  def_had_eos.cr=gr/def_had_eos.mr;
	
  def_had_eos.xi=0.0; 
  def_had_eos.zeta=0.0;
  def_had_eos.a1=0.0;
  def_had_eos.a2=0.0;
  def_had_eos.a3=0.0;
  def_had_eos.a4=0.0;
  def_had_eos.a5=0.0;
  def_had_eos.a6=0.0;
  def_had_eos.b1=0.0;
  def_had_eos.b2=0.0;
  def_had_eos.b3=0.0;

  n1=0.0;
  n0=0.16;
  surften=1.1;
  nn=0.0;
  np=0.0;
  coul_coeff=1.0;
      
  nfit=4;

  large_vals_unphys=false;
}

double nucmass_ldrop::mass_excess_d(double Z, double N) {
  double ret=0.0;
  
  ret=drip_binding_energy_d(Z,N,0.0,0.0,0.0,0.0);
      
  // Convert from binding energy to mass excess
  ret-=((N+Z)*o2scl_mks::unified_atomic_mass-Z*o2scl_mks::mass_electron-
	N*o2scl_mks::mass_neutron-Z*o2scl_mks::mass_proton)*
    o2scl_settings.get_convert_units().convert("kg","MeV",1.0);
  
  return ret;
      
}

double nucmass_ldrop::drip_binding_energy_d
(double Z, double N, double npout, double nnout, double chi, double T) {

  double ret=0.0, A=Z+N, nL;
      
  // Force inc_rest_mass to true
  n->inc_rest_mass=true;
  p->inc_rest_mass=true;

  // Determine the inner densities
  double delta=(1-2.0*Z/A);
  nL=n0+n1*delta*delta;
  np=nL*(1.0-delta)/2.0;
  nn=nL*(1.0+delta)/2.0;
  if (nn>0.20 || np>0.20) {
    if (large_vals_unphys) return 1.0e99;
    std::cout << "Densities too large 1 (n0,n1,nn,np): "
              << n0 << " " << n1 << " "
	      << nn << " " << np << std::endl;
    O2SCL_ERR2("Densities too large in ",
               "nucmass_ldrop::drip_binding_energy_d().",
               o2scl::exc_efailed);
  }

  // Determine radii
  Rn=cbrt(3.0*N/4.0/o2scl_const::pi/nn);
  Rp=cbrt(3.0*Z/4.0/o2scl_const::pi/np);
      
  // Compute bulk energy per baryon
  n->n=nn;
  p->n=np;

  // In principle, these next two lines shouldn't be needed
  // but they're here for now just in case
  n->mu=n->m;
  p->mu=p->m;

  int err=heos->calc_e(*n,*p,th);
  if (err!=0) {
    O2SCL_ERR2("Hadronic EOS failed in ",
	       "nucmass_ldrop::drip_binding_energy_d().",exc_efailed);
  }
  bulk=(th.ed-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm;
  ret+=bulk;

  // Determine surface energy per baryon
  surf=surften*cbrt(36.0*o2scl_const::pi*nL)/nL/cbrt(A);
  ret+=surf;
      
  // Add Coulomb energy per baryon
  coul=coul_coeff*0.8*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*Rp*np*np/nL;
  ret+=coul;
      
  // Convert to total binding energy
  ret*=A;
      
  return ret;
}

int nucmass_ldrop::fit_fun(size_t nv, const ubvector &x) {
  surften=x[0];
  n1=x[1];
  n0=x[2];
  coul_coeff=x[3];
  return 0;
}

int nucmass_ldrop::guess_fun(size_t nv, ubvector &x) {
  x[0]=surften;
  x[1]=n1;
  x[2]=n0;
  x[3]=coul_coeff;
  return 0;
}

nucmass_ldrop_skin::nucmass_ldrop_skin() {
  doi=0.8;
  ss=0.5;
  nfit=6;
  full_surface=true;
  new_skin_mode=false;
  rel_vacuum=true;
  
  pp=1.25;
  a0=0.935;
  a2=-5.1;
  a4=-1.1;
  Tchalf=20.085/o2scl_const::hc_mev_fm;
}

int nucmass_ldrop_skin::fit_fun(size_t nv, const ubvector &x) {
  doi=x[0];
  surften=x[1];
  ss=x[2];
  coul_coeff=x[3];
  n1=x[4];
  n0=x[5];
  return 0;
}

int nucmass_ldrop_skin::guess_fun(size_t nv, ubvector &x) {
  x[0]=doi;
  x[1]=surften;
  x[2]=ss;
  x[3]=coul_coeff;
  x[4]=n1;
  x[5]=n0;
  return 0;
}

double nucmass_ldrop_skin::drip_binding_energy_d
(double Z, double N, double npout, double nnout, double chi, double T) {
  
  int err;
  double ret=0.0, A=Z+N, nL;

  // Force inc_rest_mass to true
  n->inc_rest_mass=true;
  p->inc_rest_mass=true;

  // Determine the inner densities
  double I=1.0-2.0*Z/A, delta;
  double X=Z/A;
  nL=n0+n1*I*I;
      
  delta=I*doi;
  np=nL*(1.0-delta)/2.0;
  nn=nL*(1.0+delta)/2.0;
  if (nn>0.20 || np>0.20) {
    if (large_vals_unphys) return 1.0e99;
    std::cout << "Densities too large 2 (n0,n1,nn,np):\n  "
              << n0 << " " << n1 << " "
	      << nn << " " << np << std::endl;
    std::cout << "nL,delta,I: " << nL << " " << delta << " " << I << endl;
    O2SCL_ERR2("Densities too large in ",
               "nucmass_ldrop::drip_binding_energy_d().",
               o2scl::exc_efailed);
  }

  if (!std::isfinite(nn) || !std::isfinite(np)) {
    O2SCL_ERR2("Neutron or proton density not finite in ",
	       "nucmass_ldrop::drip_binding_energy_d().",exc_efailed);
    return 0.0;
  }

  // Determine radii

  Rn=cbrt(3.0*N/nn/4.0/o2scl_const::pi);
  Rp=cbrt(3.0*Z/np/4.0/o2scl_const::pi);
	
  // Bulk part of the free energy per baryon

  if (!new_skin_mode) {

    // If new_skin_mode is false, just compute the 
    // bulk energy once, given nn and np
    n->n=nn;
    p->n=np;
    n->mu=n->m;
    p->mu=p->m;
    
    if (n->n<0.0) n->n=1.0e-3;
    if (p->n<0.0) p->n=1.0e-3;

    if (T<=0.0) {
      err=heos->calc_e(*n,*p,th);
      bulk=(th.ed-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm;
    } else {
      err=heos->calc_temp_e(*n,*p,T,th);
      bulk=(th.ed-T*th.en-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm;
    }
    if (err!=0) {
      O2SCL_ERR2("Hadronic eos failed in ",
		 "nucmass_ldrop_skin::drip_binding_energy_d().",
		 exc_efailed);
    }
    ret+=bulk;

  } else {

    // Otherwise, try to separate out the contribution
    // from the skin. 

    // First compute the relative strength of the core
    // and skin contributions 
    double Acore, Askin;
    if (N>=Z) {
      Acore=Z*(nn+np)/np;
      Askin=A-Acore;
    } else {
      Acore=N*(nn+np)/nn;
      Askin=A-Acore;
    }
    if (Acore>A) {
      Askin=0.0;
      Acore=A;
    }

    // The "core" contribution
    n->n=nn;
    p->n=np;
    n->mu=n->m;
    p->mu=p->m;
    if (T<=0.0) {
      err=heos->calc_e(*n,*p,th);
      bulk=(th.ed-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm*(Acore/A);
    } else {
      err=heos->calc_temp_e(*n,*p,T,th);
      bulk=(th.ed-T*th.en-nn*n->m-np*p->m)/nL*
	o2scl_const::hc_mev_fm*(Acore/A);
    }
    if (err!=0) {
      O2SCL_ERR2("Hadronic eos failed in ",
		 "nucmass_ldrop_skin::drip_binding_energy_d().",
		 exc_efailed);
    }

    // Note that, for this model, Rn>Rp iff N>Z.
    if (Rn>Rp) {

      // A neutron skin
      n->n=nn;
      p->n=npout;
      n->mu=n->m;
      p->mu=p->m;
      if (T<=0.0) {
	err=heos->calc_e(*n,*p,th);
	bulk+=(th.ed-nn*n->m)/nL*o2scl_const::hc_mev_fm*(Askin/A);
      } else {
	err=heos->calc_temp_e(*n,*p,T,th);
	bulk+=(th.ed-T*th.en-nn*n->m)/nL*
	  o2scl_const::hc_mev_fm*(Askin/A);
      }
      if (err!=0) {
	O2SCL_ERR2("Hadronic eos failed in ",
		   "nucmass_ldrop_skin::drip_binding_energy_d().",
		   exc_efailed);
      }

    } else if (Rp>Rn) {

      // A proton skin
      n->n=nnout;
      p->n=np;
      n->mu=n->m;
      p->mu=p->m;
      if (T<=0.0) {
	err=heos->calc_e(*n,*p,th);
	bulk+=(th.ed-np*p->m)/nL*o2scl_const::hc_mev_fm*(Askin/A);
      } else {
	err=heos->calc_temp_e(*n,*p,T,th);
	bulk+=(th.ed-T*th.en-np*p->m)/nL*
	  o2scl_const::hc_mev_fm*(Askin/A);
      }
      if (err!=0) {
	O2SCL_ERR2("Hadronic eos failed in ",
		   "nucmass_ldrop_skin::drip_binding_energy_d().",
		   exc_efailed);
      }
    }
    
    if (!rel_vacuum && (nnout>0.0 || npout>0.0)) {
      
      double Rbig;
      if (Rn>Rp) Rbig=Rn;
      else Rbig=Rp;
      
      n->n=nnout;
      p->n=npout;
      n->mu=n->m;
      p->mu=p->m;
      if (T<=0.0) {
	err=heos->calc_e(*n,*p,th);
	bulk-=(th.ed-nnout*n->m-npout*p->m)*
	  4.0*o2scl_const::pi/3.0*pow(Rbig,3.0)/A*o2scl_const::hc_mev_fm;
      } else {
	err=heos->calc_temp_e(*n,*p,T,th);
	bulk-=(th.ed-T*th.en-nnout*n->m-npout*p->m)*
	  4.0*o2scl_const::pi/3.0*pow(Rbig,3.0)/A*o2scl_const::hc_mev_fm;
      }
    }

    ret+=bulk;

  }

  // Determine surface energy per baryon

  if (full_surface) {

    double x=np/(nn+np);
    double x3=x*x*x;
    double omx=1.0-x;
    double omx3=omx*omx*omx;
    double bcoeff;
    if (ss==0.0) bcoeff=-16.0+96.0*surften/0.5;
    else bcoeff=-16.0+96.0*surften/ss;
    double bfun=(16.0+bcoeff)/(1.0/x3+bcoeff+1.0/omx3);
    double y=0.5-x;
    double y2=y*y, y4=y2*y2;
    double a=a0+a2*y2+a4*y4;
    double arg=1.0-3.313*y2-7.362*y4;
    double Tc=Tchalf*sqrt(1.0-3.313*y2-7.362*y4);
	
    if (T<Tc) {
      bfun*=pow((1.0-T*T/Tc/Tc)/(1.0+a*T*T/Tc/Tc),pp);
    } else {
      bfun=0.0;
    }
	
    surf=surften*bfun*cbrt(36.0*o2scl_const::pi*nL)/nL/cbrt(A);
	
  } else {
	
    surf=surften*(1.0-ss*delta*delta)*
      cbrt(36.0*o2scl_const::pi*nL)/nL/cbrt(A);

  }
      
  ret+=surf;
      
  // Add Coulomb energy per baryon

  double chip;
  if (Rn>Rp) {
    chip=chi*pow(Rp/Rn,3.0);
  } else {
    chip=chi;
  }

  // This mass formula only works for d=3, but the full
  // formula is given here for reference
  double dim=3.0;
  double fdu=(2.0/(dim-2.0)*(1.0-0.5*dim*pow(chip,1.0-2.0/dim))+chip)/
    (dim+2.0);
  coul=coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*
    Rp*Rp*pow(fabs(np-npout),2.0)/nL*fdu;
  ret+=coul;
  
  // Convert to total binding energy
  ret*=A;

  return ret;
}

int nucmass_ldrop_pair::fit_fun(size_t nv, const ubvector &x) {
  doi=x[0];
  surften=x[1];
  ss=x[2];
  coul_coeff=x[3];
  n1=x[4];
  n0=x[5];
  Epair=x[6];
  return 0;
}

int nucmass_ldrop_pair::guess_fun(size_t nv, ubvector &x) {
  x[0]=doi;
  x[1]=surften;
  x[2]=ss;
  x[3]=coul_coeff;
  x[4]=n1;
  x[5]=n0;
  x[6]=Epair;
  return 0;
}

double nucmass_ldrop_pair::drip_binding_energy_d
(double Z, double N, double npout, double nnout, double chi, double T) {
  
  double A=(Z+N);
  
  pair=-Epair*(cos(Z*o2scl_const::pi)+cos(N*o2scl_const::pi))/
    2.0/pow(A,1.5);
  
  return A*pair+nucmass_ldrop_skin::drip_binding_energy_d
    (Z,N,npout,nnout,chi,T);
}
