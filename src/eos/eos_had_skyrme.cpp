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

#include <o2scl/eos_had_skyrme.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_skyrme::eos_had_skyrme() {

  parent_method=false;

  neutron->init(939.0/hc_mev_fm,2.0);
  proton->init(939.0/hc_mev_fm,2.0);

  fet=&nrf;
}

void eos_had_skyrme::hamiltonian_coeffs(double &ham1, double &ham2,
					double &ham3, double &ham4,
					double &ham5, double &ham6) {
  
  ham1=0.5*t0*(1.0+0.5*x0);
  ham2=-0.5*t0*(0.5+x0);
  ham3=a*t3/6.0*(1.0+0.5*x3);
  ham4=a*t3*pow(2.0,alpha-2.0)/6.0*(1.0-x3);
  ham5=b*t3/12.0*(1.0+0.5*x3);
  ham6=-b*t3/12.0*(0.5+x3);
  
  return;
}

int eos_had_skyrme::calc_deriv_temp_e(fermion_deriv &ne, fermion_deriv &pr,
				      double ltemper, thermo &locth,
				      thermo_np_deriv_helm &locthd) {
  
#if !O2SCL_NO_RANGE_CHECK
  check_input(ne,pr,ltemper);
#endif

  double term, term2;
  eff_mass(ne,pr,term,term2);

  if (ne.ms<0.0 || pr.ms<0.0) {
    string err=((string)"Effective masses negative (m^{*}_n=")+
      o2scl::dtos(ne.ms)+" and m^{*}_p="+o2scl::dtos(pr.ms)+
      ") at densities n_n="+o2scl::dtos(ne.n)+" and n_p="+
      o2scl::dtos(pr.n)+" in eos_had_skyrme::calc_deriv_temp_e().";
    std::cout << err << std::endl;
    O2SCL_CONV_RET(err.c_str(),exc_einval,this->err_nonconv);
  }

  // See note in class documentation about zero density
  if (ltemper>0.0 && ne.n==0.0) {
    ne.nu=-std::numeric_limits<double>::infinity();
    ne.ed=0.0;
    ne.pr=0.0;
    ne.en=0.0;
  } else {
    nrfd.calc_density(ne,ltemper);
  }
  if (ltemper>0.0 && pr.n==0.0) {
    pr.nu=-std::numeric_limits<double>::infinity();
    pr.ed=0.0;
    pr.pr=0.0;
    pr.en=0.0;
  } else {
    nrfd.calc_density(pr,ltemper);
  }

  // Compute the coefficients of different powers of density
  // in the hamiltonian
  double ham1, ham2, ham3, ham4, ham5, ham6;
  hamiltonian_coeffs(ham1,ham2,ham3,ham4,ham5,ham6);
  
  // Compute the base thermodynamic properties
  base_thermo(ne,pr,ltemper,locth,term,term2,
	      ham1,ham2,ham3,ham4,ham5,ham6);
  
  // Compute the second derivatives
  second_deriv(ne,pr,ltemper,locth,locthd,term,term2,
	       ham1,ham2,ham3,ham4,ham5,ham6);

  return success;
}

int eos_had_skyrme::calc_temp_e(fermion &ne, fermion &pr, 
				double ltemper, thermo &locth) {

#if !O2SCL_NO_RANGE_CHECK
  check_input(ne,pr,ltemper);
#endif

  double term, term2;
  eff_mass(ne,pr,term,term2);

  if (ne.ms<0.0 || pr.ms<0.0) {
    string err=((string)"Effective masses negative (m^{*}_n=")+
      o2scl::dtos(ne.ms)+" and m^{*}_p="+o2scl::dtos(pr.ms)+
      ") at densities n_n="+o2scl::dtos(ne.n)+" and n_p="+
      o2scl::dtos(pr.n)+" in eos_had_skyrme::calc_temp_e().";
    std::cout << err << std::endl;
    O2SCL_CONV_RET(err.c_str(),exc_einval,this->err_nonconv);
  }

  // See note in class documentation about zero density
  if (ltemper>0.0 && ne.n==0.0) {
    if (ne.inc_rest_mass) {
      ne.nu=ne.m;
    } else {
      ne.nu=0.0;
    }
    ne.ed=0.0;
    ne.pr=0.0;
    ne.en=0.0;
  } else {
    nrf.calc_density(ne,ltemper);
  }
  if (ltemper>0.0 && pr.n==0.0) {
    if (pr.inc_rest_mass) {
      pr.nu=pr.m;
    } else {
      pr.nu=0.0;
    }
    pr.ed=0.0;
    pr.pr=0.0;
    pr.en=0.0;
  } else {
    nrf.calc_density(pr,ltemper);
  }

  // Compute the coefficients of different powers of density
  // in the hamiltonian
  double ham1, ham2, ham3, ham4, ham5, ham6;
  hamiltonian_coeffs(ham1,ham2,ham3,ham4,ham5,ham6);
  
  // Compute the base thermodynamic properties
  base_thermo(ne,pr,ltemper,locth,term,term2,
	      ham1,ham2,ham3,ham4,ham5,ham6);
  
  return success;
}

int eos_had_skyrme::calc_e(fermion &ne, fermion &pr, thermo &locth) {
  return calc_temp_e(ne,pr,0.0,locth);
}

int eos_had_skyrme::calc_deriv_e(fermion_deriv &ne, fermion_deriv &pr,
				 thermo &locth,
				 thermo_np_deriv_helm &thd) {
  return calc_deriv_temp_e(ne,pr,0.0,locth,thd);
}

double eos_had_skyrme::feoa_symm(double nb) {
  double ret, kr23, beta, t3p;

  if (parent_method) {
    return eos_had_base::feoa(nb);
  }
  t3p=(a+b)*t3;
  kr23=0.6/(neutron->m+proton->m)*pow(1.5*pi2*nb,2.0/3.0);
  beta=(neutron->m+proton->m)/4.0*(0.25*(3.0*t1+5.0*t2)+t2*x2);
  ret=kr23*(1.0+beta*nb)+0.375*t0*nb+0.0625*t3p*pow(nb,1.0+alpha);

  return ret;
}

double eos_had_skyrme::fmsom_symm(double nb) {
  double ret, beta;

  if (parent_method) {
    return eos_had_base::fmsom(nb);
  }
  beta=(neutron->m+proton->m)/4.0*(0.25*(3.0*t1+5.0*t2)+t2*x2);
  ret=1.0/(1.0+beta*nb);

  return ret;
}

double eos_had_skyrme::fcomp_nuc(double nb) {
  double ret, kr23, beta, t3p;

  if (parent_method) {
    return eos_had_base::fcomp(nb);
  }
  t3p=(a+b)*t3;
  kr23=0.6/(neutron->m+proton->m)*pow(1.5*pi2*nb,2.0/3.0);
  beta=(neutron->m+proton->m)/4.0*(0.25*(3.0*t1+5.0*t2)+t2*x2);

  // This only works at saturation density:
  //  ret=-2.0*kr23+10.0*kr23*beta*nb+
  //    9.0/16.0*alpha*(alpha+1.0)*t3p*pow(nb,1.0+alpha);

  ret=10.0*kr23+27.0*nb*t0/4.0+40.0*kr23*beta*nb+
    9.0/16.0*alpha*(alpha+1.0)*t3p*pow(nb,1.0+alpha)+
    9.0/8.0*t3p*(1.0+alpha)*pow(nb,1.0+alpha);
  
  return ret;
}

double eos_had_skyrme::fesym(double nb, double pf) {
  double ret, kr23;

  if (pf!=0.5 || parent_method) {
    return eos_had_base::fesym(nb,pf);
  }
  kr23=0.6/(neutron->m+proton->m)*pow(1.5*pi2*nb,2.0/3.0);
  ret=5.0/9.0*kr23+10.0/6.0*(neutron->m+proton->m)*kr23*nb*
    (t2/6.0*(1.0+1.25*x2)-0.125*t1*x1)-
    b*t3/24.0*(0.5+x3)*pow(nb,1.0+alpha)-
    0.25*t0*(0.5+x0)*nb-a/96*pow(nb,1.0+alpha)*t3*
    (2.0-alpha*(3.0+alpha)+x3*(4.0+alpha*(3.0+alpha)));
  
  return ret;
}

double eos_had_skyrme::fkprime_nuc(double nb) {
  double ret, kr23, t3p, beta, lmsom;

  if (parent_method) {
    return eos_had_base::fkprime(nb);
  }
  t3p=(a+b)*t3;
  kr23=0.6/(neutron->m+proton->m)*pow(1.5*pi2*nb,2.0/3.0);
  beta=0.5*(neutron->m+proton->m)/2.0*(0.25*(3.0*t1+5.0*t2)+t2*x2);
  lmsom=1.0/(1.0+beta*nb);
  ret=2.0*kr23*(9.0-5.0/lmsom)+
    t3p*pow(nb,1.0+alpha)*27.0/16.0*alpha*(alpha*alpha-1.0);

  return ret;
}

int eos_had_skyrme::calpar(double gt0, double gt3, double galpha,
			   double gt1, double gt2) {

  ubvector x(3);

  fixn0=n0;
  fixmsom=msom;
  fixeoa=eoa;
  fixcomp=comp;
  fixesym=esym;

  x[0]=gt0;
  x[1]=gt3;
  x[2]=galpha;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_skyrme::calparfun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  mm_funct fmf2=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_skyrme::calparfun2),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  if (eos_mroot->msolve(3,x,fmf)!=0) {
    O2SCL_CONV_RET("Solution failed in calparfun().",exc_efailed,
                   this->err_nonconv);
  }
  t0=x[0];
  t3=x[1];
  alpha=x[2];
  
  x[0]=gt1;
  x[1]=gt2;
  if (eos_mroot->msolve(2,x,fmf2)!=0) {
    O2SCL_CONV_RET("Solution failed in calparfun2().",exc_efailed,
                   this->err_nonconv);
  }
  t1=x[0];
  t2=x[1];

  return success;
}

int eos_had_skyrme::calparfun(size_t nv, const ubvector &x, ubvector &y) {
  double t3p, kr23;
  double pres, beta, mnuc=(neutron->m+proton->m)/2.0;

  t0=x[0];
  t3=x[1];
  alpha=x[2];

  t3p=(a+b)*t3;
  kr23=0.3/mnuc*pow(1.5*pi2*fixn0,2.0/3.0);
  beta=(1.0/fixmsom-1.0)/fixn0;
  pres=kr23*(fixn0*fixn0*beta+2.0/3.0*fixn0/fixmsom)+0.375*t0*fixn0*fixn0
    +0.0625*t3p*(1.0+alpha)*pow(fixn0,2.0+alpha);
  eoa=kr23*(1.0+beta*fixn0)+0.375*t0*fixn0+0.0625*t3p*pow(fixn0,1.0+alpha);
  comp=-2.0*kr23+10.0*kr23*beta*fixn0+
    9.0/16.0*alpha*(alpha+1.0)*t3p*pow(fixn0,1.0+alpha);

  y[0]=pres;
  y[1]=eoa-fixeoa;
  y[2]=comp-fixcomp;
    
  return 0;
}

int eos_had_skyrme::calparfun2(size_t nv, const ubvector &x, 
			       ubvector &y) {
  double kr23, beta;
  double mnuc=(neutron->m+proton->m)/2.0;
  
  t1=x[0];
  t2=x[1];

  kr23=0.3/mnuc*pow(1.5*pi2*fixn0,2.0/3.0);
  beta=0.5*mnuc*(0.25*(3.0*t1+5.0*t2)+t2*x2);
  msom=1.0/(1.0+beta*fixn0);
  esym=5.0/9.0*kr23+10.0/3.0*mnuc*kr23*fixn0*
    (1.0/6.0*t2*(1.0+1.25*x2)-0.125*t1*x1)-
    1.0/24.0*b*t3*(0.5+x3)*pow(fixn0,1.0+alpha)-
    0.25*t0*(0.5+x0)*fixn0-a/96*pow(fixn0,1.0+alpha)*t3*
    (2.0-alpha*(3.0+alpha)+x3*(4.0+alpha*(3.0+alpha)));
  
  y[0]=msom-fixmsom;
  y[1]=esym-fixesym;
  
  return 0;
}

int eos_had_skyrme::check_landau(double nb, double m) {
  double f0, g0, f0p, g0p, f1, g1, f1p, g1p;
  
  landau_nuclear(nb,m,f0,g0,f0p,g0p,f1,g1,f1p,g1p);
  if (f0<-1.0) return 1;
  if (g0<-1.0) return 2;
  if (f0p<-1.0) return 3;
  if (g0p<-1.0) return 4;
  if (f1<-3.0) return 5;
  if (g1<-3.0) return 6;
  if (f1p<-3.0) return 7;
  if (g1p<-3.0) return 8;
  landau_neutron(nb,m,f0,g0,f1,g1);
  if (f0<-1.0) return 9;
  if (g0<-1.0) return 10;
  if (f1<-3.0) return 11;
  if (g1<-3.0) return 12;

  return 0;
}

void eos_had_skyrme::landau_nuclear
(double nb, double m, double &f0, double &g0, double &f0p,
 double &g0p, double &f1, double &g1, double &f1p, double &g1p) {
  
  double T0, T1, T2, T3, x, y, z, kf, mstar;
  
  kf=pow(3.0*pi2*nb/2.0,1.0/3.0);
  x=t1*x1;
  y=t2*x2;
  z=t3*x3;
  T0=0.125*(3.0*t1+5.0*t2+4.0*y);
  T1=0.125*(2.0*x+2.0*y-t1+t2);
  T2=0.125*(2.0*x-2.0*y+t1-t2);
  T3=0.125*(t1-t2);
  mstar=1.0/(1.0/m+T0*nb);

  f1=-3.0*T0*mstar*nb;
  g1=-3.0*T1*mstar*nb;
  f1p=3.0*T2*mstar*nb;
  g1p=3.0*T3*mstar*nb;
  f0=(0.75*t0+0.0625*(alpha+1.0)*(alpha+2.0)*t3*pow(nb,alpha))*
    2.0*mstar*kf/pi2-f1;
  g0=(0.25*t0*(2.0*x0-1.0)+1.0/24.0*t3*pow(nb,alpha)*(2.0*x3-1.0))*
    2.0*mstar*kf/pi2-g1;
  f0p=(-0.25*t0*(2.0*x0+1.0)-1.0/24.0*t3*pow(nb,alpha)*(2.0*x3+1.0))*
    2.0*mstar*kf/pi2-f1p;
  g0p=(-0.25*t0-1.0/24.0*t3*pow(nb,alpha))*2.0*mstar*kf/pi2-g1p;
  
  return;
}

void eos_had_skyrme::landau_neutron
(double nb, double m, double &f0, double &g0, double &f1, double &g1) {
  
  double T0, T1, T2, T3, x, y, z, kf, mstar;
  
  kf=pow(3.0*pi2*nb,1.0/3.0);
  x=t1*x1;
  y=t2*x2;
  z=t3*x3;
  T0=0.125*(3.0*t1+5.0*t2+4.0*y);
  T1=0.125*(2.0*x+2.0*y-t1+t2);
  T2=0.125*(2.0*x-2.0*y+t1-t2);
  T3=0.125*(t1-t2);
  mstar=1.0/(1.0/m+(T0-T2)*nb);

  f1=-3.0*(T0-T2)*mstar*nb;
  g1=-3.0*(T1-T3)*mstar*nb;
  
  f0=(0.5*t0*(1.0-x0)+1.0/24.0*(alpha+1.0)*(alpha+2.0)*t3*pow(nb,alpha)*
      (1.0-x3))*mstar*kf/pi2-f1;
  g0=(0.5*t0*(x0-1.0)+1.0/12.0*t3*pow(nb,alpha)*(x3-1.0))*mstar*kf/pi2-g1;
  
  return;
}

void eos_had_skyrme::alt_params_saturation
(double n0t, double EoA, double K, double Ms_star, double a, double L,
 double Mv_star, double CrDr0, double CrDr1, double CrnJ0, double CrnJ1) {

  // This quantity has units of fm
  double h2o2m=1.0/(def_neutron.m+def_proton.m);
  // This quantity is unitless
  double C_k=0.6*pow(1.5*pi2,2.0/3.0);
  // This quantity has units of fm^{-2}
  double tau_c=C_k*pow(n0t,2.0/3.0);
  // This quantity is unitless
  double alpha2=(tau_c*(4.0/Ms_star-3.0)*h2o2m-K-9.0*EoA)/
    (tau_c*(6.0/Ms_star-9.0)*h2o2m+9.0*EoA);
  // This quantity has units of fm^2
  double Crr00=(((2.0-3.0*alpha2)/Ms_star-3.0)*tau_c*h2o2m+
		3.0*(1.0+alpha2)*EoA)/(3.0*alpha2*n0t);
  // This quantity has units of fm^{3*alpha+2}, where
  // alpha is referred to as "gamma" in Kortelainen et al. (2010)
  double Crr0D=((3.0-2.0/Ms_star)*tau_c*h2o2m-3.0*EoA)/
    (3.0*alpha2*pow(n0t,1.0+alpha2));
  // This quantity has units of fm^4
  double Crt0=(1.0/Ms_star-1.0)/n0t*h2o2m;
  // This quantity has units of fm^4
  double Crt1=Crt0-(1.0/Mv_star-1.0)/n0t*h2o2m;
  // This quantity has units of fm^2
  double Crr10=(27.0*(1.0+alpha2)*a-9.0*L+5.0*tau_c*(2.0-3.0*alpha2)*
		(Crt0+3.0*Crt1)*n0t-5.0*tau_c*(1.0+3.0*alpha2)*h2o2m)/
    (27.0*alpha2*n0t);
  // This quantity has units of fm^{3*alpha+2}, where
  // alpha is referred to as "gamma" in Kortelainen et al. (2010)
  double Crr1D=(-27.0*a+9.0*L+5.0*(h2o2m-2.0*n0t*(Crt0+3.0*Crt1))*
		tau_c)/(27.0*alpha2*pow(n0t,1.0+alpha2));
  
  alt_params_set(Crr00,Crr10,Crr0D,Crr1D,Crt0,Crt1,CrDr0,CrDr1,
		 CrnJ0,CrnJ1,alpha2);

  return;
}

void eos_had_skyrme::alt_params_set
(double Crr00, double Crr10, double Crr0D, double Crr1D, double Crt0,
 double Crt1, double CrDr0, double CrDr1, double CrnJ0, double CrnJ1,
 double alpha2) {
  a=0.0;
  b=1.0;
  t0=8.0/3.0*Crr00;
  x0=(-Crr00-3.0*Crr10)/2.0/Crr00;
  t3=16.0*Crr0D;
  x3=(-Crr0D-3.0*Crr1D)/2.0/Crr0D;
  t1=-4.0/3.0*(4.0*CrDr0-Crt0);
  x1=(3.0*Crt1+Crt0-4.0*CrDr0-12.0*CrDr1)/2.0/(4.0*CrDr0-Crt0);
  t2=4.0/3.0*(4.0*CrDr0-8.0*CrDr1+3.0*Crt0-6.0*Crt1);
  x2=(20.0*CrDr1+15.0*Crt1-3.0*Crt0-4.0*CrDr0)/2.0/
    (4.0*CrDr0-8.0*CrDr1+3.0*Crt0-6.0*Crt1);
  b4=-CrnJ0+CrnJ1;
  b4p=-2.0*CrnJ1;
  alpha=alpha2;
  return;
}

void eos_had_skyrme::alt_params_get
(double &Crr00, double &Crr10, double &Crr0D, double &Crr1D, double &Crt0,
 double &Crt1, double &CrDr0, double &CrDr1, double &CrnJ0, double &CrnJ1,
 double &alpha2) {
  Crr00=0.375*t0;
  Crr10=-0.25*t0*(0.5+x0);
  Crr0D=0.0625*t3;
  Crr1D=-t3/24.0*(0.5+x3);
  Crt0=0.1875*t1+0.25*t2*(1.25+x2);
  Crt1=-0.125*t1*(0.5+x1)+0.125*t2*(0.5+x2);
  CrDr0=-0.140625*t1+0.0625*t2*(1.25+x2);
  CrDr1=0.09375*t1*(0.5+x1)+0.03125*(0.5+x2);
  CrnJ0=-b4-b4p/2.0;
  CrnJ1=-b4p/2.0;
  alpha2=alpha;
  return;
}

/*
  int eos_had_skyrme::calpar_new(double m) {
  double T0=(1.0/msom/m-1.0/m)*n0;
  double kf=pow(3.0*pi2*n0/2.0,1.0/3.0);
  alpha=(comp/9.0+eoa+(0.1/m-2.0/15.0/msom/m)*kf*kf)/
  (-eoa+(0.3/m-0.2/msom/m)*kf*kf);
  t3=16.0/pow(n0,1.0+alpha)/alpha*(-eoa+(0.3/m-0.2/msom/m)*kf*kf);
  t0=8.0/3.0/n0*(eoa-0.3/msom/m*kf*kf-0.0625*t3*pow(n0,1.0+alpha));
  
  return 0;
  }
*/

