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

#include <o2scl/eos_had_skyrme.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_skyrme::eos_had_skyrme() {

  parent_method=false;
  mu_at_zero_density=false;

  neutron->init(939.0/hc_mev_fm,2.0);
  proton->init(939.0/hc_mev_fm,2.0);

  fet=&nrf;
}

int eos_had_skyrme::calc_temp_e(fermion &ne, fermion &pr, 
			    double ltemper, thermo &locth) {
  
  double n, x, hamk, ham, ham1, ham2, ham3, ham4, ham5, ham6;
  double dhdnn, dhdnp, na, npa, nna, term, term2, common, gn, gp;
 
#if !O2SCL_NO_RANGE_CHECK
  if (!o2scl::is_finite(ne.n) || !o2scl::is_finite(ne.n) ||
      !o2scl::is_finite(ltemper)) {
    O2SCL_ERR2("Nucleon densities or temperature not finite in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true) {
    O2SCL_ERR2("Neutron or protons non-interacting in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
#endif

  //---------------------------------------
  // Some local variables of interest:
  //
  // hamk is just the kinetic part of the hamiltonian 
  //   hbar^2 tau / (2 m^{star})
  // ham{1-6} are remaining parts of the hamiltonian, modulo
  //   factors of density
  // dhdn{n,p} are the partial derivatives of the hamiltonian wrt the 
  //   neutron and proton densities (hold energy densities constant)

  // If the temperature is too small, just use 
  // the zero-temperature code

  if (ne.n<=0.0 && pr.n<=0.0) {
    ne.ms=ne.m;
    pr.ms=pr.m;
    ne.mu=ne.m;
    pr.mu=pr.m;
    ne.pr=0.0;
    pr.pr=0.0;
    ne.ed=0.0;
    pr.ed=0.0;
    locth.pr=0.0;
    locth.ed=0.0;
    locth.en=0.0;
    return success;
  } else if (ne.n<=0.0) {
    ne.n=0.0;
  } else if (pr.n<=0.0) {
    pr.n=0.0;
  }

  if (ltemper<=0.0) {
    calc_e(ne,pr,locth);
    return 0;
  }

  n=ne.n+pr.n;
  x=pr.n/n;

  // Landau effective masses
  term=0.25*(t1*(1.0+x1/2.0)+t2*(1.0+x2/2.0));
  term2=0.25*(t2*(0.5+x2)-t1*(0.5+x1));
  ne.ms=ne.m/(1.0+2.0*(n*term+ne.n*term2)*ne.m);
  pr.ms=pr.m/(1.0+2.0*(n*term+pr.n*term2)*pr.m);
  
  nrf.calc_density(ne,ltemper);
  nrf.calc_density(pr,ltemper);

  // Single particle potentials and energy density

  na=pow(fabs(n),alpha);
  npa=pow(fabs(pr.n),alpha);
  nna=pow(fabs(ne.n),alpha);

  hamk=ne.ed+pr.ed;
  ham1=0.5*t0*(1.0+0.5*x0);
  ham2=-0.5*t0*(0.5+x0);
  ham3=a*t3/6.0*(1.0+0.5*x3);
  ham4=a*t3*pow(2.0,alpha)/96.0*(1.0-x3);
  ham5=b*t3/12.0*(1.0+0.5*x3);
  ham6=-b*t3/12.0*(0.5+x3);

  ham=hamk+ham1*n*n+ham2*(ne.n*ne.n+pr.n*pr.n)+
    ham3*na*ne.n*pr.n+ham4*(nna*ne.n*ne.n+npa*pr.n*pr.n)+
    ham5*n*n*na+ham6*(ne.n*ne.n+pr.n*pr.n)*na;
  
  if (ne.inc_rest_mass) {
    gn=2.0*ne.ms*(ne.ed-ne.n*ne.m);
  } else {
    gn=2.0*ne.ms*ne.ed;
  }
  if (pr.inc_rest_mass) {
    gp=2.0*pr.ms*(pr.ed-pr.n*pr.m);
  } else {
    gp=2.0*pr.ms*pr.ed;
  }
  common=(gn+gp)*term+2.0*ham1*n+ham5*(2.0+alpha)*n*na;
  if (!mu_at_zero_density && ne.n<=0.0) {
    dhdnn=0.0;
  } else {
    dhdnn=common+ne.nu+gn*term2+
      2.0*ham2*ne.n+ham3*na*pr.n*(alpha*ne.n/n+1.0)+
      ham4*(nna*ne.n*(2.0+alpha))+
      ham6*(2.0*ne.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/n);
  }
  if (!mu_at_zero_density && pr.n<=0.0) {
    dhdnp=0.0;
  } else {
    dhdnp=common+pr.nu+gp*term2+
      2.0*ham2*pr.n+ham3*na*ne.n*(alpha*pr.n/n+1.0)+
      ham4*(npa*pr.n*(2.0+alpha))+
      ham6*(2.0*pr.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/n);
  }
  
  ne.mu=dhdnn;
  pr.mu=dhdnp;
  
  // Thermodynamics
  locth.ed=ham;
  locth.en=ne.en+pr.en;
  locth.pr=ltemper*locth.en+ne.mu*ne.n+pr.mu*pr.n-locth.ed;

  return success;
}

int eos_had_skyrme::calc_e(fermion &ne, fermion &pr, thermo &locth) {

#if !O2SCL_NO_RANGE_CHECK
  if (!o2scl::is_finite(ne.n) || !o2scl::is_finite(ne.n)) {
    O2SCL_ERR2("Nucleon densities not finite in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true) {
    O2SCL_ERR2("Neutron or protons non-interacting in ",
	       "eos_had_skyrme::calc_eq_temp_p().",exc_einval);
  }
#endif

  double x, n, hamk, ham, ham1, ham2, ham3, ham4, ham5, ham6;
  double dhdnn, dhdnp, na, npa, nna, term, term2, common, gn, gp;

  ne.non_interacting=false;
  pr.non_interacting=false;

  //---------------------------------------
  // Some local variables of interest:
  //
  // hamk is just the kinetic part of the hamiltonian 
  //   hbar^2 tau / (2 m^{star})
  // ham{1-6} are remaining parts of the hamiltonian, modulo
  //   factors of density
  // dhdn{n,p} are the total derivatives of the hamiltonian wrt the 
  //   neutron and proton densities (takes into account chain rule
  //   contributions from energy density)

  if (ne.n<=0.0 && pr.n<=0.0) {
    ne.ms=ne.m;
    pr.ms=pr.m;
    ne.mu=ne.m;
    pr.mu=pr.m;
    ne.pr=0.0;
    pr.pr=0.0;
    ne.ed=0.0;
    pr.ed=0.0;
    locth.pr=0.0;
    locth.ed=0.0;
    locth.en=0.0;
    return success;
  } else if (ne.n<=0.0) {
    ne.n=0.0;
  } else if (pr.n<=0.0) {
    pr.n=0.0;
  }

  n=ne.n+pr.n;
  x=pr.n/n;

  na=pow(fabs(n),alpha);
  npa=pow(fabs(pr.n),alpha);
  nna=pow(fabs(ne.n),alpha);

  term=0.25*(t1*(1.0+x1/2.0)+t2*(1.0+x2/2.0));
  term2=0.25*(t2*(0.5+x2)-t1*(0.5+x1));
  ne.ms=ne.m/(1.0+2.0*(n*term+ne.n*term2)*ne.m);
  pr.ms=pr.m/(1.0+2.0*(n*term+pr.n*term2)*pr.m);

  // We don't record error values, since these functions usually
  // always succeed
  nrf.calc_density_zerot(ne);
  nrf.calc_density_zerot(pr);
  
  hamk=ne.ed+pr.ed;

  ham1=0.5*t0*(1.0+0.5*x0);
  ham2=-0.5*t0*(0.5+x0);
  ham3=a*t3/6.0*(1.0+0.5*x3);
  ham4=a*t3*pow(2.0,alpha-2.0)/6.0*(1.0-x3);
  ham5=b*t3/12.0*(1.0+0.5*x3);
  ham6=-b*t3/12.0*(0.5+x3);

  ham=hamk+ham1*n*n+ham2*(ne.n*ne.n+pr.n*pr.n)+
    ham3*na*ne.n*pr.n+ham4*(nna*ne.n*ne.n+npa*pr.n*pr.n)+
    ham5*n*n*na+ham6*(ne.n*ne.n+pr.n*pr.n)*na;

  if (ne.inc_rest_mass) {
    gn=2.0*ne.ms*(ne.ed-ne.n*ne.m);
  } else {
    gn=2.0*ne.ms*ne.ed;
  }
  if (pr.inc_rest_mass) {
    gp=2.0*pr.ms*(pr.ed-pr.n*pr.m);
  } else {
    gp=2.0*pr.ms*pr.ed;
  }
  common=(gn+gp)*term+2.0*ham1*n+ham5*(2.0+alpha)*n*na;
  if (!mu_at_zero_density && ne.n<=0.0) {
    dhdnn=0.0;
  } else {
    dhdnn=common+ne.nu+gn*term2+
      2.0*ham2*ne.n+ham3*na*pr.n*(alpha*ne.n/n+1.0)+
      ham4*(nna*ne.n*(2.0+alpha))+
      ham6*(2.0*ne.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/n);
  }
  if (!mu_at_zero_density && pr.n<=0.0) {
    dhdnp=0.0;
  } else {
    dhdnp=common+pr.nu+gp*term2+
      2.0*ham2*pr.n+ham3*na*ne.n*(alpha*pr.n/n+1.0)+
      ham4*(npa*pr.n*(2.0+alpha))+
      ham6*(2.0*pr.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/n);
  }
    
  ne.mu=dhdnn;
  pr.mu=dhdnp;
  
  locth.ed=ham;
  locth.pr=-locth.ed+ne.mu*ne.n+pr.mu*pr.n;
  locth.en=0.0;

  if (!o2scl::is_finite(locth.pr)) {
    O2SCL_ERR("Pressure not finite in calc_e()",exc_efailed);
  }

  return success;
}

double eos_had_skyrme::feoa(double nb) {
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

double eos_had_skyrme::fmsom(double nb) {
  double ret, beta;

  if (parent_method) {
    return eos_had_base::fmsom(nb);
  }
  beta=(neutron->m+proton->m)/4.0*(0.25*(3.0*t1+5.0*t2)+t2*x2);
  ret=1.0/(1.0+beta*nb);

  return ret;
}

double eos_had_skyrme::fcomp(double nb) {
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

double eos_had_skyrme::fkprime(double nb) {
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

  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_skyrme::calparfun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  mm_funct11 fmf2=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_skyrme::calparfun2),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  if (eos_mroot->msolve(3,x,fmf)!=0) {
    O2SCL_ERR("Solution failed in calparfun().",exc_efailed);
  }
  t0=x[0];
  t3=x[1];
  alpha=x[2];
  
  x[0]=gt1;
  x[1]=gt2;
  if (eos_mroot->msolve(2,x,fmf2)!=0) {
    O2SCL_ERR("Solution failed in calparfun2().",exc_efailed);
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
(double n0, double EoA, double K, double Ms_star, double a, double L,
 double Mv_star, double CrDr0, double CrDr1, double CrnJ0, double CrnJ1) {
  
  double h2o2m=20.73553/197.33;
  double C_k=0.6*pow(1.5*pi2,2.0/3.0);
  double tau_c=C_k*pow(n0,2.0/3.0);
  double alpha2=(tau_c*(4.0/Ms_star-3.0)*h2o2m-K-9.0*EoA)/
    (tau_c*(6.0/Ms_star-9.0)*h2o2m+9.0*EoA);
  double Crr00=(((2.0-3.0*alpha2)/Ms_star-3.0)*tau_c*h2o2m+
		3.0*(1.0+alpha2)*EoA)/(3.0*alpha2*n0);
  double Crr0D=((3.0-2.0/Ms_star)*tau_c*h2o2m-3.0*EoA)/
    (3.0*alpha2*pow(n0,1.0+alpha2));
  double Crt0=(1.0/Ms_star-1.0)/n0*h2o2m;
  double Crt1=Crt0-(1.0/Mv_star-1.0)/n0*h2o2m;
  double Crr10=(27.0*(1.0+alpha2)*a-9.0*L+5.0*tau_c*(2.0-3.0*alpha2)*
		(Crt0+3.0*Crt1)*n0-5.0*tau_c*(1.0+3.0*alpha2)*h2o2m)/
    (27.0*alpha2*n0);
  double Crr1D=(-27.0*a+9.0*L+5.0*(h2o2m-2.0*n0*(Crt0+3.0*Crt1))*
		tau_c)/(27.0*alpha2*pow(n0,1.0+alpha2));
  
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

