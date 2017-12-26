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

#include <cstdlib>

#include <o2scl/boson_eff.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

// Constructor, Destructor
boson_eff::boson_eff() {
  
  density_mroot=&def_density_mroot;
  psi_root=&def_psi_root;
  
  load_coefficients(cf_bosejel34);
}

boson_eff::~boson_eff() {
}

// Load coefficients for finite-temperature approximations
int boson_eff::load_coefficients(int ctype) {
  
  if (ctype==cf_bosejel21) {
    parma=0.978;
    sizem=2;
    sizen=1;
    Pmnb.resize(sizem+1,sizen+1);
    Pmnb(0,0)=1.63146;
    Pmnb(0,1)=2.11571;
    Pmnb(1,0)=3.0*Pmnb(0,0);
    Pmnb(1,1)=3.0*Pmnb(0,1);
    Pmnb(2,0)=3.31275;
    Pmnb(2,1)=5.15372;
  } else if (ctype==cf_bosejel22) {
    parma=0.914;
    sizem=2;
    sizen=2;
    Pmnb.resize(sizem+1,sizen+1);
    Pmnb(0,0)=1.68131;
    Pmnb(0,1)=3.47558;
    Pmnb(0,2)=2.16582;
    Pmnb(1,0)=3.0*Pmnb(0,0);
    Pmnb(1,1)=3.0*Pmnb(0,1);
    Pmnb(1,2)=3.0*Pmnb(0,2);
    Pmnb(2,0)=3.25053;
    Pmnb(2,1)=7.82859;
    Pmnb(2,1)=5.19126;
  } else if (ctype==cf_bosejel34) {
    parma=1.029;
    sizem=3;
    sizen=4;
    Pmnb.resize(sizem+1,sizen+1);
    Pmnb(0,0)=1.68134;
    Pmnb(0,1)=6.85070;
    Pmnb(0,2)=10.8537;
    Pmnb(0,3)=7.81843;
    Pmnb(0,4)=2.16461;
    Pmnb(1,0)=4.0*Pmnb(0,0);
    Pmnb(1,1)=4.0*Pmnb(0,1);
    Pmnb(1,2)=4.0*Pmnb(0,2);
    Pmnb(1,3)=4.0*Pmnb(0,3);
    Pmnb(1,4)=4.0*Pmnb(0,4);
    Pmnb(2,0)=8.49651;
    Pmnb(2,1)=35.6058;
    Pmnb(2,2)=57.7134;
    Pmnb(2,3)=42.3593;
    Pmnb(2,4)=11.8199;
    Pmnb(3,0)=3.45614;
    Pmnb(3,1)=15.1152;
    Pmnb(3,2)=25.5254;
    Pmnb(3,3)=19.2745;
    Pmnb(3,4)=5.51757;
  } else if (ctype==cf_bosejel34cons) {
    parma=1.040;
    sizem=3;
    sizen=4;
    Pmnb.resize(sizem+1,sizen+1);
    Pmnb(0,0)=1.68130;
    Pmnb(0,1)=6.85060;
    Pmnb(0,2)=10.8539;
    Pmnb(0,3)=7.81762;
    Pmnb(0,4)=2.16465;
    Pmnb(1,0)=4.0*Pmnb(0,0);
    Pmnb(1,1)=4.0*Pmnb(0,1);
    Pmnb(1,2)=4.0*Pmnb(0,2);
    Pmnb(1,3)=4.0*Pmnb(0,3);
    Pmnb(1,4)=4.0*Pmnb(0,4);
    Pmnb(2,0)=8.51373;
    Pmnb(2,1)=35.6576;
    Pmnb(2,2)=57.7975;
    Pmnb(2,3)=42.4049;
    Pmnb(2,4)=11.8321;
    Pmnb(3,0)=3.47433;
    Pmnb(3,1)=15.1995;
    Pmnb(3,2)=25.6536;
    Pmnb(3,3)=19.3811;
    Pmnb(3,4)=5.54423;
  } else {
    O2SCL_ERR("Invalid type in boson_eff::load_coefficients().",
		  exc_efailed);
  }
  
  return success;
}

void boson_eff::calc_mu(boson &b, double temper) {
  int nn, mm, nvar=1, ret=0;
  double xx[2], pren, preu, prep, sumn, sumu, sump;
  double gg, opg, nc, h, oph, psi;

  if (b.non_interacting) { b.nu=b.mu; b.ms=b.m; }

  // Massless boson_effs
  if (b.ms==0.0) {
    b.massless_calc(temper);
    return;
  }

  psi=(b.nu-b.ms)/temper;
  if (psi>=0.0) {
    h=0.0;
  } else {
    if (psi>-0.05) {
      xx[0]=sqrt(-2.0*parma*psi);
    } else if (psi>-1.0) {
      xx[0]=-sqrt(parma)*(-3.0+4.0*psi+4.0*log(2.0));
    } else {
      xx[0]=sqrt(parma)*exp(1.0-psi);
    }
    
    funct mfs=std::bind(std::mem_fn<double(double,double)>
			  (&boson_eff::solve_fun),
			  this,std::placeholders::_1,psi);
    
    int psi_root_err=psi_root->solve(xx[0],mfs);
    if (psi_root_err!=0) {
      O2SCL_ERR("psi_root failed in nepn_mroot().",psi_root_err);
    }
    h=xx[0];
  }

  oph=1.0+h;
  gg=temper/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);
  
  preu=pow(gg,2.5)/pow(oph,(double)sizem+1)/pow(opg,sizen-1.5);
  prep=preu;
  pren=pow(sqrt(parma)+h,2.0)*
    pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
    
  sumn=0.0;
  sumu=0.0;
  sump=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      if (h!=0.0 || mm>=2) {
	sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	  (h*(sizem+1.0-mm)-mm);
      }
      sumu+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn)*
	(1.5+nn+gg/opg*(1.5-sizen));
      sump+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn);
    }
  }
  if (h==0.0) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(1,nn)*pow(gg,(double)nn)*sizem;
    }
  }
    
  b.n=b.g/2.0*pren*sumn*nc;
  b.pr=b.g/2.0*prep*sump*nc*b.ms;
  b.ed=b.g/2.0*preu*sumu*nc*b.ms+b.n*b.ms;
  b.en=(b.ed+b.pr-b.mu*b.n)/temper;
  
  return;
}

double boson_eff::solve_fun(double x, double psi) {
  double h, sqt, y;

  h=x;
  sqt=sqrt(parma)+h;
  y=(h/sqt-log(sqt/sqrt(parma)))/psi-1.0;

  return y;
}

void boson_eff::calc_density(boson &b, double temper) {
  double h, oph, gg, opg, nc, preu, prep, pren;
  double sumn, sumu, sump, sqt, psi;
  int ret=0, mm, nn;
  ubvector xx(2);

  if (b.non_interacting) { b.ms=b.m; b.nu=b.mu; }

  fix_density=b.n;

  psi=(b.nu-b.ms)/temper;
  
  // If psi is too small, then we won't be able to solve for the 
  // right density, so we try to repair this here. It is possible
  // that -20 isn't large enough. This should be examined.
  if (psi<-20.0) psi=-20.0;
  
  mm_funct mfd=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,boson &,double)>
     (&boson_eff::density_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(b),temper);

  if (psi>=0.0) {
    h=0.0;
  } else {
    if (psi>-0.05) {
      xx[0]=sqrt(-2.0*parma*psi);
    } else if (psi>-1.0) {
      xx[0]=-sqrt(parma)*(-3.0+4.0*psi+4.0*log(2.0));
    } else {
      xx[0]=sqrt(parma)*exp(1.0-psi);
    }
    ret=density_mroot->msolve(1,xx,mfd);
    if (ret!=0) {
      O2SCL_ERR("Solver failed in boson_eff::calc_density().",
		exc_efailed);
    }
    h=xx[0];
  }
  if (!std::isfinite(h)) {
    O2SCL_ERR("Variable h not finite in calc_density_mroot().",exc_einval);
  }
  sqt=sqrt(parma)+h;
  b.nu=(h/sqt-log(sqt/sqrt(parma)))*temper+b.ms;

  if (b.non_interacting) { b.mu=b.nu; }

  oph=1.0+h;
  gg=temper/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  preu=pow(gg,2.5)/pow(oph,(double)sizem+1)/pow(opg,sizen-1.5);
  prep=preu;
  pren=pow(sqrt(parma)+h,2.0)*
    pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
    
  sumn=0.0;
  sumu=0.0;
  sump=0.0;

  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      if (h!=0.0 || mm>=2) {
	sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	  (h*(sizem+1.0-mm)-mm);
      }
      sumu+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn)*
	(1.5+nn+gg/opg*(1.5-sizen));
      sump+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn);
    }
  }

  if (h==0.0) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(1,nn)*pow(gg,(double)nn)*sizem;
    }
  }

  b.n=b.g/2.0*pren*sumn*nc;
  b.pr=b.g/2.0*prep*sump*nc*b.ms;
  b.ed=b.g/2.0*preu*sumu*nc*b.ms+b.n*b.ms;
  b.en=(b.ed+b.pr-b.mu*b.n)/temper;

  return;
}

int boson_eff::density_fun(size_t nv, const ubvector &x, 
			   ubvector &y, boson &b, double T) {
  double h,gg,opg,nc,oph,sumn,pren;
  int mm, nn;

  h=x[0];

  oph=1.0+h;
  gg=T/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  pren=pow(sqrt(parma)+h,2.0)*
    pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
  sumn=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	(h*(sizem+1.0-mm)-mm);
    }
  }

  y[0]=b.g/2.0*pren*sumn*nc/fix_density-1.0;
  if (!std::isfinite(y[0])) {
    O2SCL_ERR("Not finite in boson_eff::density_fun().",
	      exc_efailed);
  }

  return success;
}

void boson_eff::pair_mu(boson &b, double temper) {

  if (b.non_interacting) { b.nu=b.mu; b.ms=b.m; }
  calc_mu(b,temper);

  boson antip(b.ms,b.g);
  b.anti(antip);
  
  if (b.non_interacting) { antip.nu=antip.mu; antip.ms=antip.m; }
  calc_mu(antip,temper);
  b.n-=antip.n;
  b.pr+=antip.pr;
  b.ed+=antip.ed;
  b.en+=antip.en;

  return;
}

void boson_eff::pair_density(boson &b, double temper) {
  double oph, opg, gg, h, pren, preu, prep, sumn;
  double sumu, sump, nc, sqt, psi;
  int mm, nn;
  ubvector xx(2);

  if (b.non_interacting) { b.ms=b.m; b.nu=b.mu; }

  if (b.ms==0.0) return b.massless_calc(temper);

  fix_density=b.n;

  psi=(b.nu-b.ms)/temper;
  if (psi>-0.05) {
    xx[0]=sqrt(-2.0*parma*psi);
  } else if (psi>-1.0) {
    xx[0]=-sqrt(parma)*(-3.0+4.0*psi+4.0*log(2.0));
  } else {
    xx[0]=sqrt(parma)*exp(1.0-psi);
  }
  psi=(-b.nu-b.ms)/temper;

  mm_funct mfd=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,boson &,double)>
     (&boson_eff::pair_density_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(b),temper);

  if (psi>-0.05) {
    xx[1]=sqrt(-2.0*parma*psi);
  } else if (psi>-1.0) {
    xx[1]=-sqrt(parma)*(-3.0+4.0*psi+4.0*log(2.0));
  } else {
    xx[1]=sqrt(parma)*exp(1.0-psi);
  }
  int density_mroot_err=0;
  //int density_mroot_err=density_mroot->msolve(2,xx,this,mfd);
  if (density_mroot_err!=0) {
    O2SCL_ERR("mroot failed in pair_density_mroot().",density_mroot_err);
    // We continue execution anyway
  }
  h=xx[1];
  sqt=sqrt(parma)+h;
  b.nu=(h/sqt-log(sqt/sqrt(parma)))*temper+b.ms;

  if (b.non_interacting) { b.mu=b.nu; }

  oph=1.0+h;
  gg=temper/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  preu=pow(gg,2.5)/pow(oph,(double)sizem+1)/pow(opg,sizen-1.5);
  prep=preu;
  pren=pow(sqrt(parma)+h,2.0)*
    pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
    
  sumn=0.0;
  sumu=0.0;
  sump=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	(h*(sizem+1.0-mm)-mm);
      sumu+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn)*
	(1.5+nn+gg/opg*(1.5-sizen));
      sump+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn);
    }
  }
  b.n=-b.g/2.0*pren*sumn*nc;
  b.pr=b.g/2.0*prep*sump*nc*b.ms;
  b.ed=b.g/2.0*preu*sumu*nc*b.ms+b.n*b.ms;

  h=xx[0];
  sqt=sqrt(parma)+h;
  b.nu=(h/sqt-log(sqt/sqrt(parma)))*temper+b.ms;

  if (b.non_interacting) { b.mu=b.nu; }

  oph=1.0+h;
  gg=temper/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  preu=pow(gg,2.5)/pow(oph,(double)sizem+1)/pow(opg,sizen-1.5);
  prep=preu;
  pren=pow(sqrt(parma)+h,2.0)*
    pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
    
  sumn=0.0;
  sumu=0.0;
  sump=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	(h*(sizem+1.0-mm)-mm);
      sumu+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn)*
	(1.5+nn+gg/opg*(1.5-sizen));
      sump+=Pmnb(mm,nn)*pow(h,(double)mm)*pow(gg,(double)nn);
    }
  }
  b.n+=b.g/2.0*pren*sumn*nc;
  b.pr+=b.g/2.0*prep*sump*nc*b.ms;
  b.ed+=b.g/2.0*preu*sumu*nc*b.ms+b.n*b.ms;
  b.en=(b.ed+b.pr-b.mu*b.n)/temper;

  return;
}

int boson_eff::pair_density_fun(size_t nv, const ubvector &x, 
				ubvector &y, boson &b, double T) {

  double h,gg,opg,nc,oph,sumn,pren,sqt,psi;
  int mm, nn;

  if (x[0]<0.0 || x[1]<0.0) return exc_einval;

  // Solve for particles:
  h=x[0];
  sqt=sqrt(parma)+h;
  if (sqt<0.0) return exc_einval;
  psi=h/sqt-log(sqt/sqrt(parma));
  y[0]=psi;

  oph=1.0+h;
  gg=T/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  pren=sqt*sqt*pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
  sumn=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	(h*(sizem+1.0-mm)-mm);
    }
  }

  y[1]=b.g/2.0*pren*sumn*nc;
  
  // Solve for anti-particles:
  h=x[1];
  sqt=sqrt(parma)+h;
  if (sqt<0.0) return exc_einval;
  psi=h/sqt-log(sqt/sqrt(parma));
  y[0]+=psi;

  oph=1.0+h;
  gg=T/b.ms;
  opg=1.0+gg;
  nc=1.0/pi2*pow(b.ms,3.0);

  pren=sqt*sqt*pow(gg,1.5)/pow(oph,sizem+2.0)/pow(opg,sizen-1.5);
  sumn=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnb(mm,nn)*pow(h,(double)mm-2.0)*pow(gg,(double)nn)*
	(h*(sizem+1.0-mm)-mm);
    }
  }
  
  y[1]-=b.g/2.0*pren*sumn*nc;

  y[0]=y[0]/(2.0*b.m/T)-1.0;
  y[1]=y[1]/fix_density-1.0;
  //  if (!std::isfinite(y[1]) || !std::isfinite(y[2])) return 1;

  return success;
}

