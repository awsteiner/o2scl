/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <o2scl/boson_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// boson_rel class

boson_rel::boson_rel() {
  density_root=&def_density_root;
  nit=&def_nit;
  dit=&def_dit;
}

boson_rel::~boson_rel() {
}

void boson_rel::set_inte(inte<> &l_nit, inte<> &l_dit) {
  nit=&l_nit;
  dit=&l_dit;
  return;
}

void boson_rel::calc_mu(boson &b, double temper) {

  if (temper<=0.0) {
    O2SCL_ERR2("Temperature less than or equal to zero in ",
	       "boson_rel::calc_mu().",exc_einval);
  }
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

  bool deg=true;
  double deg_limit=2.0;
  double psi;
  if (b.inc_rest_mass) {
    psi=(b.nu-b.ms)/temper;
  } else {
    psi=(b.nu+(b.m-b.ms))/temper;
  }
  if (psi<deg_limit) deg=false;
  
  if (deg) {
    
    // Compute the upper limit for degenerate integrals

    double arg, upper_limit_fac=20.0;
    if (b.inc_rest_mass) {
      arg=pow(upper_limit_fac*temper+b.nu,2.0)-b.ms*b.ms;
    } else {
      arg=pow(upper_limit_fac*temper+b.nu+b.m,2.0)-b.ms*b.ms;
    }
    double ul;
    if (arg>0.0) {
      ul=sqrt(arg);
    } else {
      O2SCL_ERR2("Zero density in degenerate limit in boson_rel::",
		 "calc_mu().",exc_efailed);
      return;
    }
    
    funct fd=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_density_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
    funct fe=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_energy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
    funct fs=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_entropy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
    
    b.n=dit->integ(fd,0.0,ul);
    b.n*=b.g/2.0/pi2;
    
    b.ed=dit->integ(fe,0.0,ul);
    b.ed*=b.g/2.0/pi2;
    
    b.en=dit->integ(fs,0.0,ul);
    b.en*=b.g/2.0/pi2;
    
  } else {
    
    // If the temperature is large enough, perform the full integral
    
    funct mfd=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::density_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
    funct mfe=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::energy_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
    funct mfs=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::entropy_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
      
    double prefac=b.g*pow(temper,3.0)/2.0/pi2;

    // Compute the number density
    
    b.n=nit->integ(mfd,0.0,0.0);
    b.n*=prefac;

    // Compute the energy density

    b.ed=nit->integ(mfe,0.0,0.0);
    b.ed*=prefac*temper;
    if (!b.inc_rest_mass) b.ed-=b.n*b.m;
    
    // Compute the entropy

    b.en=nit->integ(mfs,0.0,0.0);
    b.en*=prefac;
    
  }

  b.pr=-b.ed+temper*b.en+b.mu*b.n;

  return;
}

void boson_rel::nu_from_n(boson &b, double temper) {
  double nex;

  nex=b.nu/temper;
  funct mf=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::solve_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
  density_root->solve(nex,mf);
  b.nu=nex*temper;
  
  return;
}

void boson_rel::calc_density(boson &b, double temper) {

  if (temper<=0.0) {
    O2SCL_ERR2("Temperature less than or equal to zero in ",
	       "boson_rel::calc_density().",exc_einval);
  }
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

  nu_from_n(b,temper);

  funct fe=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_energy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
  funct fs=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_entropy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);

  b.ed=dit->integ(fe,0.0,sqrt(pow(20.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.ed*=b.g/2.0/pi2;
  b.en=dit->integ(fs,0.0,sqrt(pow(20.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.en*=b.g/2.0/pi2;

  b.pr=-b.ed+temper*b.en+b.mu*b.n;

  return;
}

double boson_rel::deg_density_fun(double k, boson &b, double T) {

  double E=sqrt(k*k+b.ms*b.ms), ret;

  ret=k*k/(exp(E/T-b.nu/T)-1.0);

  if (!std::isfinite(ret)) {
    cout << "1: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }
  
  return ret;
}
  
double boson_rel::deg_energy_fun(double k, boson &b, double T) {

  double E=sqrt(k*k+b.ms*b.ms), ret;

  ret=k*k*E/(exp(E/T-b.nu/T)-1.0);
  
  if (!std::isfinite(ret)) {
    cout << "2: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}
  
double boson_rel::deg_entropy_fun(double k, boson &b, double T) {

  double E=sqrt(k*k+b.ms*b.ms), nx, ret;
  nx=1.0/(exp(E/T-b.nu/T)-1.0);
  ret=-k*k*(nx*log(nx)-(1.0+nx)*log(1.0+nx));
  
  if (!std::isfinite(ret)) {
    cout << "3: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}
  
double boson_rel::density_fun(double u, boson &b, double T) {
  double ret, y, eta;

  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  eta=b.ms/T;

  if (y-u>200.0 && eta-u>200.0) {
    if (eta+u+y>100.0) {
      ret=0.0;
    } else {
      ret=(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1.0);
    }
  } else {
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
  }

  if (!std::isfinite(ret)) {
    cout << "4: " << u << " " << y << " " << eta << " " 
	 << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}

double boson_rel::energy_fun(double u, boson &b, double T) {
  double ret, y, eta;

  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  eta=b.ms/T;
  
  if (y-u>200.0 && eta-u>200.0) {
    if (eta+u+y>100.0) {
      ret=0.0;
    } else {
      ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1.0);
    }
  } else {
    ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
  }
  
  if (!std::isfinite(ret)) {
    cout << "5: " << u << " " << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}

double boson_rel::entropy_fun(double u, boson &b, double T) {
  double ret, y, eta, term1, term2;

  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  eta=b.ms/T;

  if (u-eta>200.0 && u-y>200.0) {
    ret=0.0;
  } else {
    term1=exp(eta+u)*log(1.0/1.0-exp(y-eta-u));
    term2=exp(y)*log(1.0/(exp(eta+u-y)-1.0));
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
      (exp(eta+u)-exp(y));
  }

  if (!std::isfinite(ret)) {
    return 0.0;
  }

  /*
  if (false) {
    cout << "6: " << u << " " << eta << " " << y << endl;
    cout << b.ms << " " << b.nu << " " << T << endl;

    u=200.0;
    term1=exp(eta+u)*log(1.0/1.0-exp(y-eta-u));
    term2=exp(y)*log(1.0/(exp(eta+u-y)-1.0));
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
      (exp(eta+u)-exp(y));
    cout << ret << endl;
    
    exit(-1);
  }
  */

  return ret;
}

double boson_rel::solve_fun(double x, boson &b, double T) {
  double nden, yy;
  
  funct fd=std::bind(std::mem_fn<double(double,boson &b,double)>
		       (&boson_rel::deg_density_fun),
		       this,std::placeholders::_1,std::ref(b),T);
  
  b.nu=T*x;
  nden=dit->integ(fd,0.0,sqrt(pow(20.0*T+b.nu,2.0)-b.ms*b.ms));
  nden*=b.g/2.0/pi2;
  yy=nden/b.n-1.0;

  return yy;
}

void boson_rel::pair_mu(boson &b, double temper) {
  
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }
  calc_mu(b,temper);
  
  boson antip(b.ms,b.g);
  b.anti(antip);
  calc_mu(antip,temper);
  b.n-=antip.n;
  b.pr+=antip.pr;
  b.ed+=antip.ed;
  b.en+=antip.en;

  return;
}

