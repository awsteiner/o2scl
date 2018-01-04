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

  funct fd=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_density_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
  funct fe=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_energy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
  funct fs=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_entropy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);

  b.n=dit->integ(fd,0.0,sqrt(pow(15.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.n*=b.g/2.0/pi2;
  b.ed=dit->integ(fe,0.0,sqrt(pow(15.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.ed*=b.g/2.0/pi2;
  b.en=dit->integ(fs,0.0,sqrt(pow(15.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.en*=b.g/2.0/pi2;

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

  return ret;
}
  
double boson_rel::deg_energy_fun(double k, boson &b, double T) {

  double E=sqrt(k*k+b.ms*b.ms), ret;

  ret=k*k*E/(exp(E/T-b.nu/T)-1.0);
  
  return ret;
}
  
double boson_rel::deg_entropy_fun(double k, boson &b, double T) {

  double E=sqrt(k*k+b.ms*b.ms), nx, ret;
  nx=1.0/(exp(E/T-b.nu/T)-1.0);
  ret=-k*k*(nx*log(nx)-(1.0+nx)*log(1.0+nx));
  
  return ret;
}
  
double boson_rel::density_fun(double u, boson &b, double T) {
  double ret, y, mx;

  y=b.nu/T;
  mx=b.ms/T;
  
  ret=(mx+u)*sqrt(u*u+2.0*mx*u)*exp(u+y)/(exp(y)-exp(mx+u));

  return ret;
}

double boson_rel::energy_fun(double u, boson &b, double T) {
  double ret, y, mx;

  y=b.nu/T;
  mx=b.ms/T;
  
  ret=(mx+u)*(mx+u)*sqrt(u*u+2.0*mx*u)*exp(u+y)/(exp(y)-exp(mx+u));
  
  return ret;
}

double boson_rel::entropy_fun(double u, boson &b, double T) {
  double ret, y, mx, term1, term2;

  y=b.mu/T;
  mx=b.ms/T;

  term1=log(exp(y-mx-u)-1.0)/(exp(y-mx-u)-1.0);
  term2=log(1.0-exp(mx+u-y))/(1.0-exp(mx+u-y));
  ret=(mx+u)*exp(u)*sqrt(u*u+2.0*mx*u)*(term1+term2);
  
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

