/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/fermion_eff.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

fermion_eff::fermion_eff() {

  err_nonconv=true;

  tlimit=0.0;

  parma=0.42;
  Pmnf.resize(4,4);
  Pmnf(0,0)=5.51219;
  Pmnf(0,1)=18.6215;
  Pmnf(0,2)=22.0078;
  Pmnf(0,3)=8.71963;
  Pmnf(1,0)=17.4343;
  Pmnf(1,1)=57.7188;
  Pmnf(1,2)=66.1098;
  Pmnf(1,3)=25.6323;
  Pmnf(2,0)=18.1239;
  Pmnf(2,1)=58.7781;
  Pmnf(2,2)=65.1429;
  Pmnf(2,3)=24.4246;
  Pmnf(3,0)=6.30952;
  Pmnf(3,1)=19.8967;
  Pmnf(3,2)=21.1375;
  Pmnf(3,3)=7.55866;
  sizem=3;
  sizen=3;

  psi_root=&def_psi_root;
  density_root=&def_density_root;
  min_psi=-4.0;

  verify_ti=false;
}

fermion_eff::~fermion_eff() {
}

void fermion_eff::load_coefficients(int ctype) {
  
  if (ctype==cf_fermilat3) {
    sizem=3;
    sizen=3;
    Pmnf.resize(4,4);
    Pmnf(0,0)=5.447327;
    Pmnf(0,1)=18.40237;
    Pmnf(0,2)=21.74882;
    Pmnf(0,3)=8.692826;
    Pmnf(1,0)=17.20333;
    Pmnf(1,1)=56.90495;
    Pmnf(1,2)=65.12704;
    Pmnf(1,3)=25.23724;
    Pmnf(2,0)=17.86753;
    Pmnf(2,1)=57.84077;
    Pmnf(2,2)=63.99390;
    Pmnf(2,3)=23.96116;
    Pmnf(3,0)=6.216856;
    Pmnf(3,1)=19.54390;
    Pmnf(3,2)=20.69990;
    Pmnf(3,3)=7.381855;
    parma=0.425;
  } else if (ctype==cf_fermijel2) {
    sizem=2;
    sizen=2;
    Pmnf.resize(3,3);
    Pmnf(0,0)=5.23810;
    Pmnf(0,1)=12.4991;
    Pmnf(0,2)=8.36058;
    Pmnf(1,0)=11.1751;
    Pmnf(1,1)=25.3687;
    Pmnf(1,2)=15.6543;
    Pmnf(2,0)=5.91800;
    Pmnf(2,1)=12.4945;
    Pmnf(2,2)=6.82530;
    parma=0.442;
  } else if (ctype==cf_fermijel3) {
    sizem=3;
    sizen=3;
    Pmnf.resize(4,4);
    Pmnf(0,0)=5.51219;
    Pmnf(0,1)=18.6215;
    Pmnf(0,2)=22.0078;
    Pmnf(0,3)=8.71963;
    Pmnf(1,0)=17.4343;
    Pmnf(1,1)=57.7188;
    Pmnf(1,2)=66.1098;
    Pmnf(1,3)=25.6323;
    Pmnf(2,0)=18.1239;
    Pmnf(2,1)=58.7781;
    Pmnf(2,2)=65.1429;
    Pmnf(2,3)=24.4246;
    Pmnf(3,0)=6.30952;
    Pmnf(3,1)=19.8967;
    Pmnf(3,2)=21.1375;
    Pmnf(3,3)=7.55866;
    parma=0.42;
  } else if (ctype==cf_fermijel3cons) {
    sizem=3;
    sizen=3;
    Pmnf.resize(4,4);
    Pmnf(0,0)=5.34599;
    Pmnf(0,1)=18.0517;
    Pmnf(0,2)=21.3422;
    Pmnf(0,3)=8.53240;
    Pmnf(1,0)=16.8441;
    Pmnf(1,1)=55.7051;
    Pmnf(1,2)=63.6901;
    Pmnf(1,3)=24.6213;
    Pmnf(2,0)=17.4708;
    Pmnf(2,1)=56.3902;
    Pmnf(2,2)=62.1319;
    Pmnf(2,3)=23.2602;
    Pmnf(3,0)=6.07364;
    Pmnf(3,1)=18.9992;
    Pmnf(3,2)=20.0285;
    Pmnf(3,3)=7.11153;
    parma=0.433;
  } else {
    O2SCL_ERR("Invalid type in fermion_eff::load_coefficients().",
		  exc_efailed);
  }

  return;
}

int fermion_eff::calc_mu(fermion &f, double temper) {
  int nn, mm;
  double xx, pren, preu, prep, sumn, sumu, sump;
  double ff, gg, opf, opg, nc;

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
  if (f.ms<0.0) {
    O2SCL_ERR2("Effective mass less than zero ",
	       "in fermion_eff::calc_mu().",exc_einval);
  }

  // Use zero temperature results if necessary
  if (temper<=tlimit) {
    calc_mu_zerot(f);
    return 0;
  }
  
  // Compute psi
  double psi;
  if (f.inc_rest_mass) {
    psi=(f.nu-f.ms)/temper;
  } else {
    psi=(f.nu+f.m-f.ms)/temper;
  }

  // Try the nondegenerate expansion if we can
  if (psi<-4.0) {
    bool acc=calc_mu_ndeg(f,temper,1.0e-7);
    if (acc) return 0;
  }

  // Try the degenerate expansion if we can
  if (psi>10.0) {
    bool acc=calc_mu_deg(f,temper,1.0e-10);
    if (acc) return 0;
  }

  // The following code assumes that the chemical potential includes
  // the rest mass, so we add it back in if necessary, and we'll
  // subtract it back out before finishing.
  if (!f.inc_rest_mass) f.nu+=f.m;

  // Create an initial guess according to the value of psi
  if (psi<0.0) {
    // small f expansion
    xx=4.0*parma*exp(-2.0+psi);
  } else {
    // large f expansion
    xx=parma/8.0*(4.0+psi*psi+psi*sqrt(8.0+psi*psi));
  }

  if (!std::isfinite(psi) || !std::isfinite(xx)) {
    O2SCL_ERR("Psi or xx not finite in fermion_eff::calc_mu().",
	      exc_efailed);
  }

  // Perform the solution
  funct mfs=std::bind(std::mem_fn<double(double,double)>
			(&fermion_eff::solve_fun),
			this,std::placeholders::_1,psi);
  bool psi_root_ec=psi_root->err_nonconv;
  psi_root->err_nonconv=false;
  int ret=psi_root->solve(xx,mfs);
  if (ret!=0) {
    o2scl::root_brent_gsl<> rbg;
    rbg.err_nonconv=false;
    int ret2=rbg.solve(xx,mfs);
    if (ret2!=0) {
      O2SCL_ERR("Psi solver failed in fermion_eff::calc_mu().",
		o2scl::exc_efailed);
    }
  }
  psi_root->err_nonconv=psi_root_ec;
  ff=xx;

  // Compute Eqs. 12 and 16 in Johns, et al. (1996)
  opf=1.0+ff;
  gg=temper/f.ms*sqrt(opf);
  if (gg<=0.0) gg=0.0;
  opg=1.0+gg;
  nc=1.0/pi2*pow(f.ms,3.0);
    
  pren=ff*pow(gg,1.5)/sqrt(1+ff/parma)/pow(opf,sizem+0.5)/pow(opg,sizen-1.5);
  preu=ff*pow(gg,2.5)*pow(opg,1.5)/pow(opf,sizem+1.0)/pow(opg,(double)sizen);
  prep=ff*pow(gg,2.5)*pow(opg,1.5)/pow(opf,sizem+1.0)/pow(opg,(double)sizen);
  
  sumn=0.0;
  sumu=0.0;
  sump=0.0;
  for(mm=0;mm<=sizem;mm++) {
    for(nn=0;nn<=sizen;nn++) {
      sumn+=Pmnf(mm,nn)*pow(ff,(double)mm)*pow(gg,(double)nn)*
	(1.0+mm+ff/opf*(0.25+nn/2.0-sizem)+ff*gg/opf/opg*(0.75-sizen/2.0));
      sumu+=Pmnf(mm,nn)*pow(ff,(double)mm)*pow(gg,(double)nn)*
	(1.5+nn+gg/opg*(1.5-sizen));
      sump+=Pmnf(mm,nn)*pow(ff,(double)mm)*pow(gg,(double)nn);
    }
  } 

  f.n=f.g/2.0*pren*sumn*nc;
  f.pr=f.g/2.0*prep*sump*nc*f.ms;

  // Note that here, we are adding the rest mass contribution with the
  // effective mass, since the EFF coefficients, even though they are
  // fully relativistic, were designed to compute the energy without
  // the rest mass term. The 'bare' rest mass energy is removed later
  // if necessary.
  f.ed=f.g/2.0*preu*sumu*nc*f.ms+f.n*f.ms;

  f.en=(f.ed+f.pr-f.n*f.nu)/temper;

  // If necessary, remove the rest mass
  if (!f.inc_rest_mass) {
    f.nu-=f.m;
    f.ed-=f.n*f.m;
  }

  return 0;
}

int fermion_eff::calc_density(fermion &f, double temper) {

  double T=temper;

  // Handle special cases

  if (temper<=tlimit) {
    calc_density_zerot(f);
    return 0;
  }

  if (f.non_interacting) { f.ms=f.m; f.nu=f.mu;}
  
  if (f.n==0.0) {
    f.nu=f.ms;
    if (!f.inc_rest_mass) f.nu=0.0;
    f.ed=0.0;
    f.pr=0.0;
    f.en=0.0;
    if (f.non_interacting) { f.mu=f.nu; }
    return 0;
  }
  
  // If the initial guess leads to a zero density, change the initial
  // guess. There's probably a better way of getting an initial guess
  // for mu which could be implemented in the future.
  double nstore=f.n;
  if (f.non_interacting) f.mu=f.nu;
  calc_mu(f,temper);
  if (f.n==0.0) {
    if (f.inc_rest_mass) f.nu=f.ms;
    else f.nu=f.ms-f.m;
  }
  f.n=nstore;
  
  // Perform the solution
  funct df2=std::bind(std::mem_fn<double(double,fermion &,double)>
			(&fermion_eff::density_fun),
			this,std::placeholders::_1,std::ref(f),temper);
  int ret=density_root->solve(f.nu,df2);
  if (ret!=0) {
    O2SCL_CONV2_RET("Function calc_density() failed in fermion_eff::",
		    "calc_density().",exc_efailed,this->err_nonconv);
  }
  df2(f.nu);

  if (f.non_interacting) f.mu=f.nu;

  return 0;
}

double fermion_eff::solve_fun(double x, double psi) {
  double ff, sqt, aarg, y;

  ff=x;

  if (ff<=0.0) {
    string str="Variable 'f' is less than or equal to zero ("+dtos(ff)+
      ") in fermion_eff::solve_fun().";
    O2SCL_ERR(str.c_str(),exc_efailed);
  }

  // For small f, we lose too much precision using
  // the usual formula, so we use the expansion
  if (ff/parma<1.0e-10) {
    sqt=1.0+ff/parma/2.0;
    aarg=ff/parma/2.0/(sqt+1.0);
  } else {
    sqt=sqrt(1.0+ff/parma);
    aarg=(sqt-1.0)/(sqt+1.0);
  }
  
  if (fabs(psi)<1.0e-10) {
    y=(2.0*sqt+log(aarg))-psi;
  } else {
    y=(2.0*sqt+log(aarg))/psi-1.0;
  }
  
  if (!std::isfinite(y)) {
    O2SCL_ERR2("Variable 'y' not finite in ",
		   "fermion_eff::solve_fun().",exc_efailed);
  }

  return y;
}

double fermion_eff::density_fun(double x, fermion &f, double temper) {
  double nold, nnew, y;

  if (f.non_interacting) f.mu=x;
  else f.nu=x;
  nold=f.n;
  calc_mu(f,temper);
  nnew=f.n;
  y=nnew/nold-1.0;
  f.n=nold;
  
  return y;
}
				 
void fermion_eff::pair_mu(fermion &f, double temper) {
  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

  if (f.ms==0.0) {

    massless_pair_mu(f,temper);

  } else {

    // Create the antiparticle with same mass and degeneracy
    fermion antip(f.ms,f.g);
    f.anti(antip);
    
    // Compute particle contribution
    calc_mu(f,temper);
    
    // Compute antiparticle contribution
    calc_mu(antip,temper);

    // Add the antiparticle contribution
    f.n-=antip.n;
    f.pr+=antip.pr;
    f.ed+=antip.ed;
    f.en+=antip.en;
  }
  
  return;
}

int fermion_eff::pair_density(fermion &f, double temper) {

  if (f.non_interacting) { f.ms=f.m; f.nu=f.mu; }
  if (f.ms==0.0) {
    massless_pair_density(f,temper);
    return success;
  }
  
  funct pdf2=std::bind(std::mem_fn<double(double,fermion &,double)>
			 (&fermion_eff::pair_density_fun),
			 this,std::placeholders::_1,std::ref(f),temper);
  bool density_root_ec=density_root->err_nonconv;
  density_root->err_nonconv=false;
  int ret=density_root->solve(f.nu,pdf2);
  if (ret!=0) {
    // Bracket the root and use root_brent_gsl
    double blow=f.nu*0.99;
    double bhigh=f.nu*1.01;
    double ylow=pdf2(blow);
    double yhigh=pdf2(bhigh);
    for(size_t j=0;j<5 && yhigh*ylow>0.0;j++) {
      bhigh+=(bhigh-blow)*10.0;
      yhigh=pdf2(bhigh);
    }
    for(size_t j=0;j<5 && yhigh*ylow>0.0;j++) {
      blow-=(bhigh-blow)*10.0;
      ylow=pdf2(blow);
    }
    if (ylow*yhigh<0.0) {
      root_brent_gsl<> rbg;
      rbg.err_nonconv=false;
      int ret2=rbg.solve_bkt(blow,bhigh,pdf2);
      if (ret2!=0) {
	O2SCL_CONV2_RET("Solvers failed in fermion_eff::",
			"pair_density().",exc_efailed,this->err_nonconv);
      }
      f.nu=blow;
    }
  }
  density_root->err_nonconv=density_root_ec;
  pdf2(f.nu);
  
  if (f.non_interacting) { f.mu=f.nu; }

  return success;
}

double fermion_eff::pair_density_fun(double x, fermion &f, double temper) {
  double nold, nnew, y;

  if (f.non_interacting) f.mu=x;
  else f.nu=x;
  nold=f.n;
  pair_mu(f,temper);
  nnew=f.n;
  y=nnew/nold-1.0;
  f.n=nold;

  return y;
}

