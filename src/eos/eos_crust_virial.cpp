/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2014, Andrew W. Steiner
  
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
#include <o2scl/eos_crust_virial.h>
#include <o2scl/fit_nonlin.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::matrix<double> ubmatrix;

eos_crust_virial::eos_crust_virial() {

  alpha.init(o2scl_settings.get_convert_units().convert
	     ("kg","1/fm",o2scl_mks::mass_alpha),1.0);
  deuteron.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_deuteron),1.0);

  // Data from the tables in Horowitz05. The rows for T<=10 MeV
  // are from Table 1 and Table 2 and the rows for T>10 MeV are
  // from Table 3.
  double arr[16][9]={
    {0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
    {1,0.288,0.032,19.4,-43.8,1.51,1.21,2.55,1.59},
    {2,0.303,0.012,6.10,-7.39,2.26,0.90,4.12,2.95},
    {3,0.306,0.005,4.01,-3.54,2.57,0.63,5.64,4.81},
    {4,0.307,0.002,3.19,-2.30,2.73,0.44,7.44,7.90},
    {5,0.308,0.002,2.74,-1.73,2.81,0.32,9.57,11.3},
    {6,0.308,0.003,2.46,-1.40,2.86,0.23,11.9,14.3},
    {7,0.308,0.004,2.26,-1.18,2.89,0.18,14.3,16.3},
    {8,0.309,0.006,2.11,-1.04,2.92,0.15,16.5,17.3},
    {9,0.310,0.008,2.00,-0.93,2.93,0.14,18.6,17.5},
    {10,0.311,0.010,1.91,-0.85,2.95,0.13,20.4,17.0},
    {12,0.313,0.013,1.76,-0.73,2.97,0.12,23.3,14.7},
    {14,0.315,0.014,1.66,-0.65,2.98,0.10,25.4,11.7},
    {16,0.317,0.014,1.57,-0.59,3.00,0.07,26.7,8.48},
    {18,0.319,0.013,1.51,-0.55,3.00,0.02,27.5,5.44},
    {20,0.320,0.011,1.45,-0.52,3.00,-0.04,28.0,2.69}};
  for(size_t i=0;i<16;i++) {
    Tv.push_back(arr[i][0]);
    bnv.push_back(arr[i][1]);
    Tbnpv.push_back(arr[i][2]);
    bpnv.push_back(arr[i][3]);
    Tbpnpv.push_back(arr[i][4]);
    banv.push_back(arr[i][5]);
    Tbanpv.push_back(arr[i][6]);
    bav.push_back(arr[i][7]);
    Tbapv.push_back(arr[i][8]);
  }
  ibn.set(16,Tv,bnv,itp_linear);
  iTbnp.set(16,Tv,Tbnpv,itp_linear);
  ibpn.set(16,Tv,bpnv,itp_linear);
  iTbpnp.set(16,Tv,Tbpnpv,itp_linear);
  iban.set(16,Tv,banv,itp_linear);
  iTbanp.set(16,Tv,Tbanpv,itp_linear);
  iba.set(16,Tv,bav,itp_linear);
  iTbap.set(16,Tv,Tbapv,itp_linear);
}

double eos_crust_virial::bn(double T) {
  return ibn.eval(T);
}

double eos_crust_virial::ban(double T) {
  return iban.eval(T);
}

double eos_crust_virial::ba(double T) {
  return iba.eval(T);
}

double eos_crust_virial::bpn(double T) {
  return ibpn.eval(T);
}

double eos_crust_virial::Tbn_prime(double T) {
  return iTbnp.eval(T);
}

double eos_crust_virial::Tban_prime(double T) {
  return iTbanp.eval(T);
}

double eos_crust_virial::Tba_prime(double T) {
  return iTbap.eval(T);
}
    
double eos_crust_virial::Tbpn_prime(double T) {
  return iTbpnp.eval(T);
}

int eos_crust_virial::calc_temp_p(fermion &n, fermion &p, double T, 
				  thermo &th) {
      
  return calc_temp_p_alpha(n,p,deuteron,alpha,T,th);
}

int eos_crust_virial::calc_temp_p_alpha
(fermion &n, fermion &p, boson &d, boson &a, double T, thermo &th) {
      
#if !O2SCL_NO_RANGE_CHECK
  if (!o2scl::is_finite(n.mu) || !o2scl::is_finite(n.mu) ||
      !o2scl::is_finite(T)) {
    O2SCL_ERR2("Chemical potentials or temperature not finite in ",
	       "eos_crust_virial::calc_eq_temp_p_alpha().",exc_einval);
  }
  if (fabs(n.g-2.0)>1.0e-10 || fabs(p.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_crust_virial::calc_eq_temp_p_alpha().",exc_einval);
  }
  if (fabs(n.m-4.5)>1.0 || fabs(p.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_crust_virial::calc_eq_temp_p_alpha().",exc_einval);
  }
#endif

  // Temperature in MeV
  double TMeV=T*o2scl_const::hc_mev_fm;

  // Alpha and deuteron chemical potentials in nuc. stat. eq.
  a.mu=2.0*p.mu+2.0*n.mu;
  d.mu=p.mu+n.mu;

  // Fugacities for neutrons, protons, alphas, and deuterons
  double zn, zp, za, zd;
  if (n.inc_rest_mass) {
    zn=exp((n.mu-n.m)/T);
  } else {
    zn=exp(n.mu/T);
  }
  if (p.inc_rest_mass) {
    zp=exp((p.mu-p.m)/T);
  } else {
    zp=exp(p.mu/T);
  }
  za=zp*zp*zn*zn*exp(-(a.m-2.0*n.m-2.0*p.m)/T);
  zd=zp*zn*exp(-(d.m-n.m-p.m)/T);
      
  // de Broglie wavelengths for nucleons, alphas, and deutrons
  double lambda=sqrt(4.0*o2scl_const::pi/(n.m+p.m)/T);
  double lambda3=lambda*lambda*lambda;
  double lambda_a=sqrt(2.0*o2scl_const::pi/a.m/T);
  double lambda_a3=lambda_a*lambda_a*lambda_a;
  double lambda_d=sqrt(2.0*o2scl_const::pi/d.m/T);
  double lambda_d3=lambda_d*lambda_d*lambda_d;
      
  // Pressure
  th.pr=T*(2.0/lambda3*(zn+zp+(zn*zn+zp*zp)*bn(TMeV)+2.0*zp*zn*bpn(TMeV))
	   +1.0/lambda_a3*(za+za*za*ba(TMeV)+2.0*za*(zn+zp)*ban(TMeV)));
      
  // Densities
  n.n=2.0/lambda3*(zn+2.0*zn*zn*bn(TMeV)+2.0*zp*zn*bpn(TMeV)+
		   8.0*za*zn*ban(TMeV));
  p.n=2.0/lambda3*(zp+2.0*zp*zp*bn(TMeV)+2.0*zp*zn*bpn(TMeV)+
		   8.0*za*zp*ban(TMeV));
  a.n=1.0/lambda_a3*(za+2.0*za*za*ba(TMeV)+2.0*za*(zn+zp)*ban(TMeV));
  d.n=1.0/lambda_d3*zd;
      
  // Entropy
  th.en=5.0*th.pr/2.0/T-n.n*log(zn)-p.n*log(zp)-a.n*log(za)+
    2.0/lambda3*((zn*zn+zp*zp)*Tbn_prime(TMeV)+
		 2.0*zp*zn*Tbpn_prime(TMeV))+
    1.0/lambda_a3*(za*za*Tba_prime(TMeV)+2.0*za*(zn+zp)*Tban_prime(TMeV));

  // Energy density
  th.ed=-th.pr+T*th.en;
  if (n.inc_rest_mass) {
    th.ed+=n.n*n.mu;
  } else {
    th.ed+=n.n*n.mu+n.n*n.m;
  }
  if (p.inc_rest_mass) {
    th.ed+=p.n*p.mu;
  } else {
    th.ed+=p.n*p.mu+p.n*p.m;
  }
  if (a.inc_rest_mass) {
    th.ed+=a.n*a.mu;
  } else {
    th.ed+=a.n*a.mu+a.n*a.m;
  }
  /*
    Alternate expression for the energy density

    double ed2=th.pr*1.5-a.n*Ealpha+2.0*T/lambda3*
    ((zn*zn+zp*zp)*Tbn_prime(TMeV)+2.0*zn*zp*Tbpn_prime(TMeV))+
    T/lambda_a3*(za*za*Tba_prime(TMeV)+2.0*za*(zn+zp)*Tban_prime(TMeV));
  */

  return 0;
}

void eos_crust_virial::fit() {

  vector<double> berr(16);

  typedef fit_funct11_strings<vector<double> > fit_func;

  // Fitter class
  fit_nonlin<chi_fit_funct<vector<double>,ubmatrix,fit_func>,
	     vector<double>,ubmatrix> fn;

  vector<double> params(3);
  ubmatrix covar(3,3);
  double chi2;

#ifdef O2SCL_NEVER_DEFINED
}{
#endif

  // --------------------------------------------
  // Fit neutron virial coefficient

  for(size_t i=0;i<16;i++) {
    berr[i]=fabs(bnv[i])/1.0e2;
  }

  // Replace the T=0 point
  double Tv0=Tv[0];
  double bnv0=bnv[0];
  Tv[0]=0.999;
  bnv[0]=0.288;
  berr[0]=bnv[0]/1.0e2;

  // Fitting function
  fit_func ffs("a-b/(T+0.5)+c*T","a,b,c","T");

  // Chi-squared and fitting data
  chi_fit_funct<vector<double>,ubmatrix,fit_func> 
    cff(16,Tv,bnv,berr,ffs);

  params[0]=0.3;
  params[1]=0.01;
  params[2]=0.0006;

  fn.fit(3,params,covar,chi2,cff);
  cout << chi2 << " " 
       << params[0] << " " << params[1] << " " << params[2] << endl;

  // --------------------------------------------
  // Fit neutron-proton virial coefficient
  
  for(size_t i=0;i<16;i++) {
    berr[i]=fabs(bpnv[i])/1.0e2;
  }

  // Replace the T=0 point
  double bpnv0=bpnv[0];
  Tv[0]=0.999;
  bpnv[0]=19.4;
  berr[0]=bpnv[0]/1.0e2;
  
  // Fitting function
  fit_func ffs2("-a+b*exp(2.099/T)-c*T","a,b,c","T");

  // Chi-squared and fitting data
  chi_fit_funct<vector<double>,ubmatrix,fit_func> 
    cff2(16,Tv,bpnv,berr,ffs2);

  params[0]=1.0;
  params[1]=2.5;
  params[2]=0.018;

  fn.fit(3,params,covar,chi2,cff2);
  cout << chi2 << " " 
       << params[0] << " " << params[1] << " " << params[2] << endl;

  // --------------------------------------------
  // Return original values

  Tv[0]=Tv0;
  bnv[0]=bnv0;
  bpnv[0]=bpnv0;

  return;
}
