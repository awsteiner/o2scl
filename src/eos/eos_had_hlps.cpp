 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2014, Andrew W. Steiner
  
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
#include <o2scl/eos_had_hlps.h>
#include <o2scl/hdf_file.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

// --------------------------------------------------------------

eos_had_hlps::eos_had_hlps() {

  // Some default values from the paper
  kprime=0.0;
  n0=0.16;
  gamma=4.0/3.0;
  alpha=5.87;
  eta=3.81;
  alphaL=1.4;
  etaL=0.86;
}

#ifdef O2SCL_NEVER_DEFINED
void eos_had_hlps::func(double K) {

  // These limits are from Thomas Krueger and are detailed in 
  // the PRL and in PRC 88 025802 (2013)
  
  // Upper limit for neutron matter from chiral effective theory
  double upper_limit[22][2]={{3.37737e-05,0.109},
			     {0.00027019,0.395997},
			     {0.000911891,0.809968},
			     {0.00216152,1.30979},
			     {0.00422172,1.88308},
			     {0.00729513,2.54298},
			     {0.0115844,3.26487},
			     {0.0172921,4.00818},
			     {0.024621,4.9124},
			     {0.0337737,5.90204},
			     {0.0449528,6.9384},
			     {0.058361,8.23317},
			     {0.0742009,9.69399},
			     {0.0926751,11.6035},
			     {0.113986,13.9669},
			     {0.125769,15.4819},
			     {0.138337,17.1633},
			     {0.151716,18.984},
			     {0.16593,20.9512},
			     {0.181006,23.1962},
			     {0.196968,25.5491},
			     {0.213843,28.1975}};

  // Lower limit for neutron matter from chiral effective theory
  double lower_limit[25][2]={{3.37737e-05,0.106},
			     {0.00027019,0.377998},
			     {0.000911891,0.765969},
			     {0.00216152,1.23579},
			     {0.00422172,1.78708},
			     {0.00729513,2.44496},
			     {0.00814476,2.59203},
			     {0.0115844,3.13234},
			     {0.0172921,3.87178},
			     {0.0214105,4.28327},
			     {0.024621,4.55192},
			     {0.0337737,5.08014},
			     {0.0449528,5.66559},
			     {0.058361,6.22853},
			     {0.0742009,6.95932},
			     {0.0926751,8.07315},
			     {0.113986,9.54407},
			     {0.125769,10.4475},
			     {0.138337,11.5665},
			     {0.145046,12.2088},
			     {0.151716,12.7498},
			     {0.16593,14.0645},
			     {0.181006,15.504},
			     {0.196968,17.2011},
			     {0.213843,18.9384}};

  std::vector<double> lower_nB, lower_E, upper_nB, upper_E;
  cout.setf(ios::scientific);
  for(size_t i=0;i<25;i++) {
    lower_nB.push_back(lower_limit[i][0]);
    lower_E.push_back(lower_limit[i][1]);
    if (i<22) {
      upper_nB.push_back(upper_limit[i][0]);
      upper_E.push_back(upper_limit[i][1]);
    }
  }
  interp<std::vector<double> > it(itp_cspline);

  cout << "Lower and upper bounds for the energy of neutron matter: " 
       << endl;
  cout << it.eval(0.16,25,lower_nB,lower_E) << endl;
  cout << it.eval(0.16,22,upper_nB,upper_E) << endl;
  cout << "Lower and upper bounds for the derivative of the "
       << "energy of neutron matter: " << endl;
  cout << it.deriv(0.16,25,lower_nB,lower_E) << endl;
  cout << it.deriv(0.16,22,upper_nB,upper_E) << endl;

  //1.349511e+01
  //2.011534e+01
  //9.534376e+01
  //1.382914e+02
  
  return;
}

void eos_had_hlps::fix_comp(double K) {
  
  if (K<210.0 || K>270.0) {
    O2SCL_ERR("K out of range in eos_had_hlps::fix_comp().",
	      exc_einval);
  }

  ubvector vcomp, valpha, vgamma, veta;
  o2scl::interp<> itp;

  // Compressibility, alpha, gamma, eta
  double arr[21][4]={
    {1.114901e+00,7.018317e+00,1.256284e+00,4.949759e+00},
    {1.125036e+00,6.842056e+00,1.265751e+00,4.773499e+00},
    {1.135172e+00,6.677136e+00,1.275261e+00,4.608578e+00},
    {1.145307e+00,6.523232e+00,1.284771e+00,4.454674e+00},
    {1.155443e+00,6.379277e+00,1.294281e+00,4.310720e+00},
    {1.165578e+00,6.244337e+00,1.303790e+00,4.175780e+00},
    {1.175714e+00,6.117591e+00,1.313300e+00,4.049033e+00},
    {1.185849e+00,5.998313e+00,1.322809e+00,3.929755e+00},
    {1.195985e+00,5.885851e+00,1.332320e+00,3.817294e+00},
    {1.206120e+00,5.779659e+00,1.341829e+00,3.711102e+00},
    {1.216255e+00,5.679206e+00,1.351339e+00,3.610649e+00},
    {1.226391e+00,5.584049e+00,1.360849e+00,3.515491e+00},
    {1.236526e+00,5.493788e+00,1.370358e+00,3.425231e+00},
    {1.246662e+00,5.408039e+00,1.379868e+00,3.339481e+00},
    {1.256797e+00,5.326478e+00,1.389378e+00,3.257920e+00},
    {1.266933e+00,5.248807e+00,1.398887e+00,3.180249e+00},
    {1.277068e+00,5.174753e+00,1.408397e+00,3.106195e+00},
    {1.287204e+00,5.104070e+00,1.417907e+00,3.035512e+00},
    {1.297339e+00,5.036533e+00,1.427416e+00,2.967975e+00},
    {1.307475e+00,4.971936e+00,1.436926e+00,2.903378e+00},
    {1.317610e+00,4.910091e+00,1.446435e+00,2.841533e+00}};
  vcomp.resize(21);
  valpha.resize(21);
  vgamma.resize(21);
  veta.resize(21);
  for(size_t i=0;i<21;i++) {
    vcomp[i]=arr[i][0];
    valpha[i]=arr[i][1];
    vgamma[i]=arr[i][2];
    veta[i]=arr[i][3];
  }
  alpha=itp.interp(21,vcomp,valpha,K/hc_mev_fm);
  gamma=itp.interp(21,vcomp,vgamma,K/hc_mev_fm);
  eta=itp.interp(21,vcomp,veta,K/hc_mev_fm);

  return;
}

#endif

void eos_had_hlps::fix_coeffs(double M, double local_n0, double B, double K) {

  double lastg;
  for(size_t i=0;i<100;i++) {

    // First fix eta and alpha from B, M, and n0
    double p=20.0*cbrt(3.0)*B*M/cbrt(local_n0)/cbrt(local_n0);
    double q=3.0*cbrt(2.0*pi2*pi2);
    double C=cbrt(4.0/pi2/pi2)/15.0/(gamma-1.0);
    eta=C*(q-p);
    alpha=0.8+C*(q-p)*gamma;
    lastg=gamma;

    // Now compute new gamma from K, M and n0
    double C2=4.0*K*M/3.0*cbrt(4.0/9.0/local_n0/local_n0/pi2/pi2);
    double a0=4.0-C2-6.0*alpha;
    double a1=3.0*eta;
    double a2=3.0*eta;
    std::complex<double> x1, x2;
    quad.solve_rc(a2,a1,a0,x1,x2);
    if (fabs(x1.real()-4.0/3.0)<fabs(x2.real()-4.0/3.0)) {
      gamma=x1.real();
    } else {
      gamma=x2.real();
    }

    // Stop early if we've converged
    if (fabs(lastg-gamma)<1.0e-6) i=100;
  }

  // Throw an exception if we failed
  if (fabs(lastg-gamma)>1.0e-3) {
    O2SCL_ERR("hlps_eos::fix_coeffs() failed.",exc_efailed);
  }

  return;
}

void eos_had_hlps::fix_neutron_matter
(double M, double local_n0, double Eneut, double dEneut) {
  double C=cbrt(4.0/local_n0/local_n0/pi2/pi2)/5.0/(gamma-1.0);
  etaL=C*(-10.0*Eneut*M+10.0*dEneut*M*local_n0+
	  cbrt(9.0*local_n0*local_n0*pi2*pi2))/cbrt(9.0);
  alphaL=C*(10.0*cbrt(3.0)*dEneut*M*local_n0-10.0*cbrt(3.0)*Eneut*M*gamma+
	    3.0*cbrt(local_n0*local_n0*pi2*pi2)*(3.0*gamma-2.0))/3.0;
  return;
}

void eos_had_hlps::fix_SL(double M, double local_n0, double S, double L) {
  double C=pow(1.5*local_n0*pi2,2.0/3.0)/4.0/M;
  etaL=(3.0*(L-3.0*S)+C*(2.0+9.0*(gamma-1.0)*eta))/(18.0*C*(gamma-1.0));
  alphaL=(C*(-4.0+9.0*alpha*(gamma-1.0)+6.0*gamma)+3.0*(L-3.0*S*gamma))/
    (18.0*C*(gamma-1.0));
  return;
}

int eos_had_hlps::calc_e(fermion &ln, fermion &lp, thermo &lth) {
  double xp, barn;
  
  barn=ln.n+lp.n;
  
  if (barn<=0.0) {
    xp=0.0;
    lth.ed=0.0;
    ln.mu=0.0;
    lp.mu=0.0;
    lth.pr=0.0;
    return 0;
  } else {
    xp=lp.n/barn;
  }

  double M=(ln.m+lp.m)/2.0;
  double T=cbrt(3.0*pi2*n0/2.0);
  T*=T/2.0/M;
  double x23=cbrt(xp)*cbrt(xp);
  double x53=x23*xp;
  double omx23=cbrt(1-xp)*cbrt(1-xp);
  double omx53=omx23*(1-xp);
  double u=barn/n0;
  double twou13=cbrt(u*2.0);
  
  // Energy per baryon
  double eoa=T*(0.6*(x53+omx53)*twou13*twou13-
		((2.0*alpha-4.0*alphaL)*xp*(1.0-xp)+alphaL)*u+
		((2.0*eta-4.0*etaL)*xp*(1.0-xp)+etaL)*pow(u,gamma));
  // Energy density
  lth.ed=ln.m*ln.n+lp.m*lp.n+barn*eoa;
  
  // Derivatives of energy per baryon with respect to u and x
  double deoadu=T*(0.6*(x53+omx53)*4.0/3.0/twou13-
		   ((2.0*alpha-4.0*alphaL)*xp*(1.0-xp)+alphaL)+
		   gamma*((2.0*eta-4.0*etaL)*xp*(1.0-xp)+etaL)*
		   pow(u,gamma-1.0));
  double deoadx=T*(twou13*twou13*(x23-omx23)-
		   (2.0*alpha-4.0*alphaL)*u*(1.0-2.0*xp)+
		   (2.0*eta-4.0*etaL)*pow(u,gamma)*(1.0-2.0*xp));
  
  // For chain rule
  double dxdnn=-lp.n/barn/barn;
  double dxdnp=-lp.n/barn/barn+1.0/barn;
  double dudnn=1.0/n0;
  double dudnp=1.0/n0;

  // Chemical potentials
  ln.mu=ln.m+(deoadu*barn+eoa*n0)*dudnn+deoadx*barn*dxdnn;
  lp.mu=lp.m+(deoadu*barn+eoa*n0)*dudnp+deoadx*barn*dxdnp;

  // Pressure
  lth.pr=-lth.ed+ln.mu*ln.n+lp.mu*lp.n;

  return 0;
}

