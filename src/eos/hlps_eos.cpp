#include <o2scl/hlps_eos.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

// --------------------------------------------------------------

hlps_eos::hlps_eos() {

#ifdef O2SCL_NEVER_DEFINED
  
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

  cout << it.eval(0.16,25,lower_nB,lower_E) << endl;
  cout << it.eval(0.16,22,upper_nB,upper_E) << endl;
  cout << it.deriv(0.16,25,lower_nB,lower_E) << endl;
  cout << it.deriv(0.16,22,upper_nB,upper_E) << endl;

  //1.349511e+01
  //2.011534e+01
  //9.534376e+01
  //1.382914e+02
    
  exit(-1);

#endif

  // Some default values from the paper
  kprime=0.0;
  n0=0.16;
  gamma=4.0/3.0;
  alpha=5.87;
  eta=3.81;
  alphaL=1.4;
  etaL=0.86;
}

void hlps_eos::fix_coeffs(double M, double B, double K) {

  double lastg;
  for(size_t i=0;i<100;i++) {

    // First fix eta and alpha from B, M, and n0
    double p=20.0*cbrt(3.0)*B*M/cbrt(n0)/cbrt(n0);
    double q=3.0*cbrt(2.0*pi2*pi2);
    double C=cbrt(4.0/pi2/pi2)/15.0/(gamma-1.0);
    eta=C*(q-p);
    alpha=0.8+C*(q-p)*gamma;
    lastg=gamma;

    // Now compute new gamma from K, M and n0
    double C2=4.0*K*M/3.0*cbrt(4.0/9.0/n0/n0/pi2/pi2);
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
    O2SCL_ERR("hlpe_eos::fix_coeffs() failed.",exc_efailed);
  }

  return;
}

void hlps_eos::fix_neutron_matter(double M, double Eneut, 
				  double dEneut) {
  double C=cbrt(4.0/n0/n0/pi2/pi2)/5.0/(gamma-1.0);
  etaL=C*(-10.0*Eneut*M+10.0*dEneut*M*n0+cbrt(9.0*n0*n0*pi2*pi2))/cbrt(9.0);
  alphaL=C*(10.0*cbrt(3.0)*dEneut*M*n0-10.0*cbrt(3.0)*Eneut*M*gamma+
	    3.0*cbrt(n0*n0*pi2*pi2)*(3.0*gamma-2.0))/3.0;
  return;
}

int hlps_eos::calc_e(fermion &ln, fermion &lp, thermo &lth) {
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

