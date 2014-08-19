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

#include <o2scl/fermion.h>
// For gsl_hypot()
#include <gsl/gsl_sys.h>

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

fermion::fermion(double mass, double dof) : part(mass,dof) {
  kf=0.0;
  del=0.0;
}

void fermion_zerot::kf_from_density(fermion &f) {
  f.kf=cbrt(6.0*pi2/f.g*f.n);
  return;
}

void fermion_zerot::energy_density_zerot(fermion &f) {
  double r,efs;
  if (f.kf>0.0) {
    if (f.ms<=0.0) {
      f.ed=f.g*(pow(f.kf,4.0)/8.0/pi2);
    } else {
      efs=gsl_hypot(f.kf,f.ms);
      r=(f.kf+efs)/f.ms;
      f.ed=f.g/16.0/pi2*(2.0*f.kf*pow(efs,3.0)-f.kf*efs*f.ms*f.ms
			 -pow(f.ms,4.0)*log(r));
    }
  } else {
    f.ed=0.0;
  }
  return;
}

void fermion_zerot::pressure_zerot(fermion &f) {
  double r,efs;
  if (f.kf>0.0) {
    if (f.ms<=0.0) {
      f.pr=f.g*(pow(f.kf,4.0)/24.0/pi2);
    } else {
      efs=gsl_hypot(f.kf,f.ms);
      r=(f.kf+efs)/f.ms;
      f.pr=f.g/48.0/pi2*(2.0*efs*pow(f.kf,3.0)-3.0*f.kf*efs*f.ms*f.ms
			 +3.0*pow(f.ms,4.0)*log(r));
    }
  } else {
    f.pr=0.0;
  }
  return;
}

fermion_eval_thermo::fermion_eval_thermo() {
  massless_root=&def_massless_root;
}

void fermion_zerot::calc_mu_zerot(fermion &f) {
  bool nulessthan0;
  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
  if (f.inc_rest_mass) {
    if (f.nu>f.ms) {
      nulessthan0=false;
      f.kf=sqrt(f.nu*f.nu-f.ms*f.ms);
    } else {
      nulessthan0=false;
      f.kf=0.0;
    }
  } else {
    double nupm=f.nu+f.m;
    if ((nupm)>f.ms) {
      nulessthan0=false;
      f.kf=sqrt(nupm*nupm-f.ms*f.ms);
    } else {
      nulessthan0=false;
      f.kf=0.0;
    }
  }
  f.n=f.g/6.0/pi2*pow(f.kf,3.0);
  energy_density_zerot(f);
  pressure_zerot(f);
  f.en=0.0;
  if (!f.inc_rest_mass) f.ed-=f.n*f.m;
  if (nulessthan0==true) {
    f.n*=-1.0;
    f.kf*=-1.0;
  }

  return;
}

void fermion_zerot::calc_density_zerot(fermion &f) {
  if (f.non_interacting) { f.ms=f.m; }

  f.kf=cbrt(6.0*pi2/f.g*f.n);
  f.nu=gsl_hypot(f.kf,f.ms);
  f.n=f.g/6.0/pi2*pow(f.kf,3.0);
  energy_density_zerot(f);
  pressure_zerot(f);
  f.en=0.0;

  if (!f.inc_rest_mass) {
    f.nu-=f.m;
    f.ed-=f.n*f.m;
  }

  if (f.non_interacting) { f.mu=f.nu; }

  return;
}

void fermion_eval_thermo::massless_calc_mu(fermion &f, double temper) {
  
  double fm2, fm3;

  if (f.non_interacting) { f.nu=f.mu; }

  fm2=gsl_sf_fermi_dirac_int(2,f.nu/temper);
  fm3=gsl_sf_fermi_dirac_int(3,f.nu/temper);
  
  f.n=f.g/pi2*pow(temper,3.0)*fm2;
  f.ed=f.g*3.0/pi2*pow(temper,4.0)*fm3;
  f.pr=f.ed/3.0;
  f.en=(f.ed+f.pr-f.n*f.nu)/temper;

  return;
}

double fermion_eval_thermo::massless_fun::operator()(double x) {
  double fm2;
  
  fm2=gsl_sf_fermi_dirac_int(2,x/(temper_));
  double ret=f_.g*pow(temper_,3.0)*fm2/pi2/f_.n-1.0;
  return ret;
}

void fermion_eval_thermo::massless_calc_density(fermion &f, double temper) {
  double x, T=temper;
  
  x=f.ms+temper;
  massless_fun mf(f,temper);
  funct11 mf2=std::bind(std::mem_fn<double(double)>
			(&massless_fun::operator()),
			&mf,std::placeholders::_1);
  massless_root->solve(x,mf2);
  f.nu=x;

  massless_calc_mu(f,temper);

  // If the particle is non-interacting, then need to set
  // mu=nu to get the entropy right
  if (f.non_interacting) { f.mu=f.nu; }

  return;
}

void fermion_eval_thermo::massless_pair_mu(fermion &f, double temper) {
  double pitmu, pitmu2, nu2;

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  if (f.nu==0.0) {
    f.n=0.0;
    f.ed=f.g/8.0/pi2*7.0/15.0*pow(pi*temper,4.0);
    f.pr=f.ed/3.0;
    f.en=(f.ed+f.pr-f.n*f.mu)/temper;
  } else {
    nu2=f.nu*f.nu;
    pitmu=pi*temper/f.nu;
    pitmu2=pitmu*pitmu;
    f.n=f.g*f.nu*nu2/6.0/pi2*(1.0+pitmu2);
    f.ed=f.g*nu2*nu2/8.0/pi2*(1.0+2.0*pitmu2+7.0/15.0*pitmu2*pitmu2);
    f.pr=f.ed/3.0;
    f.en=(f.ed+f.pr-f.n*f.mu)/temper;
    
    // Might the following work better for the energy density?
    // pit=pi*temper;
    // pit2=pit*pit;
    // ed=g/8.0/pi2*(nu2*nu2+2.0*pit2*nu2+7.0/15.0*pit2*pit2);
    
  }

  return;
}

void fermion_eval_thermo::massless_pair_density(fermion &f, double temper) {

  double t2=temper*temper,pitmu,pitmu2,nu2;
  double cbt, alpha, two13, alpha16;

  if (f.non_interacting) { f.ms=f.m; }
  if (f.n<=0.0) {
    f.nu=0.0;
    f.ed=f.g/8.0/pi2*7.0/15.0*pow(pi*temper,4.0);
    f.pr=f.ed/3.0;
  } else {
    alpha=f.g*f.g*pi2*t2*t2*t2/243.0/f.n/f.n;
    if (alpha>1.0e4) {
      f.nu=(2.0/3.0/sqrt(alpha)-8.0/81.0/pow(alpha,1.5)+
	    32.0/729.0/pow(alpha,2.5))*pi*temper/sqrt(3.0);
    } else if (alpha<3.0e-4) {
      two13=cbrt(2.0);
      alpha16=pow(alpha,1.0/6.0);
      f.nu=(two13/alpha16-alpha16/two13+alpha/alpha16/6.0/two13/two13
	    +alpha*alpha16/12.0/two13-alpha*alpha/alpha16/18.0/two13/two13-
	    5.0*alpha*alpha*alpha16/144.0/two13+
	    77.0/2592.0*alpha*alpha*alpha/alpha16/two13/two13)*
	pi*temper/sqrt(3.0);
    } else {
      cbt=pow(-1.0+sqrt(1.0+alpha),1.0/3.0)/pow(alpha,1.0/6.0);
      f.nu=pi*temper/sqrt(3.0)*(1.0/cbt-cbt);
    }
    pitmu=pi*temper/f.nu;
    pitmu2=pitmu*pitmu;
    nu2=f.nu*f.nu;
    f.ed=f.g*nu2*nu2/8.0/pi2*(1.0+2.0*pitmu2+7.0/15.0*pitmu2*pitmu2);
    f.pr=f.ed/3.0;

    if (!o2scl::is_finite(f.nu)) {
      string str="Chemical potential not finite ("+dtos(f.nu)+
	") in fermion::massless_pair_density().";
      O2SCL_ERR(str.c_str(),exc_efailed);
    }
  }

  if (f.non_interacting) { f.mu=f.nu; }
  f.en=(f.ed+f.pr-f.n*f.nu)/temper;

  return;
}

double fermion_eval_thermo::calibrate
(fermion &f, int verbose, std::string fname) {

  double ret=0;
  
  // ----------------------------------------------------------------
  // Will return to these original values afterwards

  part orig;
  orig.mu=f.mu;
  orig.m=f.m;
  orig.ms=f.ms;
  orig.g=f.g;
  orig.non_interacting=f.non_interacting;
  orig.inc_rest_mass=f.inc_rest_mass;

  // ----------------------------------------------------------------
  // Read data file

  if (fname=="") {
    fname=o2scl_settings.get_data_dir()+"fermion_cal.o2";
  }

  if (verbose>1) {
    cout << "In fermion_T::calibrate(), loading file named\n\t" 
	 << fname << "\n" << endl;
  }
  table<> tab;
  hdf_file hf;
  hf.open(fname);
  string name;
  hdf_input(hf,tab,name);
  hf.close();
  
  if (tab.get_nlines()==0) {
    string str="Failed to load data from file '"+fname+
      "' in fermion_T::calibrate(). Bad filename?";
    O2SCL_ERR(str.c_str(),exc_efilenotfound);
  }
  
  // ----------------------------------------------------------------

  f.g=2.0;
  
  size_t cnt=0;
  part bad, dev, exact;
  double m_bad=0.0, mu_bad=0.0, T_bad=0.0, mot_bad=0.0, psi_bad=0.0;
  f.non_interacting=true;
  
  // ----------------------------------------------------------------
  // First pass, test calc_mu() 

  // k=0 is with rest mass, k=1 is without
  for(size_t k=0;k<2;k++) {

    // Initialize storage
    dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
    bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
    // Temperature loop
    for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {

      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("mot",i);
	double psi=tab.get("psi",i);
	exact.n=tab.get("n",i);
	exact.ed=tab.get("ed",i);
	exact.pr=tab.get("pr",i);
	exact.en=tab.get("en",i);
      
	if (k==0) {
	  
	  f.inc_rest_mass=true;
          
	  f.m=mot*T;
	  f.mu=f.m+T*psi;
	  
	} else {
	  
	  f.inc_rest_mass=false;
	  
	  f.m=mot*T;
	  f.mu=T*psi;
	  
	}
	
	calc_mu(f,T);
	
	exact.n*=pow(T,3.0);
	if (k==0) {
	  exact.ed*=pow(T,4.0);
	} else {
	  exact.ed=exact.ed*pow(T,4.0)-exact.n*f.m;
	}
	exact.pr*=pow(T,4.0);
	exact.en*=pow(T,3.0);
	
	if (verbose>1) {
	  cout << "T,m,mu: " << T << " " << f.m << " " << f.mu << endl;
	  cout << "n,ed,pr,en: " << endl;
	  cout << "approx: " << f.n << " " << f.ed << " " << f.pr << " " 
	       << f.en << endl;
	  cout << "exact : " << exact.n << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "bad   : " << bad.mu << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << endl;
	  if (verbose>2) {
	    char ch;
	    cin >> ch;
	  }
	}

	dev.n+=fabs((f.n-exact.n)/exact.n);
	dev.ed+=fabs((f.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((f.pr-exact.pr)/exact.pr);
	dev.en+=fabs((f.en-exact.en)/exact.en);
	
	cnt++;
	if (fabs((f.n-exact.n)/exact.n)>bad.n) {
	  bad.n=fabs((f.n-exact.n)/exact.n);
	  if (bad.n>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.n;
	  }
	}
	if (fabs((f.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((f.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((f.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((f.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((f.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((f.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}

	// End of loop over points in data file
      }
      // End of temperature loop
    }

    dev.n/=cnt;
    dev.ed/=cnt;
    dev.pr/=cnt;
    dev.en/=cnt;

    if (verbose>0) {
      if (k==0) {
	cout << "Function calc_mu(), include rest mass:" << endl;
      } else {
	cout << "Function calc_mu(), without rest mass:" << endl;
      }

      cout << "Average performance: " << endl;
      cout << "n: " << dev.n << " ed: " << dev.ed << " pr: " 
	   << dev.pr << " en: " << dev.en << endl;
      cout << "Worst case: " << endl;
      cout << "n: " << bad.n << " ed: " << bad.ed << " pr: " 
	   << bad.pr << " en: " << bad.en << endl;
      cout << "mu: " << mu_bad << " m: " << m_bad << " T: " << T_bad 
	   << " mot: " << mot_bad << "\n\tpsi: " << psi_bad << endl;
      cout << endl;
      if (verbose>2) {
	char ch;
	cin >> ch;
      }
    }

    // End of k loop
  }

  // ----------------------------------------------------------------
  // Second pass, test calc_density()

  // k=0 is with rest mass, k=1 is without
  for(size_t k=0;k<2;k++) {

    // Initialize storage
    dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
    bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
    // Temperature loop
    for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
      // Loop over each point in the data file
      for(size_t i=0;i<tab.get_nlines();i++) {
	
	double mot=tab.get("mot",i);
	double psi=tab.get("psi",i);
	f.n=tab.get("n",i);	
	exact.ed=tab.get("ed",i);
	exact.pr=tab.get("pr",i);
	exact.en=tab.get("en",i);

	f.m=mot*T;
	if (k==0) {
	  f.inc_rest_mass=true;
	  exact.mu=f.m+T*psi;
	} else {
	  f.inc_rest_mass=false;
	  exact.mu=T*psi;
	}

	f.n*=pow(T,3.0);
	if (k==0) {
	  exact.ed*=pow(T,4.0);
	} else {
	  exact.ed=exact.ed*pow(T,4.0)-f.n*f.m;
	}
	exact.pr*=pow(T,4.0);
	exact.en*=pow(T,3.0);

	// Give it a guess for the chemical potential
	f.mu=f.m;

	calc_density(f,T);
	
	if (verbose>1) {
	  cout << "T, m, n: " << T << " " << f.m << " " << f.n << endl;
	  cout << "mu,ed,pr,en: " << endl;
	  cout << "approx: " << f.mu << " " << f.ed << " " << f.pr << " " 
	       << f.en << endl;
	  cout << "exact : " << exact.mu << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "bad   : " << bad.mu << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << endl;
	  if (verbose>2) {
	    char ch;
	    cin >> ch;
	  }
	}

	dev.mu+=fabs((f.mu-exact.mu)/exact.mu);
	dev.ed+=fabs((f.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((f.pr-exact.pr)/exact.pr);
	dev.en+=fabs((f.en-exact.en)/exact.en);
	
	cnt++;
	if (fabs((f.mu-exact.mu)/exact.mu)>bad.mu) {
	  bad.mu=fabs((f.mu-exact.mu)/exact.mu);
	  if (bad.n>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.n;
	  }
	}
	if (fabs((f.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((f.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((f.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((f.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((f.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((f.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    mu_bad=f.mu;
	    m_bad=f.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}

	// End of loop over points in data file
      }
      // End of temperature loop
    }

    dev.mu/=cnt;
    dev.ed/=cnt;
    dev.pr/=cnt;
    dev.en/=cnt;

    if (verbose>0) {
      if (k==0) {
	cout << "Function calc_density(), include rest mass:" << endl;
      } else {
	cout << "Function calc_density(), without rest mass:" << endl;
      }

      cout << "Average performance: " << endl;
      cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
	   << dev.pr << " en: " << dev.en << endl;
      cout << "Worst case: " << endl;
      cout << "mu: " << mu_bad << " m: " << m_bad << " T: " << T_bad 
	   << " mot: " << mot_bad << "\n\tpsi: " << psi_bad << endl;
      cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
	   << bad.pr << " en: " << bad.en << endl;
      cout << endl;
      if (verbose>2) {
	char ch;
	cin >> ch;
      }
    }

    // End of k loop
  }

  // ----------------------------------------------------------------
  // Return to the original values 

  f.mu=orig.mu;
  f.m=orig.m;
  f.ms=orig.ms;
  f.g=orig.g;
  f.non_interacting=orig.non_interacting;
  f.inc_rest_mass=orig.inc_rest_mass;

  cout << "Calibration result: " << ret << "\n" << endl;
  cout << endl;
  
  return ret;
}

bool fermion_eval_thermo::calc_mu_deg(fermion &f, double temper, 
				      double prec) {
  
  // Handle the zero-temperature limit
  if (temper==0.0) {
    calc_mu_zerot(f);
    return true;
  }

  // Double check to ensure T and mass are positive
  if (temper<0.0 || f.ms<0.0) {
    O2SCL_ERR2("Temperature or mass negative in fermion_eval_thermo",
	       "::calc_mu_deg().",exc_einval);
  }
  
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
  // Compute psi and tt
  double psi;
  if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
  else psi=(f.nu+f.m-f.ms)/temper;
  double tt=temper/f.ms;

  // Return false immediately psi<0 where the expressions below
  // don't work because of the square roots
  if (psi<0.0) return false;
  
  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);
  
  // Define x = psi * t = (mu/m - 1) and related values
  double x=psi*tt;
  double sx=sqrt(x);
  double s2x=sqrt(2.0+x);
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  
  // Evaluate the first and last term for the pressure
  double pterm1;
  if (x>1.0e-5) {
    pterm1=(x*(1.0+x)*(2.0+x)*(-3.0+2.0*x*(2.0+x))+6.0*sx*s2x*
	    log((sx+s2x)/sqrt(2.0)))/24.0/sx/s2x;
  } else {
    pterm1=x2*sx*(29568.0+15840.0*x+1540.0*x2-105.0*x3)/55440.0/sqrt(2.0);
  }
  double pterm4=-31.0*pow(pi*tt,6.0)/1008.0*(1.0+x)*
    sx*s2x/pow(x*(2.0+x),4.0);

  // Check if we're going to succeed
  if (fabs(pterm4)/fabs(pterm1)>prec) {
    return false;
  }
  
  // First order density term (first order entropy term is zero)
  double nterm1=sx*s2x*x*(2.0+x)/3.0/f.ms;
  
  // Second order terms
  double pterm2=tt*tt*pi2/6.0*(1.0+x)*sx*s2x;
  double nterm2=tt*tt*pi2/6.0*(1.0+4.0*x+2.0*x2)/
    f.ms/sx/s2x;
  double enterm2=tt*pi2/3.0*(1.0+x)*sx*s2x/f.ms;

  // Third order terms
  double pterm3=7.0*pow(pi*tt,4.0)/360.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/pow(x*(2.0+x),1.5);
  double nterm3=7.0*pow(pi*tt,4.0)/120.0/sx/s2x/
    x2/(x+2.0)/(x+2.0)/f.ms;
  double enterm3=7.0*pow(pi*tt,4.0)/tt/90.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/f.ms/sx/s2x/x/(x+2.0);

  // Fourth order terms for density and entropy
  double nterm4=31.0*pow(pi*tt,6.0)/1008.0*sx*s2x*
    (7.0+12.0*x+6.0*x2)/f.ms/pow(x*(2.0+x),5.0);
  double enterm4=-31.0*pow(pi*tt,6.0)/tt/168.0*sx*s2x*
    (1.0+x)/pow(x*(2.0+x),4.0);

  // Add up all the terms
  f.pr=prefac*(pterm1+pterm2+pterm3+pterm4);
  f.n=prefac*(nterm1+nterm2+nterm3+nterm4);
  f.en=prefac*(enterm2+enterm3+enterm4);
  f.ed=-f.pr+f.nu*f.n+temper*f.en;

  return true;
}

bool fermion_eval_thermo::calc_mu_ndeg(fermion &f, double temper, 
				       double prec) {

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

  // Compute psi and tt
  double psi;
  if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
  else psi=(f.nu+f.m-f.ms)/temper;
  double tt=temper/f.ms;

  // Return false immediately if we're degenerate
  if (psi>-1.0) return false;

  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);

  // One term is always used, so only values of max_term greater than
  // 0 are useful.
  static const size_t max_term=20;
  
  // Maximum argument for exponential
  // double log_dbl_max=709.78;

  // Return zero if psi+1/t is too small
  if (psi+1.0/tt<-700.0) {
    f.n=0.0;
    f.ed=0.0;
    f.pr=0.0;
    f.en=0.0;
    return true;
  }
  
  if (max_term/tt<700.0) {
    // Ratio of last term to first term
    double rat=exp((max_term-1.0)*(psi+1.0/tt))*
      gsl_sf_bessel_Kn(2,max_term/tt)/
      gsl_sf_bessel_Kn(2,1.0/tt)/pow(((double)max_term),2.0);
    
    // If the ratio between the last term and the first term is 
    // not small enough, return false
    if (rat>prec) {
      return false;
    }
  }
    
  double first_term=0.0;
  f.pr=0.0;
  f.n=0.0;
  f.en=0.0;

  for(size_t j=1;j<=max_term;j++) {
    double dj=((double)j);
    double pterm, nterm, enterm;

    double novert=dj/tt;

    // If n/t is large, we have to use an expansion for
    // the bessel functions and part of the exponential
    if (novert>700.0) {
      double novert2=novert*novert;
      double novert3=novert2*novert;
      double novert4=novert2*novert2;
      double novert5=novert3*novert2;
      double novert6=novert3*novert3;
      double novert7=novert4*novert3;
      double novert8=novert4*novert4;
      double novert9=novert5*novert4;
      double novert10=novert5*novert5;
      double novert11=novert6*novert5;
      double novert12=novert6*novert6;
      pterm=exp(dj*psi)/novert/novert*sqrt(pi/2.0/novert);
      pterm*=(1.0+15.0/8.0/novert+105.0/128.0/novert2-
	      315.0/1024.0/novert3+10395.0/32768.0/novert4-
	      135135.0/262144.0/novert5+4729725.0/4194304.0/novert6-
	      103378275.0/33554432.0/novert7+
	      21606059475.0/2147483648.0/novert8-
	      655383804075.0/17179869184.0/novert9+
	      45221482481175.0/274877906944.0/novert10-
	      1747193641318125.0/2199023255552.0/novert11+
	      298770112665399375.0/70368744177664.0/novert12);
      if (j%2==0) pterm*=-1.0;
      nterm=pterm*dj/temper;
      double ebk1=1.0+3.0/8.0/novert-15.0/128.0/novert2+
	105.0/1024.0/novert3-4725.0/32768.0/novert4+
	72765.0/262144.0/novert5-2837835.0/4194304.0/novert6+
	66891825.0/33554432.0/novert7-14783093325.0/2147483648.0/novert8+
	468131288625.0/17179869184.0/novert9-
	33424574007825.0/274877906944.0/novert10+
	1327867167401775.0/2199023255552.0/novert11-
	232376754295310625.0/70368744177664.0/novert12;
      double ebk3=1.0+35.0/8.0/novert+945.0/128.0/novert2+
	3465.0/1024.0/novert3-45045.0/32768.0/novert4+
	405405.0/262144.0/novert5-11486475.0/4194304.0/novert6+
	218243025.0/33554432.0/novert7-41247931725.0/2147483648.0/novert8+
	1159525191825.0/17179869184.0/novert9-
	75369137468625.0/274877906944.0/novert10+
	2774954606799375.0/2199023255552.0/novert11-
	456017540384030625.0/70368744177664.0/novert12;
      if (j%2==0) {
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj-
		exp(dj*psi)*sqrt(pi/2.0/novert)*
		(ebk1+ebk3)/2.0/dj)/f.ms;
      } else {
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj+
		exp(dj*psi)*sqrt(pi/2.0/novert)*
		(ebk1+ebk3)/2.0/dj)/f.ms;
      }
    } else {
      pterm=exp(dj*(psi+1.0/tt))/novert/novert*gsl_sf_bessel_Kn(2,novert);
      if (j%2==0) pterm*=-1.0;
      nterm=pterm*dj/temper;
      if (j%2==0) {
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj-
		exp(dj*(psi+1.0/tt))/2.0/dj*
		(gsl_sf_bessel_K1(dj/tt)+
		 gsl_sf_bessel_Kn(3,dj/tt)))/f.ms;
      } else {
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj+
		exp(dj*(psi+1.0/tt))/2.0/dj*
		(gsl_sf_bessel_K1(dj/tt)+
		 gsl_sf_bessel_Kn(3,dj/tt)))/f.ms;
      }
    }

    if (j==1) first_term=pterm;
    f.pr+=pterm;
    f.n+=nterm;
    f.en+=enterm;

    // If the first term is zero, then the rest of the terms
    // will be zero so just return early
    if (first_term==0.0) {
      f.pr=0.0;
      f.n=0.0;
      f.ed=0.0;
      f.en=0.0;
      return true;
    }

    // Stop early if we can
    if (j>1 && fabs(pterm)<prec*fabs(first_term)) {
      f.pr*=prefac;
      f.n*=prefac;
      f.en*=prefac;
      f.ed=-f.pr+f.nu*f.n+temper*f.en;
      return true;
    } else if (j==max_term) {
      // The last iteration failed to converge
      return false;
    }
  }

  O2SCL_ERR2("Loop exited in fermion_eval_thermo::",
	     "calc_mu_ndeg().",exc_esanity);
  return false;
}
