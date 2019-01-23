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

#include <o2scl/boson.h>

// For reading HDF5 files in calibrate()
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

//--------------------------------------------
// boson class

boson::boson(double mass, double dof) : part(mass,dof) {
  co=0.0;
}

void boson::massless_calc(double temper) {
  
  n=g*pow(temper,3.0)*zeta3/pi2;
  ed=g*pi2*pow(temper,4.0)/30.0;
  pr=ed/3.0;
  en=g*pi2*pow(temper,3.0)/22.5;

  return;
}

double boson_thermo::calibrate(boson &b, int verbose, bool test_pair,
			std::string fname) {
  
  double ret=0;

  // ----------------------------------------------------------------
  // Will return to these original values afterwards

  boson orig=b;

  // ----------------------------------------------------------------
  // Read data file

  if (fname=="") {
    fname=o2scl_settings.get_data_dir()+"boson_cal2.o2";
  }

  if (verbose>1) {
    cout << "In boson_T::calibrate(), loading file named\n\t" 
	 << fname << "\n" << endl;
  }
  table<> tab;
  hdf_file hf;
  hf.open(fname);
  string name;
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input(hf,tab,name);
#endif
  hf.close();

  tab.function_column("ed_mot*mot","ed");
  tab.function_column("pair_ed_mot*mot","pair_ed");
  
  table<> tab2;
  hf.open(fname);
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input(hf,tab2,name);
#endif
  hf.close();
  
  tab2.function_column("ed_mot*mot","ed");
  tab2.function_column("pair_ed_mot*mot","pair_ed");
  
  if (tab.get_nlines()==0) {
    string str="Failed to load data from file '"+fname+
      "' in boson_thermo::calibrate(). Bad filename?";
    O2SCL_ERR(str.c_str(),exc_efilenotfound);
  }
  
  // ----------------------------------------------------------------

  b.g=2.0;
  
  size_t cnt=0;
  part bad, dev, exact;
  double m_bad=0.0, mu_bad=0.0, T_bad=0.0, mot_bad=0.0, psi_bad=0.0;
  b.non_interacting=true;
  
  // ----------------------------------------------------------------
  // First pass, test calc_mu() 

  // k=0,2 are with rest mass, k=1,3 are without
  // k=0,1 are non-interacting, k=2,3 are interacting
  for(size_t k=0;k<4;k++) {

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
      
	if (k%2==0) {
	  
	  b.inc_rest_mass=true;

	  if (k>=2) {
	    b.non_interacting=false;
	    b.ms=mot*T;
	    b.m=b.ms*1.5;
	    b.nu=b.ms+T*psi;
	    b.mu=0.0;
	  } else {
	    b.non_interacting=true;
	    b.m=mot*T;
	    b.mu=b.m+T*psi;
	    b.nu=0.0;
	    b.ms=0.0;
	  }
	  
	} else {
	  
	  b.inc_rest_mass=false;
	  
	  if (k>=2) {
	    b.non_interacting=false;
	    b.ms=mot*T;
	    b.m=b.ms*1.5;
	    b.nu=T*psi-b.m+b.ms;
	    b.mu=0.0;
	  } else {
	    b.non_interacting=true;
	    b.m=mot*T;
	    b.mu=T*psi;
	    b.nu=0.0;
	    b.ms=0.0;
	  }
	  
	}
	
	calc_mu(b,T);
	
	exact.n*=pow(T,3.0);
	if (k%2==0) {
	  exact.ed*=pow(T,4.0);
	} else {
	  exact.ed=exact.ed*pow(T,4.0)-exact.n*b.m;
	}
	exact.pr*=pow(T,4.0);
	exact.en*=pow(T,3.0);
	
	dev.n+=fabs((b.n-exact.n)/exact.n);
	dev.ed+=fabs((b.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((b.pr-exact.pr)/exact.pr);
	dev.en+=fabs((b.en-exact.en)/exact.en);
	
	cnt++;
	if (fabs((b.n-exact.n)/exact.n)>bad.n) {
	  bad.n=fabs((b.n-exact.n)/exact.n);
	  if (bad.n>ret) {
	    if (k>=2) {
	      mu_bad=b.nu;
	      m_bad=b.ms;
	    } else {
	      mu_bad=b.mu;
	      m_bad=b.m;
	    }
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.n;
	  }
	}
	if (fabs((b.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((b.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    if (k>=2) {
	      mu_bad=b.nu;
	      m_bad=b.ms;
	    } else {
	      mu_bad=b.mu;
	      m_bad=b.m;
	    }
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((b.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((b.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    if (k>=2) {
	      mu_bad=b.nu;
	      m_bad=b.ms;
	    } else {
	      mu_bad=b.mu;
	      m_bad=b.m;
	    }
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((b.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((b.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    if (k>=2) {
	      mu_bad=b.nu;
	      m_bad=b.ms;
	    } else {
	      mu_bad=b.mu;
	      m_bad=b.m;
	    }
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}

	if (verbose>1) {
	  cout.precision(5);
	  if (k>=2) {
	    cout << "T,ms,nu,psi,mot: " << T << " " << b.ms << " " << b.nu
		 << " " << psi << " " << mot << endl;
	  } else {
	    cout << "T,m,mu,psi,mot: " << T << " " << b.m << " " << b.mu
		 << " " << psi << " " << mot << endl;
	  }
	  cout.precision(5);
	  cout << "n,ed,pr,en: " << endl;
	  cout << "approx: " << b.n << " " << b.ed << " " << b.pr << " " 
	       << b.en << endl;
	  cout << "exact : " << exact.n << " " << exact.ed << " " 
	       << exact.pr << " " << exact.en << endl;
	  cout << "bad   : " << bad.n << " " << bad.ed << " " 
	       << bad.pr << " " << bad.en << endl;
	  cout << endl;
	  if (verbose>2) {
	    char ch;
	    cin >> ch;
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
      } else if (k==1) {
	cout << "Function calc_mu(), without rest mass:" << endl;
      } else if (k==2) {
	cout << "Function calc_mu(), include rest mass, interacting:" << endl;
      } else {
	cout << "Function calc_mu(), without rest mass, interacting:" << endl;
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

    // Reset b.non_interacting
    b.non_interacting=true;
    
    // End of k loop
  }

#ifdef O2SCL_NEVER_DEFINED
  
  // ----------------------------------------------------------------
  // Second pass, test calc_density()

  // k=0,2 are with rest mass, k=1,3 are without
  // k=0,1 are non-interacting, k=2,3 are interacting
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
	b.n=tab.get("n",i);	
	exact.ed=tab.get("ed",i);
	exact.pr=tab.get("pr",i);
	exact.en=tab.get("en",i);

	if (k>=2) {
	  b.non_interacting=false;
	  b.ms=mot*T;
	  b.m=b.ms*1.5;
	} else {
	  b.non_interacting=true;
	  b.m=mot*T;
	  b.ms=0.0;
	}
	if (k%2==0) {
	  b.inc_rest_mass=true;
	  if (k>=2) {
	    exact.nu=b.m+T*psi;
	    exact.mu=0.0;
	  } else {
	    exact.mu=b.m+T*psi;
	    exact.nu=0.0;
	  }
	} else {
	  b.inc_rest_mass=false;
	  if (k>=2) {
	    exact.nu=T*psi-b.m+b.ms;
	    exact.mu=0.0;
	  } else {
	    exact.mu=T*psi;
	    exact.nu=0.0;
	  }
	}

	b.n*=pow(T,3.0);
	if (k==0) {
	  exact.ed*=pow(T,4.0);
	} else {
	  exact.ed=exact.ed*pow(T,4.0)-b.n*b.m;
	}
	exact.pr*=pow(T,4.0);
	exact.en*=pow(T,3.0);

	// Give it a guess for the chemical potential
	if (k>=2) {
	  b.nu=b.m;
	} else {
	  b.mu=b.m;
	}

	calc_density(f,T);
	
	if (k>=2) {
	  dev.nu+=fabs((b.nu-exact.nu)/exact.nu);
	} else {
	  dev.mu+=fabs((b.mu-exact.mu)/exact.mu);
	}
	dev.ed+=fabs((b.ed-exact.ed)/exact.ed);
	dev.pr+=fabs((b.pr-exact.pr)/exact.pr);
	dev.en+=fabs((b.en-exact.en)/exact.en);
	
	cnt++;
	if (k>=2) {
	  if (fabs((b.nu-exact.nu)/exact.nu)>bad.mu) {
	    bad.mu=fabs((b.nu-exact.nu)/exact.nu);
	    if (bad.mu>ret) {
	      mu_bad=b.nu;
	      m_bad=b.ms;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.n;
	    }
	  }
	} else {
	  if (fabs((b.mu-exact.mu)/exact.mu)>bad.mu) {
	    bad.mu=fabs((b.mu-exact.mu)/exact.mu);
	    if (bad.mu>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.n;
	    }
	  }
	}
	if (fabs((b.ed-exact.ed)/exact.ed)>bad.ed) {
	  bad.ed=fabs((b.ed-exact.ed)/exact.ed);
	  if (bad.ed>ret) {
	    mu_bad=b.mu;
	    m_bad=b.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.ed;
	  }
	}
	if (fabs((b.pr-exact.pr)/exact.pr)>bad.pr) {
	  bad.pr=fabs((b.pr-exact.pr)/exact.pr);
	  if (bad.pr>ret) {
	    mu_bad=b.mu;
	    m_bad=b.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.pr;
	  }
	}
	if (fabs((b.en-exact.en)/exact.en)>bad.en) {
	  bad.en=fabs((b.en-exact.en)/exact.en);
	  if (bad.en>ret) {
	    mu_bad=b.mu;
	    m_bad=b.m;
	    T_bad=T;
	    mot_bad=mot;
	    psi_bad=psi;
	    ret=bad.en;
	  }
	}

	if (verbose>1) {
	  cout.precision(5);
	  if (k>=2) {
	    cout << "T,ms,n,psi,mot: " << T << " " << b.ms << " " << b.n
		 << " " << psi << " " << mot << endl;
	  } else {
	    cout << "T,m,n,psi,mot: " << T << " " << b.m << " " << b.n
		 << " " << psi << " " << mot << endl;
	  }
	  cout.precision(6);
	  cout << "mu,ed,pr,en: " << endl;
	  cout << "approx: " << b.mu << " " << b.ed << " " << b.pr << " " 
	       << b.en << endl;
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
      } else if (k==1) {
	cout << "Function calc_density(), without rest mass:" << endl;
      } else if (k==2) {
	cout << "Function calc_density(), include rest mass, interacting:"
	     << endl;
      } else {
	cout << "Function calc_density(), without rest mass, interacting:"
	     << endl;
      }

      cout << "Average performance: " << endl;
      cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
	   << dev.pr << " en: " << dev.en << endl;
      cout << "Worst case: " << endl;
      cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
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

  if (test_pair) {
  
    // ----------------------------------------------------------------
    // Third pass, test pair_mu() 

    // k=0,2 are with rest mass, k=1,3 are without
    // k=0,1 are non-interacting, k=2,3 are interacting
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
	  exact.n=tab2.get("pair_n",i);
	  exact.ed=tab2.get("pair_ed",i);
	  exact.pr=tab2.get("pair_pr",i);
	  exact.en=tab2.get("pair_en",i);
      
	  if (k%2==0) {

	    b.inc_rest_mass=true;

	    if (k>=2) {
	      b.non_interacting=false;
	      b.ms=mot*T;
	      b.m=b.ms*1.5;
	      b.nu=b.m+T*psi;
	      b.mu=0.0;
	    } else {
	      b.non_interacting=true;
	      b.m=mot*T;
	      b.mu=b.m+T*psi;
	      b.nu=0.0;
	      b.ms=0.0;
	    }
	  
	  } else {
	  
	    b.inc_rest_mass=false;
	  
	    if (k>=2) {
	      b.non_interacting=false;
	      b.ms=mot*T;
	      b.m=b.ms*1.5;
	      b.nu=b.ms+T*psi-b.m;
	      b.mu=0.0;
	    } else {
	      b.non_interacting=true;
	      b.m=mot*T;
	      b.mu=T*psi;
	      b.nu=0.0;
	      b.ms=0.0;
	    }
	  
	  }
	
	  pair_mu(f,T);
	
	  exact.n*=pow(T,3.0);
	  if (k==0) {
	    exact.ed*=pow(T,4.0);
	  } else {
	    exact.ed=exact.ed*pow(T,4.0)-exact.n*b.m;
	  }
	  exact.pr*=pow(T,4.0);
	  exact.en*=pow(T,3.0);
	
	  dev.n+=fabs((b.n-exact.n)/exact.n);
	  dev.ed+=fabs((b.ed-exact.ed)/exact.ed);
	  dev.pr+=fabs((b.pr-exact.pr)/exact.pr);
	  dev.en+=fabs((b.en-exact.en)/exact.en);
	
	  cnt++;
	  if (fabs((b.n-exact.n)/exact.n)>bad.n) {
	    bad.n=fabs((b.n-exact.n)/exact.n);
	    if (bad.n>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.n;
	    }
	  }
	  if (fabs((b.ed-exact.ed)/exact.ed)>bad.ed) {
	    bad.ed=fabs((b.ed-exact.ed)/exact.ed);
	    if (bad.ed>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.ed;
	    }
	  }
	  if (fabs((b.pr-exact.pr)/exact.pr)>bad.pr) {
	    bad.pr=fabs((b.pr-exact.pr)/exact.pr);
	    if (bad.pr>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.pr;
	    }
	  }
	  if (fabs((b.en-exact.en)/exact.en)>bad.en) {
	    bad.en=fabs((b.en-exact.en)/exact.en);
	    if (bad.en>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.en;
	    }
	  }

	  if (verbose>1) {
	    cout.precision(5);
	    cout << "T,m,mu,psi,mot: " << T << " " << b.m << " " << b.mu
		 << " " << psi << " " << mot << endl;
	    cout.precision(6);
	    cout << "n,ed,pr,en: " << endl;
	    cout << "approx: " << b.n << " " << b.ed << " " << b.pr << " " 
		 << b.en << endl;
	    cout << "exact : " << exact.n << " " << exact.ed << " " 
		 << exact.pr << " " << exact.en << endl;
	    cout << "bad   : " << bad.n << " " << bad.ed << " " 
		 << bad.pr << " " << bad.en << endl;
	    cout << endl;
	    if (verbose>2) {
	      char ch;
	      cin >> ch;
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
	  cout << "Function pair_mu(), include rest mass:" << endl;
	} else {
	  cout << "Function pair_mu(), without rest mass:" << endl;
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
    // Fourth pass, test pair_density()

    // k=0,2 are with rest mass, k=1,3 are without
    // k=0,1 are non-interacting, k=2,3 are interacting
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
	  b.n=tab2.get("pair_n",i);	
	  exact.ed=tab2.get("pair_ed",i);
	  exact.pr=tab2.get("pair_pr",i);
	  exact.en=tab2.get("pair_en",i);
	  
	  if (k>=2) {
	    b.non_interacting=false;
	    b.ms=mot*T;
	    b.m=b.ms*1.5;
	  } else {
	    b.non_interacting=true;
	    b.m=mot*T;
	    b.ms=0.0;
	  }
	  if (k%2==0) {
	    b.inc_rest_mass=true;
	    if (k>=2) {
	      exact.nu=b.m+T*psi;
	      exact.mu=0.0;
	    } else {
	      exact.mu=b.m+T*psi;
	      exact.nu=0.0;
	    }
	  } else {
	    b.inc_rest_mass=false;
	    if (k>=2) {
	      exact.nu=T*psi-b.m+b.ms;
	      exact.mu=0.0;
	    } else {
	      exact.mu=T*psi;
	      exact.nu=0.0;
	    }
	  }

	  b.n*=pow(T,3.0);
	  if (k==0) {
	    exact.ed*=pow(T,4.0);
	  } else {
	    exact.ed=exact.ed*pow(T,4.0)-b.n*b.m;
	  }
	  exact.pr*=pow(T,4.0);
	  exact.en*=pow(T,3.0);

	  // Give it a guess for the chemical potential
	  b.mu=b.m;

	  pair_density(f,T);
	
	  dev.mu+=fabs((b.mu-exact.mu)/exact.mu);
	  dev.ed+=fabs((b.ed-exact.ed)/exact.ed);
	  dev.pr+=fabs((b.pr-exact.pr)/exact.pr);
	  dev.en+=fabs((b.en-exact.en)/exact.en);
	
	  cnt++;
	  if (fabs((b.mu-exact.mu)/exact.mu)>bad.mu) {
	    bad.mu=fabs((b.mu-exact.mu)/exact.mu);
	    if (bad.n>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.n;
	    }
	  }
	  if (fabs((b.ed-exact.ed)/exact.ed)>bad.ed) {
	    bad.ed=fabs((b.ed-exact.ed)/exact.ed);
	    if (bad.ed>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.ed;
	    }
	  }
	  if (fabs((b.pr-exact.pr)/exact.pr)>bad.pr) {
	    bad.pr=fabs((b.pr-exact.pr)/exact.pr);
	    if (bad.pr>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.pr;
	    }
	  }
	  if (fabs((b.en-exact.en)/exact.en)>bad.en) {
	    bad.en=fabs((b.en-exact.en)/exact.en);
	    if (bad.en>ret) {
	      mu_bad=b.mu;
	      m_bad=b.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret=bad.en;
	    }
	  }

	  if (verbose>1) {
	    cout.precision(5);
	    cout << "T,m,n,psi,mot: " << T << " " << b.m << " " << b.n
		 << " " << psi << " " << mot << endl;
	    cout.precision(6);
	    cout << "mu,ed,pr,en: " << endl;
	    cout << "approx: " << b.mu << " " << b.ed << " " << b.pr << " " 
		 << b.en << endl;
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
	  cout << "Function pair_density(), include rest mass:" << endl;
	} else {
	  cout << "Function pair_density(), without rest mass:" << endl;
	}

	cout << "Average performance: " << endl;
	cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
	     << dev.pr << " en: " << dev.en << endl;
	cout << "Worst case: " << endl;
	cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
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

    // End of 'if (test_pair)'
  }

  // ----------------------------------------------------------------
  // Return to the original values 

  f=orig;

  cout << "Calibration result: " << ret << "\n" << endl;

#endif
  
  return ret;
}
