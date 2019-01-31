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
#ifndef O2SCL_PART_H
#define O2SCL_PART_H

#include <string>
#include <iostream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/inte.h>
#include <o2scl/funct.h>
#include <o2scl/mroot.h>

// To get directories for calibrate function
#include <o2scl/lib_settings.h>
// To read tables in calibrate function
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

/** \file part.h
    \brief File defining \ref o2scl::thermo and \ref o2scl::part 
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A class containing three thermodynamical variables (energy
      density, pressure, entropy density)
  */
  class thermo {

  public:

    /// pressure
    double pr;
    /// energy density
    double ed;
    /// entropy density
    double en;

    /// Return string denoting type ("thermo")
    const char *type() { return "thermo"; }

    // Default constructor
    thermo() {
    }

    /// Copy constructor
    thermo(const thermo &t) {
      ed=t.ed;
      pr=t.pr;
      en=t.en;
    }

    /// Copy construction with operator=()
    thermo &operator=(const thermo &t) {
      if (this!=&t) {
	ed=t.ed;
	pr=t.pr;
	en=t.en;
      }
      return *this;
    }

  };

  /** \brief Addition operator
   */
  extern thermo operator+(const thermo &left, const thermo &right);

  /** \brief Subtraction operator
   */
  extern thermo operator-(const thermo &left, const thermo &right);

  /** \brief Particle base class 
   */
  class part {
    
  public:

    /// Degeneracy (e.g. spin and color if applicable)
    double g;
    /// Mass
    double m;
    /// Number density
    double n;
    /// Energy density
    double ed;
    /// Pressure
    double pr;
    /// Chemical potential
    double mu;
    /// Entropy density
    double en;
    /// Effective mass (Dirac unless otherwise specified)
    double ms;
    /// Effective chemical potential
    double nu;
    /** \brief If true, include the mass in the energy 
	density and chemical potential (default true) 
    */
    bool inc_rest_mass;
    /// True if the particle is non-interacting (default true)
    bool non_interacting;

    /// Copy constructor
    part(const part &p) {
      g=p.g;
      m=p.m;
      ms=p.ms;
      n=p.n;
      ed=p.ed;
      pr=p.pr;
      mu=p.mu;
      en=p.en;
      nu=p.nu;
      inc_rest_mass=p.inc_rest_mass;
      non_interacting=p.non_interacting;
    }

    /// Copy construction with operator=()
    part &operator=(const part &p) {
      if (this!=&p) {
	g=p.g;
	m=p.m;
	ms=p.ms;
	n=p.n;
	ed=p.ed;
	pr=p.pr;
	mu=p.mu;
	en=p.en;
	nu=p.nu;
	inc_rest_mass=p.inc_rest_mass;
	non_interacting=p.non_interacting;
      }
      return *this;
    }
    
    /// Make a particle of mass \c mass and degeneracy \c dof.
    part(double mass=0.0, double dof=0.0);

    virtual ~part();
  
    /// Set the mass \c mass and degeneracy \c dof.
    virtual void init(double mass, double dof);

    /** \brief Make \c ap an anti-particle with the same mass
	and degeneracy

	This sets the \ref m, \ref g, \ref ms, \ref inc_rest_mass
	and \ref non_interacting fields of \c ap equal to that
	of the current object. If \ref inc_rest_mass is true,
	then it sets 
	\f[
	\mu_{\mathrm{anti}} = - \mu
	\qquad\mathrm{and}\qquad
	\nu_{\mathrm{anti}} = - \nu
	\f]
	and if \ref inc_rest_mass is false, it sets
	\f[
	\mu_{\mathrm{anti}} = - \mu - 2 m
	\qquad\mathrm{and}\qquad
	\nu_{\mathrm{anti}} = - \nu - 2 m
	\f]
    */
    virtual void anti(part &ap);

    /// Return string denoting type ("part")
    virtual const char *type() { return "part"; }
    
  };

  /** \brief Addition operator
   */
  extern thermo operator+(const thermo &left, const part &right);

  /** \brief Subtraction operator
   */
  extern thermo operator-(const thermo &left, const part &right);

  /** \brief Calibrate a particle thermodynamics class with results
      stored in a table

      This compares the approximation to the exact results using
      calc_density(), calc_mu(), pair_density() and pair_mu(). It
      tries each function twelve times. It tries three different
      temperatures, setting both <tt>inc_rest_mass</tt> and
      <tt>non_interacting</tt> equal to <tt>true</tt> and
      <tt>false</tt>.
      
      The <tt>verbose</tt> parameter controls the amount of output.
      
      \future Also calibrate massless fermions?
  */
  template<class part_t, class thermo_t>
    double part_calibrate(part_t &p, thermo_t &th, bool test_pair,
			  std::string file, bool nr_mode=false,
			  int verbose=0, bool external=false) {
			  
    double ret=0;
  
    // ----------------------------------------------------------------
    // Will return to these original values afterwards

    part_t orig=p;

    // ----------------------------------------------------------------
    // Read data file

    std::string fname;
    if (external==false) {
      fname=o2scl_settings.get_data_dir()+file;
    } else {
      fname=file;
    }

    if (verbose>1) {
      std::cout << "In part_calibrate(), loading file named\n\t'" 
		<< fname << "'.\n" << std::endl;
    }
    o2scl::table<> tab;
    o2scl_hdf::hdf_file hf;
    hf.open(fname);
    std::string name;
#ifndef O2SCL_NO_HDF_INPUT  
    hdf_input(hf,tab,name);
#endif
    hf.close();
  
    if (tab.get_nlines()==0) {
      std::string str="Failed to load data from file '"+fname+
	"' in part_calibrate(). Bad filename?";
      O2SCL_ERR(str.c_str(),exc_efilenotfound);
    }

    if (!tab.is_column("ed")) {
      tab.function_column("ed_mot*mot","ed");
      if (test_pair) {
	tab.function_column("pair_ed_mot*mot","pair_ed");
      }
    }
    
    // ----------------------------------------------------------------

    p.g=2.0;
  
    size_t cnt=0;
    part bad, dev, exact;
    double m_bad=0.0, mu_bad=0.0, T_bad=0.0, mot_bad=0.0, psi_bad=0.0;
    p.non_interacting=true;
  
    // ----------------------------------------------------------------
    // First pass, test calc_mu() 

    // k=0,2 are with rest mass, k=1,3 are without
    // k=0,1 are non-interacting, k=2,3 are interacting
    for(size_t k=0;k<4;k++) {

      double ret_local=0.0;
      
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
	  
	    p.inc_rest_mass=true;

	    if (k>=2) {
	      p.non_interacting=false;
	      p.ms=mot*T;
	      p.m=p.ms*1.5;
	      p.nu=p.ms+T*psi;
	      p.mu=0.0;
	    } else {
	      p.non_interacting=true;
	      p.ms=0.0;
	      p.m=mot*T;
	      p.nu=0.0;
	      p.mu=p.m+T*psi;
	    }
	  
	  } else {
	  
	    p.inc_rest_mass=false;
	  
	    if (k>=2) {
	      p.non_interacting=false;
	      p.ms=mot*T;
	      p.m=p.ms*1.5;
	      p.nu=T*psi-p.m+p.ms;
	      p.mu=0.0;
	    } else {
	      p.non_interacting=true;
	      p.ms=0.0;
	      p.m=mot*T;
	      p.nu=0.0;
	      p.mu=T*psi;
	    }
	  
	  }
	
	  th.calc_mu(p,T);
	
	  exact.n*=pow(T,3.0);
	  if (k%2==0) {
	    exact.ed*=pow(T,4.0);
	  } else {
	    exact.ed=exact.ed*pow(T,4.0)-exact.n*p.m;
	  }
	  exact.pr*=pow(T,4.0);
	  exact.en*=pow(T,3.0);
	
	  dev.n+=fabs((p.n-exact.n)/exact.n);
	  dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	  dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	  dev.en+=fabs((p.en-exact.en)/exact.en);
	  if (fabs((p.pr-exact.pr)/exact.pr)>1.0e-2) {
	    std::cout << "ni,icm: " << p.non_interacting << " "
		      << p.inc_rest_mass << std::endl;
	    std::cout << "psi,mot,T: "
		      << psi << " " << mot << " " << T << std::endl;
	    std::cout << "m,ms,mu,nu: " << p.m << " " << p.ms << " "
		      << p.mu << " " << p.nu << std::endl;
	    std::cout << "approx n,ed,pr,en: " << p.n << " " << p.ed << " "
		      << p.pr << " " << p.en << std::endl;
	    std::cout << "exact  n,ed,pr,en: "
		      << exact.n << " " << exact.ed << " "
		      << exact.pr << " " << exact.en
		      << std::endl;
	    exit(-1);
	  }
	
	  cnt++;
	  if (fabs((p.n-exact.n)/exact.n)>bad.n) {
	    bad.n=fabs((p.n-exact.n)/exact.n);
	    if (bad.n>ret_local) {
	      if (k>=2) {
		mu_bad=p.nu;
		m_bad=p.ms;
	      } else {
		mu_bad=p.mu;
		m_bad=p.m;
	      }
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.n;
	    }
	  }
	  if (fabs((p.ed-exact.ed)/exact.ed)>bad.ed) {
	    bad.ed=fabs((p.ed-exact.ed)/exact.ed);
	    if (bad.ed>ret_local) {
	      if (k>=2) {
		mu_bad=p.nu;
		m_bad=p.ms;
	      } else {
		mu_bad=p.mu;
		m_bad=p.m;
	      }
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.ed;
	    }
	  }
	  if (fabs((p.pr-exact.pr)/exact.pr)>bad.pr) {
	    bad.pr=fabs((p.pr-exact.pr)/exact.pr);
	    if (bad.pr>ret_local) {
	      if (k>=2) {
		mu_bad=p.nu;
		m_bad=p.ms;
	      } else {
		mu_bad=p.mu;
		m_bad=p.m;
	      }
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.pr;
	    }
	  }
	  if (fabs((p.en-exact.en)/exact.en)>bad.en) {
	    bad.en=fabs((p.en-exact.en)/exact.en);
	    if (bad.en>ret_local) {
	      if (k>=2) {
		mu_bad=p.nu;
		m_bad=p.ms;
	      } else {
		mu_bad=p.mu;
		m_bad=p.m;
	      }
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.en;
	    }
	  }

	  if (verbose>1) {
	    std::cout.precision(5);
	    if (k>=2) {
	      std::cout << "T,ms,nu,psi,mot: " << T << " "
			<< p.ms << " " << p.nu << " "
			<< psi << " " << mot << std::endl;
	    } else {
	      std::cout << "T,m,mu,psi,mot: " << T << " "
			<< p.m << " " << p.mu << " " 
			<< psi << " " << mot << std::endl;
	    }
	    std::cout.precision(5);
	    std::cout << "n,ed,pr,en: " << std::endl;
	    std::cout << "approx: " << p.n << " " << p.ed << " "
		      << p.pr << " " << p.en << std::endl;
	    std::cout << "exact : " << exact.n << " " << exact.ed << " " 
		      << exact.pr << " " << exact.en << std::endl;
	    std::cout << "bad   : " << bad.n << " " << bad.ed << " " 
		      << bad.pr << " " << bad.en << std::endl;
	    std::cout << "ret_local,ret: " << ret_local << " "
		      << ret << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
	      char ch;
	      std::cin >> ch;
	    }
	  }

	  if (ret_local>ret) {
	    ret=ret_local;
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
	  std::cout << "Function calc_mu(), include rest mass"
		    << std::endl;
	} else if (k==1) {
	  std::cout << "Function calc_mu(), without rest mass"
		    << std::endl;
	} else if (k==2) {
	  std::cout << "Function calc_mu(), include rest mass, "
		    << "interacting" << std::endl;
	} else {
	  std::cout << "Function calc_mu(), without rest mass, "
		    << "interacting" << std::endl;
	}

	std::cout << "Average performance: " << std::endl;
	std::cout << "n: " << dev.n << " ed: " << dev.ed << " pr: " 
		  << dev.pr << " en: " << dev.en << std::endl;
	std::cout << "Worst case: " << std::endl;
	std::cout << "n: " << bad.n << " ed: " << bad.ed << " pr: " 
		  << bad.pr << " en: " << bad.en << std::endl;
	std::cout << "mu: " << mu_bad << " m: " << m_bad 
		  << " T: " << T_bad << " mot: " << mot_bad
		  << "\n\tpsi: " << psi_bad << std::endl;
	std::cout << "ret_local,ret: " << ret_local << " "
		  << ret << std::endl;
	std::cout << std::endl;
	if (verbose>2) {
	  char ch;
	  std::cin >> ch;
	}
      }

      // Reset p.non_interacting
      p.non_interacting=true;

      // End of k loop
    }
    
    // ----------------------------------------------------------------
    // Second pass, test calc_density()

    // k=0,2 are with rest mass, k=1,3 are without
    // k=0,1 are non-interacting, k=2,3 are interacting
    for(size_t k=0;k<4;k++) {

      double ret_local=0.0;
      
      // Initialize storage
      dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
      bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
      // Temperature loop
      for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	// Loop over each point in the data file
	for(size_t i=0;i<tab.get_nlines();i++) {
	
	  double mot=tab.get("mot",i);
	  double psi=tab.get("psi",i);
	  p.n=tab.get("n",i)*pow(T,3.0);
	  exact.ed=tab.get("ed",i);
	  exact.pr=tab.get("pr",i);
	  exact.en=tab.get("en",i);

	  if (k>=2) {
	    p.non_interacting=false;
	    p.ms=mot*T;
	    p.m=p.ms*1.5;
	  } else {
	    p.non_interacting=true;
	    p.m=mot*T;
	    p.ms=0.0;
	  }
	  if (k%2==0) {
	    p.inc_rest_mass=true;
	    if (k>=2) {
	      if (nr_mode) {
		p.nu=T*psi+p.m;
	      } else {
		p.nu=T*psi+p.ms;
	      }
	      p.mu=0.0;
	    } else {
	      p.nu=0.0;
	      p.mu=T*psi+p.m;
	    }
	  } else {
	    p.inc_rest_mass=false;
	    if (k>=2) {
	      if (nr_mode) {
		p.nu=T*psi;
	      } else {
		p.nu=T*psi-p.m+p.ms;
	      }
	      p.mu=0.0;
	    } else {
	      p.nu=0.0;
	      p.mu=T*psi;
	    }
	  }
	  
	  exact.n=p.n;
	  if (k%2==0) {
	    exact.ed*=pow(T,4.0);
	  } else {
	    exact.ed=exact.ed*pow(T,4.0)-exact.n*p.m;
	  }
	  exact.pr*=pow(T,4.0);
	  exact.en*=pow(T,3.0);

	  // Give it a guess for the chemical potential
	  if (k>=2) {
	    p.nu=p.m;
	  } else {
	    p.mu=p.m;
	  }

	  th.calc_density(p,T);
	
	  if (k>=2) {
	    dev.nu+=fabs((p.nu-exact.nu)/exact.nu);
	  } else {
	    dev.mu+=fabs((p.mu-exact.mu)/exact.mu);
	  }
	  dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	  dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	  dev.en+=fabs((p.en-exact.en)/exact.en);

	  cnt++;
	  if (k>=2) {
	    if (fabs((p.nu-exact.nu)/exact.nu)>bad.mu) {
	      bad.mu=fabs((p.nu-exact.nu)/exact.nu);
	      if (bad.mu>ret_local) {
		mu_bad=p.nu;
		m_bad=p.ms;
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.n;
	      }
	    }
	  } else {
	    if (fabs((p.mu-exact.mu)/exact.mu)>bad.mu) {
	      bad.mu=fabs((p.mu-exact.mu)/exact.mu);
	      if (bad.mu>ret_local) {
		mu_bad=p.mu;
		m_bad=p.m;
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.n;
	      }
	    }
	  }
	  if (fabs((p.ed-exact.ed)/exact.ed)>bad.ed) {
	    bad.ed=fabs((p.ed-exact.ed)/exact.ed);
	    if (bad.ed>ret_local) {
	      mu_bad=p.mu;
	      m_bad=p.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.ed;
	    }
	  }
	  if (fabs((p.pr-exact.pr)/exact.pr)>bad.pr) {
	    bad.pr=fabs((p.pr-exact.pr)/exact.pr);
	    if (bad.pr>ret_local) {
	      mu_bad=p.mu;
	      m_bad=p.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.pr;
	    }
	  }
	  if (fabs((p.en-exact.en)/exact.en)>bad.en) {
	    bad.en=fabs((p.en-exact.en)/exact.en);
	    if (bad.en>ret_local) {
	      mu_bad=p.mu;
	      m_bad=p.m;
	      T_bad=T;
	      mot_bad=mot;
	      psi_bad=psi;
	      ret_local=bad.en;
	    }
	  }

	  if (verbose>1) {
	    std::cout.precision(5);
	    if (k>=2) {
	      std::cout << "T,ms,n,psi,mot: " << T << " "
			<< p.ms << " " << p.n << " "
			<< psi << " " << mot << std::endl;
	    } else {
	      std::cout << "T,m,n,psi,mot: " << T << " "
			<< p.m << " " << p.n << " "
			<< psi << " " << mot << std::endl;
	    }
	    std::cout.precision(6);
	    if (k>=2) {
	      std::cout << "nu,ed,pr,en: " << std::endl;
	      std::cout << "approx: " << p.nu << " " << p.ed << " "
			<< p.pr << " " << p.en << std::endl;
	      std::cout << "exact : " << exact.nu << " " << exact.ed << " " 
			<< exact.pr << " " << exact.en << std::endl;
	    } else {
	      std::cout << "mu,ed,pr,en: " << std::endl;
	      std::cout << "approx: " << p.mu << " " << p.ed << " "
			<< p.pr << " " << p.en << std::endl;
	      std::cout << "exact : " << exact.mu << " " << exact.ed << " " 
			<< exact.pr << " " << exact.en << std::endl;
	    }
	    std::cout << "bad   : " << bad.mu << " " << bad.ed << " " 
		      << bad.pr << " " << bad.en << std::endl;
	    std::cout << "ret_local,ret: " << ret_local << " "
		      << ret << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
	      char ch;
	      std::cin >> ch;
	    }
	  }

	  if (ret_local>ret) {
	    ret=ret_local;
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
	  std::cout << "Function calc_density(), include rest mass"
		    << std::endl;
	} else if (k==1) {
	  std::cout << "Function calc_density(), without rest mass"
		    << std::endl;
	} else if (k==2) {
	  std::cout << "Function calc_density(), include rest mass, "
		    << "interacting" << std::endl;
	} else {
	  std::cout << "Function calc_density(), without rest mass, "
		    << "interacting" << std::endl;
	}
	
	std::cout << "Average performance: " << std::endl;
	std::cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
		  << dev.pr << " en: " << dev.en << std::endl;
	std::cout << "Worst case: " << std::endl;
	std::cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
		  << bad.pr << " en: " << bad.en << std::endl;
	std::cout << "mu: " << mu_bad << " m: " << m_bad
		  << " T: " << T_bad << " mot: " << mot_bad
		  << "\n\tpsi: " << psi_bad << std::endl;
	std::cout << "ret_local,ret: " << ret_local << " "
		  << ret << std::endl;
	std::cout << std::endl;
	if (verbose>2) {
	  char ch;
	  std::cin >> ch;
	}
      }

      // End of k loop
    }

    if (test_pair) {
  
      // ----------------------------------------------------------------
      // Third pass, test pair_mu() 

      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	double ret_local=0.0;

	// Initialize storage
	dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	// Temperature loop
	for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    double mot=tab.get("mot",i);
	    double psi=tab.get("psi",i);
	    exact.n=tab.get("pair_n",i);
	    exact.ed=tab.get("pair_ed",i);
	    exact.pr=tab.get("pair_pr",i);
	    exact.en=tab.get("pair_en",i);
      
	    if (k>=2) {
	      p.non_interacting=false;
	      p.ms=mot*T;
	      p.m=p.ms*1.5;
	    } else {
	      p.non_interacting=true;
	      p.m=mot*T;
	      p.ms=0.0;
	    }
	    if (k%2==0) {
	      p.inc_rest_mass=true;
	      if (k>=2) {
		if (nr_mode) {
		  p.nu=T*psi+p.m;
		} else {
		  p.nu=T*psi+p.ms;
		}
		p.mu=0.0;
	      } else {
		p.nu=0.0;
		p.mu=T*psi+p.m;
	      }
	    } else {
	      p.inc_rest_mass=false;
	      if (k>=2) {
		if (nr_mode) {
		  p.nu=T*psi;
		} else {
		  p.nu=T*psi-p.m+p.ms;
		}
		p.mu=0.0;
	      } else {
		p.nu=0.0;
		p.mu=T*psi;
	      }
	    }
	
	    th.pair_mu(p,T);
	
	    exact.n*=pow(T,3.0);
	    if (k%2==0) {
	      exact.ed*=pow(T,4.0);
	    } else {
	      exact.ed=exact.ed*pow(T,4.0)-exact.n*p.m;
	    }
	    exact.pr*=pow(T,4.0);
	    exact.en*=pow(T,3.0);
	
	    dev.n+=fabs((p.n-exact.n)/exact.n);
	    dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	    dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	    dev.en+=fabs((p.en-exact.en)/exact.en);

	    cnt++;
	    if (fabs((p.n-exact.n)/exact.n)>bad.n) {
	      bad.n=fabs((p.n-exact.n)/exact.n);
	      if (bad.n>ret_local) {
		if (k>=2) {
		  mu_bad=p.mu;
		  m_bad=p.ms;
		} else {
		  mu_bad=p.mu;
		  m_bad=p.m;
		}
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.n;
	      }
	    }
	    if (fabs((p.ed-exact.ed)/exact.ed)>bad.ed) {
	      bad.ed=fabs((p.ed-exact.ed)/exact.ed);
	      if (bad.ed>ret_local) {
		if (k>=2) {
		  mu_bad=p.mu;
		  m_bad=p.ms;
		} else {
		  mu_bad=p.mu;
		  m_bad=p.m;
		}
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.ed;
	      }
	    }
	    if (fabs((p.pr-exact.pr)/exact.pr)>bad.pr) {
	      bad.pr=fabs((p.pr-exact.pr)/exact.pr);
	      if (bad.pr>ret_local) {
		if (k>=2) {
		  mu_bad=p.mu;
		  m_bad=p.ms;
		} else {
		  mu_bad=p.mu;
		  m_bad=p.m;
		}
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.pr;
	      }
	    }
	    if (fabs((p.en-exact.en)/exact.en)>bad.en) {
	      bad.en=fabs((p.en-exact.en)/exact.en);
	      if (bad.en>ret_local) {
		if (k>=2) {
		  mu_bad=p.mu;
		  m_bad=p.ms;
		} else {
		  mu_bad=p.mu;
		  m_bad=p.m;
		}
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.en;
	      }
	    }

	    if (verbose>1) {
	      std::cout.precision(5);
	      std::cout << "T,m,mu,psi,mot: " << T << " " << p.m << " "
			<< p.mu << " " << psi << " " << mot << std::endl;
	      std::cout.precision(6);
	      std::cout << i << " " << exact.m << " " << exact.ms << " "
			<< exact.mu << " " << exact.nu << std::endl;
	      std::cout << "n,ed,pr,en: " << std::endl;
	      std::cout << "approx: " << p.n << " " << p.ed << " "
			<< p.pr << " " << p.en << std::endl;
	      std::cout << "exact : " << exact.n << " " << exact.ed << " " 
			<< exact.pr << " " << exact.en << std::endl;
	      std::cout << "bad   : " << bad.n << " " << bad.ed << " " 
			<< bad.pr << " " << bad.en << std::endl;
	      std::cout << "ret_local,ret: " << ret_local << " "
			<< ret << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
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
	    std::cout << "Function pair_mu(), include rest mass"
		      << std::endl;
	  } else if (k==1) {
	    std::cout << "Function pair_mu(), without rest mass"
		      << std::endl;
	  } else if (k==2) {
	    std::cout << "Function pair_mu(), include rest mass, "
		      << "interacting" << std::endl;
	  } else {
	    std::cout << "Function pair_mu(), without rest mass, "
		      << "interacting" << std::endl;
	  }

	  std::cout << "Average performance: " << std::endl;
	  std::cout << "n: " << dev.n << " ed: " << dev.ed << " pr: " 
		    << dev.pr << " en: " << dev.en << std::endl;
	  std::cout << "Worst case: " << std::endl;
	  std::cout << "n: " << bad.n << " ed: " << bad.ed << " pr: " 
		    << bad.pr << " en: " << bad.en << std::endl;
	  std::cout << "mu: " << mu_bad << " m: " << m_bad
		    << " T: " << T_bad << " mot: " << mot_bad
		    << "\n\tpsi: " << psi_bad << std::endl;
	  std::cout << "ret_local,ret: " << ret_local << " "
		    << ret << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
	    char ch;
	    std::cin >> ch;
	  }
	}

	// End of k loop
      }

      // ----------------------------------------------------------------
      // Fourth pass, test pair_density()

      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	double ret_local=0.0;
	
	// Initialize storage
	dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	// Temperature loop
	for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    double mot=tab.get("mot",i);
	    double psi=tab.get("psi",i);
	    p.n=tab.get("pair_n",i)*pow(T,3.0);	
	    exact.ed=tab.get("pair_ed",i);
	    exact.pr=tab.get("pair_pr",i);
	    exact.en=tab.get("pair_en",i);
	  
	    if (k>=2) {
	      p.non_interacting=false;
	      p.ms=mot*T;
	      p.m=p.ms*1.5;
	    } else {
	      p.non_interacting=true;
	      p.m=mot*T;
	      p.ms=0.0;
	    }
	    if (k%2==0) {
	      p.inc_rest_mass=true;
	      if (k>=2) {
		if (nr_mode) {
		  exact.nu=T*psi+p.m;
		} else {
		  exact.nu=T*psi+p.ms;
		}
		exact.mu=0.0;
	      } else {
		exact.nu=0.0;
		exact.mu=T*psi+p.m;
	      }
	    } else {
	      p.inc_rest_mass=false;
	      if (k>=2) {
		if (nr_mode) {
		  exact.nu=T*psi;
		} else {
		  exact.nu=T*psi-p.m+p.ms;
		}
		exact.mu=0.0;
	      } else {
		exact.nu=0.0;
		exact.mu=T*psi;
	      }
	    }

	    exact.n=p.n;
	    if (k%2==0) {
	      exact.ed*=pow(T,4.0);
	    } else {
	      exact.ed=exact.ed*pow(T,4.0)-exact.n*p.m;
	    }
	    exact.pr*=pow(T,4.0);
	    exact.en*=pow(T,3.0);

	    // Give it a guess for the chemical potential
	    p.mu=p.m;

	    th.pair_density(p,T);

	    if (k>=2) {
	      dev.nu+=fabs((p.nu-exact.nu)/exact.nu);
	    } else {
	      dev.mu+=fabs((p.mu-exact.mu)/exact.mu);
	    }
	    dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	    dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	    dev.en+=fabs((p.en-exact.en)/exact.en);
	
	    cnt++;
	    if (k>=2) {
	      if (fabs((p.nu-exact.nu)/exact.nu)>bad.mu) {
		bad.mu=fabs((p.nu-exact.nu)/exact.nu);
		if (bad.mu>ret_local) {
		  mu_bad=p.nu;
		  m_bad=p.ms;
		  T_bad=T;
		  mot_bad=mot;
		  psi_bad=psi;
		  ret_local=bad.n;
		}
	      }
	    } else {
	      if (fabs((p.mu-exact.mu)/exact.mu)>bad.mu) {
		bad.mu=fabs((p.mu-exact.mu)/exact.mu);
		if (bad.mu>ret_local) {
		  mu_bad=p.mu;
		  m_bad=p.m;
		  T_bad=T;
		  mot_bad=mot;
		  psi_bad=psi;
		  ret_local=bad.n;
		}
	      }
	    }
	    if (fabs((p.ed-exact.ed)/exact.ed)>bad.ed) {
	      bad.ed=fabs((p.ed-exact.ed)/exact.ed);
	      if (bad.ed>ret_local) {
		mu_bad=p.mu;
		m_bad=p.m;
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.ed;
	      }
	    }
	    if (fabs((p.pr-exact.pr)/exact.pr)>bad.pr) {
	      bad.pr=fabs((p.pr-exact.pr)/exact.pr);
	      if (bad.pr>ret_local) {
		mu_bad=p.mu;
		m_bad=p.m;
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.pr;
	      }
	    }
	    if (fabs((p.en-exact.en)/exact.en)>bad.en) {
	      bad.en=fabs((p.en-exact.en)/exact.en);
	      if (bad.en>ret_local) {
		mu_bad=p.mu;
		m_bad=p.m;
		T_bad=T;
		mot_bad=mot;
		psi_bad=psi;
		ret_local=bad.en;
	      }
	    }
	    
	    if (verbose>1) {
	      std::cout.precision(5);
	      if (k>=2) {
		std::cout << "T,ms,n,psi,mot: " << T << " " << p.ms << " " 
			  << p.n << " " << psi << " " << mot << std::endl;
	      } else {
		std::cout << "T,m,n,psi,mot: " << T << " " << p.m << " " 
			  << p.n << " " << psi << " " << mot << std::endl;
	      }
	      std::cout.precision(6);
	      if (k>=2) {
		std::cout << "nu,ed,pr,en: " << std::endl;
		std::cout << "approx: " << p.nu << " " << p.ed << " "
			  << p.pr << " " << p.en << std::endl;
		std::cout << "exact : " << exact.nu << " " << exact.ed << " " 
			  << exact.pr << " " << exact.en << std::endl;
	      } else {
		std::cout << "mu,ed,pr,en: " << std::endl;
		std::cout << "approx: " << p.mu << " " << p.ed << " "
			  << p.pr << " " << p.en << std::endl;
		std::cout << "exact : " << exact.mu << " " << exact.ed << " " 
			  << exact.pr << " " << exact.en << std::endl;
	      }
	      std::cout << "bad   : " << bad.mu << " " << bad.ed << " " 
			<< bad.pr << " " << bad.en << std::endl;
	      std::cout << "ret_local,ret: " << ret_local << " "
			<< ret << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
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
	    std::cout << "Function pair_density(), include rest mass"
		      << std::endl;
	  } else if (k==1) {
	    std::cout << "Function pair_density(), without rest mass"
		      << std::endl;
	  } else if (k==2) {
	    std::cout << "Function pair_density(), include rest mass, "
		      << "interacting" << std::endl;
	  } else {
	    std::cout << "Function pair_density(), without rest mass, "
		      << "interacting" << std::endl;
	  }

	  std::cout << "Average performance: " << std::endl;
	  std::cout << "mu: " << dev.mu << " ed: " << dev.ed << " pr: " 
		    << dev.pr << " en: " << dev.en << std::endl;
	  std::cout << "Worst case: " << std::endl;
	  std::cout << "mu: " << bad.mu << " ed: " << bad.ed << " pr: " 
		    << bad.pr << " en: " << bad.en << std::endl;
	  std::cout << "mu: " << mu_bad << " m: " << m_bad
		    << " T: " << T_bad << " mot: " << mot_bad
		    << "\n\tpsi: " << psi_bad << std::endl;
	  std::cout << "ret_local,ret: " << ret_local << " "
		    << ret << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
	    char ch;
	    std::cin >> ch;
	  }
	}

	if (ret_local>ret) {
	  ret=ret_local;
	}
    
	// End of k loop
      }

      // End of 'if (test_pair)'
    }

    // ----------------------------------------------------------------
    // Return to the original values 

    p=orig;
  
    return ret;
  }
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
