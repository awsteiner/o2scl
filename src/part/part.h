/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
    \brief File defining \ref o2scl::thermo_tl and \ref o2scl::part_tl 
*/

namespace o2scl {

  /** \brief A class containing three thermodynamical variables (energy
      density, pressure, entropy density)
  */
  template<class fp_t=double> class thermo_tl {

  public:

    /// Pressure
    fp_t pr;
    /// Energy density
    fp_t ed;
    /// Entropy density
    fp_t en;

    /// Return string denoting type ("thermo")
    const char *type() { return "thermo"; }

    // Default constructor
    thermo_tl() {
    }

    /// Copy constructor
    thermo_tl(const thermo_tl &t) {
      ed=t.ed;
      pr=t.pr;
      en=t.en;
    }

    /// Copy construction with operator=()
    thermo_tl &operator=(const thermo_tl &t) {
      if (this!=&t) {
        ed=t.ed;
        pr=t.pr;
        en=t.en;
      }
      return *this;
    }

  };

  /** \brief Double-precision thermodynamics object
   */
  typedef thermo_tl<double> thermo;

  /** \brief Addition operator
   */
  extern thermo operator+(const thermo &left, const thermo &right);

  /** \brief Subtraction operator
   */
  extern thermo operator-(const thermo &left, const thermo &right);

  /** \brief Particle base class 
   */
  template<class fp_t=double> class part_tl {
    
  public:

    /// Degeneracy (e.g. spin and color if applicable)
    fp_t g;
    /// Mass
    fp_t m;
    /// Number density
    fp_t n;
    /// Energy density
    fp_t ed;
    /// Pressure
    fp_t pr;
    /// Chemical potential
    fp_t mu;
    /// Entropy density
    fp_t en;
    /// Effective mass (Dirac unless otherwise specified)
    fp_t ms;
    /// Effective chemical potential
    fp_t nu;
    /** \brief If true, include the mass in the energy 
        density and chemical potential (default true) 
    */
    bool inc_rest_mass;
    /// True if the particle is non-interacting (default true)
    bool non_interacting;

    /// Copy constructor
    part_tl(const part_tl &p) {
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
    part_tl &operator=(const part_tl &p) {
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
    part_tl(fp_t mass=0.0, fp_t dof=0.0) {
      m=mass; 
      ms=mass; 
      g=dof;
    
      non_interacting=true;
      inc_rest_mass=true;
    }    
  
    virtual ~part_tl() {
    }
  
    /// Set the mass \c mass and degeneracy \c dof.
    virtual void init(fp_t mass, fp_t dof) {
      m=mass; 
      ms=mass; 
      g=dof;
      return;
    }

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
    virtual void anti(part_tl &ax) {
      ax.g=g;
      ax.m=m;
      ax.ms=ms;
      ax.inc_rest_mass=inc_rest_mass;
      ax.non_interacting=non_interacting;
      if (inc_rest_mass) {
        ax.nu=-nu;
        ax.mu=-mu;
      } else {
        ax.nu=-nu-2.0*m;
        ax.mu=-mu-2.0*m;
      }
      return;
    }    

    /// Return string denoting type ("part_tl")
    virtual const char *type() { return "part_tl"; }
    
  };

  /** \brief Double-precision thermodynamics object
   */
  typedef part_tl<double> part;

  /** \brief Addition operator
   */
  extern thermo operator+(const thermo &left, const part &right);

  /** \brief Subtraction operator
   */
  extern thermo operator-(const thermo &left, const part &right);

  /** \brief Object to organize calibration of particle classes
      to results stored in a table
  */
  template <class fp_t=double> class part_calibrate_class_tl {
    
  public:

    /** \brief Set mass and flags from mot, T, and the index k
     */
    template<class part_t>
    void set_mass_flags(part_t &p, fp_t mot, fp_t T, size_t k) {
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
      } else {
	p.inc_rest_mass=false;
      }
      return;
    }
    
    /** \brief Set chemical potential from psi, T, the index k,
	and the flag nr_mode
    */
    template<class part_t>
    void set_chem_pot(part_t &p, fp_t psi, fp_t T, size_t k,
                      bool nr_mode) {
      if (k%2==0) {
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
      return;
    }

    /** \brief Check the density against the exact result 
	and update 
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_density(part1_t &p, part2_t &exact, part3_t &bad, size_t k,
                       fp_t T, fp_t mot, fp_t psi,
                       fp_t &mu_bad, fp_t &m_bad,
                       fp_t &T_bad, fp_t &mot_bad, fp_t &psi_bad,
                       fp_t &ret_local) {
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
      return;
    }      
    
    /** \brief Check the chemical potential against the exact result 
	and update 
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_chem_pot(part1_t &p, part2_t &exact, part3_t &bad, size_t k,
                        fp_t T, fp_t mot, fp_t psi,
                        fp_t &mu_bad, fp_t &m_bad,
                        fp_t &T_bad, fp_t &mot_bad, fp_t &psi_bad,
                        fp_t &ret_local) {
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
      return;
    }
      
    /** \brief Check the energy density, pressure, and entropy against
	the exact result and update
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_eps(part1_t &p, part2_t &exact, part3_t &bad, size_t k,
                   fp_t T, fp_t mot, fp_t psi,
                   fp_t &mu_bad, fp_t &m_bad,
                   fp_t &T_bad, fp_t &mot_bad, fp_t &psi_bad,
                   fp_t &ret_local) {
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
      return;
    }      
    
    /** \brief Check the three derivatives against
	the exact result and update
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_derivs(part1_t &p, part2_t &exact, part3_t &bad, size_t k,
                      fp_t T, fp_t mot, fp_t psi,
                      fp_t &mu_bad, fp_t &m_bad,
                      fp_t &T_bad, fp_t &mot_bad, fp_t &psi_bad,
                      fp_t &ret_local) {
      if (fabs((p.dndT-exact.dndT)/exact.dndT)>bad.dndT) {
	bad.dndT=fabs((p.dndT-exact.dndT)/exact.dndT);
	if (bad.dndT>ret_local) {
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
	  ret_local=bad.dndT;
	}
      }
      if (fabs((p.dndmu-exact.dndmu)/exact.dndmu)>bad.dndmu) {
	bad.dndmu=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	if (bad.dndmu>ret_local) {
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
	  ret_local=bad.dndmu;
	}
      }
      if (fabs((p.dsdT-exact.dsdT)/exact.dsdT)>bad.dsdT) {
	bad.dsdT=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	if (bad.dsdT>ret_local) {
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
	  ret_local=bad.dsdT;
	}
      }
      return;
    }      
    
    /** \brief Calibrate a particle thermodynamics class with results
	stored in a table
	
	This compares the approximation to the exact results using
	calc_density(), calc_mu(), pair_density() and pair_mu(). It
	tries each function twelve times. It tries three different
	temperatures, setting both <tt>inc_rest_mass</tt> and
	<tt>non_interacting</tt> equal to <tt>true</tt> and
	<tt>false</tt>.
      
	The <tt>verbose</tt> parameter controls the amount of output.
      
        \verbatim embed:rst

        .. todo::

           In function pair_calibrate()

           - Future: Also calibrate massless fermions?

        \endverbatim
    */
    template<class part_t, class thermo_t>
    fp_t part_calibrate(part_t &p, thermo_t &th, bool test_pair,
                        std::string file, bool nr_mode=false,
                        int verbose=0, bool external=false) {
      
      fp_t ret=0;
  
      // ----------------------------------------------------------------
      // Will return to these original values afterwards
      
      part_t orig=p;

      // ----------------------------------------------------------------
      // Read data file

      fp_t ti_test=0;
      fp_t ti_local=0;
      fp_t ti_max=0;
      size_t ti_count=0;
      
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
      hdf_input_n(hf,tab,name);
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
      part_t bad, dev, exact;
      fp_t m_bad=0.0, mu_bad=0.0, T_bad=0.0, mot_bad=0.0, psi_bad=0.0;
      p.non_interacting=true;
  
      // ----------------------------------------------------------------
      // First pass, test calc_mu()
      
      if (verbose>1) {
        std::cout << "First pass, testing calc_mu().\n" << std::endl;
      }

      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	fp_t ret_local=0;
        ti_local=0;
      
	// Initialize storage
	dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	// Temperature loop
	for(fp_t T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    fp_t mot=tab.get("mot",i);
	    fp_t psi=tab.get("psi",i);
	    exact.n=tab.get("n",i);
	    exact.ed=tab.get("ed",i);
	    exact.pr=tab.get("pr",i);
	    exact.en=tab.get("en",i);

	    set_mass_flags(p,mot,T,k);
	    set_chem_pot(p,psi,T,k,nr_mode);

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
            }

	    th.calc_mu(p,T);
	
	    exact.n*=pow(T,3.0);
	    if (nr_mode) {
              exact.ed*=pow(T,4.0);
	      if (k%2==0) {
                exact.ed+=exact.n*p.m;
	      }
	    } else {
              exact.ed*=pow(T,4.0);
	      if (k%2!=0) {
		exact.ed-=exact.n*p.m;
	      }
	    }
	    exact.pr*=pow(T,4.0);
	    exact.en*=pow(T,3.0);
	
            if (th.verify_ti) {
              double val;
              if (p.pr==0.0) {
                val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/
                  abs(-p.ed+p.n*p.nu+p.en*T);
              } else {
                val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/abs(p.pr);
              }
              ti_test+=val;
              ti_count++;
              if (ti_test>ti_local) ti_local=ti_test;
              
              if (val>0.1 || !isfinite(ti_test)) {
                std::cout << "     T,m*,m: " << T << " " << p.ms
                          << " " << p.m << std::endl;
                std::cout << "     n1,n2: " << p.n << " "
                          << exact.n << std::endl;
                std::cout << "     ed1,ed2: " << p.ed << " "
                          << exact.ed << std::endl;
                std::cout << "     en1,en2: " << p.en << " "
                          << exact.en << std::endl;
                std::cout << "     pr1,pr2,pr3: " << p.pr << " "
                          << exact.pr << " "
                          << -p.ed+p.n*p.nu+p.en*T << std::endl;
                std::cout << "     ed,en,n,mu: " << p.ed << " "
                          << p.en << " " << p.n << " " << p.mu << std::endl;

                O2SCL_ERR("Thermodynamic identity violated.",
                          o2scl::exc_esanity);
              }
            }
	    dev.n+=fabs((p.n-exact.n)/exact.n);
	    dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	    dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	    dev.en+=fabs((p.en-exact.en)/exact.en);
	
	    cnt++;
	  
	    check_density<part_t,part_t,part_t>
              (p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
               mot_bad,psi_bad,ret_local);
	  
	    check_eps<part_t,part_t,part_t>
              (p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
               mot_bad,psi_bad,ret_local);

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
              if (ti_count>0) {
                std::cout << "check ti: " << ti_test/ti_count << " "
                          << ti_local << std::endl;
              }
	      std::cout << "ret_local,ret: " << ret_local << " "
			<< ret << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
	    }
	    if (ti_local>ti_max) {
	      ti_max=ret_local;
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
          if (ti_count>0) {
            std::cout << "check ti: " << ti_test/ti_count << " "
                      << ti_local << std::endl;
          }
	  std::cout << "ret_local,ret: " << ret_local << " "
		    << ret << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
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

      if (verbose>1) {
        std::cout << "Second pass, testing calc_density().\n" << std::endl;
      }
      
      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	fp_t ret_local=0.0;
      
	// Initialize storage
	dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	// Temperature loop
	for(fp_t T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    fp_t mot=tab.get("mot",i);
	    fp_t psi=tab.get("psi",i);
	    p.n=tab.get("n",i)*pow(T,3.0);
	    exact.ed=tab.get("ed",i);
	    exact.pr=tab.get("pr",i);
	    exact.en=tab.get("en",i);

	    set_mass_flags(p,mot,T,k);
	    set_chem_pot(exact,psi,T,k,nr_mode);
	    
	    exact.n=p.n;
	    if (nr_mode) {
	      if (k%2==0) {
		exact.ed=exact.ed*pow(T,4.0)+exact.n*p.m;
	      } else {
		exact.ed=exact.ed*pow(T,4.0);
	      }
	    } else {
	      if (k%2==0) {
		exact.ed*=pow(T,4.0);
	      } else {
		exact.ed=exact.ed*pow(T,4.0)-exact.n*p.m;
	      }
	    }
	    exact.pr*=pow(T,4.0);
	    exact.en*=pow(T,3.0);

	    // Give it a guess for the chemical potential
	    if (k>=2) {
	      p.nu=p.m;
	    } else {
	      p.mu=p.m;
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
            }
            
	    th.calc_density(p,T);
            if (th.verify_ti) {
              double val;
              if (p.pr==0.0) {
                val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/
                  abs(-p.ed+p.n*p.nu+p.en*T);
              } else {
                val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/abs(p.pr);
              }
              ti_test+=val;
              ti_count++;
              if (ti_test>ti_local) ti_local=ti_test;
              if (val>0.1 || !isfinite(ti_test)) {
                std::cout << "     " << p.ed << " " << p.en << " "
                          << p.n << " " << p.nu << std::endl;
                O2SCL_ERR("Thermodynamic identity violated.",
                          o2scl::exc_esanity);
              }
            }
	
	    if (k>=2) {
	      dev.nu+=fabs((p.nu-exact.nu)/exact.nu);
	    } else {
	      dev.mu+=fabs((p.mu-exact.mu)/exact.mu);
	    }
	    dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	    dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	    dev.en+=fabs((p.en-exact.en)/exact.en);

	    cnt++;

	    check_chem_pot<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,
                                   m_bad,T_bad,mot_bad,psi_bad,ret_local);
				   
	  
	    check_eps<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
			      mot_bad,psi_bad,ret_local);

	    if (verbose>1) {
	      std::cout.precision(6);
	      if (k>=2) {
		std::cout << "nu,ed,pr,en: " << std::endl;
		std::cout << "approx: " << p.nu << " " << p.ed << " "
			  << p.pr << " " << p.en << std::endl;
		std::cout << "exact : " << exact.nu << " " << exact.ed 
			  << " " << exact.pr << " " << exact.en << std::endl;
	      } else {
		std::cout << "mu,ed,pr,en: " << std::endl;
		std::cout << "approx: " << p.mu << " " << p.ed << " "
			  << p.pr << " " << p.en << std::endl;
		std::cout << "exact : " << exact.mu << " " << exact.ed 
			  << " " << exact.pr << " " << exact.en << std::endl;
	      }
	      std::cout << "bad   : " << bad.mu << " " << bad.ed << " " 
			<< bad.pr << " " << bad.en << std::endl;
              if (ti_count>0) {
                std::cout << "check ti: " << ti_test/ti_count << " "
                          << ti_local << std::endl;
              }
	      std::cout << "ret_local,ret: " << ret_local << " "
			<< ret << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
	    }
	    if (ti_local>ti_max) {
	      ti_max=ret_local;
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
          if (ti_count>0) {
            std::cout << "check ti: " << ti_test/ti_count << " "
                      << ti_local << std::endl;
          }
	  std::cout << "ret_local,ret: " << ret_local << " "
		    << ret << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
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

	  fp_t ret_local=0.0;

	  // Initialize storage
	  dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	  bad.n=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	  // Temperature loop
	  for(fp_t T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	    // Loop over each point in the data file
	    for(size_t i=0;i<tab.get_nlines();i++) {
	
	      fp_t mot=tab.get("mot",i);
	      fp_t psi=tab.get("psi",i);
	      exact.n=tab.get("pair_n",i);
	      exact.ed=tab.get("pair_ed",i);
	      exact.pr=tab.get("pair_pr",i);
	      exact.en=tab.get("pair_en",i);
	    
	      set_mass_flags(p,mot,T,k);
	      set_chem_pot(p,psi,T,k,nr_mode);
	
	      if (verbose>1) {
		std::cout.precision(5);
		std::cout << "T,m,mu,psi,mot: " << T << " " << p.m << " "
			  << p.mu << " " << psi << " " << mot << std::endl;
              }
              
	      th.pair_mu(p,T);
              if (th.verify_ti) {
                double val;
                if (p.pr==0.0) {
                  val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/
                    abs(-p.ed+p.n*p.nu+p.en*T);
                } else {
                  val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/abs(p.pr);
                }
                ti_test+=val;
                ti_count++;
              if (ti_test>ti_local) ti_local=ti_test;
                if (val>0.1 || !isfinite(ti_test)) {
                  std::cout << "     " << p.ed << " " << p.en << " "
                            << p.n << " " << p.nu << std::endl;
                  O2SCL_ERR("Thermodynamic identity violated.",
                            o2scl::exc_esanity);
                }
              }
	
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

	      check_density<part_t,part_t,part_t>
                (p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
                 mot_bad,psi_bad,ret_local);
	    
	      check_eps<part_t,part_t,part_t>
                (p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
                 mot_bad,psi_bad,ret_local);
              
	      if (verbose>1) {
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
                if (ti_count>0) {
                  std::cout << "check ti: " << ti_test/ti_count << " "
                            << ti_local << std::endl;
                }
		std::cout << "ret_local,ret: " << ret_local << " "
			  << ret << std::endl;
		std::cout << std::endl;
		if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
		  char ch;
		  std::cin >> ch;
		}
	      }

	      if (ret_local>ret) {
		ret=ret_local;
	      }
	    if (ti_local>ti_max) {
	      ti_max=ret_local;
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
            if (ti_count>0) {
              std::cout << "check ti: " << ti_test/ti_count << " "
                        << ti_local << std::endl;
            }
            std::cout << "ret_local,ret: " << ret_local << " "
		      << ret << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
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

	  fp_t ret_local=0.0;
	
	  // Initialize storage
	  dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	  bad.mu=0.0; bad.ed=0.0; bad.pr=0.0; bad.en=0.0;
    
	  // Temperature loop
	  for(fp_t T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	    // Loop over each point in the data file
	    for(size_t i=0;i<tab.get_nlines();i++) {
	
	      fp_t mot=tab.get("mot",i);
	      fp_t psi=tab.get("psi",i);
	      p.n=tab.get("pair_n",i)*pow(T,3.0);	
	      exact.ed=tab.get("pair_ed",i);
	      exact.pr=tab.get("pair_pr",i);
	      exact.en=tab.get("pair_en",i);
	  
	      set_mass_flags(p,mot,T,k);
	      set_chem_pot(exact,psi,T,k,nr_mode);

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

	      if (verbose>1) {
		std::cout.precision(5);
		if (k>=2) {
		  std::cout << "T,ms,n,psi,mot: " << T << " " << p.ms << " " 
			    << p.n << " " << psi << " " << mot << std::endl;
		} else {
		  std::cout << "T,m,n,psi,mot: " << T << " " << p.m << " " 
			    << p.n << " " << psi << " " << mot << std::endl;
		}
              }
              
	      th.pair_density(p,T);
              if (th.verify_ti) {
                double val;
                if (p.pr==0.0) {
                  val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/
                    abs(-p.ed+p.n*p.nu+p.en*T);
                } else {
                  val=abs(p.pr+p.ed-p.n*p.nu-p.en*T)/abs(p.pr);
                }
                ti_test+=val;
                ti_count++;
              if (ti_test>ti_local) ti_local=ti_test;
                if (val>0.1 || !isfinite(ti_test)) {
                  std::cout << "     " << p.ed << " " << p.en << " "
                            << p.n << " "
                            << p.nu << std::endl;
                  O2SCL_ERR("Thermodynamic identity violated.",
                            o2scl::exc_esanity);
                }
              }
              
	      if (k>=2) {
		dev.nu+=fabs((p.nu-exact.nu)/exact.nu);
	      } else {
		dev.mu+=fabs((p.mu-exact.mu)/exact.mu);
	      }
	      dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
	      dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
	      dev.en+=fabs((p.en-exact.en)/exact.en);
	
	      cnt++;

	      check_chem_pot<part_t>(p,exact,bad,k,T,mot,psi,
                                     mu_bad,m_bad,T_bad,
				     mot_bad,psi_bad,ret_local);

	      check_eps<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,
                                m_bad,T_bad,
				mot_bad,psi_bad,ret_local);
	    
	    
	      if (verbose>1) {
		std::cout.precision(6);
		if (k>=2) {
		  std::cout << "nu,ed,pr,en: " << std::endl;
		  std::cout << "approx: " << p.nu << " " << p.ed << " "
			    << p.pr << " " << p.en << std::endl;
		  std::cout << "exact : " << exact.nu << " "
                            << exact.ed << " " << exact.pr << " "
                            << exact.en << std::endl;
		} else {
		  std::cout << "mu,ed,pr,en: " << std::endl;
		  std::cout << "approx: " << p.mu << " " << p.ed << " "
			    << p.pr << " " << p.en << std::endl;
		  std::cout << "exact : " << exact.mu << " "
                            << exact.ed << " " << exact.pr << " "
                            << exact.en << std::endl;
		}
		std::cout << "bad   : " << bad.mu << " " << bad.ed << " " 
			  << bad.pr << " " << bad.en << std::endl;
                if (ti_count>0) {
                  std::cout << "check ti: " << ti_test/ti_count
                            << " " << ti_local << std::endl;
                }
		std::cout << "ret_local,ret: " << ret_local << " "
			  << ret << std::endl;
		std::cout << std::endl;
		if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
		  char ch;
		  std::cin >> ch;
		}
	      }

	      if (ret_local>ret) {
		ret=ret_local;
	      }
	    if (ti_local>ti_max) {
	      ti_max=ret_local;
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
            if (ti_count>0) {
              std::cout << "check ti: " << ti_test/ti_count << " "
                        << ti_local << std::endl;
            }
	    std::cout << "ret_local,ret: " << ret_local << " "
		      << ret << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
                std::cout << "Waiting for character: " << std::flush;
	      char ch;
	      std::cin >> ch;
	    }
	  }

	  if (ret_local>ret) {
	    ret=ret_local;
	  }
          if (ti_local>ti_max) {
            ti_max=ret_local;
          }
    
	  // End of k loop
	}

	// End of 'if (test_pair)'
      }

      // ----------------------------------------------------------------
      // Return to the original values 

      p=orig;

      if (th.verify_ti) {
        std::cout << "check ti, ti_max: " << ti_test/ti_count << " "
                  << ti_max << std::endl;
      }
  
      return ret;
    }

  };

  typedef part_calibrate_class_tl<double> part_calibrate_class;

}

#endif
