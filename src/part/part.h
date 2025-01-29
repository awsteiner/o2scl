/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_PART_H
#define O2SCL_PART_H

#include <string>
#include <iostream>
#include <cmath>

#ifdef O2SCL_MULTIP
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_SET_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif
#endif

#include <o2scl/constants.h>
#include <o2scl/inte.h>
#include <o2scl/funct.h>
#include <o2scl/mroot.h>
#include <o2scl/test_mgr.h>

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

  /** \brief Long double-precision thermodynamics object
   */
  typedef thermo_tl<long double> thermo_ld;
  
#ifdef O2SCL_MULTIP
  
  /** \brief 25-digit precision thermodynamics object
   */
  typedef thermo_tl<boost::multiprecision::number<
                       boost::multiprecision::cpp_dec_float<25> > >
  thermo_cdf25;
  
#endif
  
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

  /** \brief Calibrate particle classes by comparing double to
      multiprecision
  */
#ifdef O2SCL_MULTIP
  template <class fp1_t=double,
            class fp2_t=long double,
            class fp3_t=cpp_dec_float_25> class part_cal_new
#else
  template <class fp1_t=double,
            class fp2_t=long double,
            class fp3_t=long double> class part_cal_new
#endif
    {
    
  public:
    
    /** \brief Test the ``calc_mu`` function
     */
    template<class part1_t, class part2_t, class part3_t,
             class thermo1_t, class thermo2_t, class thermo3_t>
    void test_calc_mu(part1_t &f, part2_t &fld, part3_t &f25,
                      thermo1_t &fr, thermo2_t &frld,
                      thermo3_t &fr25, o2scl::test_mgr &t,
                      int &count, int first_test,
                      int cmu_n_max, int cmu_en_max,
                      int cmu_ld_n_max, int cmu_ld_en_max,
                      int cmu_ti_max, int cmu_ld_ti_max,
                      int cmu_25_ti_max) {

      std::cout.precision(4);
      
      // Sums of calc_mu() comparisons between fp types
      int cmu_n=0, cmu_en=0, cmu_ld_n=0, cmu_ld_en=0;
      // Sums of calc_mu() accuracy via. thermodynamic identity
      int cmu_ti=0, cmu_ld_ti=0, cmu_25_ti=0;
      
      std::cout << " cnt m          T           mu/n       "
                << "d-ld  ld-25 ti verify" << std::endl;
      
      for(int im=-2;im<=1;im++) {
        
        f.m=im;
        f.m=pow(10,f.m);
        fld.m=im;
        fld.m=pow(10,fld.m);
        f25.m=im;
        f25.m=pow(10,f25.m);
        
        for(int iT=-2;iT<=1;iT++) {
          
          fp1_t T=iT;
          T=pow(10,T);
          fp2_t Tld=iT;
          Tld=pow(10,Tld);
          fp3_t T25=iT;
          T25=pow(10,T25);

          for(int imu=0;imu<8;imu++) {

            if (imu<4) {
              f.mu=imu-2;
              f.mu=pow(10,f.mu);
              fld.mu=imu-2;
              fld.mu=pow(10,fld.mu);
              f25.mu=imu-2;
              f25.mu=pow(10,f25.mu);
            } else {
              f.mu=imu-6;
              f.mu=-pow(10,f.mu);
              fld.mu=imu-6;
              fld.mu=-pow(10,fld.mu);
              f25.mu=imu-6;
              f25.mu=-pow(10,f25.mu);
            }
            
            std::cout.width(4);
            std::cout << count << " ";
            
            if (count>=first_test) {
              
              std::cout << f.m << " " << T << " ";
              std::cout.setf(std::ios::showpos);
              std::cout << f.mu << " ";
              std::cout.unsetf(std::ios::showpos);
              
              int ret=fr.calc_mu(f,T);
              int retld=frld.calc_mu(fld,Tld);
              int ret25=fr25.calc_mu(f25,T25);
              
              if (ret==0) {
                
                // Test with inc_rest_mass=false
                part1_t f2=f;
                f2.inc_rest_mass=false;
                f2.mu-=f2.m;
                fr.calc_mu(f2,T);
                t.test_rel(f.ed,f2.ed+f2.m*f2.n,1.0e-11,"irm false ed");
                t.test_rel(f.en,f2.en,1.0e-14,"irm false en");
                
                // Test with non_interacting=false
                part1_t f3=f;
                f3.non_interacting=false;
                f3.nu=f3.mu;
                f3.ms=f3.m;
                f3.m*=sqrt(2.0);
                fr.calc_mu(f3,T);
                t.test_rel(f.n,f3.n,1.0e-14,"ni false n");
                t.test_rel(f.en,f3.en,1.0e-14,"ni false en");
                
                // Test with both equal to false
                part1_t f4=f3;
                f4.inc_rest_mass=false;
                f4.non_interacting=false;
                f4.nu-=f4.m;
                fr.calc_mu(f4,T);
                t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed");
                t.test_rel(f3.en,f4.en,1.0e-13,"both false en");
                
                if (t.get_success()==false) {
                  std::cout << "1" << std::endl;
                  exit(-1);
                }
              }
              
              int idn=-2, iden=-2, ildn=-2, ilden=-2;
              if (ret==0 && retld==0) {
                idn=count_digits_same(f.n,fld.n);
                iden=count_digits_same(f.en,fld.en);
              }
              cmu_n+=idn;
              cmu_en+=iden;
              if (retld==0 && ret25==0) {
                ildn=count_digits_same(fld.n,f25.n);
                ilden=count_digits_same(fld.en,f25.en);
              }
              cmu_ld_en+=ilden;
              cmu_ld_n+=ildn;
              
              std::cout.width(2);
              std::cout << idn << " ";
              std::cout.width(2);
              std::cout << iden << " ";
              std::cout.width(2);
              std::cout << ildn << " ";
              std::cout.width(2);
              std::cout << ilden << " ";
              
              int x=-2, xld=-2, x25=-2;
              
              if (ret==0) {
                fp1_t pr2=-f.ed+f.n*f.mu+T*f.en;
                x=count_digits_same(f.pr,pr2);
              }
              cmu_ti+=x;
              std::cout.width(2);
              std::cout << x << " ";
              
              if (retld==0) {
                fp2_t pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
                xld=count_digits_same(fld.pr,pr2ld);
              }
              cmu_ld_ti+=xld;
              std::cout.width(2);
              std::cout << xld << " ";
              
              if (ret25==0) {
                fp3_t pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
                x25=count_digits_same(f25.pr,pr225);
              }
              cmu_25_ti+=x25;
              std::cout.width(2);
              std::cout << x25 << " cmu" << std::endl;
              
            } else {
              
              std::cout << std::endl;
              
            }
            
            count++;
            
          }
        }
      }
      
      std::cout << "calc_mu density (double <-> long double): "
           << cmu_n << std::endl;
      std::cout << "calc_mu entropy (double <-> long double): "
           << cmu_en << std::endl;
      std::cout << "calc_mu density (long double <-> cdf_25): "
           << cmu_ld_n << std::endl;
      std::cout << "calc_mu entropy (long double <-> cdf_25): "
           << cmu_ld_en << std::endl;
      std::cout << "calc_mu ti: " << cmu_ti << std::endl;
      std::cout << "calc_mu long double ti: " << cmu_ld_ti << std::endl;
      std::cout << "calc_mu cpp_dec_float_25 ti: " << cmu_25_ti << std::endl;
      std::cout << std::endl;

      t.test_gen(cmu_n>=cmu_n_max,"cmu_n");
      t.test_gen(cmu_en>=cmu_en_max,"cmu_en");
      t.test_gen(cmu_ld_n>=cmu_ld_n_max,"cmu_ld_n");
      t.test_gen(cmu_ld_en>=cmu_ld_en_max,"cmu_ld_en");
      t.test_gen(cmu_ti>=cmu_ti_max,"cmu_ti");
      t.test_gen(cmu_ld_ti>=cmu_ld_ti_max,"cmu_ld_ti");
      t.test_gen(cmu_25_ti>=cmu_25_ti_max,"cmu_25_ti");
      
      std::cout.precision(6);
      
      return;
    }
    
    /** \brief Test the ``pair_mu`` function
     */
    template<class part1_t, class part2_t, class part3_t,
             class thermo1_t, class thermo2_t, class thermo3_t>
    void test_pair_mu(part1_t &f, part2_t &fld, part3_t &f25,
                      thermo1_t &fr, thermo2_t &frld,
                      thermo3_t &fr25, o2scl::test_mgr &t,
                      int &count, int first_test,
                      int pmu_n_max, int pmu_en_max,
                      int pmu_ld_n_max, int pmu_ld_en_max,
                      int pmu_ti_max, int pmu_ld_ti_max,
                      int pmu_25_ti_max) {

      std::cout.precision(4);
      
      // Sums of pair_mu() comparisons between fp types
      int pmu_n=0, pmu_en=0, pmu_ld_n=0, pmu_ld_en=0;
      // Sums of pair_mu() accuracy via. thermodynamic identity
      int pmu_ti=0, pmu_ld_ti=0, pmu_25_ti=0;
      
      std::cout << " cnt m          T           mu/n       "
                << "d-ld  ld-25 ti verify" << std::endl;
      
      for(int im=-2;im<=1;im++) {
        
        f.m=im;
        f.m=pow(10,f.m);
        fld.m=im;
        fld.m=pow(10,fld.m);
        f25.m=im;
        f25.m=pow(10,f25.m);
        
        for(int iT=-2;iT<=1;iT++) {
          
          fp1_t T=iT;
          T=pow(10,T);
          fp2_t Tld=iT;
          Tld=pow(10,Tld);
          fp3_t T25=iT;
          T25=pow(10,T25);

          for(int imu=0;imu<8;imu++) {

            if (imu<4) {
              f.mu=imu-2;
              f.mu=pow(10,f.mu);
              fld.mu=imu-2;
              fld.mu=pow(10,fld.mu);
              f25.mu=imu-2;
              f25.mu=pow(10,f25.mu);
            } else {
              f.mu=imu-6;
              f.mu=-pow(10,f.mu);
              fld.mu=imu-6;
              fld.mu=-pow(10,fld.mu);
              f25.mu=imu-6;
              f25.mu=-pow(10,f25.mu);
            }
            
            std::cout.width(4);
            std::cout << count << " ";
            
            if (count>=first_test) {
              
              std::cout << f.m << " " << T << " ";
              std::cout.setf(std::ios::showpos);
              std::cout << f.mu << " ";
              std::cout.unsetf(std::ios::showpos);
              
              int ret=fr.pair_mu(f,T);
              int retld=frld.pair_mu(fld,Tld);
              int ret25=fr25.pair_mu(f25,T25);
              
              if (ret==0) {
                
                // Test with inc_rest_mass=false
                part1_t f2=f;
                f2.inc_rest_mass=false;
                f2.mu-=f2.m;
                fr.pair_mu(f2,T);
                t.test_rel(f.ed,f2.ed+f2.m*f2.n,1.0e-11,"irm false ed");
                t.test_rel(f.en,f2.en,1.0e-13,"irm false en");
                
                // Test with non_interacting=false
                part1_t f3=f;
                f3.non_interacting=false;
                f3.nu=f3.mu;
                f3.ms=f3.m;
                f3.m*=sqrt(2.0);
                fr.pair_mu(f3,T);
                t.test_rel(f.n,f3.n,1.0e-14,"ni false n");
                t.test_rel(f.en,f3.en,1.0e-14,"ni false en");
                
                // Test with both equal to false
                part1_t f4=f3;
                f4.inc_rest_mass=false;
                f4.non_interacting=false;
                f4.nu-=f4.m;
                fr.pair_mu(f4,T);
                t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed");
                t.test_rel(f3.en,f4.en,1.0e-13,"both false en");
                
                if (t.get_success()==false) {
                  std::cout << "1" << std::endl;
                  exit(-1);
                }
              }
              
              int idn=-2, iden=-2, ildn=-2, ilden=-2;
              if (ret==0 && retld==0) {
                idn=count_digits_same(f.n,fld.n);
                iden=count_digits_same(f.en,fld.en);
              }
              pmu_n+=idn;
              pmu_en+=iden;
              if (retld==0 && ret25==0) {
                ildn=count_digits_same(fld.n,f25.n);
                ilden=count_digits_same(fld.en,f25.en);
              }
              pmu_ld_en+=ilden;
              pmu_ld_n+=ildn;
              
              std::cout.width(2);
              std::cout << idn << " ";
              std::cout.width(2);
              std::cout << iden << " ";
              std::cout.width(2);
              std::cout << ildn << " ";
              std::cout.width(2);
              std::cout << ilden << " ";
              
              int x=-2, xld=-2, x25=-2;
              
              if (ret==0) {
                fp1_t pr2=-f.ed+f.n*f.mu+T*f.en;
                x=count_digits_same(f.pr,pr2);
              }
              pmu_ti+=x;
              std::cout.width(2);
              std::cout << x << " ";
              
              if (retld==0) {
                fp2_t pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
                xld=count_digits_same(fld.pr,pr2ld);
              }
              pmu_ld_ti+=xld;
              std::cout.width(2);
              std::cout << xld << " ";
              
              if (ret25==0) {
                fp3_t pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
                x25=count_digits_same(f25.pr,pr225);
              }
              pmu_25_ti+=x25;
              std::cout.width(2);
              std::cout << x25 << " pmu" << std::endl;
              
            } else {
              
              std::cout << std::endl;
              
            }
            
            count++;
            
          }
        }
      }
      
      std::cout << "pair_mu density (double <-> long double): "
           << pmu_n << std::endl;
      std::cout << "pair_mu entropy (double <-> long double): "
           << pmu_en << std::endl;
      std::cout << "pair_mu density (long double <-> cdf_25): "
           << pmu_ld_n << std::endl;
      std::cout << "pair_mu entropy (long double <-> cdf_25): "
           << pmu_ld_en << std::endl;
      std::cout << "pair_mu ti: " << pmu_ti << std::endl;
      std::cout << "pair_mu long double ti: " << pmu_ld_ti << std::endl;
      std::cout << "pair_mu cpp_dec_float_25 ti: " << pmu_25_ti << std::endl;
      std::cout << std::endl;

      t.test_gen(pmu_n>=pmu_n_max,"pmu_n");
      t.test_gen(pmu_en>=pmu_en_max,"pmu_en");
      t.test_gen(pmu_ld_n>=pmu_ld_n_max,"pmu_ld_n");
      t.test_gen(pmu_ld_en>=pmu_ld_en_max,"pmu_ld_en");
      t.test_gen(pmu_ti>=pmu_ti_max,"pmu_ti");
      t.test_gen(pmu_ld_ti>=pmu_ld_ti_max,"pmu_ld_ti");
      t.test_gen(pmu_25_ti>=pmu_25_ti_max,"pmu_25_ti");
      
      std::cout.precision(6);
      
      return;
    }
    
  };
  
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
        // It's important that m<m* for bosons to work
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

        This function sets either \c mu or \c nu, depending on whether
        or not we're comparing interacting or noninteracting
        particles. The other chemical potential is set to zero to
        ensure that the correct one is being used in the calculation.
    */
    template<class part_t>
    void set_chem_pot(part_t &p, fp_t psi, fp_t T, size_t k,
                      bool nr_mode) {
      /*
        if (k==0) {
        std::cout << "Function calc_mu(), include rest mass:"
        << std::endl;
        } else if (k==1) {
        std::cout << "Function calc_mu(), without rest mass:"
        << std::endl;
        } else if (k==2) {
        std::cout << "Function calc_mu(), include rest mass, "
        << "interacting:" << std::endl;
        } else {
        std::cout << "Function calc_mu(), without rest mass, "
        << "interacting:" << std::endl;
        }
      */
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

        This function compares the number densities in \c p and \c
        exact. If the difference is greater than that recorded in \c
        max, and if so then all the "max" quantities are updated.
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_density(part1_t &p, part2_t &exact, part3_t &max, size_t k,
                       fp_t T, fp_t mot, fp_t psi,
                       fp_t &mu_max, fp_t &m_max,
                       fp_t &T_max, fp_t &mot_max, fp_t &psi_max,
                       fp_t &ret_local) {
      if (fabs((p.n-exact.n)/exact.n)>max.n) {
	max.n=fabs((p.n-exact.n)/exact.n);
	if (max.n>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.n;
	}
      }
      return;
    }      
    
    /** \brief Check the chemical potential against the exact result 
	and update 
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_chem_pot(part1_t &p, part2_t &exact, part3_t &max, size_t k,
                        fp_t T, fp_t mot, fp_t psi,
                        fp_t &mu_max, fp_t &m_max,
                        fp_t &T_max, fp_t &mot_max, fp_t &psi_max,
                        fp_t &ret_local) {
      if (k>=2) {
	if (fabs((p.nu-exact.nu)/exact.nu)>max.mu) {
	  max.mu=fabs((p.nu-exact.nu)/exact.nu);
	  if (max.mu>ret_local) {
	    mu_max=p.nu;
	    m_max=p.ms;
	    T_max=T;
	    mot_max=mot;
	    psi_max=psi;
	    ret_local=max.n;
	  }
	}
      } else {
	if (fabs((p.mu-exact.mu)/exact.mu)>max.mu) {
	  max.mu=fabs((p.mu-exact.mu)/exact.mu);
	  if (max.mu>ret_local) {
	    mu_max=p.mu;
	    m_max=p.m;
	    T_max=T;
	    mot_max=mot;
	    psi_max=psi;
	    ret_local=max.n;
	  }
	}
      }
      return;
    }
      
    /** \brief Check the energy density, pressure, and entropy against
	the exact result and update
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_eps(part1_t &p, part2_t &exact, part3_t &max, size_t k,
                   fp_t T, fp_t mot, fp_t psi,
                   fp_t &mu_max, fp_t &m_max,
                   fp_t &T_max, fp_t &mot_max, fp_t &psi_max,
                   fp_t &ret_local) {
      if (fabs((p.ed-exact.ed)/exact.ed)>max.ed) {
	max.ed=fabs((p.ed-exact.ed)/exact.ed);
	if (max.ed>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.ed;
	}
      }
      if (fabs((p.pr-exact.pr)/exact.pr)>max.pr) {
	max.pr=fabs((p.pr-exact.pr)/exact.pr);
	if (max.pr>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.pr;
	}
      }
      if (fabs((p.en-exact.en)/exact.en)>max.en) {
	max.en=fabs((p.en-exact.en)/exact.en);
	if (max.en>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.en;
	}
      }
      return;
    }      
    
    /** \brief Check the three derivatives against
	the exact result and update
    */
    template<class part1_t, class part2_t, class part3_t>
    void check_derivs(part1_t &p, part2_t &exact, part3_t &max, size_t k,
                      fp_t T, fp_t mot, fp_t psi,
                      fp_t &mu_max, fp_t &m_max,
                      fp_t &T_max, fp_t &mot_max, fp_t &psi_max,
                      fp_t &ret_local) {
      if (fabs((p.dndT-exact.dndT)/exact.dndT)>max.dndT) {
	max.dndT=fabs((p.dndT-exact.dndT)/exact.dndT);
	if (max.dndT>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.dndT;
	}
      }
      if (fabs((p.dndmu-exact.dndmu)/exact.dndmu)>max.dndmu) {
	max.dndmu=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	if (max.dndmu>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.dndmu;
	}
      }
      if (fabs((p.dsdT-exact.dsdT)/exact.dsdT)>max.dsdT) {
	max.dsdT=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	if (max.dsdT>ret_local) {
	  if (k>=2) {
	    mu_max=p.nu;
	    m_max=p.ms;
	  } else {
	    mu_max=p.mu;
	    m_max=p.m;
	  }
	  T_max=T;
	  mot_max=mot;
	  psi_max=psi;
	  ret_local=max.dsdT;
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
    fp_t part_calibrate(part_t &p, thermo_t &th,
                        std::string file, 
                        bool test_pair=true,
                        bool test_density=true,
                        bool nr_mode=false,
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
	  "' in part_calibrate(). Max filename?";
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
      part_t max, dev, exact;
      fp_t m_max=0.0, mu_max=0.0, T_max=0.0, mot_max=0.0, psi_max=0.0;
      p.non_interacting=true;

      // This counts tests, 4*4*3 times the number of lines in the
      // calibrate data file
      int count=0;
      
      // ----------------------------------------------------------------
      // First pass, test calc_mu()
      
      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {
        
        if (verbose>1) {
          if (k==0) {
            std::cout << "Function calc_mu(), include rest mass:\n"
                      << std::endl;
          } else if (k==1) {
            std::cout << "Function calc_mu(), without rest mass:\n"
                      << std::endl;
          } else if (k==2) {
            std::cout << "Function calc_mu(), include rest mass, "
                      << "interacting:\n" << std::endl;
          } else {
            std::cout << "Function calc_mu(), without rest mass, "
                      << "interacting:\n" << std::endl;
          }
        }

	fp_t ret_local=0;
        ti_local=0;
      
	// Initialize storage
	dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	max.n=0.0; max.ed=0.0; max.pr=0.0; max.en=0.0;
    
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
            count++;

	    if (verbose>1) {
	      std::cout.precision(4);
	      if (k>=2) {
		std::cout << "T,ms,nu,psi,mot,cnt: " << T << " "
			  << p.ms << " " << p.nu << " "
			  << psi << " " << mot << " " << count << std::endl;
	      } else {
		std::cout << "T,m,mu,psi,mot,cnt: " << T << " "
			  << p.m << " " << p.mu << " " 
			  << psi << " " << mot << " " << count << std::endl;
	      }
	      std::cout.precision(6);
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
              fp_t val;
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
              (p,exact,max,k,T,mot,psi,mu_max,m_max,T_max,
               mot_max,psi_max,ret_local);
	  
	    check_eps<part_t,part_t,part_t>
              (p,exact,max,k,T,mot,psi,mu_max,m_max,T_max,
               mot_max,psi_max,ret_local);

	    if (ret_local>ret) {
	      ret=ret_local;
	    }
	    if (ti_local>ti_max) {
	      ti_max=ret_local;
	    }
	  
	    if (verbose>1) {
	      std::cout.precision(5);
              /*
                AWS, 6/21/22: took this out because it appeared
                to be duplicated by the output above
                
                if (k>=2) {
		std::cout << "T,ms,nu,psi,mot: " << T << " "
                << p.ms << " " << p.nu << " "
                << psi << " " << mot << std::endl;
                } else {
		std::cout << "T,m,mu,psi,mot: " << T << " "
                << p.m << " " << p.mu << " " 
                << psi << " " << mot << std::endl;
                }
              */
	      std::cout.precision(5);
	      std::cout << "n,ed,pr,en: " << std::endl;
	      std::cout << "comput: " << p.n << " " << p.ed << " "
			<< p.pr << " " << p.en << std::endl;
	      std::cout << "exact : " << exact.n << " " << exact.ed << " " 
			<< exact.pr << " " << exact.en << std::endl;
	      std::cout << "maxdev: " << max.n << " " << max.ed << " " 
			<< max.pr << " " << max.en << std::endl;
              if (ti_count>0) {
                std::cout << "check ti (avg, curr, max): "
                          << ti_test/ti_count << " "
                          << ti_local << " " << ti_max << std::endl;
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
	    std::cout << "Function calc_mu(), include rest mass:"
		      << std::endl;
	  } else if (k==1) {
	    std::cout << "Function calc_mu(), without rest mass:"
		      << std::endl;
	  } else if (k==2) {
	    std::cout << "Function calc_mu(), include rest mass, "
		      << "interacting:" << std::endl;
	  } else {
	    std::cout << "Function calc_mu(), without rest mass, "
		      << "interacting:" << std::endl;
	  }

	  std::cout << "Average performance: " << std::endl;
	  std::cout << "n: " << dev.n << " ed: " << dev.ed << " pr: " 
		    << dev.pr << " en: " << dev.en << std::endl;
	  std::cout << "Worst case: " << std::endl;
	  std::cout << "n: " << max.n << " ed: " << max.ed << " pr: " 
		    << max.pr << " en: " << max.en << std::endl;
          std::cout << "Worst case occurred at:" << std::endl;
	  std::cout << "mu: " << mu_max << " m: " << m_max 
		    << " T: " << T_max << " mot: " << mot_max
		    << "\n\tpsi: " << psi_max << std::endl;
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

      if (test_density) {
      
        // k=0,2 are with rest mass, k=1,3 are without
        // k=0,1 are non-interacting, k=2,3 are interacting
        for(size_t k=0;k<4;k++) {

          if (verbose>1) {
            if (k==0) {
              std::cout << "Function calc_density(), include rest mass, "
                        << "noninteracting:" << std::endl;
            } else if (k==1) {
              std::cout << "Function calc_density(), without rest mass, "
                        << "noninteracting:" << std::endl;
            } else if (k==2) {
              std::cout << "Function calc_density(), include rest mass, "
                        << "interacting:" << std::endl;
            } else {
              std::cout << "Function calc_density(), without rest mass, "
                        << "interacting:" << std::endl;
            }
          }
        
          fp_t ret_local=0.0;
          ti_local=0;
      
          // Initialize storage
          dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
          max.mu=0.0; max.ed=0.0; max.pr=0.0; max.en=0.0;
    
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
              set_mass_flags(exact,mot,T,k);
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
            
              count++;

              if (verbose>1) {
                std::cout.precision(5);
                if (k>=2) {
                  std::cout << "T,ms,n,psi,mot,count: " << T << " "
                            << p.ms << " " << p.n << " "
                            << psi << " " << mot << " " << count << std::endl;
                } else {
                  std::cout << "T,m,n,psi,mot,count: " << T << " "
                            << p.m << " " << p.n << " "
                            << psi << " " << mot << " " << count << std::endl;
                }
              }

              th.calc_density(p,T);
              if (th.verify_ti) {
                fp_t val;
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
                if (false && dev.mu>0.1) {
                  std::cout << "p.mu,exact.mu: "
                            << p.mu << " " << exact.mu << std::endl;
                  std::cout << "T,m,n,psi,mot,count: " << T << " "
                            << p.m << " " << p.n << " "
                            << psi << " " << mot << " " << count << std::endl;
                  O2SCL_ERR("Chemical potential match failed.",
                            o2scl::exc_einval);
                }
              }
              dev.ed+=fabs((p.ed-exact.ed)/exact.ed);
              dev.pr+=fabs((p.pr-exact.pr)/exact.pr);
              dev.en+=fabs((p.en-exact.en)/exact.en);

              cnt++;

              check_chem_pot<part_t>(p,exact,max,k,T,mot,psi,mu_max,
                                     m_max,T_max,mot_max,psi_max,ret_local);
              check_eps<part_t>(p,exact,max,k,T,mot,psi,mu_max,m_max,T_max,
                                mot_max,psi_max,ret_local);
            
              if (verbose>1) {
                std::cout.precision(6);
                if (k>=2) {
                  std::cout << "nu,ed,pr,en: " << std::endl;
                  std::cout << "comput: " << p.nu << " " << p.ed << " "
                            << p.pr << " " << p.en << std::endl;
                  std::cout << "exact : " << exact.nu << " " << exact.ed 
                            << " " << exact.pr << " " << exact.en << std::endl;
                } else {
                  std::cout << "mu,ed,pr,en: " << std::endl;
                  std::cout << "comput: " << p.mu << " " << p.ed << " "
                            << p.pr << " " << p.en << std::endl;
                  std::cout << "exact : " << exact.mu << " " << exact.ed 
                            << " " << exact.pr << " " << exact.en << std::endl;
                }
                std::cout << "maxdev: " << max.mu << " " << max.ed << " " 
                          << max.pr << " " << max.en << std::endl;
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
            std::cout << "mu: " << max.mu << " ed: " << max.ed << " pr: " 
                      << max.pr << " en: " << max.en << std::endl;
            std::cout << "mu: " << mu_max << " m: " << m_max
                      << " T: " << T_max << " mot: " << mot_max
                      << "\n\tpsi: " << psi_max << std::endl;
            if (ti_count>0) {
              std::cout << "check ti: " << ti_test/ti_count << " "
                        << ti_local << std::endl;
            }
            std::cout << "ret_local,ret x: " << ret_local << " "
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

      }

      if (test_pair) {
  
	// ----------------------------------------------------------------
	// Third pass, test pair_mu() 

	// k=0,2 are with rest mass, k=1,3 are without
	// k=0,1 are non-interacting, k=2,3 are interacting
	for(size_t k=0;k<4;k++) {

	  fp_t ret_local=0.0;
          ti_local=0;

	  // Initialize storage
	  dev.n=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
	  max.n=0.0; max.ed=0.0; max.pr=0.0; max.en=0.0;
    
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
              count++;
            
	      if (verbose>1) {
		std::cout.precision(5);
		std::cout << "T,m,mu,psi,mot,count: " << T << " "
                          << p.m << " "
			  << p.mu << " " << psi << " " << mot << " "
                          << count << std::endl;
              }
              
	      th.pair_mu(p,T);
              if (th.verify_ti) {
                fp_t val;
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
                (p,exact,max,k,T,mot,psi,mu_max,m_max,T_max,
                 mot_max,psi_max,ret_local);
	    
	      check_eps<part_t,part_t,part_t>
                (p,exact,max,k,T,mot,psi,mu_max,m_max,T_max,
                 mot_max,psi_max,ret_local);
              
	      if (verbose>1) {
		std::cout.precision(6);
		std::cout << "n,ed,pr,en: " << std::endl;
		std::cout << "comput: " << p.n << " " << p.ed << " "
			  << p.pr << " " << p.en << std::endl;
		std::cout << "exact : " << exact.n << " " << exact.ed << " " 
			  << exact.pr << " " << exact.en << std::endl;
		std::cout << "maxdev: " << max.n << " " << max.ed << " " 
			  << max.pr << " " << max.en << std::endl;
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
	    std::cout << "n: " << max.n << " ed: " << max.ed << " pr: " 
		      << max.pr << " en: " << max.en << std::endl;
	    std::cout << "mu: " << mu_max << " m: " << m_max
		      << " T: " << T_max << " mot: " << mot_max
		      << "\n\tpsi: " << psi_max << std::endl;
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

        if (test_density) {
        
          // k=0,2 are with rest mass, k=1,3 are without
          // k=0,1 are non-interacting, k=2,3 are interacting
          for(size_t k=0;k<4;k++) {

            fp_t ret_local=0.0;
            ti_local=0;
	
            // Initialize storage
            dev.mu=0.0; dev.ed=0.0; dev.pr=0.0; dev.en=0.0;
            max.mu=0.0; max.ed=0.0; max.pr=0.0; max.en=0.0;
    
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
                count++;

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
                    std::cout << "T,ms,n,psi,mot,count: "
                              << T << " " << p.ms << " " 
                              << p.n << " " << psi << " " << mot << " "
                              << count << std::endl;
                  } else {
                    std::cout << "T,m,n,psi,mot,count: "
                              << T << " " << p.m << " " 
                              << p.n << " " << psi << " " << mot << " "
                              << count << std::endl;
                  }
                }
              
                th.pair_density(p,T);
                if (th.verify_ti) {
                  fp_t val;
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

                check_chem_pot<part_t>(p,exact,max,k,T,mot,psi,
                                       mu_max,m_max,T_max,
                                       mot_max,psi_max,ret_local);

                check_eps<part_t>(p,exact,max,k,T,mot,psi,mu_max,
                                  m_max,T_max,
                                  mot_max,psi_max,ret_local);
	    
	    
                if (verbose>1) {
                  std::cout.precision(6);
                  if (k>=2) {
                    std::cout << "nu,ed,pr,en: " << std::endl;
                    std::cout << "comput: " << p.nu << " " << p.ed << " "
                              << p.pr << " " << p.en << std::endl;
                    std::cout << "exact : " << exact.nu << " "
                              << exact.ed << " " << exact.pr << " "
                              << exact.en << std::endl;
                  } else {
                    std::cout << "mu,ed,pr,en: " << std::endl;
                    std::cout << "comput: " << p.mu << " " << p.ed << " "
                              << p.pr << " " << p.en << std::endl;
                    std::cout << "exact : " << exact.mu << " "
                              << exact.ed << " " << exact.pr << " "
                              << exact.en << std::endl;
                  }
                  std::cout << "maxdev: " << max.mu << " " << max.ed << " " 
                            << max.pr << " " << max.en << std::endl;
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
              std::cout << "mu: " << max.mu << " ed: " << max.ed << " pr: " 
                        << max.pr << " en: " << max.en << std::endl;
              std::cout << "mu: " << mu_max << " m: " << m_max
                        << " T: " << T_max << " mot: " << mot_max
                        << "\n\tpsi: " << psi_max << std::endl;
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

          // End of 'if (test_density)'
        }
        
	// End of 'if (test_pair)'
      }

      // ----------------------------------------------------------------
      // Return to the original values 

      p=orig;

      if (th.verify_ti && verbose>0) {
        std::cout << "check ti, ti_max: " << ti_test/ti_count << " "
                  << ti_max << std::endl;
      }
  
      return ret;
    }

  };

  typedef part_calibrate_class_tl<double> part_calibrate_class;

}

#endif
