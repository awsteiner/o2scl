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
#ifndef O2SCL_BOSON_H
#define O2SCL_BOSON_H

/** \file boson.h
    \brief File defining \ref o2scl::boson_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/polylog.h>

#include <o2scl/part.h>

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_SET_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif
#endif

namespace o2scl {

  /** \brief Boson class
   */
  template<class fp_t=double> class boson_tl : public part_tl<fp_t> {
    
  public:

    /// Create a boson with mass \c mass and degeneracy \c dof 
    boson_tl(fp_t mass=0.0, fp_t dof=0.0) {
      co=0.0;
      this->m=mass;
      this->g=dof;
    }      

    /** \brief The condensate
	
	The condensate variable is provided principally for
	user storage and is provided mostly for the user
    */
    fp_t co;

    /** \brief Calculate properties of massless bosons 
	
	The expressions used are exact. The chemical potentials are
	ignored.
    */
    virtual void massless_calc(fp_t temper) {

      fp_t zeta3_loc=boost::math::constants::zeta_three<fp_t>();
      this->n=this->g*pow(temper,3.0)*zeta3_loc/o2scl_const::pi2;
      this->ed=this->g*o2scl_const::pi2*pow(temper,4.0)/30;
      this->pr=this->ed/3;
      this->en=this->g*o2scl_const::pi2*pow(temper,3.0)/45*2;
      
      return;
    }      
    
    /// Return string denoting type ("boson")
    virtual const char *type() { return "boson"; }

  };

  typedef boson_tl<double> boson;

  typedef boson_tl<long double> boson_ld;
  
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  typedef boson_tl<boost::multiprecision::number<
                     boost::multiprecision::cpp_dec_float<25> > > boson_cdf25;
  
#endif

  /** \brief Compute the thermodynamic properties of a boson 
      [abstract base]
   */
  template<class be_inte_t=bessel_K_exp_integ_gsl,
           class fp_t=double> class boson_thermo_tl {
    
  public:

    /// Object for Bessel-exp integrals
    be_inte_t be_integ;
    
    /** \brief Compute a term in the nondegenerate expansion
     */
    void ndeg_terms(size_t j, fp_t tt,
                    fp_t xx, fp_t m, bool inc_rest_mass,
                    bool inc_antip, fp_t &pterm, fp_t &nterm,
                    fp_t &enterm, fp_t &edterm) {
      
      fp_t dj=((fp_t)j);
      fp_t jot=dj/tt;

      fp_t K2j=be_integ.K2exp(jot);
      fp_t K1j=be_integ.K1exp(jot);
      if (inc_antip==false) {
        pterm=exp(jot*xx)/jot/jot*K2j;
        nterm=pterm*jot/m;
        fp_t enterm1=(4*tt-dj*xx-dj)/dj/tt*nterm;
        fp_t enterm2=exp(jot*xx)/dj*K1j/m;
        enterm=enterm1+enterm2;
        edterm=(K1j*dj+3*K2j*tt)/jot/dj*exp(xx*jot);
      } else {
        fp_t co=cosh(jot*(xx+1));
        fp_t si=sinh(jot*(xx+1));
        pterm=exp(-jot)*2*co/jot/jot*K2j;
        nterm=pterm*tanh(jot*(xx+1))*jot/m;
        enterm=(-2*K2j*(xx+1)*si/dj/m+
                2*(K1j*dj+4*K2j*tt)*co/dj/dj/m)*exp(-jot);
        edterm=2*exp(-jot)/jot*(K1j*co-2*K2j*xx*si-2*K2j*si+
                                3*K2j*tt*co/dj);
      }
                    
      return;
    }
    
    /** \brief Non-degenerate expansion for bosons
        
    */
    bool calc_mu_ndeg(boson_tl<fp_t> &b, fp_t temper, 
                      fp_t prec=1.0e-18, bool inc_antip=false) {
      
      if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }
      
      // Compute psi and tt
      fp_t psi, psi_num;
      if (b.inc_rest_mass) {
        psi_num=b.nu-b.ms;
      } else {
        if (b.non_interacting) {
          psi_num=b.nu;
        } else {
          psi_num=b.nu+b.m-b.ms;
        }
      }
      psi=psi_num/temper;
      fp_t tt=temper/b.ms;
      
      // Prefactor 'd' in Johns96
      fp_t prefac=b.g/2/o2scl_const::pi2*pow(b.ms,4);
      
      // One term is always used, so only values of max_term greater than
      // 0 are useful.
      static const size_t max_term=200;
      
      // ─────────────────────────────────────────────────────────────────
      // Return early if the last term is going to be too large.
      
      // Ratio of last term to first term in the pressure expansion
      fp_t rat, pterm_1, pterm_max;
      fp_t pterm, nterm, enterm, edterm;

      ndeg_terms(1,tt,psi*tt,b.ms,b.inc_rest_mass,inc_antip,
                 pterm_1,nterm,enterm,edterm);
      ndeg_terms(max_term,tt,psi*tt,b.ms,b.inc_rest_mass,inc_antip,
                 pterm_max,nterm,enterm,edterm);
      rat=pterm_max/pterm_1;
      
      // If the ratio between the last term and the first term is 
      // not small enough, return false
      if (isfinite(rat) && rat>prec) {
        return false;
      }
      
      // ─────────────────────────────────────────────────────────────────
      // Go through term by term and see if we obtain the requested
      // precision

      fp_t first_term=0;
      b.pr=0;
      b.n=0;
      b.en=0;
      
      for(size_t j=1;j<=max_term;j++) {
        
        ndeg_terms(j,tt,psi*tt,b.ms,b.inc_rest_mass,inc_antip,
                   pterm,nterm,enterm,edterm);
        
        if (j==1) first_term=pterm;
        b.pr+=pterm;
        b.n+=nterm;
        b.en+=enterm;
        
        // If the first term is zero, then the rest of the terms
        // will be zero so just return early
        if (first_term==0) {
          b.pr=0;
          b.n=0;
          b.ed=0;
          b.en=0;
          return true;
        }
        
        // Stop if the last term is sufficiently small compared to
        // the first term
        if (j>1 && fabs(pterm)<prec*fabs(first_term)) {
          b.pr*=prefac;
          b.n*=prefac;
          b.en*=prefac;
          b.ed=-b.pr+b.nu*b.n+temper*b.en;
          return true;
        }
        
        // End of 'for(size_t j=1;j<=max_term;j++)'
      }

      // We failed to add enough terms, so return false
      return false;
    }
    
    /** \brief Calculate thermodynamic properties as function of
	chemical potential
    */
    virtual void calc_mu(boson_tl<fp_t> &b, fp_t temper)=0;

    /** \brief Calculate thermodynamic properties as function of
	density
    */
    virtual void calc_density(boson_tl<fp_t> &b, fp_t temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of chemical potential
    */
    virtual void pair_mu(boson_tl<fp_t> &b, fp_t temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of density
    */
    virtual void pair_density(boson_tl<fp_t> &b, fp_t temper)=0;

  };

  typedef boson_thermo_tl<> boson_thermo;
  
}

#endif
