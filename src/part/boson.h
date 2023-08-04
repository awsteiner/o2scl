/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
    \brief File defining \ref o2scl::boson
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

namespace o2scl {

  /** \brief Boson class
   */
  template<class fp_t=double> class boson_tl : public part_tl<fp_t> {
    
  public:

    /// Create a boson with mass \c mass and degeneracy \c dof 
    boson_tl(fp_t mass=0.0, fp_t dof=0.0);

    /** \brief The condensate
	
	The condensate variable is provided principally for
	user storage and is mostly ignored by \o2p classes. 
    */
    fp_t co;

    /** \brief Calculate properties of massless bosons 
	
	The expressions used are exact. The chemical potentials are
	ignored.
    */
    virtual void massless_calc(fp_t temper);
    
    /// Return string denoting type ("boson")
    virtual const char *type() { return "boson"; }

  };

  typedef boson_tl<double> boson;

  /** \brief Compute the thermodynamic properties of a boson 
      [abstract base]
   */
  class boson_thermo {
    
  public:

    /// Object for Bessel-exp integrals
    bessel_K_exp_integ_gsl be_integ;
    
    /** \brief Compute a term in the nondegenerate expansion
     */
    void ndeg_terms(size_t j, double tt,
                    double xx, double m, bool inc_rest_mass,
                    bool inc_antip, double &pterm, double &nterm,
                    double &enterm, double &edterm) {
      
      double dj=((double)j);
      double jot=dj/tt;

      double K2j=be_integ.K2exp(jot);
      double K1j=be_integ.K1exp(jot);
      if (inc_antip==false) {
        pterm=exp(jot*xx)/jot/jot*K2j;
        nterm=pterm*jot/m;
        double enterm1=(4*tt-dj*xx-dj)/dj/tt*nterm;
        double enterm2=exp(jot*xx)/dj*K1j/m;
        enterm=enterm1+enterm2;
        edterm=(K1j*dj+3.0*K2j*tt)/jot/dj*exp(xx*jot);
      } else {
        double co=cosh(jot*(xx+1));
        double si=sinh(jot*(xx+1));
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
    bool calc_mu_ndeg(boson &b, double temper, 
                      double prec=1.0e-18, bool inc_antip=false) {
      
      if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }
      
      // Compute psi and tt
      double psi, psi_num;
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
      double tt=temper/b.ms;
      
      // Prefactor 'd' in Johns96
      double prefac=b.g/2.0/o2scl_const::pi2*pow(b.ms,4.0);
      
      // One term is always used, so only values of max_term greater than
      // 0 are useful.
      static const size_t max_term=200;
      
      // ─────────────────────────────────────────────────────────────────
      // Return early if the last term is going to be too large.
      
      // Ratio of last term to first term in the pressure expansion
      double rat, pterm_1, pterm_max;
      double pterm, nterm, enterm, edterm;

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

      double first_term=0.0;
      b.pr=0.0;
      b.n=0.0;
      b.en=0.0;
      
      for(size_t j=1;j<=max_term;j++) {
        
        ndeg_terms(j,tt,psi*tt,b.ms,b.inc_rest_mass,inc_antip,
                   pterm,nterm,enterm,edterm);
        
        if (j==1) first_term=pterm;
        b.pr+=pterm;
        b.n+=nterm;
        b.en+=enterm;
        
        // If the first term is zero, then the rest of the terms
        // will be zero so just return early
        if (first_term==0.0) {
          b.pr=0.0;
          b.n=0.0;
          b.ed=0.0;
          b.en=0.0;
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
    virtual void calc_mu(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties as function of
	density
    */
    virtual void calc_density(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of chemical potential
    */
    virtual void pair_mu(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of density
    */
    virtual void pair_density(boson &b, double temper)=0;

  };
  
}

#endif
