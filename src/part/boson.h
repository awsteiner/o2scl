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
  class boson : public part {
    
  public:

    /// Create a boson with mass \c mass and degeneracy \c dof 
    boson(double mass=0.0, double dof=0.0);

    /** \brief The condensate
	
	The condensate variable is provided principally for
	user storage and is mostly ignored by \o2p classes. 
    */
    double co;

    /** \brief Calculate properties of massless bosons 
	
	The expressions used are exact. The chemical potentials are
	ignored.
    */
    virtual void massless_calc(double temper);
    
    /// Return string denoting type ("boson")
    virtual const char *type() { return "boson"; }

  };

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
      if (inc_antip==false) {
        double K1j=be_integ.K1exp(jot);
        pterm=exp(jot*xx)/jot/jot*K2j;
        nterm=pterm*jot/m;
        double enterm1=(4*tt-dj*xx-dj)/dj/tt*nterm;
        double enterm2=exp(jot*xx)/dj*K1j/m;
        enterm=enterm1+enterm2;
        edterm=(K1j*dj+3.0*K2j*tt)/jot/dj*exp(xx*jot);
      } else {
        double K3j=be_integ.K3exp(jot);
        pterm=exp(-jot)*2*cosh(jot*(xx+1))/jot/jot*K2j;
        nterm=pterm*tanh(jot*(xx+1))*jot/m;
        // entropy and energy density terms not right yet
        double enterm1=-(1+xx)/tt*nterm/m;
        double enterm2=2*exp(-jot)/dj*cosh(jot*(xx+1))*K3j/m;
          enterm=enterm1-enterm2;
        edterm=2/jot/dj*exp(-jot)*(K3j*dj*cosh(jot*(xx+1))-
                                     2*K2j*dj*xx*sinh(jot*(xx+1))-
                                     2*K2j*dj*sinh(jot*(xx+1))-
                                     K2j*tt*cosh(jot*(xx+1)));
      }
                    
      return;
    }
    
    /** \brief Non-degenerate expansion for boson
        
        Attempts to evaluate thermodynamics of a non-degenerate
        fermion. If the result is accurate to within the requested
        precision, this function returns <tt>true</tt>, and otherwise
        this function returns <tt>false</tt> and the values in stored
        in the <tt>pr</tt>, <tt>n</tt>, <tt>en</tt>, and <tt>ed</tt>
        field are meaningless.
        
        If \f$ \mu \f$ is negative and sufficiently far from zero,
        then the thermodynamic quantities are smaller than the smallest
        representable double-precision number. In this case,
        this function will return <tt>true</tt> and report all
        quantities as zero.
        
        \verbatim embed:rst
        The following uses the notation of [Johns96]_.
        \endverbatim

        Defining \f$ \psi \equiv (\mu-m)/T \f$, \f$ t \equiv T/m \f$,
        and \f$ d \equiv g~m^4/(2 \pi^2) \f$ the pressure 
        in the non-degenerate limit (\f$ \psi \rightarrow - \infty \f$)
        is
        \f[
        P = d \sum_{n=1}^{\infty} P_n
        \f]
        where 
        \f[
        P_n \equiv \left(-1\right)^{n+1} \left(\frac{t^2}{n^2}\right)
        e^{n \left(\psi+1/t\right)} K_2 \left( \frac{n}{t} \right)
        \f]
        The density is then
        \f[
        n = d \sum_{n=1}^{\infty} \frac{n P_n}{T}
        \f]
        and the entropy density is
        \f[
        s = \frac{d}{m} \sum_{n=1}^{\infty} \left\{ \frac{2 P_n}{t}
        -\frac{n P_n}{t^2}+ 
        \frac{\left(-1\right)^{n+1}}{2 n} 
        e^{n \left(\psi+1/t\right)} \left[ K_1 \left( \frac{n}{t} 
        \right)+K_3 \left( \frac{n}{t} \right) \right]
        \right\}
        \f]
        
        This function is accurate over a wide range of conditions
        when \f$ \psi < -4 \f$.
        
        The ratio of the nth term to the first term in the pressure
        series is
        \f[
        R_n \equiv \frac{P_{n}}{P_{1}} = \frac{(-1)^{n+1} 
        e^{(n-1)(\psi+1/t)} K_2(n/t) }{n^2 K_2(1/t)}
        \f]
        This function currently uses 20 terms in the series and
        immediately returns <tt>false</tt> if \f$ |R_{20}| \f$
        is greater than <tt>prec</tt>
        
        In the nondegenerate and nonrelativistic (\f$ t \rightarrow 0
        \f$) limit, the argument to the Bessel functions and the
        exponential becomes too large. In this case, it's better to
        use the expansions, e.g. for \f$ x \equiv n/t \rightarrow
        \infty \f$,
        \f[
        \sqrt{\frac{2 x}{\pi}} e^{x} K_2(x) \approx
        1 + \frac{3}{8 x} - \frac{15}{128 x^2} + ...
        \f]

        \comment
        AWS, 6/28/21: This comment doesn't make sense to me, 
        so I'm taking it out. 

        The current code currently goes up to \f$ x^{-12} \f$ in the
        expansion, which is enough for the default precision of \f$
        10^{-18} \f$ since \f$ (20/700)^{12} \sim 10^{-19} \f$.
        \endcomment
    */
    bool calc_mu_ndeg(boson &f, double temper, 
                      double prec=1.0e-18, bool inc_antip=false) {
      
      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
      
      // Compute psi and tt
      double psi, psi_num;
      if (f.inc_rest_mass) {
        psi_num=f.nu-f.ms;
      } else {
        if (f.non_interacting) {
          psi_num=f.nu;
        } else {
          psi_num=f.nu+f.m-f.ms;
        }
      }
      psi=psi_num/temper;
      double tt=temper/f.ms;
      
      // Return false immediately if we're degenerate
      if (inc_antip==false && psi>0.0) return false;
      
      // Prefactor 'd' in Johns96
      double prefac=f.g/2.0/o2scl_const::pi2*pow(f.ms,4.0);
      
      // One term is always used, so only values of max_term greater than
      // 0 are useful.
      static const size_t max_term=200;
      
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
      
      // -----------------------------------------------------
      // Return early if the last term is going to be too large.
      
      // Ratio of last term to first term in the pressure expansion
      double rat;
      double dj1=((double)max_term), jot1=max_term/tt;
      double dj2=1.0, jot2=1.0/tt;

      /*
      
      if (inc_antip==false) {
        rat=exp(dj1*psi)/jot1/jot1*be_integ.K2exp(jot1);
        rat/=exp(dj2*psi)/jot2/jot2*be_integ.K2exp(jot2);
        //std::cout << "rat: " << rat << std::endl;
      } else {
        if (f.inc_rest_mass) {
          rat=exp(-jot1)*2.0*cosh(dj1*f.nu/temper)/jot1/jot1*
            be_integ.K2exp(jot1);
          rat/=exp(-jot2)*2.0*cosh(dj2*f.nu/temper)/jot2/jot2*
            be_integ.K2exp(jot2);
        } else {
          rat=exp(-jot1)*2.0*cosh(dj1*(f.nu+f.m)/temper)/jot1/jot1*
            be_integ.K2exp(jot1);
          rat/=exp(-jot2)*2.0*cosh(dj2*(f.nu+f.m)/temper)/jot2/jot2*
            be_integ.K2exp(jot2);
        }
      }
      
      // If the ratio between the last term and the first term is 
      // not small enough, return false
      if (o2scl::o2isfinite(rat) && rat>prec) {
        return false;
      }
      
      double first_term=0.0;
      f.pr=0.0;
      f.n=0.0;
      f.en=0.0;
      
      for(size_t j=1;j<=max_term;j++) {
        
        double pterm, nterm, enterm, edterm;
        
        ndeg_terms(j,tt,psi*tt,f.ms,f.inc_rest_mass,inc_antip,
                   pterm,nterm,enterm,edterm);
        
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
        
        // Stop if the last term is sufficiently small compared to
        // the first term
        if (j>1 && fabs(pterm)<prec*fabs(first_term)) {
          f.pr*=prefac;
          f.n*=prefac;
          f.en*=prefac;
          f.ed=-f.pr+f.nu*f.n+temper*f.en;
          return true;
        }
        
        // End of 'for(size_t j=1;j<=max_term;j++)'
      }

      */
      
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
