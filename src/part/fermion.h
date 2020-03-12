/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_H
#define O2SCL_FERMION_H

/** \file fermion.h
    \brief File defining \ref o2scl::fermion_tl
*/
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>

// For gsl_sf_fermi_dirac_int()
#include <gsl/gsl_specfunc.h>

#include <o2scl/constants.h>
#include <o2scl/funct.h>
#include <o2scl/root.h>
#include <o2scl/root_cern.h>
#include <o2scl/part.h>
#include <o2scl/polylog.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Fermion class

      This class adds two member data variables, \ref kf and \ref
      del, for the Fermi momentum and the gap, respectively.
  */
  template<class fp_t=double> class fermion_tl : public part_tl<fp_t> {

  public:

  /// Fermi momentum
  fp_t kf;
  /// Gap
  fp_t del;

  /// Create a fermion with mass \c mass and degeneracy \c dof.
  fermion_tl(fp_t mass=0, fp_t dof=0) : part_tl<fp_t>(mass,dof) {
    kf=0.0;
    del=0.0;
  }

  virtual ~fermion_tl() {
  }

  /// Return string denoting type ("fermion_tl")
  virtual const char *type() { return "fermion_tl"; }
    
  /// Copy constructor
  fermion_tl(const fermion_tl &f) {
    this->g=f.g;
    this->m=f.m;
    this->ms=f.ms;
    this->n=f.n;
    this->ed=f.ed;
    this->pr=f.pr;
    this->mu=f.mu;
    this->en=f.en;
    this->nu=f.nu;
    this->inc_rest_mass=f.inc_rest_mass;
    this->non_interacting=f.non_interacting;
    kf=f.kf;
    del=f.del;
  }

  /// Copy construction with operator=()
  fermion_tl &operator=(const fermion_tl &f) {
    if (this!=&f) {
      this->g=f.g;
      this->m=f.m;
      this->ms=f.ms;
      this->n=f.n;
      this->ed=f.ed;
      this->pr=f.pr;
      this->mu=f.mu;
      this->en=f.en;
      this->nu=f.nu;
      this->inc_rest_mass=f.inc_rest_mass;
      this->non_interacting=f.non_interacting;
      kf=f.kf;
      del=f.del;
    }
    return *this;
  }

  };

  typedef fermion_tl<double> fermion;

  /** \brief Fermion properties at zero temperature

      This is a base class for the computation of fermionic statistics
      at zero temperature. The more general case of finite temperature
      is taken care of by \ref fermion_thermo_tl objects. The
      primary functions are calc_mu_zerot() and calc_density_zerot()
      which compute all the thermodynamic quantities as a function of
      the chemical potential, or the density, respectively.
      
      \future Use hypot() and other more accurate functions for the
      analytic expressions for the zero temperature integrals. [Progress
      has been made, but there are probably other functions which may
      break down for small but finite masses and temperatures]
  */
  template<class fp_t=double> class fermion_zerot_tl {

  protected:

  /// Desc
  fp_t pi;
    
  /// Desc
  fp_t pi2;
    
  public:

  fermion_zerot_tl() {
    pi=boost::math::constants::pi<fp_t>();
    pi2=boost::math::constants::pi_sqr<fp_t>();
  };

  virtual ~fermion_zerot_tl() {
  }

  /// \name Zero-temperature fermions 
  //@{
  /** \brief Calculate the Fermi momentum from the density
	
      Uses the relation \f$ k_F = ( 6 \pi^2 n /g )^{1/3} \f$
  */
  void kf_from_density(fermion_tl<fp_t> &f) {
    f.kf=cbrt(6.0*pi2/f.g*f.n);
    return;
  }    

  /** \brief Energy density at T=0 from \ref o2scl::fermion_tl::kf and 
      \ref o2scl::part_tl::ms

      Calculates the integral 
      \f[
      \varepsilon = \frac{g}{2 \pi^2} \int_0^{k_F} k^2 
      \sqrt{k^2+m^{* 2}} d k
      \f]
  */
  void energy_density_zerot(fermion_tl<fp_t> &f) {
    fp_t r,efs;
    if (f.kf>0.0) {
      if (f.ms<=0.0) {
	f.ed=f.g*(pow(f.kf,4.0)/8.0/pi2);
      } else {
	efs=o2hypot(f.kf,f.ms);
	r=(f.kf+efs)/f.ms;
	f.ed=f.g/16.0/pi2*(2.0*f.kf*pow(efs,3.0)-f.kf*efs*f.ms*f.ms
			   -pow(f.ms,4.0)*log(r));
      }
    } else {
      f.ed=0.0;
    }
    return;
  }

  /** \brief Pressure at T=0 from \ref o2scl::fermion_tl::kf and 
      \ref o2scl::part_tl::ms

      Calculates the integral 
      \f[
      P=\frac{g}{6 \pi^2} \int_0^{k_F} \frac{k^4}{\sqrt{k^2+m^{* 2}}} d k
      \f]
  */
  void pressure_zerot(fermion_tl<fp_t> &f) {
    fp_t r,efs;
    if (f.kf>0.0) {
      if (f.ms<=0.0) {
	f.pr=f.g*(pow(f.kf,4.0)/24.0/pi2);
      } else {
	efs=o2hypot(f.kf,f.ms);
	r=(f.kf+efs)/f.ms;
	f.pr=f.g/48.0/pi2*(2.0*efs*pow(f.kf,3.0)-
			   3.0*f.kf*efs*f.ms*f.ms
			   +3.0*pow(f.ms,4.0)*log(r));
      }
    } else {
      f.pr=0.0;
    }
    return;
  }

  /** \brief Zero temperature fermions from \ref o2scl::part_tl::mu or 
      \ref o2scl::part_tl::nu and \ref o2scl::part_tl::ms
  */
  virtual void calc_mu_zerot(fermion_tl<fp_t> &f) {
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
      fp_t nupm=f.nu+f.m;
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
    
  /** \brief Zero temperature fermions from \ref o2scl::part_tl::n 
      and \ref o2scl::part_tl::ms
  */
  virtual void calc_density_zerot(fermion_tl<fp_t> &f) {
    if (f.non_interacting) { f.ms=f.m; }

    f.kf=cbrt(6.0*pi2/f.g*f.n);
    f.nu=o2hypot(f.kf,f.ms);
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
  //@}

  };

  /** \brief Double-precision version of \ref o2scl::fermion_zerot_tl 
   */
  typedef fermion_zerot_tl<double> fermion_zerot;
  
  /** \brief Fermion with finite-temperature thermodynamics
      [abstract base]

      This is an abstract base for the computation of
      finite-temperature fermionic statistics. Different children
      (e.g. \ref fermion_eff and \ref fermion_rel_tl) use different
      techniques to computing the momentum integrations.

      Because massless fermions at finite temperature are much
      simpler, there are separate member functions included in this
      class to handle them. The functions massless_calc_density() and
      massless_calc_mu() compute the thermodynamics of massless
      fermions at finite temperature given the density or the chemical
      potentials. The functions massless_pair_density() and
      massless_pair_mu() perform the same task, but automatically
      include antiparticles.

      The function massless_calc_density() uses a \ref root object to
      solve for the chemical potential as a function of the density.
      The default is an object of type root_cern. The function
      massless_pair_density() does not need to use the \ref root
      object because of the simplification afforded by the inclusion
      of antiparticles.

      \future Create a Chebyshev approximation for inverting the 
      the Fermi functions for massless_calc_density() functions?
  */
  template<class fd_inte_t=fermi_dirac_integ_gsl,
    class be_inte_t=bessel_K_exp_integ_gsl, class fp_t=double>
    class fermion_thermo_tl : public fermion_zerot_tl<fp_t> {
    
  public:
  
  fermion_thermo_tl() {
    massless_root=&def_massless_root;
  }
  
  virtual ~fermion_thermo_tl() {
  }
  
  /// Object for Fermi-Dirac integrals
  fd_inte_t fd_integ;
  
  /// Object for Bessel-exp integrals
  be_inte_t be_integ;
  
  /** \brief Calculate thermodynamic properties from the chemical
      potential using a nondegenerate expansion
  */
  template<class fermion_t>
  bool calc_mu_ndeg_tlate(fermion_t &f, fp_t temper, 
			  fp_t prec, bool inc_antip) {
  
    if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

    // Compute psi and tt
    fp_t psi, psi_num;
    if (f.inc_rest_mass) {
      psi_num=f.nu-f.ms;
    } else {
      psi_num=f.nu+f.m-f.ms;
    }
    psi=psi_num/temper;
    fp_t tt=temper/f.ms;

    // Return false immediately if we're degenerate
    if (inc_antip==false && psi>-1.0) return false;

    // Prefactor 'd' in Johns96
    fp_t prefac=f.g/2.0/this->pi2*pow(f.ms,4.0);

    // One term is always used, so only values of max_term greater than
    // 0 are useful.
    static const size_t max_term=200;
  
    // Maximum argument for exponential
    // fp_t log_dbl_max=709.78;

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
    fp_t rat;
    fp_t dj1=((fp_t)max_term), jot1=max_term/tt;
    fp_t dj2=1.0, jot2=1.0/tt;

    if (inc_antip==false) {
      rat=exp(dj1*psi)/jot1/jot1*be_integ.K2exp(jot1);
      rat/=exp(dj2*psi)/jot2/jot2*be_integ.K2exp(jot2);
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
    if (std::isfinite(rat) && rat>prec) {
      return false;
    }
  
    fp_t first_term=0.0;
    f.pr=0.0;
    f.n=0.0;
    f.en=0.0;

    for(size_t j=1;j<=max_term;j++) {

      fp_t pterm, nterm, enterm;
	
      ndeg_terms(j,tt,psi*tt,f.ms,f.inc_rest_mass,inc_antip,
		 pterm,nterm,enterm);

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

    // We failed to add enough terms, so return false
    return false;
  }
    
  /** \brief Calculate thermodynamic properties from the chemical
      potential using a degenerate expansion
  */
  template<class fermion_t>
  bool calc_mu_deg_tlate(fermion_t &f, fp_t temper, 
			 fp_t prec) {

    // Handle the zero-temperature limit
    if (temper==0.0) {
      this->calc_mu_zerot(f);
      return true;
    }

    // Double check to ensure T and mass are positive
    if (temper<0.0 || f.ms<0.0) {
      O2SCL_ERR2("Temperature or mass negative in fermion_thermo",
		 "::calc_mu_deg_tlate().",exc_einval);
    }
  
    if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
    // Compute psi and tt
    fp_t psi;
    if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
    else psi=(f.nu+f.m-f.ms)/temper;
    fp_t tt=temper/f.ms;

    // Return false immediately psi<0 where the expressions below
    // don't work because of the square roots
    if (psi<0.0) return false;
  
    // Prefactor 'd' in Johns96
    fp_t prefac=f.g/2.0/this->pi2*pow(f.ms,4.0);
  
    // Define x = psi * t = (mu/m - 1) and related values
    fp_t x=psi*tt;
    fp_t sx=sqrt(x);
    fp_t s2x=sqrt(2.0+x);
    fp_t x2=x*x;
    fp_t x3=x2*x;
    fp_t x4=x2*x2;
  
    // Evaluate the first and last term for the pressure
    fp_t pterm1;
    if (x>1.0e-5) {
      pterm1=(x*(1.0+x)*(2.0+x)*(-3.0+2.0*x*(2.0+x))+6.0*sx*s2x*
	      log((sx+s2x)/sqrt(2.0)))/24.0/sx/s2x;
    } else {
      pterm1=x2*sx*(29568.0+15840.0*x+1540.0*x2-105.0*x3)/55440.0/sqrt(2.0);
    }
    fp_t pterm4=-31.0*pow(this->pi*tt,6.0)/1008.0*(1.0+x)*
    sx*s2x/pow(x*(2.0+x),4.0);

    // Check if we're going to succeed
    if (fabs(pterm4)/fabs(pterm1)>prec) {
      return false;
    }
  
    // First order density term (first order entropy term is zero)
    fp_t nterm1=sx*s2x*x*(2.0+x)/3.0/f.ms;
  
    // Second order terms
    fp_t pterm2=tt*tt*this->pi2/6.0*(1.0+x)*sx*s2x;
    fp_t nterm2=tt*tt*this->pi2/6.0*(1.0+4.0*x+2.0*x2)/
    f.ms/sx/s2x;
    fp_t enterm2=tt*this->pi2/3.0*(1.0+x)*sx*s2x/f.ms;

    // Third order terms
    fp_t pterm3=7.0*pow(this->pi*tt,4.0)/360.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/pow(x*(2.0+x),1.5);
    fp_t nterm3=7.0*pow(this->pi*tt,4.0)/120.0/sx/s2x/
    x2/(x+2.0)/(x+2.0)/f.ms;
    fp_t enterm3=7.0*pow(this->pi*tt,4.0)/tt/90.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/f.ms/sx/s2x/x/(x+2.0);

    // Fourth order terms for density and entropy
    fp_t nterm4=31.0*pow(this->pi*tt,6.0)/1008.0*sx*s2x*
    (7.0+12.0*x+6.0*x2)/f.ms/pow(x*(2.0+x),5.0);
    fp_t enterm4=-31.0*pow(this->pi*tt,6.0)/tt/168.0*sx*s2x*
    (1.0+x)/pow(x*(2.0+x),4.0);

    // Add up all the terms
    f.pr=prefac*(pterm1+pterm2+pterm3+pterm4);
    f.n=prefac*(nterm1+nterm2+nterm3+nterm4);
    f.en=prefac*(enterm2+enterm3+enterm4);
    f.ed=-f.pr+f.nu*f.n+temper*f.en;

    return true;
  }

  /** \brief Non-degenerate expansion for fermions

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
      
      Defining \f$ \psi \equiv (\mu-m)/T \f$, \f$ t \equiv T/m \f$,
      and \f$ d \equiv g~m^4/(2 \pi^2) \f$ the pressure 
      in the non-degenerate limit (\f$ \psi \rightarrow - \infty \f$)
      is (\ref Johns96)
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
      The current code currently goes up to \f$ x^{-12} \f$ in the
      expansion, which is enough for the default precision of \f$
      10^{-18} \f$ since \f$ (20/700)^{12} \sim 10^{-19} \f$.
  */
  virtual bool calc_mu_ndeg(fermion &f, fp_t temper, 
			    fp_t prec=1.0e-18, bool inc_antip=false) {
    return calc_mu_ndeg_tlate<fermion>(f,temper,prec,inc_antip);
  }

  /** \brief Degenerate expansion for fermions

      Attempts to evaulate thermodynamics of a degenerate fermion.
      If the result is accurate to within the requested precision,
      this function returns <tt>true</tt>, and otherwise this
      function returns <tt>false</tt> and the values in stored in
      the <tt>pr</tt>, <tt>n</tt>, <tt>en</tt>, and <tt>ed</tt>
      field are meaningless.

      The pressure, density, and energy density, should be accurate
      to the requested precision, but the first term in the series
      expansion for the entropy is zero, so the entropy is one order
      lower in accuracy.

      \future Make a function like this for dndm, dsdT, etc. 
      for fermion_deriv .
  */
  virtual bool calc_mu_deg(fermion &f, fp_t temper, 
			   fp_t prec=1.0e-18) {
    return calc_mu_deg_tlate<fermion>(f,temper,prec);
  }      

  /** \brief Calculate properties as function of chemical potential
   */
  virtual void calc_mu(fermion &f, fp_t temper)=0;

  /** \brief Calculate properties as function of density

      \note This function returns an integer value, in contrast to
      \ref calc_mu(), because of the potential for non-convergence.
  */
  virtual int calc_density(fermion &f, fp_t temper)=0;

  /** \brief Calculate properties with antiparticles as function of
      chemical potential
  */
  virtual void pair_mu(fermion &f, fp_t temper)=0;

  /** \brief Calculate properties with antiparticles as function of
      density

      \note This function returns an integer value, in contrast to
      \ref pair_mu(), because of the potential for non-convergence.
  */
  virtual int pair_density(fermion &f, fp_t temper)=0;

  /// \name Massless fermions
  //@{
  /// Finite temperature massless fermions
  virtual void massless_calc_mu(fermion &f, fp_t temper) {
  
    fp_t fm2, fm3;

    if (f.non_interacting) { f.nu=f.mu; }

    fm2=this->fd_integ.calc_2(f.nu/temper)/2.0;
    fm3=this->fd_integ.calc_3(f.nu/temper)/6.0;
  
    f.n=f.g/this->pi2*pow(temper,3.0)*fm2;
    f.ed=f.g*3.0/this->pi2*pow(temper,4.0)*fm3;
    f.pr=f.ed/3.0;
    f.en=(f.ed+f.pr-f.n*f.nu)/temper;

    return;
  }    

  /// Finite temperature massless fermions
  virtual void massless_calc_density(fermion &f, fp_t temper) {
    fp_t x, T=temper;
  
    x=f.ms+temper;
    funct mf2=std::bind(std::mem_fn<fp_t(fp_t,fermion &,fp_t)>
			(&fermion_thermo_tl<fd_inte_t,be_inte_t,fp_t>::massless_solve_fun),
			this,std::placeholders::_1,std::ref(f),temper);
    massless_root->solve(x,mf2);
    f.nu=x;

    massless_calc_mu(f,temper);

    // If the particle is non-interacting, then need to set
    // mu=nu to get the entropy right
    if (f.non_interacting) { f.mu=f.nu; }

    return;
  }

  /** \brief Finite temperature massless fermions and antifermions 
   */
  virtual void massless_pair_mu(fermion &f, fp_t temper) {
    fp_t pitmu, pitmu2, nu2;

    if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
    if (f.nu==0.0) {
      f.n=0.0;
      f.ed=f.g/8.0/this->pi2*7.0/15.0*
      pow(this->pi*temper,4.0);
      f.pr=f.ed/3.0;
      f.en=(f.ed+f.pr-f.n*f.mu)/temper;
    } else {
      nu2=f.nu*f.nu;
      pitmu=this->pi*temper/f.nu;
      pitmu2=pitmu*pitmu;
      f.n=f.g*f.nu*nu2/6.0/this->pi2*(1.0+pitmu2);
      f.ed=f.g*nu2*nu2/8.0/this->pi2*(1.0+2.0*pitmu2+
				      7.0/15.0*pitmu2*pitmu2);
      f.pr=f.ed/3.0;
      f.en=(f.ed+f.pr-f.n*f.mu)/temper;
    
      // Might the following work better for the energy density?
      // pit=pi*temper;
      // pit2=pit*pit;
      // ed=g/8.0/pi2*(nu2*nu2+2.0*pit2*nu2+7.0/15.0*pit2*pit2);
    
    }

    return;
  }

  /** \brief Finite temperature massless fermions and antifermions 

      In the cases \f$ n^3 \gg T \f$ and \f$ T \gg n^3 \f$ ,
      expansions are used instead of the exact formulas to avoid
      loss of precision.

      In particular, using the parameter
      \f[
      \alpha = \frac{g^2 \pi^2 T^6}{243 n^2}
      \f]
      and defining the expression
      \f[
      \mathrm{cbt} = \alpha^{-1/6} \left( -1 + \sqrt{1+\alpha}\right)^{1/3}
      \f]
      we can write the chemical potential as
      \f[
      \mu = \frac{\pi T}{\sqrt{3}} \left(\frac{1}{\mathrm{cbt}} -
      \mathrm{cbt} \right)
      \f]
	
      These expressions, however, do not work well when \f$ \alpha
      \f$ is very large or very small, so series expansions are
      used whenever \f$ \alpha > 10^{4} \f$ or 
      \f$ \alpha < 3 \times 10^{-4} \f$. For small \f$ \alpha \f$, 
      \f[
      \left(\frac{1}{\mathrm{cbt}} -
      \mathrm{cbt} \right) \approx
      \frac{2^{1/3}}{\alpha^{1/6}} - 
      \frac{\alpha^{1/6}}{2^{1/3}} +
      \frac{\alpha^{5/6}}{6{\cdot}2^{2/3}} +
      \frac{\alpha^{7/6}}{12{\cdot}2^{1/3}} -
      \frac{\alpha^{11/6}}{18{\cdot}2^{2/3}} -
      \frac{5 \alpha^{13/6}}{144{\cdot}2^{1/3}} +
      \frac{77 \alpha^{17/6}}{2592{\cdot}2^{2/3}}
      \f]
      and for large \f$ \alpha \f$, 
      \f[
      \left(\frac{1}{\mathrm{cbt}} -
      \mathrm{cbt} \right) \approx
      \frac{2}{3} \sqrt{\frac{1}{\alpha}} - 
      \frac{8}{81} \left(\frac{1}{\alpha}\right)^{3/2} +
      \frac{32}{729} \left(\frac{1}{\alpha}\right)^{5/2}
      \f]

      This approach works to within about 1 \part in \f$ 10^{12} \f$,
      and is tested in <tt>fermion_ts.cpp</tt>.
	
      \future This could be improved by including more terms
      in the expansions.
  */
  virtual void massless_pair_density(fermion &f, fp_t temper) {

    fp_t t2=temper*temper,pitmu,pitmu2,nu2;
    fp_t cbt, alpha, two13, alpha16;

    if (f.non_interacting) { f.ms=f.m; }
    if (f.n<=0.0) {
      f.nu=0.0;
      f.ed=f.g/8.0/this->pi2*7.0/15.0*pow(this->pi*temper,4.0);
      f.pr=f.ed/3.0;
    } else {
      alpha=f.g*f.g*this->pi2*t2*t2*t2/243.0/f.n/f.n;
      if (alpha>1.0e4) {
	f.nu=(2.0/3.0/sqrt(alpha)-8.0/81.0/pow(alpha,1.5)+
	      32.0/729.0/pow(alpha,2.5))*this->pi*temper/sqrt(3.0);
      } else if (alpha<3.0e-4) {
	two13=cbrt(2.0);
	alpha16=pow(alpha,1.0/6.0);
	f.nu=(two13/alpha16-alpha16/two13+alpha/alpha16/6.0/two13/two13
	      +alpha*alpha16/12.0/two13-alpha*alpha/alpha16/18.0/two13/two13-
	      5.0*alpha*alpha*alpha16/144.0/two13+
	      77.0/2592.0*alpha*alpha*alpha/alpha16/two13/two13)*
	this->pi*temper/sqrt(3.0);
      } else {
	cbt=pow(-1.0+sqrt(1.0+alpha),1.0/3.0)/pow(alpha,1.0/6.0);
	f.nu=this->pi*temper/sqrt(3.0)*(1.0/cbt-cbt);
      }
      pitmu=this->pi*temper/f.nu;
      pitmu2=pitmu*pitmu;
      nu2=f.nu*f.nu;
      f.ed=f.g*nu2*nu2/8.0/this->pi2*
      (1.0+2.0*pitmu2+7.0/15.0*pitmu2*pitmu2);
      f.pr=f.ed/3.0;

      if (!std::isfinite(f.nu)) {
	std::string str="Chemical potential not finite ("+dtos(f.nu)+
	  ") in fermion::massless_pair_density().";
	O2SCL_ERR(str.c_str(),exc_efailed);
      }
    }

    if (f.non_interacting) { f.mu=f.nu; }
    f.en=(f.ed+f.pr-f.n*f.nu)/temper;

    return;
  }
  //@}
    
  /** \brief Set the solver for use in massless_calc_density() */ 
  void set_massless_root(root<> &rp) {
    massless_root=&rp;
    return;
  }

  /** \brief The default solver for massless_calc_density()
	
      We default to a solver of type root_cern here since we
      don't have a bracket or a derivative.
  */
  root_cern<> def_massless_root;

  /// Return string denoting type ("fermion_thermo")
  virtual const char *type() { return "fermion_thermo"; }

  /** \brief Compute a term in the nondegenerate expansion
   */
  void ndeg_terms(size_t j, fp_t tt,
		  fp_t xx, fp_t m, bool inc_rest_mass,
		  bool inc_antip, fp_t &pterm, fp_t &nterm,
		  fp_t &enterm) {
      
    fp_t dj=((fp_t)j);
    fp_t jot=dj/tt;

    if (inc_antip==false) {
      pterm=exp(jot*xx)/jot/jot*be_integ.K2exp(jot);
      if (j%2==0) pterm*=-1.0;
      nterm=pterm*jot/m;
      fp_t enterm1=(4.0*tt-dj*xx-dj)/dj/tt*nterm;
      fp_t enterm2=exp(jot*xx)/dj*be_integ.K1exp(jot)/m;
      if (j%2==0) {
	enterm=enterm1-enterm2;
      } else {
	enterm=enterm1+enterm2;
      }
    } else {
      pterm=exp(-jot)*2.0*cosh(jot*(xx+1.0)/tt)/jot/jot*
      be_integ.K2exp(jot);
      if (j%2==0) pterm*=-1.0;
      nterm=pterm*tanh(jot*(xx+1.0))*jot;
      fp_t enterm1=-(1.0+xx)/tt*nterm/m;
      fp_t enterm2=2.0*exp(-jot*xx)/dj*cosh(jot*(xx+1.0))*
      be_integ.K3exp(jot)/m;
      if (j%2==0) {
	enterm=enterm1-enterm2;
      } else {
	enterm=enterm1+enterm2;
      }
    }
		    
    return;
  }
    
#ifndef DOXYGEN_NO_O2NS

  protected:
    
  /// A pointer to the solver for massless fermions
  root<> *massless_root;

  /// Solve for the chemical potential for massless fermions
  fp_t massless_solve_fun(fp_t x, fermion &f, fp_t temper) {
    fp_t fm2=this->fd_integ.calc_2(x/temper)/2.0;
    return f.g*pow(temper,3.0)*fm2/this->pi2/f.n-1.0;
  }    
  
#endif
    
  };

  /** \brief Double-precision version of \ref o2scl::fermion_thermo_tl 
   */
  typedef fermion_thermo_tl<> fermion_thermo;
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
