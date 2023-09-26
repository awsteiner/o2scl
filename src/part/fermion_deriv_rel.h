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
#ifndef O2SCL_FERMION_DERIV_REL_H
#define O2SCL_FERMION_DERIV_REL_H

/** \file fermion_deriv_rel.h
    \brief File defining \ref o2scl::fermion_deriv_rel_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/root_cern.h>
#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/part_deriv.h>
#include <o2scl/fermion_rel.h>

namespace o2scl {

  /** \brief Desc 
   */
  template<class fp_t> class fermion_deriv_rel_integ {

  public:

    fermion_deriv_rel_integ() {
      method=automatic;
      intl_method=by_parts;
    }
    
    /// The internal integration method
    int intl_method;
    
    /** \name Method of computing derivatives
     */
    //@{
    /// Method (default is \ref automatic)
    int method;
    /// Automatically choose method
    static const int automatic=0;
    /// In the form containing \f$ f(1-f) \f$ .
    static const int direct=1;
    /// Integrate by parts
    static const int by_parts=2;
    //@}

    /** \name The integrands, as a function of \f$ u=k/T \f$, for 
	non-degenerate integrals
    */
    //@{
    /** \brief Integrand for derivative of density with respect to
        temperature for non-degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t density_T_fun(internal_fp_t u, fp_t m2, fp_t ms2,
                                fp_t nu2, fp_t T2, bool inc_rest_mass) {
      
      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      
      internal_fp_t k=u*T, E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=k*k*(E-nu)/T*ff*(1-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(nu)/T-k*k*(nu)/T/E)*
	    T*fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function(x2);
	  ret=k*k*(E-nu)/T*ff*(1-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(nu+m)/T-k*k*(nu+m)/T/E)*
	    T*fermi_function(x2);
	}
      }
      return ret;
    }

    /** \brief Integrand for derivative of density with respect to
        chemical potential for non-degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t density_mu_fun(internal_fp_t u, fp_t m2, fp_t ms2,
                                 fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t k=u*T, E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=k*k*ff*(1-ff);
	} else {
	  ret=T*(E*E+k*k)/E*fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x2);
	  ret=k*k*ff*(1-ff);
	} else {
	  ret=T*(E*E+k*k)/E*fermi_function(x2);
	}
      }
      return ret;
    }
    
    /** \brief Integrand for derivative of entropy with respect to
        temperature for non-degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t entropy_T_fun(internal_fp_t u, fp_t m2, fp_t ms2,
                                fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t k=u*T, E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=T*k*k*ff*(1-ff)*pow(E-nu,2.0)/pow(T,3);
	} else {
	  ret=(E-nu)/E/T*
	    (pow(E,3)+3*E*k*k-(E*E+k*k)*(nu))*
	    fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function(x2);
	  ret=T*k*k*ff*(1-ff)*pow(E-nu,2.0)/pow(T,3);
	} else {
	  ret=(E-m-nu)/E/T*
	    (pow(E,3)+3*E*k*k-(E*E+k*k)*(nu+m))*
	    fermi_function(x2);
	}
      }
      return ret;
    }

    /** \brief Integrand for derivative of density with respect to
        effective mass for non-degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t density_ms_fun(internal_fp_t u, fp_t m2, fp_t ms2,
                                 fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t k=u*T, E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function((E-nu)/T);
	  ret=-k*k*ms/(E)/T*ff*(1-ff);
	} else {
	  ret=-ms*fermi_function((E-nu)/T);
	}
      } else {
	E=hypot(k,ms);
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function((E-nu)/T);
	  ret=-k*k*ms/(E+m)/T*ff*(1-ff);
	} else {
	  ret=-ms*fermi_function((E-m-nu)/T);
	}
      }
      ret*=T;
      return ret;
    }
    //@}

    /** \name The integrands, as a function of momentum, for the 
	degenerate integrals
    */
    //@{
    /** \brief Integrand for derivative of density with respect to
        temperature for degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t deg_density_T_fun(internal_fp_t k, fp_t m2, fp_t ms2,
                                fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=k*k*(E-nu)/T/T*ff*(1-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(nu)/T-k*k*(nu)/T/E)*
	    fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function(x2);
	  ret=k*k*(E-nu)/T/T*ff*(1-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(nu+m)/T-k*k*(nu+m)/T/E)*
	    fermi_function(x2);
	}
      }
      return ret;
    }

    /** \brief Integrand for derivative of density with respect to
        chemical potential for degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t deg_density_mu_fun(internal_fp_t k, fp_t m2, fp_t ms2,
                                fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=k*k/T*ff*(1-ff);
	} else {
	  ret=(E*E+k*k)/E*fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function(x2);
	  ret=k*k/T*ff*(1-ff);
	} else {
	  ret=(E*E+k*k)/E*fermi_function(x2);
	}
      }
      return ret;
    }

    /** \brief Integrand for derivative of entropy with respect to
        temperature for degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t deg_entropy_T_fun(internal_fp_t k, fp_t m2, fp_t ms2,
                                fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);

      internal_fp_t E, ret;
      E=hypot(k,ms);
      if (inc_rest_mass) {
        internal_fp_t x=(E-nu)/T;
	internal_fp_t ff=fermi_function(x);
	if (intl_method==direct) {
	  ret=k*k*ff*(1-ff)*pow(E-nu,2.0)/pow(T,3);
	} else {
	  ret=(E-nu)/E/T/T*
	    (pow(E,3)+3*E*k*k-(E*E+k*k)*nu)*ff;
	}
      } else {
        internal_fp_t x2=(E-m-nu)/T;
	internal_fp_t ff=fermi_function(x2);
	if (intl_method==direct) {
	  ret=k*k*ff*(1-ff)*pow(E-nu-m,2.0)/pow(T,3);
	} else {
	  ret=(E-m-nu)/E/T/T*
	    (pow(E,3)+3*E*k*k-(E*E+k*k)*(nu+m))*ff;
	}
      }
      return ret;
    }

    /** \brief Integrand for derivative of density with respect to
        effective mass for degenerate particles
     */
    template<class internal_fp_t>
    internal_fp_t deg_density_ms_fun(internal_fp_t k, fp_t m2, fp_t ms2,
                                     fp_t nu2, fp_t T2, bool inc_rest_mass) {

      internal_fp_t m=static_cast<internal_fp_t>(m2);
      internal_fp_t ms=static_cast<internal_fp_t>(ms2);
      internal_fp_t nu=static_cast<internal_fp_t>(nu2);
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      
      internal_fp_t E, ret;
      if (inc_rest_mass) {
	E=hypot(k,ms);
        internal_fp_t x=(E-nu)/T;
	if (intl_method==direct) {
	  internal_fp_t ff=fermi_function(x);
	  ret=-k*k*ms/E/T*ff*(1-ff);
	} else {
	  ret=-ms*fermi_function(x);
	}
      } else {
	E=hypot(k,ms);
        internal_fp_t x2=(E-m-nu)/T;
	if (intl_method==direct) {
	  E-=m;
	  internal_fp_t ff=fermi_function(x2);
	  ret=-k*k*ms/(E+m)/T*ff*(1-ff);
	} else {
	  ret=-ms*fermi_function(x2);
	}
      }
      return ret;
    }
    
  };

  /** \brief Equation of state for a relativistic fermion

      \note This class only has preliminary support for
      inc_rest_mass=true (more testing should be done, particularly
      for the "pair" functions)

      This implements an equation of state for a relativistic fermion
      using direct integration. After subtracting the rest mass from
      the chemical potentials, the distribution function is
      \f[
      \left\{1+\exp[(\sqrt{k^2+m^{* 2}}-m-\nu)/T]\right\}^{-1}
      \f]
      where \f$ k \f$ is the momentum, \f$ \nu \f$ is the effective
      chemical potential, \f$ m \f$ is the rest mass, and \f$ m^{*}
      \f$ is the effective mass.  For later use, we define \f$ E^{*} =
      \sqrt{k^2 + m^{*2}} \f$ . The degeneracy parameter is
      \f[
      \psi=(\nu+(m-m^{*}))/T 
      \f] 
      For \f$ \psi \f$ greater than \ref deg_limit (degenerate
      regime), a finite interval integrator is used and for \f$ \psi
      \f$ less than \ref deg_limit (non-degenerate regime), an
      integrator over the interval from \f$ [0,\infty) \f$ is
      used. The upper limit on the degenerate integration is given by
      the value of the momentum \f$ k \f$ which is the solution of
      \f[
      (\sqrt{k^2+m^{*,2}}-m-\nu)/T=\mathrm{f{l}imit}
      \f]
      which is
      \f[
      \sqrt{(m+{\cal L})^2-m^{*2}}
      \f]
      where \f$ {\cal L}\equiv\mathrm{f{l}imit}\times T+\nu \f$ .

      For the entropy integration, we set the lower limit
      to
      \f[
      2 \sqrt{\nu^2+2 \nu m} - \mathrm{upper~limit}
      \f]
      since the only contribution to the entropy is at the Fermi surface.
      \comment
      I'm not sure, but I think this is an expression determined
      from a small T taylor expansion of the argument of the 
      exponential.
      \endcomment

      In the non-degenerate regime, we make the substitution \f$ u=k/T
      \f$ to help ensure that the variable of integration scales
      properly.

      Uncertainties are given in \ref unc.

      \b Evaluation \b of \b the \b derivatives

      The relevant
      derivatives of the distribution function are
      \f[
      \frac{\partial f}{\partial T}=
      f(1-f)\frac{E^{*}-m-\nu}{T^2}
      \f]
      \f[
      \frac{\partial f}{\partial \nu}=
      f(1-f)\frac{1}{T}
      \f]
      \f[
      \frac{\partial f}{\partial k}=
      -f(1-f)\frac{k}{E^{*} T}
      \f]
      \f[
      \frac{\partial f}{\partial m^{*}}=
      -f(1-f)\frac{m^{*}}{E^{*} T}
      \f]

      We also need the derivative of the entropy integrand w.r.t. the 
      distribution function, which is
      \f[
      {\cal S}\equiv f \ln f +(1-f) \ln (1-f) \qquad
      \frac{\partial {\cal S}}{\partial f} = \ln 
      \left(\frac{f}{1-f}\right) = 
      \left(\frac{\nu-E^{*}+m}{T}\right)
      \f]
      where the entropy density is
      \f[
      s = - \frac{g}{2 \pi^2} \int_0^{\infty} {\cal S} k^2 d k
      \f]

      The derivatives can be integrated directly (\ref method = \ref
      direct) or they may be converted to integrals over the
      distribution function through an integration by parts (\ref
      method = \ref by_parts)
      \f[
      \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) 
      g(k)\right|_{k=a}^{k=b}
      - \int_a^b g(k) \frac{d f(k)}{dk} dk 
      \f]
      using the distribution function for \f$ f(k) \f$ and 0 and 
      \f$ \infty \f$ as the limits, we have
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{d g(k)}{dk} f dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} g(k) f (1-f) \frac{k}{E^{*} T} dk 
      \f]
      as long as \f$ g(k) \f$ vanishes at \f$ k=0 \f$ .
      Rewriting,
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} h(k) f (1-f) dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} f \frac{T}{k} 
      \left[ h^{\prime} E^{*}-\frac{h E^{*}}{k}+\frac{h k}{E^{*}} \right] dk
      \f]
      as long as \f$ h(k)/k \f$ vanishes at \f$ k=0 \f$ .

      \b Explicit \b forms

      1) The derivative of the density wrt the chemical potential
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2}{T} f (1-f) dk
      \f]
      Using \f$ h(k)=k^2/T \f$ we get
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      \left(\frac{k^2+E^{*2}}{E^{*}}\right) f dk
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-m-\nu)}{T^2} 
      f (1-f) dk
      \f]
      Using \f$ h(k)=k^2(E^{*}-\nu)/T^2 \f$ we get
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
      \left[2 k^2+E^{*2}-E^{*}\left(\nu+m\right)-
      k^2 \left(\frac{\nu+m}{E^{*}}\right)\right] dk
      \f]

      3) The derivative of the entropy wrt the chemical potential
      \f[
      \left(\frac{d s}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-m-\nu)}{T^2} dk
      \f]
      This verifies the Maxwell relation
      \f[
      \left(\frac{d s}{d \mu}\right)_T =
      \left(\frac{d n}{d T}\right)_{\mu}
      \f]

      4) The derivative of the entropy wrt the temperature
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-m-\nu)^2}{T^3} dk
      \f]
      Using \f$ h(k)=k^2 (E^{*}-\nu)^2/T^3 \f$ 
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f(E^{*}-m-\nu)}{E^{*}T^2} 
      \left[E^{* 3}+3 E^{*} k^2- (E^{* 2}+k^2)(\nu+m)\right] d k
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      -\frac{g}{2 \pi^2} \int_0^{\infty} 
      \frac{k^2 m^{*}}{E^{*} T} f (1-f) dk
      \f]
      Using \f$ h(k)=-(k^2 m^{*})/(E^{*} T) \f$ we get
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      -\frac{g}{2 \pi^2} \int_0^{\infty} 
      m^{*} f dk
      \f]
      \comment
      This derivative may be written in terms of the 
      others
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = \frac{3 n}{m^{*}}
      - \frac{T}{m^{*}}\left[ \left(\frac{d n}{d T}\right)_{\mu}
      +\frac{\mu}{T} \left(\frac{d n}{d \mu}\right)_{T}
      \right] - \left(\frac{d n}{d \mu}\right)_{T}
      \f]
      \endcomment

      \note The dsdT integration may fail if the system is
      very degenerate. When method is byparts, the integral involves a
      large cancellation between the regions from \f$ k \in (0,
      \mathrm{ulimit/2}) \f$ and \f$ k \in (\mathrm{ulimit/2},
      \mathrm{ulimit}) \f$. Switching to method=direct and setting the
      lower limit to \f$ \mathrm{llimit} \f$, may help, but recent
      testing on this gave negative values for dsdT. For very
      degenerate systems, an expansion may be better than trying
      to perform the integration. The value of the integrand
      at k=0 also looks like it might be causing difficulties.
      
      \verbatim embed:rst

      .. todo::

         In class fermion_deriv_rel_tl:

         - Future: The option err_nonconv=false is not really
           implemented yet.

         - Future: The \ref pair_density() function is a bit slow
           because it computes the non-derivative thermodynamic
           quantities twice, and this could be improved.

      \endverbatim
  */
  template<class fermion_deriv_t=fermion_deriv_tl<double>,
           class fermion_rel_t=fermion_rel_tl<fermion_deriv_t>,
	   class fp_t=double>
  class fermion_deriv_rel_tl : public fermion_deriv_thermo_tl<fp_t>,
                               public fermion_deriv_rel_integ<fp_t> {
    
  public:

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_deriv_rel_tl() {
  
      deg_limit=2.0;
      upper_limit_fac=20;

      nit=&def_nit;
      dit=&def_dit;
  
      this->method=this->automatic;
      this->intl_method=this->by_parts;

      err_nonconv=true;

      last_method=0;
      last_method_s="";
      verify_ti=false;

      verbose=0;
      multip=false;
      tol_expan=1e-14;
    }
  
    virtual ~fermion_deriv_rel_tl() {
    }
    
    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2.0)
    */
    fp_t deg_limit;
    
    /** \brief The limit for the Fermi functions (default 20)
	
	fermion_deriv_rel will ignore corrections smaller than about
	\f$ \exp(-\mathrm{f{l}imit}) \f$ . 
    */
    fp_t upper_limit_fac;
    
    /// Storage for the most recently calculated uncertainties 
    fermion_deriv_t unc;

    /// Object for computing non-derivative quantities
    fermion_rel_t fr;

    /// Verbosity parameter (default 0)
    int verbose;

    /// If true, use multiprecision to improve the integrations
    bool multip;
    
    /// If true, verify the thermodynamic identity (default false)
    bool verify_ti;
    
    /// Tolerance for expansions (default \f$ 10^{-14} \f$)
    fp_t tol_expan;
    
    /** \brief An integer indicating the last numerical method used

	The function \ref calc_mu() sets this integer to a two-digit
	number. It is equal to 10 times the value reported by \ref
	o2scl::fermion_rel::calc_mu() plus a value from the list
	below corresponding to the method used for the derivatives
	- 1: nondegenerate expansion
	- 2: degenerate expansion
	- 3: nondegenerate integrand, using \ref by_parts for
	\ref method
	- 4: nondegenerate integrand, using user-specified value
	for \ref method
	- 5: degenerate integrand, using \ref direct
	- 6: degenerate integrand, using \ref by_parts
	- 7: degenerate integrand, using user-specified 
	value for \ref method

	The function \ref nu_from_n() sets this value equal to
	100 times the value reported by 
	\ref o2scl::fermion_rel_tl::nu_from_n() .

	The function \ref calc_density() sets this value equal to the
	value from \ref o2scl::fermion_deriv_rel_tl::nu_from_n() plus the
	value from \ref o2scl::fermion_deriv_rel_tl::calc_mu() .

    */
    int last_method;

    /// String detailing last method used
    std::string last_method_s;
    
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief Calculate properties as function of chemical potential
     */
    virtual int calc_mu(fermion_deriv_t &f, fp_t temper) {

      last_method=0;
      last_method_s="";
      
      fr.calc_mu(f,temper);
      last_method=fr.last_method*10;
      last_method_s=((std::string)"base: ")+fr.last_method_s+
	" deriv: ";
  
      int iret;

      if (temper<=0) {
	O2SCL_ERR("T=0 not implemented in fermion_deriv_rel().",
		  exc_eunimpl);
      }

      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      fp_t prefac=f.g/2/this->pi2;

      // Compute the degeneracy parameter

      bool deg=false;
      fp_t psi;
      if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
      else psi=(f.nu+f.m-f.ms)/temper;
      if (psi>deg_limit) deg=true;

      // Try the non-degenerate expansion if psi is small enough
      if (psi<-4) {
	bool acc=this->calc_mu_ndeg(f,temper,tol_expan);
	if (acc) {
          if (verbose>1) {
            std::cout << "fermion_deriv_rel::calc_mu() using nondegenerate "
                      << "expansion." << std::endl;
          }
	  unc.n=f.n*tol_expan;
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  unc.dndT=f.dndT*tol_expan;
	  unc.dsdT=f.dsdT*tol_expan;
	  unc.dndmu=f.dndmu*tol_expan;
	  last_method+=1;
          if (last_method_s.length()>200) {
            O2SCL_ERR("Last method problem (1)",o2scl::exc_esanity);
          } else {
            last_method_s+="nondeg. exp.";
          }
	  return 0;
	}
      }

      // Try the degenerate expansion if psi is large enough
      if (psi>20) {
	bool acc=this->calc_mu_deg(f,temper,tol_expan);
	if (acc) {
          if (verbose>1) {
            std::cout << "fermion_deriv_rel::calc_mu() using degenerate "
                      << "expansion." << std::endl;
          }
	  unc.n=f.n*tol_expan;
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  unc.dndT=f.dndT*tol_expan;
	  unc.dsdT=f.dsdT*tol_expan;
	  unc.dndmu=f.dndmu*tol_expan;
	  last_method+=2;
          if (last_method_s.length()>200) {
            O2SCL_ERR("Last method problem (2)",o2scl::exc_esanity);
          } else {
            last_method_s+="deg. exp.";
          }
	  return 0;
	}
      }

      if (deg==false) {
    
	// Set integration method
	if (this->method==this->automatic) {
	  this->intl_method=this->by_parts;
	  last_method+=3;
          if (last_method_s.length()>200) {
            O2SCL_ERR("Last method problem (3)",o2scl::exc_esanity);
          } else {
            last_method_s+="nondeg. integ. (automatic)";
          }
	} else {
	  this->intl_method=this->method;
	  last_method+=4;
          if (last_method_s.length()>200) {
            O2SCL_ERR("Last method problem (4)",o2scl::exc_esanity);
          } else {
            if (this->method==this->by_parts) {
              last_method_s+="nondeg. integ. (by parts)";
            } else {
              last_method_s+="nondeg. integ. (direct)";
            }
          }
	}
        
        if (multip==true) {
          
          fp_t zero=0, tol_rel=0;
          int ix=it_multip.integ_iu_err_multip
            ([this,f,temper](auto &&k) mutable {
              return this->density_T_fun(k,f.m,f.ms,f.nu,temper,
                                         f.inc_rest_mass); },
              zero,f.dndT,unc.dndT,tol_rel);
          if (ix!=0) {
            O2SCL_ERR2("dndT integration (ndeg, multip) failed in ",
                       "fermion_deriv_rel::calc_mu().",
                       exc_efailed);
          }
          
        } else {
        
          // The non-degenerate case
          
          funct density_T_fun_f=[this,f,temper](double k) -> double
          { return this->density_T_fun(k,f.m,f.ms,f.nu,temper,
                                        f.inc_rest_mass); };
          
          iret=nit->integ_iu_err(density_T_fun_f,0,f.dndT,unc.dndT);
          if (iret!=0) {
            O2SCL_ERR2("dndT integration (ndeg) failed in ",
                       "fermion_deriv_rel::calc_mu().",
                       exc_efailed);
          }

        }
        f.dndT*=prefac;
        unc.dndT*=prefac;

        if (multip==true) {
          
          fp_t zero=0, tol_rel=0;
          int ix=it_multip.integ_iu_err_multip
            ([this,f,temper](auto &&k) mutable {
              return this->density_mu_fun(k,f.m,f.ms,f.nu,temper,
                                           f.inc_rest_mass); },
              zero,f.dndmu,unc.dndmu,tol_rel);
          if (ix!=0) {
            O2SCL_ERR2("dndmu integration (ndeg, multip) failed in ",
                       "fermion_deriv_rel::calc_mu().",
                       exc_efailed);
          }
          
        } else {
          
          funct density_mu_fun_f=[this,f,temper](double k) -> double
          { return this->density_mu_fun(k,f.m,f.ms,f.nu,temper,
                                         f.inc_rest_mass); };

          iret=nit->integ_iu_err(density_mu_fun_f,0,f.dndmu,unc.dndmu);
          if (iret!=0) {
            O2SCL_ERR2("dndmu integration (ndeg) failed in ",
                       "fermion_deriv_rel::calc_mu().",
                       exc_efailed);
          }

        }
        
	f.dndmu*=prefac;
	unc.dndmu*=prefac;
    
        if (multip==true) {
          
          fp_t zero=0, tol_rel=0;
          int ix=it_multip.integ_iu_err_multip
            ([this,f,temper](auto &&k) mutable {
              return this->entropy_T_fun(k,f.m,f.ms,f.nu,temper,
                                          f.inc_rest_mass); },
              zero,f.dsdT,unc.dsdT,tol_rel);
          if (ix!=0) {
            O2SCL_ERR2("dsdT integration (ndeg, multip) failed in ",
                       "fermion_deriv_rel::calc_mu().",
                       exc_efailed);
          }
          
        } else {
          
          funct entropy_T_fun_f=[this,f,temper](double k) -> double
          { return this->entropy_T_fun(k,f.m,f.ms,f.nu,temper,
                                        f.inc_rest_mass); };
          
          iret=nit->integ_iu_err(entropy_T_fun_f,0,f.dsdT,unc.dsdT);
          if (iret!=0) {
            O2SCL_ERR2("dsdT integration (ndeg) failed in ",
                       "fermion_deriv_rel_tl<fp_t>::calc_mu().",exc_efailed);
          }

        }
        
	f.dsdT*=prefac;
	unc.dsdT*=prefac;

      } else {

	// Compute the upper limit for degenerate integrals

	fp_t arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	}
	// Set to zero to avoid uninit'ed var. warnings
	fp_t ul=0;
	if (arg>0) {
	  ul=sqrt(arg);
	} else {
	  O2SCL_ERR2("Zero density in degenerate limit in fermion_deriv_rel::",
		     "calc_mu(). Variable deg_limit set improperly?",
		     exc_efailed);
	}
    
	// Compute the lower limit for the entropy and derivative
	// integrations
    
	fp_t ll;
	if (f.inc_rest_mass) {
	  arg=pow(-upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	  if (arg>0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1;
	  }
	} else {
	  arg=pow(-upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	  if (arg>0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1;
	  }
	}

	// Set integration method
	if (this->method==this->automatic) {
	  if ((!f.inc_rest_mass && (f.nu+f.m-f.ms)/temper>1e3) ||
	      (f.inc_rest_mass && (f.nu-f.ms)/temper>1e3)) {
	    this->intl_method=this->direct;
	    last_method+=5;
            if (last_method_s.length()>200) {
              O2SCL_ERR("Last method problem (5)",o2scl::exc_esanity);
            } else {
              last_method_s+="deg. integ. (automatic: direct)";
            }
	  } else {
	    this->intl_method=this->by_parts;
	    last_method+=6;
            if (last_method_s.length()>200) {
              O2SCL_ERR("Last method problem (6)",o2scl::exc_esanity);
            } else {
              last_method_s+="deg. integ. (by_parts)";
            }
	  }
	} else {
	  this->intl_method=this->method;
	  last_method+=7;
	  if (last_method_s.length()>200) {
            O2SCL_ERR("Last method problem (7)",o2scl::exc_esanity);
	  } else {
            if (this->method==this->by_parts) {
	      last_method_s+="deg. integ. (by_parts)";
            } else {
	      last_method_s+="deg. integ. (direct)";
            }
	  }
	}
        
        if (multip==true) {
          
          fp_t tol_rel=0;

          if (this->intl_method==this->direct && ll>0) {

            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_density_mu_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                ll,ul,f.dndmu,unc.dndmu,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dndmu integration (deg, multip, llf) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }

          } else {
            
            fp_t zero=0;
            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_density_mu_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                zero,ul,f.dndmu,unc.dndmu,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dndmu integration (deg, multip, ll0) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }

          }
          
        } else {
	  
          funct deg_density_mu_fun_f=[this,f,temper](double k) -> double
          { return this->deg_density_mu_fun(k,f.m,f.ms,f.nu,temper,
                                             f.inc_rest_mass); };
	  
          if (this->intl_method==this->direct && ll>0) {
            iret=dit->integ_err(deg_density_mu_fun_f,ll,ul,
                                f.dndmu,unc.dndmu);
          } else {
            iret=dit->integ_err(deg_density_mu_fun_f,0,ul,
                                f.dndmu,unc.dndmu);
          }
          if (iret!=0) {
            O2SCL_ERR2("dndmu integration (deg) failed in fermion_",
                       "deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
                       exc_efailed);
          }

        }
        
	f.dndmu*=prefac;
	unc.dndmu*=prefac;
    
        if (multip==true) {
          
          fp_t tol_rel=0;

          if (this->intl_method==this->direct && ll>0) {

            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_density_T_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                ll,ul,f.dndT,unc.dndT,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dndT integration (deg, multip, llf) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }

          } else {
            
            fp_t zero=0;
            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_density_T_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                zero,ul,f.dndT,unc.dndT,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dndT integration (deg, multip, ll0) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }

          }

        } else {
          
          funct deg_density_T_fun_f=[this,f,temper](double k) -> double
          { return this->deg_density_T_fun(k,f.m,f.ms,f.nu,temper,
                                            f.inc_rest_mass); };
	  
          if (this->intl_method==this->direct && ll>0) {
            iret=dit->integ_err(deg_density_T_fun_f,ll,ul,f.dndT,
				unc.dndT);
          } else {
            iret=dit->integ_err(deg_density_T_fun_f,0,ul,f.dndT,
				unc.dndT);
          }
          if (iret!=0) {
            O2SCL_ERR2("dndT integration (deg) failed in fermion_",
                       "deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
                       exc_efailed);
          }

        }
        
	f.dndT*=prefac;
	unc.dndT*=prefac;

        if (multip==true) {
          
          fp_t tol_rel=0;

          if (this->intl_method==this->direct && ll>0) {

            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_entropy_T_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                ll,ul,f.dndT,unc.dsdT,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dsdT integration (deg, multip, llf) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }

          } else {
            
            fp_t zero=0;
            int ix=it_multip.integ_err_multip
              ([this,f,temper](auto &&k) mutable {
                return this->deg_entropy_T_fun(k,f.m,f.ms,f.nu,temper,
                                                 f.inc_rest_mass); },
                zero,ul,f.dsdT,unc.dsdT,tol_rel);
            if (ix!=0) {
              O2SCL_ERR2("dsdT integration (deg, multip, ll0) failed in ",
                         "fermion_deriv_rel::calc_mu().",
                         exc_efailed);
            }


          }

        } else {

          funct deg_entropy_T_fun_f=[this,f,temper](double k) -> double
          { return this->deg_entropy_T_fun(k,f.m,f.ms,f.nu,temper,
                                            f.inc_rest_mass); };

          if (this->intl_method==this->direct && ll>0) {
            iret=dit->integ_err(deg_entropy_T_fun_f,ll,ul,f.dsdT,
				unc.dsdT);
          } else {
            iret=dit->integ_err(deg_entropy_T_fun_f,0,ul,f.dsdT,
				unc.dsdT);
          }
          if (iret!=0) {
            O2SCL_ERR2("dsdT integration (deg) failed in fermion_",
                       "deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
                       exc_efailed);
          }
          
        }
        
	f.dsdT*=prefac;
	unc.dsdT*=prefac;

      }
  
      if (!std::isfinite(f.en)) {
	O2SCL_ERR2("Entropy not finite in ",
		   "fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
		   exc_efailed);
      }
      f.pr=-f.ed+temper*f.en+f.nu*f.n;
      
      // Pressure uncertainties are not computed
      unc.pr=0;

      return 0;
    }
  
    /** \brief Calculate properties as function of density
     */
    virtual int calc_density(fermion_deriv_t &f, fp_t temper) {
  
      if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
      nu_from_n(f,temper);
      int lm=last_method;
      std::string stmp=last_method_s;
  
      if (f.non_interacting) { f.mu=f.nu; }

      calc_mu(f,temper);
      last_method+=lm;
      if (last_method_s.length()>200) {
        O2SCL_ERR("Last method problem (8)",o2scl::exc_esanity);
      } else {
        last_method_s=((std::string)"deriv nu_from_n: ")+stmp+
	  " deriv calc_mu: "+last_method_s;
      }
      
      return 0;
    }

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual int pair_mu(fermion_deriv_t &f, fp_t temper) {
      if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
      fermion_deriv antip(f.ms,f.g);
      f.anti(antip);

      calc_mu(f,temper);
      int lm=last_method*100;
      std::string stmp=last_method_s;
      
      calc_mu(antip,temper);
      last_method+=lm;
      if (last_method_s.length()>200) {
        std::cout << "9: " << last_method_s << std::endl;
        O2SCL_ERR("Last method problem (9)",o2scl::exc_esanity);
      } else {
        last_method_s="part. "+stmp+" : antipart. "+last_method_s;
      }

      f.n-=antip.n;
      f.pr+=antip.pr;
      if (f.inc_rest_mass) {
	f.ed+=antip.ed;
      } else {
	f.ed=f.ed+antip.ed+2*antip.n*f.m;
      }
      f.en+=antip.en;
      f.dsdT+=antip.dsdT;
      f.dndT-=antip.dndT;
      f.dndmu+=antip.dndmu;
  
      return 0;
    }

    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual int pair_density(fermion_deriv_t &f, fp_t temper) {
      int ret=fr.pair_density(f,temper);
      pair_mu(f,temper);
      return 0;
    }

    /// Calculate effective chemical potential from density
    virtual int nu_from_n(fermion_deriv_t &f, fp_t temper) {
      int ret=fr.nu_from_n(f,temper);
      last_method=fr.last_method*100;
      if (last_method_s.length()>200) {
        O2SCL_ERR("Last method problem (10)",o2scl::exc_esanity);
      } else {
        last_method_s=((std::string)"base: ")+fr.last_method_s;
      }
      return ret;
    }

    /// \name Integration objects
    //@{
    /** \brief Set inte objects
	
	The first integrator is used for non-degenerate integration
	and should integrate from 0 to \f$ \infty \f$ (like \ref
	o2scl::inte_qagiu_gsl). The second integrator is for the
	degenerate case, and should integrate between two finite
	values.
    */
    void set_inte(inte<> &unit, inte<> &udit) {
      nit=&unit;
      dit=&udit;
      return;
    }    
    
    /// The default integrator for the non-degenerate regime
    inte_qagiu_gsl<> def_nit;

    /// The default integrator for the degenerate regime
    inte_qag_gsl<> def_dit;

    /// Multiprecision integrator
    inte_double_exp_boost<> it_multip;
    //@}
    
    /// Return string denoting type ("fermion_deriv_rel")
    virtual const char *type() { return "fermion_deriv_rel"; };

  protected:

    /// The integrator for non-degenerate fermions
    inte<> *nit;

    /// The integrator for degenerate fermions
    inte<> *dit;

  };

  /** \brief Double-precision version of 
      \ref o2scl::fermion_deriv_rel_tl 
  */
  typedef fermion_deriv_rel_tl<> fermion_deriv_rel;

  /** \brief Long double version of 
      \ref o2scl::fermion_deriv_rel_tl 
  */
  typedef fermion_deriv_rel_tl<fermion_deriv_tl<long double>,
                               fermion_rel_ld,long double>
  fermion_deriv_rel_ld;

  
}

#endif
