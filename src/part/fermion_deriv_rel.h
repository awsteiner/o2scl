/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic fermion

      \note This class only has preliminary support for
      inc_rest_mass=true (more testing should be done, particularly
      for the "pair" functions)

      \note The testing of this class is apparently sensitive to 
      the difference between gsl_hypot and std::hypot in o2hypot 
      in misc.cpp. Further testing needs to be done to verify 
      which is more accurate. This further testing will probably 
      need to wait until the full multiprecision fermion classes
      are done.

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
      \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) g(k)\right|_{k=a}^{k=b}
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
      
      \future The option err_nonconv=false is not really implemented
      yet.

      \future The \ref pair_density() function is a bit slow because
      it computes the non-derivative thermodynamic quantities 
      twice, and this could be improved.
  */
  template<class fermion_deriv_t=fermion_deriv_tl<double>,
	   class fp_t=double>
  class fermion_deriv_rel_tl : public fermion_deriv_thermo_tl<fp_t> {
    
  public:

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_deriv_rel_tl() {
  
      deg_limit=2.0;
      upper_limit_fac=20.0;

      density_root=&def_density_root;
      nit=&def_nit;
      dit=&def_dit;
  
      method=automatic;
      intl_method=by_parts;

      exp_limit=200.0;

      err_nonconv=true;

      last_method=0;
    }
  
    virtual ~fermion_deriv_rel_tl() {
    }
    
    /** \brief Limit of arguments of exponentials for Fermi functions 
	(default 200.0)
    */
    fp_t exp_limit;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2.0)
    */
    fp_t deg_limit;
    
    /** \brief The limit for the Fermi functions (default 20.0)
	
	fermion_deriv_rel will ignore corrections smaller than about
	\f$ \exp(-\mathrm{f{l}imit}) \f$ . 
    */
    fp_t upper_limit_fac;
    
    /// Storage for the most recently calculated uncertainties 
    fermion_deriv_t unc;

    /// Object for computing non-derivative quantities
    fermion_rel_tl<fermion_deriv_t> fr;

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
    
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief Calculate properties as function of chemical potential
     */
    virtual int calc_mu(fermion_deriv_t &f, fp_t temper) {

      fr.calc_mu(f,temper);
      last_method=fr.last_method*10;
  
      int iret;

      if (temper<=0.0) {
	O2SCL_ERR("T=0 not implemented in fermion_deriv_rel().",exc_eunimpl);
      }

      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      fp_t prefac=f.g/2.0/this->pi2;

      // Compute the degeneracy parameter

      bool deg=false;
      fp_t psi;
      if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
      else psi=(f.nu+f.m-f.ms)/temper;
      if (psi>deg_limit) deg=true;

      // Try the non-degenerate expansion if psi is small enough
      if (psi<-4.0) {
	bool acc=this->calc_mu_ndeg(f,temper,1.0e-14);
	if (acc) {
	  unc.n=f.n*1.0e-14;
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  unc.dndT=f.dndT*1.0e-14;
	  unc.dsdT=f.dsdT*1.0e-14;
	  unc.dndmu=f.dndmu*1.0e-14;
	  last_method+=1;
	  return 0;
	}
      }

      // Try the degenerate expansion if psi is large enough
      if (psi>20.0) {
	bool acc=this->calc_mu_deg(f,temper,1.0e-14);
	if (acc) {
	  unc.n=f.n*1.0e-14;
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  unc.dndT=f.dndT*1.0e-14;
	  unc.dsdT=f.dsdT*1.0e-14;
	  unc.dndmu=f.dndmu*1.0e-14;
	  last_method+=2;
	  return 0;
	}
      }

      if (deg==false) {
    
	// Set integration method
	if (method==automatic) {
	  intl_method=by_parts;
	  last_method+=3;
	} else {
	  intl_method=method;
	  last_method+=4;
	}

	// The non-degenerate case

	funct density_T_fun_f=
	  std::bind(std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
		    (&fermion_deriv_rel_tl<fermion_deriv_t,
		     fp_t>::density_T_fun),
		    this,std::placeholders::_1,std::ref(f),temper);
	iret=nit->integ_iu_err(density_T_fun_f,0.0,f.dndT,unc.dndT);
	if (iret!=0) {
	  O2SCL_ERR2("dndT integration (ndeg) failed in ",
		     "fermion_deriv_rel::calc_mu().",
		     exc_efailed);
	}
	f.dndT*=prefac;
	unc.dndT*=prefac;

	funct density_mu_fun_f=
	  std::bind(std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
		    (&fermion_deriv_rel_tl<fermion_deriv_t,
		     fp_t>::density_mu_fun),
		    this,std::placeholders::_1,std::ref(f),temper);
	iret=nit->integ_iu_err(density_mu_fun_f,0.0,f.dndmu,unc.dndmu);
	if (iret!=0) {
	  O2SCL_ERR2("dndmu integration (ndeg) failed in ",
		     "fermion_deriv_rel::calc_mu().",
		     exc_efailed);
	}
	f.dndmu*=prefac;
	unc.dndmu*=prefac;
    
	funct entropy_T_fun_f=
	  std::bind(std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
		    (&fermion_deriv_rel_tl<fermion_deriv_t,
		     fp_t>::entropy_T_fun),
		    this,std::placeholders::_1,std::ref(f),temper);
	iret=nit->integ_iu_err(entropy_T_fun_f,0.0,f.dsdT,unc.dsdT);
	if (iret!=0) {
	  O2SCL_ERR2("dsdT integration (ndeg) failed in ",
		     "fermion_deriv_rel_tl<fp_t>::calc_mu().",exc_efailed);
	}
	f.dsdT*=prefac;
	unc.dsdT*=prefac;

      } else {

	// Compute the upper limit for degenerate integrals

	fp_t arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	}
	// Set to zero to avoid uninit'ed var. warnings
	fp_t ul=0.0;
	if (arg>0.0) {
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
	  arg=pow(-upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	} else {
	  arg=pow(-upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	}

	// Set integration method
	if (method==automatic) {
	  if ((!f.inc_rest_mass && (f.nu+f.m-f.ms)/temper>1.0e3) ||
	      (f.inc_rest_mass && (f.nu-f.ms)/temper>1.0e3)) {
	    intl_method=direct;
	    last_method+=5;
	  } else {
	    intl_method=by_parts;
	    last_method+=6;
	  }
	} else {
	  intl_method=method;
	  last_method+=7;
	}

	funct deg_density_mu_fun_f=
	  std::bind(std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
		    (&fermion_deriv_rel_tl<fermion_deriv_t,
		     fp_t>::deg_density_mu_fun),
		    this,std::placeholders::_1,std::ref(f),temper);
	if (intl_method==direct && ll>0.0) {
	  iret=dit->integ_err(deg_density_mu_fun_f,ll,ul,
			      f.dndmu,unc.dndmu);
	} else {
	  iret=dit->integ_err(deg_density_mu_fun_f,0.0,ul,
			      f.dndmu,unc.dndmu);
	}
	if (iret!=0) {
	  O2SCL_ERR2("dndmu integration (deg) failed in ",
		     "fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
		     exc_efailed);
	}
	f.dndmu*=prefac;
	unc.dndmu*=prefac;
    
	funct deg_density_T_fun_f=std::bind
	  (std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
	   (&fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::deg_density_T_fun),
	   this,std::placeholders::_1,std::ref(f),temper);
	if (intl_method==direct && ll>0.0) {
	  iret=dit->integ_err(deg_density_T_fun_f,ll,ul,f.dndT,unc.dndT);
	} else {
	  iret=dit->integ_err(deg_density_T_fun_f,0.0,ul,f.dndT,unc.dndT);
	}
	if (iret!=0) {
	  O2SCL_ERR2("dndT integration (deg) failed in ",
		     "fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
		     exc_efailed);
	}
	f.dndT*=prefac;
	unc.dndT*=prefac;

	funct deg_entropy_T_fun_f=std::bind
	  (std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
	   (&fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::deg_entropy_T_fun),
	   this,std::placeholders::_1,std::ref(f),temper);
	if (intl_method==direct && ll>0.0) {
	  iret=dit->integ_err(deg_entropy_T_fun_f,ll,ul,f.dsdT,unc.dsdT);
	} else {
	  iret=dit->integ_err(deg_entropy_T_fun_f,0.0,ul,f.dsdT,unc.dsdT);
	}
	if (iret!=0) {
	  O2SCL_ERR2("dsdT integration (deg) failed in ",
		     "fermion_deriv_rel_tl<fermion_deriv_t,fp_t>::calc_mu().",
		     exc_efailed);
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
      unc.pr=0.0;

      return 0;
    }
  
    /** \brief Calculate properties as function of density
     */
    virtual int calc_density(fermion_deriv_t &f, fp_t temper) {
  
      if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }
  
      nu_from_n(f,temper);
      int lm=last_method;
  
      if (f.non_interacting) { f.mu=f.nu; }

      calc_mu(f,temper);
      last_method+=lm;

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
      calc_mu(antip,temper);
      last_method+=lm;

      f.n-=antip.n;
      f.pr+=antip.pr;
      if (f.inc_rest_mass) {
	f.ed+=antip.ed;
      } else {
	f.ed=f.ed+antip.ed+2.0*antip.n*f.m;
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
      return ret;
    }
  
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
    
    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_root(root<> &rp) {
      density_root=&rp;
      return;
    }

    /// The default integrator for the non-degenerate regime
    inte_qagiu_gsl<> def_nit;

    /// The default integrator for the degenerate regime
    inte_qag_gsl<> def_dit;

    /// The default solver for npen_density() and pair_density()
    root_cern<> def_density_root;
    
    /// Return string denoting type ("fermion_deriv_rel")
    virtual const char *type() { return "fermion_deriv_rel"; };

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// The internal integration method
    int intl_method;

    /// The integrator for non-degenerate fermions
    inte<> *nit;

    /// The integrator for degenerate fermions
    inte<> *dit;

    /// The solver for calc_density() and pair_density()
    root<> *density_root;

    /** \name The integrands, as a function of \f$ u=k/T \f$, for 
	non-degenerate integrals
    */
    //@{
    /** \brief Desc
     */
    fp_t density_T_fun(fp_t u, fermion_deriv_t &f, fp_t T) {
      fp_t k=u*(T), E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*(E-f.nu)/T*ff*(1.0-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(f.nu)/T-k*k*(f.nu)/T/E)*
	    T*fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*(E-f.nu)/T*ff*(1.0-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(f.nu+f.m)/T-k*k*(f.nu+f.m)/T/E)*
	    T*fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    fp_t density_mu_fun(fp_t u, fermion_deriv_t &f, fp_t T) {
      fp_t k=u*(T), E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*ff*(1.0-ff);
	} else {
	  ret=T*(E*E+k*k)/E*fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*ff*(1.0-ff);
	} else {
	  ret=T*(E*E+k*k)/E*fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }
    
    /** \brief Desc
     */
    fp_t entropy_T_fun(fp_t u, fermion_deriv_t &f, fp_t T) {
      fp_t k=u*T, E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=T*k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
	} else {
	  ret=(E-f.nu)/E/T*
	    (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu))*
	    fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=T*k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
	} else {
	  ret=(E-f.m-f.nu)/E/T*
	    (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu+f.m))*
	    fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    fp_t density_ms_fun(fp_t u, fermion_deriv_t &f, fp_t T) {
      fp_t k=u*T, E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=-k*k*f.ms/(E)/T*ff*(1.0-ff);
	} else {
	  ret=-f.ms*fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=-k*k*f.ms/(E+f.m)/T*ff*(1.0-ff);
	} else {
	  ret=-f.ms*fermi_function(E-f.m,f.nu,T,exp_limit);
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
    /** \brief Desc
     */
    fp_t deg_density_T_fun(fp_t k, fermion_deriv_t &f, fp_t T) {
      fp_t E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*(E-f.nu)/T/T*ff*(1.0-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(f.nu)/T-k*k*(f.nu)/T/E)*
	    fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k*(E-f.nu)/T/T*ff*(1.0-ff);
	} else {
	  ret=(2.0*k*k/T+E*E/T-E*(f.nu+f.m)/T-k*k*(f.nu+f.m)/T/E)*
	    fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    fp_t deg_density_mu_fun(fp_t k, fermion_deriv_t &f, fp_t T) {
      fp_t E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k/T*ff*(1.0-ff);
	} else {
	  ret=(E*E+k*k)/E*fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=k*k/T*ff*(1.0-ff);
	} else {
	  ret=(E*E+k*k)/E*fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    fp_t deg_entropy_T_fun(fp_t k, fermion_deriv_t &f, fp_t T) {
      fp_t E, ret;
      E=o2hypot(k,f.ms);
      if (f.inc_rest_mass) {
	fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	if (intl_method==direct) {
	  ret=k*k*ff*(1.0-ff)*pow(E-f.nu,2.0)/pow(T,3.0);
	} else {
	  ret=(E-f.nu)/E/T/T*
	    (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*f.nu)*ff;
	}
      } else {
	fp_t ff=fermi_function(E-f.m,f.nu,T,exp_limit);
	if (intl_method==direct) {
	  ret=k*k*ff*(1.0-ff)*pow(E-f.nu-f.m,2.0)/pow(T,3.0);
	} else {
	  ret=(E-f.m-f.nu)/E/T/T*
	    (pow(E,3.0)+3.0*E*k*k-(E*E+k*k)*(f.nu+f.m))*ff;
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    fp_t deg_density_ms_fun(fp_t k, fermion_deriv_t &f, fp_t T) {
      fp_t E, ret;
      if (f.inc_rest_mass) {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=-k*k*f.ms/E/T*ff*(1.0-ff);
	} else {
	  ret=-f.ms*fermi_function(E,f.nu,T,exp_limit);
	}
      } else {
	E=o2hypot(k,f.ms);
	if (intl_method==direct) {
	  E-=f.m;
	  fp_t ff=fermi_function(E,f.nu,T,exp_limit);
	  ret=-k*k*f.ms/(E+f.m)/T*ff*(1.0-ff);
	} else {
	  ret=-f.ms*fermi_function(E-f.m,f.nu,T,exp_limit);
	}
      }
      return ret;
    }
    
    //@}

#endif

  };

  /** \brief Double-precision version of 
      \ref o2scl::fermion_deriv_rel_tl 
  */
  typedef fermion_deriv_rel_tl<> fermion_deriv_rel;
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
