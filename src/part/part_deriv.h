/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_PART_DERIV_H
#define O2SCL_PART_DERIV_H

/** \file part_deriv.h
    \brief File defining \ref o2scl::part_deriv
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/part.h>
#include <o2scl/fermion.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Object to store second derivatives of 
      \f$ P(\mu_n,\mu_p,T) \f$
   */
  class thermo_np_deriv_press {
    
  public:
    
    /// The quantity \f$ (\partial^2 P)/(\partial T^2) \f$
    double dsdT;
    /// The quantity \f$ (\partial^2 P)/(\partial T \partial \mu_n) \f$
    double dnndT;
    /// The quantity \f$ (\partial^2 P)/(\partial T \partial \mu_p) \f$
    double dnpdT;
    /// The quantity \f$ (\partial^2 P)/(\partial \mu_n^2) \f$
    double dnndmun;
    /// The quantity \f$ (\partial^2 P)/(\partial \mu_n \partial \mu_p) \f$
    double dndmu_mixed;
    /// The quantity \f$ (\partial^2 P)/(\partial \mu_p^2) \f$
    double dnpdmup;
    
  };
  
  /** \brief Object to store second derivatives of 
      \f$ f(n_n,n_p,T) \f$
  */
  class thermo_np_deriv_helm {
    
  public:
    
    /// The quantity \f$ (\partial^2 P)/(\partial T^2) \f$
    double dsdT;
    /// The quantity \f$ (\partial^2 P)/(\partial T \partial n_n) \f$
    double dmundT;
    /// The quantity \f$ (\partial^2 P)/(\partial T \partial n_p) \f$
    double dmupdT;
    /// The quantity \f$ (\partial^2 P)/(\partial n_n^2) \f$
    double dmundnn;
    /// The quantity \f$ (\partial^2 P)/(\partial n_n \partial n_p) \f$
    double dmudn_mixed;
    /// The quantity \f$ (\partial^2 P)/(\partial n_p^2) \f$
    double dmupdnp;
  };
  
  /** \brief Particle derivatives in the pressure representation

      This class adds the derivatives \ref dndmu, \ref dndT, and
      \ref dsdT, which correspond to
      \f[
      \left(\frac{d n}{d \mu}\right)_{T,V}, \quad
      \left(\frac{d n}{d T}\right)_{\mu,V}, \quad \mathrm{and} \quad
      \left(\frac{d s}{d T}\right)_{\mu,V}
      \f]
      respectively. All other first-order thermodynamic derivatives
      can be expressed in terms of the first three derivatives. 

      \comment

      (This is no longer required)

      In the
      case that the particle is interacting (i.e. \ref
      part::non_interacting is \c false), then the derivatives which
      are computed are
      \f[
      \left(\frac{d n}{d \nu}\right)_{T,V}, \quad
      \left(\frac{d n}{d T}\right)_{\nu,V}, \quad
      \left(\frac{d s}{d T}\right)_{\nu,V}, \quad \mathrm{and} \quad
      \left(\frac{d n}{d m^{*}}\right)_{T,\nu,V},
      \f]
      If the particles are interacting, no derivative with respect to
      the bare mass is given, since classes cannot know how to relate
      the effective mass to the bare mass.

      \endcomment
  */
  class part_deriv_press {
    
  public:
    
    /// Derivative of number density with respect to chemical potential
    double dndmu;
    
    /// Derivative of number density with respect to temperature
    double dndT;

    /// Derivative of entropy density with respect to temperature
    double dsdT;

    part_deriv_press() {
    }
    
    /// Copy constructor
    part_deriv_press(const part_deriv_press &p) {
      dndmu=p.dndmu;
      dndT=p.dndT;
      dsdT=p.dsdT;
    }

    /// Copy construction with operator=()
    part_deriv_press &operator=(const part_deriv_press &p) {
      if (this!=&p) {
	dndmu=p.dndmu;
	dndT=p.dndT;
	dsdT=p.dsdT;
      }
      return *this;
    }

    /** \brief Compute derivatives in the Helmholtz free energy
	representation from the derivatives in the pressure
	representation
    */
    void deriv_f(double &dmudn, double &dmudT, double &dsdT) {
      dmudn=1.0/dndmu;
      dmudT=-dndT/dndmu;
      dsdT=dndT*dndT/dndmu-dsdT;
      return;
    }
  };
  
  /** \brief A fermion with derivative information
   */
  class fermion_deriv : public fermion, public part_deriv_press {
    
  public:

    /// Make a particle of mass \c mass and degeneracy \c dof.
  fermion_deriv(double mass=0.0, double dof=0.0) : fermion(mass,dof) {
    }
    
    /// Copy constructor
    fermion_deriv(const fermion_deriv &p) {
      g=p.g;
      m=p.m;
      ms=p.ms;
      n=p.n;
      ed=p.ed;
      pr=p.pr;
      mu=p.mu;
      en=p.en;
      nu=p.nu;
      dndmu=p.dndmu;
      dndT=p.dndT;
      dsdT=p.dsdT;
      inc_rest_mass=p.inc_rest_mass;
      non_interacting=p.non_interacting;
    }

    /// Copy constructor
    fermion_deriv(const fermion &p) {
      g=p.g;
      m=p.m;
      ms=p.ms;
      n=p.n;
      ed=p.ed;
      pr=p.pr;
      mu=p.mu;
      en=p.en;
      nu=p.nu;
      dndmu=0.0;
      dndT=0.0;
      dsdT=0.0;
      inc_rest_mass=p.inc_rest_mass;
      non_interacting=p.non_interacting;
    }

    /// Copy construction with operator=()
    fermion_deriv &operator=(const fermion_deriv &p) {
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
	dndmu=p.dndmu;
	dndT=p.dndT;
	dsdT=p.dsdT;
	inc_rest_mass=p.inc_rest_mass;
	non_interacting=p.non_interacting;
      }
      return *this;
    }
    
    /// Copy construction with operator=()
    fermion_deriv &operator=(const fermion &p) {
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
	dndmu=0.0;
	dndT=0.0;
	dsdT=0.0;
	inc_rest_mass=p.inc_rest_mass;
	non_interacting=p.non_interacting;
      }
      return *this;
    }
    
  };
  
  /** \brief A part with derivative information
   */
  class part_deriv : public part, public part_deriv_press {
    
  public:

    /// Make a particle of mass \c mass and degeneracy \c dof.
  part_deriv(double mass=0.0, double dof=0.0) : part(mass,dof) {
    }
    
    /// Copy constructor
    part_deriv(const part_deriv &p) {
      g=p.g;
      m=p.m;
      ms=p.ms;
      n=p.n;
      ed=p.ed;
      pr=p.pr;
      mu=p.mu;
      en=p.en;
      nu=p.nu;
      dndmu=p.dndmu;
      dndT=p.dndT;
      dsdT=p.dsdT;
      inc_rest_mass=p.inc_rest_mass;
      non_interacting=p.non_interacting;
    }

    /// Copy construction with operator=()
    part_deriv &operator=(const part_deriv &p) {
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
	dndmu=p.dndmu;
	dndT=p.dndT;
	dsdT=p.dsdT;
	inc_rest_mass=p.inc_rest_mass;
	non_interacting=p.non_interacting;
      }
      return *this;
    }
    
  };
  
  /** \brief Base quantities for thermodynamic derivatives
   */
  class deriv_thermo_base {
    
  public:
    
    /** \brief The heat capacity per particle at 
	constant volume (unitless)

	This function returns 
	\f[
	c_V = \frac{T}{N} \frac{\partial S}{\partial T}_{V,N} =
	\frac{T}{n} \frac{\partial s}{\partial T}_{V,n} =
	\frac{1}{N} \frac{\partial E}{\partial T}_{V,N} 
	\f]

	This is \f$ 3/2 \f$ for an ideal gas.
    */
    template<class part_deriv_t> 
    double heat_cap_ppart_const_vol(part_deriv_t &p, double temper) {
      return (p.dsdT-p.dndT*p.dndT/p.dndmu)*temper/p.n;
    }
    
    /** \brief The heat capacity per particle 
	at constant pressure (unitless)

	This function computes
	\f[
	c_P \equiv \frac{\partial H}{\partial T}_{P,N}
	\f]
	
	This is \f$ 5/2 \f$ for an ideal gas.
    */
    template<class part_deriv_t> 
    double heat_cap_ppart_const_press(part_deriv_t &p, double temper) {
      return temper/p.n*p.dsdT+p.en*p.en*temper/p.n/p.n/p.n*p.dndmu-
	2.0*p.en*temper/p.n/p.n*p.dndT;
    }

    /** \brief The adiabatic compressibility

	This function computes
	\f[
	\beta_S \equiv - \frac{1}{V} \frac{\partial V}{\partial P}_{S,N}
	\f]
    */
    template<class part_deriv_t> 
    double compress_adiabatic(part_deriv_t &p, double temper) {
      return 0.0;
    }
    
    /** \brief The isothermal compressibility

	This function computes
	\f[
	\beta_T \equiv - \frac{1}{V} \frac{\partial V}{\partial P}_{T,N}
	\f]
     */
    template<class part_deriv_t> 
    double compress_const_tptr(part_deriv_t &p, double temper) {
      return 0.0;
    }

    /** \brief The coefficient of thermal expansion

	This function computes
	\f[
	\frac{1}{V} \frac{\partial V}{\partial T}_{P,N}
	\f]
	In units of an inverse length. 
     */
    template<class part_deriv_t> 
      double coeff_thermal_exp(part_deriv_t &p, double temper) {
      return 0.0;
    }

    /** \brief The squared sound speed (unitless)

	This function computes
	\f[
	c_s^2 = \frac{\partial P}{\partial \varepsilon}_{S,N}
	\f]
    */
    template<class part_deriv_t> 
      double squared_sound_speed(part_deriv_t &p, double temper) {
      return 0.0;
    }
    
  };
  
  /** \brief Compute properties of a fermion including derivatives
      [abstract base]

      \future Include explicit zero-temperature calculation, maybe
      by making this a child of fermion_zerot or by making a 
      new fermion_deriv_zerot? 
      \comment
      dn/dmu is just g*mu*kf/2/pi^2
      \endcomment
      \future There is also a closed form for the derivatives
      of massless fermions with pairs at finite temperature
      in Constantinou et al. 2014 which could be implemented here.
  */
  class fermion_deriv_thermo : public deriv_thermo_base {

  public:

    virtual ~fermion_deriv_thermo() {
    }

    /** \brief Calculate properties as function of chemical potential
     */
    virtual int calc_mu(fermion_deriv &f, double temper)=0;

    /** \brief Calculate properties as function of density
     */
    virtual int calc_density(fermion_deriv &f, double temper)=0;

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual int pair_mu(fermion_deriv &f, double temper)=0;

    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual int pair_density(fermion_deriv &f, double temper)=0;

    /// Calculate effective chemical potential from density
    virtual int nu_from_n(fermion_deriv &f, double temper)=0;

    /** \brief Calculate properties as a function of chemical 
	potential using a degenerate expansion

	\future There is some repetition of the code
	for this function and the function
	\ref o2scl::fermion_thermo::calc_mu_deg() .
	which could be avoided.
    */
    virtual bool calc_mu_deg(fermion_deriv &f, double temper,
			     double prec);
    
    /** \brief Calculate properties as a function of chemical 
	potential using a nondegenerate expansion

	\future There is some repetition of the code
	for this function and the function
	\ref o2scl::fermion_thermo::calc_mu_ndeg() .
	which could be avoided.
    */
    virtual bool calc_mu_ndeg(fermion_deriv &f, double temper,
			      double prec, bool inc_antip=false);

  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
