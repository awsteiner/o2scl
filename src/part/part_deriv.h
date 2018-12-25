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

      The quantities \f$ c_P \f$ computed by 
      \ref heat_cap_ppart_const_press(), \f$ c_V \f$
      computed by \ref heat_cap_ppart_const_vol(), 
      \f$ \beta_T \f$ computed by \ref compress_const_tptr(),
      \f$ \beta_S \f$ computed by \ref compress_adiabatic(),
      and \f$ \alpha_V \f$ computed by \ref coeff_thermal_exp
      are related by 
      \f[
      c_P - c_V = \frac{T \alpha_V^2}{n \beta_T}
      \f]
      and 
      \f[
      \beta_T - \beta_S = \frac{T \alpha_V^2}{n c_P}
      \f]

      For the derivatives below, the following
      Jacobian is useful
      \f{eqnarray*}
      \frac{\partial (P,S,N)}{\partial (V,\mu,T)}
      &=& -n \left[ 
      \left(\frac{\partial S}{\partial V}\right)_{\mu,T}
      \left(\frac{\partial N}{\partial T}\right)_{\mu,V}
      - \left(\frac{\partial S}{\partial T}\right)_{\mu,V}
      \left(\frac{\partial N}{\partial V}\right)_{\mu,T}
      \right] + 
      s \left[ 
      \left(\frac{\partial S}{\partial V}\right)_{\mu,T}
      \left(\frac{\partial N}{\partial \mu}\right)_{V,T}
      - \left(\frac{\partial S}{\partial \mu}\right)_{V,T}
      \left(\frac{\partial N}{\partial V}\right)_{\mu,T}
      \right] 
      \nonumber \\
      &=& - V n \left[ 
      s 
      \left(\frac{\partial n}{\partial T}\right)_{\mu}
      - n \left(\frac{\partial s}{\partial T}\right)_{\mu}
      \right] + V s \left[ s
      \left(\frac{\partial n}{\partial \mu}\right)_{T}
      - n \left(\frac{\partial n}{\partial T}\right)_{T}
      \right] 
      = V n^2 \left(\frac{\partial s}{\partial T}\right)_{\mu}
      - 2 V n s \left(\frac{\partial n}{\partial T}\right)_{\mu}
      + V s^2 \left(\frac{\partial n}{\partial \mu}\right)_{T}
      \f}
      For convenience, we define the quantity
      \f[
      X \equiv \frac{1}{V}
      \frac{\partial (P,S,N)}{\partial (V,\mu,T)}
      \f]
      Another common combination of derivatives is
      \f[
      Y \equiv \left(\frac{\partial n}{\partial T}\right)_{\mu}^2 -
      \left(\frac{\partial s}{\partial T}\right)_{\mu}
      \left(\frac{\partial n}{\partial \mu}\right)_{T}
      \f]
      
   */
  class deriv_thermo_base {
    
  public:
    
    /** \brief The heat capacity per particle at 
	constant volume (unitless)

	This function returns 
	\f[
	c_V = \frac{T}{N} 
	\left(\frac{\partial S}{\partial T}\right)_{V,N} =
	\frac{T}{n} \left(\frac{\partial s}{\partial T}\right)_{V,n} =
	\frac{1}{N} \left(\frac{\partial E}{\partial T}\right)_{V,N} 
	\f]

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f[
	\frac{T}{n} \left(\frac{\partial s}{\partial T}\right)_{V,n}
	= \frac{T}{n} \frac{\partial(s,n,V)}{\partial(T,n,V)} = 
	\frac{T}{n} \left[\frac{\partial(s,n,V)}{\partial(T,\mu,V)}\right]
	\left[\frac{\partial(T,n,V)}{\partial(T,\mu,V)}\right]^{-1}
	\f]
	\f[
	= \frac{T}{n} 
	\left[
	\left(\frac{\partial s}{\partial T}\right)_{\mu} -
	\left(\frac{\partial n}{\partial T}\right)_{\mu}^2
	\left(\frac{\partial n}{\partial \mu}\right)_{T}^{-1}  
	\right]
	\f]

	This is \f$ 3/2 \f$ for an ideal gas.
    */
    template<class part_deriv_t> 
    double heat_cap_ppart_const_vol(part_deriv_t &p, double temper) {
      return (p.dsdT-p.dndT*p.dndT/p.dndmu)*temper/p.n;
    }
    
    /** \brief The heat capacity per particle 
	at constant pressure (unitless)

	This function returns 
	\f[
	c_P = \frac{T}{N} 
	\left(\frac{\partial S}{\partial T}\right)_{P,N} =
	\frac{1}{N} \left(\frac{\partial H}{\partial T}\right)_{P,N} 
	\f]

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f[
	\frac{T}{N} \left(\frac{\partial S}{\partial T}\right)_{P,N}
	= \frac{T}{N} \frac{\partial(S,N,P)}{\partial(T,N,P)} = 
	\frac{T}{N} \left[\frac{\partial(S,N,P)}{\partial(T,\mu,V)}\right]
	\left[\frac{\partial(T,N,P)}{\partial(T,\mu,V)}\right]^{-1}
	\f]
	The first Jacobian was computed above since
	\f[
	\frac{\partial(S,N,P)}{\partial(T,\mu,V)} = -
	\frac{\partial(P,S,N)}{\partial(V,\mu,T)}
	\f]
	The second is
	\f[
	\frac{\partial(T,N,P)}{\partial(T,\mu,V)}
	=
	\left[
	\left(\frac{\partial N}{\partial \mu}\right)_{T,V}
	\left(\frac{\partial P}{\partial V}\right)_{\mu,T}
	- \left(\frac{\partial N}{\partial V}\right)_{\mu,T}
	\left(\frac{\partial P}{\partial \mu}\right)_{T,V}
	\right] = - n^2
	\f]
	The final result is
	\f[
	c_P = \frac{T}{n} \left(\frac{\partial s}{\partial T}\right)_{\mu}
	+ \frac{s^2 T}{n^3} \left(\frac{\partial n}{\partial \mu}\right)_{T}
	- \frac{2 s T}{n^2} \left(\frac{\partial n}{\partial T}\right)_{\mu}
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
	\beta_S \equiv - \frac{1}{V} 
	\left(\frac{\partial V}{\partial P}\right)_{S,N}
	\f]
	(sometimes referred to as \f$ \kappa_S \f$ or 
	\f$ \chi_S \f$)

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f[
	\left(\frac{\partial V}{\partial P}\right)_{S,N} = 
	\frac{\partial (V,S,N)}{\partial (P,S,N)} =
	\frac{\partial (V,S,N)}{\partial (V,\mu,T)}
	\left[ \frac{\partial (P,S,N)}{\partial (V,\mu,T)}\right]^{-1}
	\f]
	The first Jacobian
	\f[
	\frac{\partial (V,S,N)}{\partial (V,\mu,T)} = V^2
	\left[
	\left(\frac{\partial s}{\partial T}\right)_{\mu,V}
	\left(\frac{\partial n}{\partial \mu}\right)_{T,V}
	- \left(\frac{\partial n}{\partial T}\right)_{\mu,V}^2
	\right]
	\f]
	and the second Jacobian was computed above.
	The result is
	\f[
	\beta_S = \left[
	\left(\frac{\partial n}{\partial T}\right)_{\mu}^2 -
	\left(\frac{\partial s}{\partial T}\right)_{\mu}
	\left(\frac{\partial n}{\partial \mu}\right)_{T}
	\right]
	\left[
	n^2 \left(\frac{\partial s}{\partial T}\right)_{\mu,V}
	- 2 n s \left(\frac{\partial n}{\partial T}\right)_{\mu,V}
	+ s^2 \left(\frac{\partial n}{\partial \mu}\right)_{T,V}
	\right]^{-1}
	\f]
    */
    template<class part_deriv_t> 
    double compress_adiabatic(part_deriv_t &p, double temper) {
      return (p.dndT*p.dndT-p.dndmu*p.dsdT)/
	(p.n*p.n*p.dsdT-2.0*p.n*p.en*p.dndT+p.en*p.en*p.dndmu);
    }
    
    /** \brief The isothermal compressibility

	This function computes
	\f[
	\beta_T \equiv - \frac{1}{V} 
	\left(\frac{\partial V}{\partial P}\right)_{T,N}
	\f]
	(sometimes referred to as \f$ \kappa_T \f$ or 
	\f$ \chi_T \f$) in units of inverse length to the fourth 
	power.

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f{eqnarray*}
	- \frac{1}{V} \left(\frac{\partial V}{\partial P}\right)_{T,N} &=& 
	\frac{\partial (V,T,N)}{\partial (P,T,N)} =
	\frac{1}{V}
	\frac{\partial (V,T,N)}{\partial (V,T,\mu)} 
	\left[\frac{\partial (N,P,T)}{\partial (V,\mu,T)}\right]^{-1}
	\nonumber \\ 
	&=& \left(\frac{\partial n}{\partial \mu}\right)_{T,V} 
	\left[
	\left(\frac{\partial N}{\partial V}\right)_{\mu,T} 
	\left(\frac{\partial P}{\partial \mu}\right)_{V,T} 
	- \left(\frac{\partial P}{\partial V}\right)_{\mu,T} 
	\left(\frac{\partial N}{\partial \mu}\right)_{V,T} 
	\right]^{-1} = 
	\frac{1}{n^2} \left(\frac{\partial n}{\partial \mu}\right)_{T} 
	\f}
     */
    template<class part_deriv_t> 
    double compress_const_tptr(part_deriv_t &p, double temper) {
      return p.dndmu/p.n/p.n;
    }

    /** \brief The coefficient of thermal expansion

	This function computes
	\f[
	\alpha_V = 
	\frac{1}{V} \left(\frac{\partial V}{\partial T}\right)_{P,N}
	\f]
	in units of length. 

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f{eqnarray*}
	\left(\frac{\partial V}{\partial T}\right)_{P,N} &=&
	\frac{\partial (V,P,N)}{\partial (T,P,N)} =
	-\frac{\partial (V,P,N)}{\partial (V,T,\mu)} 
	\left[ \frac{\partial (T,P,N)}{\partial (T,V,\mu)} \right]^{-1}
	\nonumber \\
	& = & 
	- \left[ 
	\left(\frac{\partial P}{\partial T}\right)_{\mu,V} 
	\left(\frac{\partial N}{\partial \mu}\right)_{T,V} -
	\left(\frac{\partial N}{\partial T}\right)_{\mu,V} 
	\left(\frac{\partial P}{\partial \mu}\right)_{T,V} 
	\right]
	\left[ 
	\left(\frac{\partial P}{\partial V}\right)_{\mu,T} 
	\left(\frac{\partial N}{\partial \mu}\right)_{V,T} -
	\left(\frac{\partial P}{\partial \mu}\right)_{V,T} 
	\left(\frac{\partial N}{\partial V}\right)_{\mu,T} 
	\right]^{-1}
	\nonumber \\
	&=& \frac{s}{n^2} 
	\left(\frac{\partial n}{\partial \mu}\right)_{T} -
	\frac{1}{n} \left(\frac{\partial n}{\partial T}\right)_{\mu}
	\f}
     */
    template<class part_deriv_t> 
      double coeff_thermal_exp(part_deriv_t &p, double temper) {
      return p.en/p.n/p.n*p.dndmu-p.dndT/p.n;
    }

    /** \brief The squared sound speed (unitless)

	This function computes the sqaured sound speed
	(including relativistic effects)
	\f[
	c_s^2 = \left(\frac{\partial P}
	{\partial \varepsilon}\right)_{S,N}
	\f]
	The result is unitless. To get the units of a squared velocity, 
	one must multiply by \f$ c^2 \f$ . To get the 
	nonrelativistic squared sound speed, you can use
	\f$ c_{s,\mathrm{NR}}^2 = 1/(n \beta_S) \f$ where \f$ \beta_S \f$
	is computed in \ref compress_adiabatic() .

	To write this in terms of the three derivatives in 
	\ref o2scl::part_deriv_press, 
	\f[
	\left(\frac{\partial P}
	{\partial \varepsilon}\right)_{S,N} =
	\frac{\partial (P,S,N)}{\partial (\varepsilon,S,N)} =
	\frac{\partial (P,S,N)}{\partial (V,T,\mu)} 
	\left[ \frac{\partial (\varepsilon,S,N)}
	{\partial (V,T,\mu)} \right]^{-1}
	\f]
	The first Jacobian was computed above (up to a sign).
	The second is the determinant of
	\f[
	\left(
	\begin{array}{ccc}
	0 
	& \frac{\partial \varepsilon}{\partial T} 
	& \frac{\partial \varepsilon}{\partial \mu} \\
	s & V \frac{\partial s}{\partial T} 
	& V \frac{\partial n}{\partial T} \\
	n & V \frac{\partial n}{\partial T} 
	& V \frac{\partial n}{\partial \mu} 
	\end{array}
	\right)
	\f]
	with					
	\f[
	\frac{\partial \varepsilon}{\partial T} =
	-s + T \frac{\partial s}{\partial T}  
	+ \mu \frac{\partial n}{\partial T} 
	\quad \mathrm{and} \quad
	\frac{\partial \varepsilon}{\partial \mu} =
	-n + T \frac{\partial n}{\partial T}  
	+ \mu \frac{\partial n}{\partial \mu} 
	\f]
	giving 
	\f[
	\frac{\partial (\varepsilon,S,N)}
	{\partial (V,T,\mu)} = V 
	\left\{ (P + \varepsilon)
	\left[ \left(\frac{\partial n}{\partial T}\right)^2
	- \left(\frac{\partial n}{\partial \mu}\right)
	\left(\frac{\partial s}{\partial T}\right)
	\right] + n^2
	\left(\frac{\partial s}{\partial T}\right)
	- 2 n s \left(\frac{\partial n}{\partial T}\right)
	+ s^2 \left(\frac{\partial n}{\partial \mu}\right)
	\right\}
	\f]
	The final result is 
	\f[
	c_s^2 = \left\{ \frac{
	\left(P + \varepsilon\right)
	\left[ 
	\left(\frac{\partial n}{\partial \mu}\right)
	\left(\frac{\partial s}{\partial T}\right) -
	\left(\frac{\partial n}{\partial T}\right)^2
	\right]
	}{n^2 \left(\frac{\partial s}{\partial T}\right)
	- 2 n s \left(\frac{\partial n}{\partial T}\right)
	+ s^2 \left(\frac{\partial n}{\partial \mu}\right)}
	-1 \right\}^{-1}
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
