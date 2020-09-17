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
#ifndef O2SCL_PART_DERIV_H
#define O2SCL_PART_DERIV_H

/** \file part_deriv.h
    \brief File defining \ref o2scl::part_deriv_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/part.h>
#include <o2scl/fermion.h>
#include <o2scl/fermion_rel.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Object to store second derivatives of 
      \f$ P(\mu_n,\mu_p,T) \f$
  */
  template<class fp_t=double>
    class thermo_np_deriv_press_tl {
    
  public:
    
  /// The quantity \f$ (\partial^2 P)/(\partial T^2) \f$
  fp_t dsdT;
  /// The quantity \f$ (\partial^2 P)/(\partial T \partial \mu_n) \f$
  fp_t dnndT;
  /// The quantity \f$ (\partial^2 P)/(\partial T \partial \mu_p) \f$
  fp_t dnpdT;
  /// The quantity \f$ (\partial^2 P)/(\partial \mu_n^2) \f$
  fp_t dnndmun;
  /// The quantity \f$ (\partial^2 P)/(\partial \mu_n \partial \mu_p) \f$
  fp_t dndmu_mixed;
  /// The quantity \f$ (\partial^2 P)/(\partial \mu_p^2) \f$
  fp_t dnpdmup;
    
  };

  typedef thermo_np_deriv_press_tl<double> thermo_np_deriv_press;
  
  /** \brief Object to store second derivatives of 
      \f$ f(n_n,n_p,T) \f$
  */
  template<class fp_t=double>
    class thermo_np_deriv_helm_tl {
    
  public:
    
  /// The quantity \f$ (\partial^2 P)/(\partial T^2) \f$
  fp_t dsdT;
  /// The quantity \f$ (\partial^2 P)/(\partial T \partial n_n) \f$
  fp_t dmundT;
  /// The quantity \f$ (\partial^2 P)/(\partial T \partial n_p) \f$
  fp_t dmupdT;
  /// The quantity \f$ (\partial^2 P)/(\partial n_n^2) \f$
  fp_t dmundnn;
  /// The quantity \f$ (\partial^2 P)/(\partial n_n \partial n_p) \f$
  fp_t dmudn_mixed;
  /// The quantity \f$ (\partial^2 P)/(\partial n_p^2) \f$
  fp_t dmupdnp;
  };

  typedef thermo_np_deriv_helm_tl<double> thermo_np_deriv_helm;
  
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
  template<class fp_t=double>
    class part_deriv_press_tl {
    
  public:
    
  /// Derivative of number density with respect to chemical potential
  fp_t dndmu;
    
  /// Derivative of number density with respect to temperature
  fp_t dndT;

  /// Derivative of entropy density with respect to temperature
  fp_t dsdT;

  part_deriv_press_tl() {
  }
    
  /// Copy constructor
  part_deriv_press_tl(const part_deriv_press_tl &p) {
    dndmu=p.dndmu;
    dndT=p.dndT;
    dsdT=p.dsdT;
  }

  /// Copy construction with operator=()
  part_deriv_press_tl &operator=(const part_deriv_press_tl &p) {
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
  void deriv_f(fp_t &dmudn, fp_t &dmudT, fp_t &dsdT_n) {
    dmudn=1.0/dndmu;
    dmudT=-dndT/dndmu;
    dsdT_n=dsdT-dndT*dndT/dndmu;
    return;
  }
  };

  /** \brief Double-precision version of \ref o2scl::part_deriv_press_tl 
   */
  typedef part_deriv_press_tl<double> part_deriv_press;
  
  /** \brief A fermion with derivative information
   */
  template<class fp_t=double>
    class fermion_deriv_tl : public fermion_tl<fp_t>,
    public part_deriv_press_tl<fp_t> {
    
  public:

  /// Make a particle of mass \c mass and degeneracy \c dof.
  fermion_deriv_tl(fp_t mass=0.0, fp_t dof=0.0) : fermion(mass,dof) {
    }
    
  /// Copy constructor
  fermion_deriv_tl(const fermion_deriv_tl &p) {
    this->g=p.g;
    this->m=p.m;
    this->ms=p.ms;
    this->n=p.n;
    this->ed=p.ed;
    this->pr=p.pr;
    this->mu=p.mu;
    this->en=p.en;
    this->nu=p.nu;
    this->dndmu=p.dndmu;
    this->dndT=p.dndT;
    this->dsdT=p.dsdT;
    this->inc_rest_mass=p.inc_rest_mass;
    this->non_interacting=p.non_interacting;
  }

  /// Copy constructor
  fermion_deriv_tl(const fermion &p) {
    this->g=p.g;
    this->m=p.m;
    this->ms=p.ms;
    this->n=p.n;
    this->ed=p.ed;
    this->pr=p.pr;
    this->mu=p.mu;
    this->en=p.en;
    this->nu=p.nu;
    this->dndmu=0.0;
    this->dndT=0.0;
    this->dsdT=0.0;
    this->inc_rest_mass=p.inc_rest_mass;
    this->non_interacting=p.non_interacting;
  }

  /// Copy construction with operator=()
  fermion_deriv_tl &operator=(const fermion_deriv_tl &p) {
    if (this!=&p) {
      this->g=p.g;
      this->m=p.m;
      this->ms=p.ms;
      this->n=p.n;
      this->ed=p.ed;
      this->pr=p.pr;
      this->mu=p.mu;
      this->en=p.en;
      this->nu=p.nu;
      this->dndmu=p.dndmu;
      this->dndT=p.dndT;
      this->dsdT=p.dsdT;
      this->inc_rest_mass=p.inc_rest_mass;
      this->non_interacting=p.non_interacting;
    }
    return *this;
  }
    
  /// Copy construction with operator=()
  fermion_deriv_tl &operator=(const fermion &p) {
    if (this!=&p) {
      this->g=p.g;
      this->m=p.m;
      this->ms=p.ms;
      this->n=p.n;
      this->ed=p.ed;
      this->pr=p.pr;
      this->mu=p.mu;
      this->en=p.en;
      this->nu=p.nu;
      this->dndmu=0.0;
      this->dndT=0.0;
      this->dsdT=0.0;
      this->inc_rest_mass=p.inc_rest_mass;
      this->non_interacting=p.non_interacting;
    }
    return *this;
  }
    
  };

  typedef fermion_deriv_tl<double> fermion_deriv;
  
  /** \brief A part with derivative information
   */
  template<class fp_t=double>
    class part_deriv_tl : public part_tl<fp_t>,
    public part_deriv_press_tl<fp_t> {
    
  public:

  /// Make a particle of mass \c mass and degeneracy \c dof.
  part_deriv_tl(fp_t mass=0.0, fp_t dof=0.0) : part(mass,dof) {
    }
    
  /// Copy constructor
  part_deriv_tl(const part_deriv_tl &p) {
    this->g=p.g;
    this->m=p.m;
    this->ms=p.ms;
    this->n=p.n;
    this->ed=p.ed;
    this->pr=p.pr;
    this->mu=p.mu;
    this->en=p.en;
    this->nu=p.nu;
    this->dndmu=p.dndmu;
    this->dndT=p.dndT;
    this->dsdT=p.dsdT;
    this->inc_rest_mass=p.inc_rest_mass;
    this->non_interacting=p.non_interacting;
  }

  /// Copy construction with operator=()
  part_deriv_tl &operator=(const part_deriv_tl &p) {
    if (this!=&p) {
      this->g=p.g;
      this->m=p.m;
      this->ms=p.ms;
      this->n=p.n;
      this->ed=p.ed;
      this->pr=p.pr;
      this->mu=p.mu;
      this->en=p.en;
      this->nu=p.nu;
      this->dndmu=p.dndmu;
      this->dndT=p.dndT;
      this->dsdT=p.dsdT;
      this->inc_rest_mass=p.inc_rest_mass;
      this->non_interacting=p.non_interacting;
    }
    return *this;
  }
    
  };

  typedef part_deriv_tl<double> part_deriv;
  
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
      \left[ \frac{\partial (P,S,N)}{\partial (V,\mu,T)} \right]
      \f]
      Another common combination of derivatives is
      \f[
      Y \equiv \left(\frac{\partial n}{\partial T}\right)_{\mu}^2 -
      \left(\frac{\partial s}{\partial T}\right)_{\mu}
      \left(\frac{\partial n}{\partial \mu}\right)_{T}
      \f]
      
  */
  template<class fp_t=double> class deriv_thermo_base_tl {
    
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
      \ref o2scl::part_deriv_press_tl, 
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
  fp_t heat_cap_ppart_const_vol(part_deriv_t &p, fp_t temper) {
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
      \ref o2scl::part_deriv_press_tl, 
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
      c_P = \frac{T X}{n^3} = 
      \frac{T}{n} \left(\frac{\partial s}{\partial T}\right)_{\mu}
      + \frac{s^2 T}{n^3} \left(\frac{\partial n}{\partial \mu}\right)_{T}
      - \frac{2 s T}{n^2} \left(\frac{\partial n}{\partial T}\right)_{\mu}
      \f]
	
      This is \f$ 5/2 \f$ for an ideal gas.
  */
  template<class part_deriv_t> 
  fp_t heat_cap_ppart_const_press(part_deriv_t &p, fp_t temper) {
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
      \ref o2scl::part_deriv_press_tl, 
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
      \beta_S = Y/X = \left[
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
  fp_t compress_adiabatic(part_deriv_t &p, fp_t temper) {
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
      \ref o2scl::part_deriv_press_tl, 
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
  fp_t compress_const_tptr(part_deriv_t &p, fp_t temper) {
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
      \ref o2scl::part_deriv_press_tl, 
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
  fp_t coeff_thermal_exp(part_deriv_t &p, fp_t temper) {
    return p.en/p.n/p.n*p.dndmu-p.dndT/p.n;
  }

  /** \brief The squared sound speed (unitless)

      This function computes the squared sound speed
      (including relativistic effects)
      \f[
      c_s^2 = \left(\frac{\partial P}
      {\partial \varepsilon}\right)_{S,N}
      \f]
      The result is unitless. To get the units of a squared velocity, 
      one must multiply by \f$ c^2 \f$ . 

      The 
      nonrelativistic squared sound speed
      is 
      \f[
      c_{s,\mathrm{NR}}^2 = \left[\frac{\partial P}
      {\partial (N/V)}\right]_{S,N} = 
      - \frac{V^2}{N} \left(\frac{\partial P}
      {\partial V}\right)_{S,N} = \frac{1}{n \beta_S}
      \f]
      where \f$ \beta_S \f$
      is computed in \ref compress_adiabatic() .

      To write \f$ c_s^2 \f$ in terms of the three derivatives in 
      \ref o2scl::part_deriv_press_tl, 
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
      -T \frac{\partial s}{\partial T}  
      + \mu \frac{\partial n}{\partial T} 
      \quad \mathrm{and} \quad
      \frac{\partial \varepsilon}{\partial \mu} =
      T \frac{\partial n}{\partial T}  
      + \mu \frac{\partial n}{\partial \mu} 
      \f]
      giving 
      \f[
      \frac{\partial (\varepsilon,S,N)}
      {\partial (V,T,\mu)} = V 
      (P + \varepsilon)
      \left[ \left(\frac{\partial n}{\partial T}\right)^2
      - \left(\frac{\partial n}{\partial \mu}\right)
      \left(\frac{\partial s}{\partial T}\right)
      \right] = V Y \left(P+\varepsilon\right)
      \f]
      The final result is 
      \f[
      c_s^2 = 
      - \frac{X}{(P+\varepsilon)Y}
      = 
      \frac{
      n^2 \left(\frac{\partial s}{\partial T}\right)
      - 2 n s \left(\frac{\partial n}{\partial T}\right)
      + s^2 \left(\frac{\partial n}{\partial \mu}\right)
      }{
      \left(P + \varepsilon\right)
      \left[ 
      \left(\frac{\partial n}{\partial \mu}\right)
      \left(\frac{\partial s}{\partial T}\right) -
      \left(\frac{\partial n}{\partial T}\right)^2
      \right]
      }
      \f]

  */
  template<class part_deriv_t> 
  fp_t squared_sound_speed(part_deriv_t &p, fp_t temper) {
    fp_t edt;
    if (p.inc_rest_mass) {
      edt=p.ed;
    } else {
      edt=p.ed+p.n*p.m;
    }
    return (p.n*p.n*p.dsdT-2.0*p.n*p.en*p.dndT+p.en*p.en*p.dndmu)/
    (edt+p.pr)/(p.dndmu*p.dsdT-p.dndT*p.dndT);
  }
    
  };
  
  typedef deriv_thermo_base_tl<double> deriv_thermo_base;
  
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
  template<class fp_t=double>
    class fermion_deriv_thermo_tl : public deriv_thermo_base_tl<fp_t> {

  protected:

  /** \brief A fermion_thermo object 

      This is for access to fermion_thermo::ndeg_terms().
  */
    fermion_rel_tl<fermion_deriv> fr;

  /// Store \f$ \pi \f$ for convenience
  fp_t pi;
  
  /// Store \f$ \pi^2 \f$ for convenience
  fp_t pi2;
  
  public:

  fermion_deriv_thermo_tl() {
    pi=boost::math::constants::pi<fp_t>();
    pi2=boost::math::constants::pi_sqr<fp_t>();
  }    

  virtual ~fermion_deriv_thermo_tl() {
  }

  /** \brief Calculate properties as function of chemical potential
   */
  virtual int calc_mu(fermion_deriv &f, fp_t temper)=0;

  /** \brief Calculate properties as function of density
   */
  virtual int calc_density(fermion_deriv &f, fp_t temper)=0;

  /** \brief Calculate properties with antiparticles as function of
      chemical potential
  */
  virtual int pair_mu(fermion_deriv &f, fp_t temper)=0;

  /** \brief Calculate properties with antiparticles as function of
      density
  */
  virtual int pair_density(fermion_deriv &f, fp_t temper)=0;

  /// Calculate effective chemical potential from density
  virtual int nu_from_n(fermion_deriv &f, fp_t temper)=0;

  /** \brief Calculate properties as a function of chemical 
      potential using a degenerate expansion

      \future There is some repetition of the code
      for this function and the function
      \ref o2scl::fermion_thermo_tl::calc_mu_deg() .
      which could be avoided.
  */
  virtual bool calc_mu_deg(fermion_deriv &f, fp_t temper,
			   fp_t prec) {
      
    if (fr.calc_mu_deg_tlate(f,temper,prec)==false) {
      return false;
    }
      
    // Compute psi and tt
    fp_t psi;
    if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
    else psi=(f.nu+f.m-f.ms)/temper;
    fp_t tt=temper/f.ms;
      
    // Prefactor 'd' in Johns96
    fp_t prefac=f.g/2.0/o2scl_const::pi2*pow(f.ms,4.0);
      
    // Define x = psi * t = (mu/m - 1) and related values
    fp_t x=psi*tt;
    fp_t sx=sqrt(x);
    fp_t s2x=sqrt(2.0+x);
    fp_t x2=x*x;
    fp_t x3=x2*x;
    fp_t x4=x2*x2;
      
    // First order density term (first order entropy term is zero)
    fp_t dndmu_term1=sx*s2x*(1.0+x)/f.ms/f.ms;
      
    // Second order terms
    fp_t dndmu_term2=tt*tt*o2scl_const::pi2/6.0*(1.0+x)*
      (-1.0+2.0*x*(2.0+x))/
      f.ms/f.ms/sx/s2x/x/(2.0+x);
    fp_t dndT_term2=tt*o2scl_const::pi2/3.0*(1.0+2.0*x*(2.0+x))/
      f.ms/f.ms/sx/s2x;
    fp_t dsdT_term2=o2scl_const::pi2/3.0*(1.0+x)*sx*s2x/
      f.ms/f.ms;
      
    // Third order terms
    fp_t dndmu_term3=-7.0*pow(o2scl_const::pi*tt,4.0)/24.0*
      (1.0+x)/sx/s2x/x3/(x+2.0)/(x+2.0)/(x+2.0)/f.ms/f.ms;
    fp_t dndT_term3=7.0*pow(o2scl_const::pi*tt,4.0)/tt/30.0/
      f.ms/f.ms/pow(x*(2.0+x),2.5);
    fp_t dsdT_term3=7.0*pow(o2scl_const::pi*tt,2.0)*
      o2scl_const::pi2/30.0/f.ms/f.ms*(1.0+x)*(-1.0+2.0*x*(2.0+x))/
      x/(2.0+x)/sx/s2x;
      
    // Fourth order terms for density and entropy
    fp_t dndmu_term4=-31.0*pow(o2scl_const::pi*tt,6.0)/48.0*(1.0+x)*
      (3.0+2.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),5.5);
    fp_t dndT_term4=31.0*pow(o2scl_const::pi*tt,6.0)/tt/168.0*
      (7.0+6.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),4.5);
    fp_t dsdT_term4=-155.0*pow(o2scl_const::pi*tt,4.0)*
      o2scl_const::pi2/168.0*
      (1.0+x)/f.ms/f.ms/pow(x*(2.0+x),3.5);
      
    // Add up all the terms
    f.dndmu=prefac*(dndmu_term1+dndmu_term2+dndmu_term3+dndmu_term4);
    f.dndT=prefac*(dndT_term2+dndT_term3+dndT_term4);
    f.dsdT=prefac*(dsdT_term2+dsdT_term3+dsdT_term4);
      
    return true;
  }
    
  /** \brief Calculate properties as a function of chemical 
      potential using a nondegenerate expansion

      \future There is some repetition of the code
      for this function and the function
      \ref o2scl::fermion_thermo_tl::calc_mu_ndeg() .
      which could be avoided.
  */
  virtual bool calc_mu_ndeg(fermion_deriv &f, fp_t temper,
			    fp_t prec, bool inc_antip=false) {
      
    if (fr.calc_mu_ndeg_tlate(f,temper,prec,
			      inc_antip)==false) {
      return false;
    }
      
    // Compute psi and tt
    fp_t psi, psi_num;
    if (f.inc_rest_mass) {
      psi_num=f.nu-f.ms;
    } else {
      psi_num=f.nu+f.m-f.ms;
    }
    psi=psi_num/temper;
    fp_t tt=temper/f.ms;
    fp_t xx=psi*tt;
      
    // Prefactor 'd' in Johns96
    fp_t prefac=f.g/2.0/o2scl_const::pi2*pow(f.ms,4.0);
      
    // One term is always used, so only values of max_term greater than
    // 0 are useful.
    static const size_t max_term=200;
      
    fp_t first_dndT=0.0;
    fp_t first_dsdT=0.0;
    fp_t first_dndmu=0.0;
      
    fp_t nu2=f.nu;
    if (f.inc_rest_mass==false) nu2+=f.m;
      
    f.dndmu=0.0;
    f.dndT=0.0;
    f.dsdT=0.0;
      
    for(size_t j=1;j<=max_term;j++) {
	
      fp_t dj=((fp_t)j);
      fp_t jot=dj/tt;
	
      // Here, we are only computing the derivatives, but we need to
      // compute the terms in the pressure, density, and entropy because
      // they are used in computing the terms for the derivatives.
      fp_t pterm, nterm, enterm;
      fp_t dndmu_term, dndT_term, dsdT_term;
	
      fr.ndeg_terms(j,tt,psi*tt,f.ms,f.inc_rest_mass,inc_antip,
		    pterm,nterm,enterm);
	
      if (inc_antip==false) {
	dndmu_term=nterm*jot;
	dndT_term=jot*enterm-nterm/tt;
	dsdT_term=(3.0*tt-2.0*dj*xx-2.0*dj)/tt/tt*enterm+
	  (5.0*dj*tt-2.0*dj*dj*xx+5.0*dj*tt*xx-dj*dj*xx*xx)/
	  dj/tt/tt/tt*nterm;
      } else {
	dndmu_term=nterm*jot;
	dndT_term=jot*enterm*tanh(jot*(xx+1.0))-
	  (tt+2.0*dj*(1.0+xx))/sinh(jot*(xx+1.0))*nterm/tt/tt;
	dsdT_term=(2.0*dj*(1.0+xx)*tanh(jot*(xx+1.0))-3.0*tt)*enterm/tt/tt+
	  (2.0*pow(dj*1.0+xx,2.0)*tanh(jot*(xx+1.0))-
	   dj*dj*(2.0+2.0*xx+xx*xx)*cosh(jot*(xx+1.0))-
	   5.0*dj*(1.0+xx)*tt)*nterm/dj/tt/tt/tt;
      }
      dndmu_term/=f.ms;
      dndT_term/=f.ms;
      dsdT_term/=f.ms;
	
      if (j==1) {
	first_dndT=dndT_term;
	first_dsdT=dsdT_term;
	first_dndmu=dndmu_term;
      }
      f.dndmu+=dndmu_term;
      f.dndT+=dndT_term;
      f.dsdT+=dsdT_term;
	
      /*
	cout << j << " " << dj << " " << tt << " " << xx << " "
	<< f.ms << " " << pterm << " " << nterm << " "
	<< enterm << " "
	<< dndmu_term << " " << dndT_term << endl;
      */
	
      // If the first terms are zero, then the rest of the terms
      // will be zero so just return early
      if (first_dndT==0.0 && first_dndmu==0.0 && first_dsdT==0.0) {
	f.dndmu=0.0;
	f.dndT=0.0;
	f.dsdT=0.0;
	return true;
      }
	
      // Stop if the last term is sufficiently small compared to
      // the first term
      if (j>1 &&
	  fabs(dndT_term)<prec*fabs(first_dndT) &&
	  fabs(dndmu_term)<prec*fabs(first_dndmu) &&
	  fabs(dsdT_term)<prec*fabs(first_dsdT)) {
	f.dndT*=prefac;
	f.dndmu*=prefac;
	f.dsdT*=prefac;
	return true;
      }
	
      // End of 'for(size_t j=1;j<=max_term;j++)'
    }
      
    // We failed to add enough terms, so return false
    return false;
  }

  };

  typedef fermion_deriv_thermo_tl<double> fermion_deriv_thermo;

  /** \brief Object to organize calibration of derivative quantities
      in particle classes to results stored in a table
  */
  class part_deriv_calibrate_class : public part_calibrate_class {
    
  public:
  
    /** \brief Calibrate a particle thermodynamics class with
	derivatives with results stored in a table

	This compares the approximation to the exact results using
	calc_density(), calc_mu(), pair_density() and pair_mu(). It
	tries each function twelve times. It tries three different
	temperatures, setting both <tt>inc_rest_mass</tt> and
	<tt>non_interacting</tt> equal to <tt>true</tt> and
	<tt>false</tt>.
      
	The <tt>verbose</tt> parameter controls the amount of output.
    */
    template<class part_t, class thermo_t>
      double part_deriv_calibrate(part_t &p, thermo_t &th, bool test_pair,
				  std::string file, bool nr_mode=false,
				  int verbose=0, bool external=false) {
				
      double ret=0;
  
      // ----------------------------------------------------------------
      // Will return to these original values afterwards

      part_t orig=p;

      // ----------------------------------------------------------------
      // Read data file

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
      hdf_input(hf,tab,name);
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
      double m_bad=0.0, mu_bad=0.0, T_bad=0.0, mot_bad=0.0, psi_bad=0.0;
      p.non_interacting=true;
  
      // ----------------------------------------------------------------
      // First pass, test calc_mu() 

      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	double ret_local=0.0;

	// Initialize storage
	dev.n=0.0; dev.dndmu=0.0; dev.dndT=0.0; dev.dsdT=0.0;
	bad.n=0.0; bad.dndmu=0.0; bad.dndT=0.0; bad.dsdT=0.0;
    
	// Temperature loop
	for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    double mot=tab.get("mot",i);
	    double psi=tab.get("psi",i);
	    exact.n=tab.get("n",i);
	    exact.dndmu=tab.get("dndmu",i);
	    exact.dndT=tab.get("dndT",i);
	    exact.dsdT=tab.get("dsdT",i);
      
	    set_mass_flags(p,mot,T,k);
	    set_chem_pot(p,psi,T,k,nr_mode);
	
	    th.calc_mu(p,T);
	
	    exact.n*=pow(T,3.0);
	    exact.dndmu*=pow(T,2.0);
	    exact.dndT*=pow(T,2.0);
	    exact.dsdT*=pow(T,2.0);
	
	    dev.dndmu+=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	    dev.dndT+=fabs((p.dndT-exact.dndT)/exact.dndT);
	    dev.dsdT+=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	
	    cnt++;
	    
	    check_derivs<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,m_bad,T_bad,
				 mot_bad,psi_bad,ret_local);

	    if (verbose>1) {
	      std::cout.precision(5);
	      if (k>=2) {
		std::cout << "T,ms,nu,psi,mot: " << T << " "
			  << p.ms << " " << p.nu
			  << " " << psi << " " << mot << std::endl;
	      } else {
		std::cout << "T,m,mu,psi,mot: " << T << " "
			  << p.m << " " << p.mu
			  << " " << psi << " " << mot << std::endl;
	      }
	      std::cout.precision(5);
	      std::cout << "n,dndmu,dndT,dsdT: " << std::endl;
	      std::cout << "approx: " << p.n << " " << p.dndmu << " "
			<< p.dndT << " " 
			<< p.dsdT << std::endl;
	      std::cout << "exact : " << exact.n << " " << exact.dndmu << " " 
			<< exact.dndT << " " << exact.dsdT << std::endl;
	      std::cout << "bad   : " << bad.n << " " << bad.dndmu << " " 
			<< bad.dndT << " " << bad.dsdT << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
	    }
	    
	    // End of loop over points in data file
	  }
	  // End of temperature loop
	}

	dev.n/=cnt;
	dev.dndmu/=cnt;
	dev.dndT/=cnt;
	dev.dsdT/=cnt;

	if (verbose>0) {
	  if (k==0) {
	    std::cout << "Function calc_mu(), include rest mass:" << std::endl;
	  } else if (k==1) {
	    std::cout << "Function calc_mu(), without rest mass:" << std::endl;
	  } else if (k==2) {
	    std::cout << "Function calc_mu(), include rest mass, "
		      << "interacting:" << std::endl;
	  } else {
	    std::cout << "Function calc_mu(), without rest mass, "
		      << "interacting:" << std::endl;
	  }

	  std::cout << "Average performance: " << std::endl;
	  std::cout << "dndmu: " << dev.dndmu << " dndT: " 
		    << dev.dndT << " dsdT: " << dev.dsdT << std::endl;
	  std::cout << "Worst case: " << std::endl;
	  std::cout << "dndmu: " << bad.dndmu << " dndT: " 
		    << bad.dndT << " dsdT: " << bad.dsdT << std::endl;
	  std::cout << "mu: " << mu_bad << " m: " << m_bad << " T: " << T_bad 
		    << " mot: " << mot_bad
		    << "\n\tpsi: " << psi_bad << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
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

      // k=0,2 are with rest mass, k=1,3 are without
      // k=0,1 are non-interacting, k=2,3 are interacting
      for(size_t k=0;k<4;k++) {

	double ret_local=0.0;

	// Initialize storage
	dev.mu=0.0; dev.dndmu=0.0; dev.dndT=0.0; dev.dsdT=0.0;
	bad.mu=0.0; bad.dndmu=0.0; bad.dndT=0.0; bad.dsdT=0.0;
    
	// Temperature loop
	for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	  // Loop over each point in the data file
	  for(size_t i=0;i<tab.get_nlines();i++) {
	
	    double mot=tab.get("mot",i);
	    double psi=tab.get("psi",i);
	    p.n=tab.get("n",i);	
	    exact.dndmu=tab.get("dndmu",i);
	    exact.dndT=tab.get("dndT",i);
	    exact.dsdT=tab.get("dsdT",i);

	    set_mass_flags(p,mot,T,k);
	    set_chem_pot(exact,psi,T,k,nr_mode);

	    p.n*=pow(T,3.0);
	    exact.dndmu*=pow(T,2.0);
	    exact.dndT*=pow(T,2.0);
	    exact.dsdT*=pow(T,2.0);

	    // Give it a guess for the chemical potential
	    if (k>=2) {
	      p.nu=p.m;
	    } else {
	      p.mu=p.m;
	    }

	    th.calc_density(p,T);
	
	    dev.dndmu+=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	    dev.dndT+=fabs((p.dndT-exact.dndT)/exact.dndT);
	    dev.dsdT+=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	
	    cnt++;
	    
	    check_derivs<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,m_bad,
				 T_bad,mot_bad,psi_bad,ret_local);

	    if (verbose>1) {
	      std::cout.precision(5);
	      if (k>=2) {
		std::cout << "T,ms,n,psi,mot: " << T << " "
			  << p.ms << " " << p.n << " " 
			  << psi << " " << mot << std::endl;
	      } else {
		std::cout << "T,m,n,psi,mot: " << T << " " << p.m << " "
			  << p.n << " " << psi << " " << mot << std::endl;
	      }
	      std::cout.precision(6);
	      std::cout << "mu,dndmu,dndT,dsdT: " << std::endl;
	      std::cout << "approx: " << p.mu << " "
			<< p.dndmu << " " << p.dndT << " " 
			<< p.dsdT << std::endl;
	      std::cout << "exact : " << exact.mu << " "
			<< exact.dndmu << " " << exact.dndT << " "
			<< exact.dsdT << std::endl;
	      std::cout << "bad   : " << bad.mu << " " << bad.dndmu << " " 
			<< bad.dndT << " " << bad.dsdT << std::endl;
	      std::cout << std::endl;
	      if (verbose>2) {
		char ch;
		std::cin >> ch;
	      }
	    }

	    if (ret_local>ret) {
	      ret=ret_local;
	    }

	    // End of loop over points in data file
	  }
	  // End of temperature loop
	}

	dev.mu/=cnt;
	dev.dndmu/=cnt;
	dev.dndT/=cnt;
	dev.dsdT/=cnt;

	if (verbose>0) {
	  if (k==0) {
	    std::cout << "Function calc_density(), include rest mass:"
		      << std::endl;
	  } else if (k==1) {
	    std::cout << "Function calc_density(), without rest mass:"
		      << std::endl;
	  } else if (k==2) {
	    std::cout << "Function calc_density(), include "
		      << "rest mass, interacting:"
		      << std::endl;
	  } else {
	    std::cout << "Function calc_density(), without rest mass, "
		      << "interacting:"
		      << std::endl;
	  }

	  std::cout << "Average performance: " << std::endl;
	  std::cout << "dndmu: " << dev.dndmu << " dndT: " 
		    << dev.dndT << " dsdT: " << dev.dsdT << std::endl;
	  std::cout << "Worst case: " << std::endl;
	  std::cout << "dndmu: " << bad.dndmu << " dndT: " 
		    << bad.dndT << " dsdT: " << bad.dsdT << std::endl;
	  std::cout << "mu: " << mu_bad << " m: " << m_bad
		    << " T: " << T_bad << " mot: " << mot_bad
		    << "\n\tpsi: " << psi_bad << std::endl;
	  std::cout << std::endl;
	  if (verbose>2) {
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

	  double ret_local=0.0;

	  // Initialize storage
	  dev.n=0.0; dev.dndmu=0.0; dev.dndT=0.0; dev.dsdT=0.0;
	  bad.n=0.0; bad.dndmu=0.0; bad.dndT=0.0; bad.dsdT=0.0;
    
	  // Temperature loop
	  for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	    // Loop over each point in the data file
	    for(size_t i=0;i<tab.get_nlines();i++) {
	
	      double mot=tab.get("mot",i);
	      double psi=tab.get("psi",i);
	      exact.n=tab.get("pair_n",i);
	      exact.dndmu=tab.get("pair_dndmu",i);
	      exact.dndT=tab.get("pair_dndT",i);
	      exact.dsdT=tab.get("pair_dsdT",i);
      
	      set_mass_flags(p,mot,T,k);
	      set_chem_pot(p,psi,T,k,nr_mode);
	
	      th.pair_mu(p,T);
	
	      exact.n*=pow(T,3.0);
	      exact.dndmu*=pow(T,2.0);
	      exact.dndT*=pow(T,2.0);
	      exact.dsdT*=pow(T,2.0);
	
	      dev.dndmu+=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	      dev.dndT+=fabs((p.dndT-exact.dndT)/exact.dndT);
	      dev.dsdT+=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	
	      cnt++;
	      
	      check_derivs<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,m_bad,
				   T_bad,mot_bad,psi_bad,ret_local);

	      if (verbose>1) {
		std::cout.precision(5);
		std::cout << "T,m,mu,psi,mot: " << T << " " << p.m
			  << " " << p.mu
			  << " " << psi << " " << mot << std::endl;
		std::cout.precision(6);
		std::cout << "n,dndmu,dndT,dsdT: " << std::endl;
		std::cout << "approx: " << p.n << " " << p.dndmu << " "
			  << p.dndT << " " 
			  << p.dsdT << std::endl;
		std::cout << "exact : " << exact.n << " "
			  << exact.dndmu << " " << exact.dndT << " "
			  << exact.dsdT << std::endl;
		std::cout << "bad   : " << bad.n << " " << bad.dndmu << " " 
			  << bad.dndT << " " << bad.dsdT << std::endl;
		std::cout << std::endl;
		if (verbose>2) {
		  char ch;
		  std::cin >> ch;
		}
	      }

	      if (ret_local>ret) {
		ret=ret_local;
	      }

	      // End of loop over points in data file
	    }
	    // End of temperature loop
	  }

	  dev.n/=cnt;
	  dev.dndmu/=cnt;
	  dev.dndT/=cnt;
	  dev.dsdT/=cnt;

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
	    std::cout << "dndmu: " << dev.dndmu << " dndT: " 
		      << dev.dndT << " dsdT: " << dev.dsdT << std::endl;
	    std::cout << "Worst case: " << std::endl;
	    std::cout << "dndmu: " << bad.dndmu << " dndT: " 
		      << bad.dndT << " dsdT: " << bad.dsdT << std::endl;
	    std::cout << "mu: " << mu_bad << " m: " << m_bad
		      << " T: " << T_bad << " mot: " << mot_bad
		      << "\n\tpsi: " << psi_bad << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
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

	  double ret_local=0.0;

	  // Initialize storage
	  dev.mu=0.0; dev.dndmu=0.0; dev.dndT=0.0; dev.dsdT=0.0;
	  bad.mu=0.0; bad.dndmu=0.0; bad.dndT=0.0; bad.dsdT=0.0;
    
	  // Temperature loop
	  for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {
      
	    // Loop over each point in the data file
	    for(size_t i=0;i<tab.get_nlines();i++) {
	
	      double mot=tab.get("mot",i);
	      double psi=tab.get("psi",i);
	      p.n=tab.get("pair_n",i);	
	      exact.dndmu=tab.get("pair_dndmu",i);
	      exact.dndT=tab.get("pair_dndT",i);
	      exact.dsdT=tab.get("pair_dsdT",i);
	      
	      set_mass_flags(p,mot,T,k);
	      set_chem_pot(exact,psi,T,k,nr_mode);

	      p.n*=pow(T,3.0);
	      exact.dndmu*=pow(T,2.0);
	      exact.dndT*=pow(T,2.0);
	      exact.dsdT*=pow(T,2.0);

	      // Give it a guess for the chemical potential
	      p.mu=p.m;

	      th.pair_density(p,T);

	      dev.dndmu+=fabs((p.dndmu-exact.dndmu)/exact.dndmu);
	      dev.dndT+=fabs((p.dndT-exact.dndT)/exact.dndT);
	      dev.dsdT+=fabs((p.dsdT-exact.dsdT)/exact.dsdT);
	
	      cnt++;
	      
	      check_derivs<part_t>(p,exact,bad,k,T,mot,psi,mu_bad,m_bad,
				   T_bad,mot_bad,psi_bad,ret_local);

	      if (verbose>1) {
		std::cout.precision(5);
		std::cout << "T,m,n,psi,mot: " << T << " " << p.m << " "
			  << p.n << " " << psi << " " << mot << std::endl;
		std::cout.precision(6);
		std::cout << "mu,dndmu,dndT,dsdT: " << std::endl;
		std::cout << "approx: " << p.mu << " " << p.dndmu << " "
			  << p.dndT << " " << p.dsdT << std::endl;
		std::cout << "exact : " << exact.mu << " "
			  << exact.dndmu << " " << exact.dndT << " "
			  << exact.dsdT << std::endl;
		std::cout << "bad   : " << bad.mu << " " << bad.dndmu << " " 
			  << bad.dndT << " " << bad.dsdT << std::endl;
		std::cout << std::endl;
		if (verbose>2) {
		  char ch;
		  std::cin >> ch;
		}
	      }

	      if (ret_local>ret) {
		ret=ret_local;
	      }

	      // End of loop over points in data file
	    }
	    // End of temperature loop
	  }

	  dev.mu/=cnt;
	  dev.dndmu/=cnt;
	  dev.dndT/=cnt;
	  dev.dsdT/=cnt;

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
	    std::cout << "dndmu: "
		      << dev.dndmu << " dndT: " 
		      << dev.dndT << " dsdT: " << dev.dsdT << std::endl;
	    std::cout << "Worst case: " << std::endl;
	    std::cout << "dndmu: " << bad.dndmu
		      << " dndT: " << bad.dndT
		      << " dsdT: " << bad.dsdT << std::endl;
	    std::cout << "mu: " << mu_bad << " m: " << m_bad
		      << " T: " << T_bad << " mot: " << mot_bad
		      << "\n\tpsi: " << psi_bad << std::endl;
	    std::cout << std::endl;
	    if (verbose>2) {
	      char ch;
	      std::cin >> ch;
	    }
	  }

	  // End of k loop
	}

	// End of 'if (test_pair)'
      }

      // ----------------------------------------------------------------
      // Return to the original values 

      p=orig;
  
      return ret;
    }

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
