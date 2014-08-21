/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

  /** \brief A particle data class with derivatives

      This class adds the derivatives \ref dndmu, \ref dndT, 
      \ref dsdT, and \ref dndm, which correspond to
      \f[
      \left(\frac{d n}{d \mu}\right)_{T,V}, \quad
      \left(\frac{d n}{d T}\right)_{\mu,V}, \quad
      \left(\frac{d s}{d T}\right)_{\mu,V}, \quad \mathrm{and} \quad
      \left(\frac{d n}{d m}\right)_{T,\mu,V},
      \f]
      respectively. All other first-order thermodynamic derivatives
      can be expressed in terms of the first three derivatives. In the
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

      \hline

      <b>Other derivatives with respect 
      to chemical potential and temperature:</b>

      There is a Maxwell relation
      \f[
      \left(\frac{d s}{d \mu}\right)_{T,V} =
      \left(\frac{d n}{d T}\right)_{\mu,V}
      \f]
      The pressure derivatives are trivial
      \f[
      \left(\frac{d P}{d \mu}\right)_{T,V}=n, \quad
      \left(\frac{d P}{d T}\right)_{\mu,V}=s
      \f]
      The energy density derivatives are related through the 
      thermodynamic identity:
      \f[
      \left(\frac{d \varepsilon}{d \mu}\right)_{T,V}=
      \mu \left(\frac{d n}{d \mu}\right)_{T,V}+
      T \left(\frac{d s}{d \mu}\right)_{T,V}
      \f]
      \f[
      \left(\frac{d \varepsilon}{d T}\right)_{\mu,V}=
      \mu \left(\frac{d n}{d T}\right)_{\mu,V}+
      T \left(\frac{d s}{d T}\right)_{\mu,V}
      \f]

      \hline

      <b>Other derivatives:</b>

      Note that the derivative of the entropy with respect to the
      temperature above is not the specific heat per particle, \f$ c_V \f$.
      The specific heat per particle is
      \f[
      c_V = \frac{T}{N} \left( \frac{\partial S}{\partial T} \right)_{V,N}
      = \frac{T}{n} \left( \frac{\partial s}{\partial T} \right)_{V,n}
      \f] 
      As noted in \ref part_section in the User's Guide for \o2p, we
      work in units so that \f$ \hbar = c = k_B = 1 \f$. In this case,
      \f$ c_V \f$ is unitless as defined here. To compute \f$ c_V \f$
      in terms of the derivatives above, note that the
      descendants of part_deriv provide all of the thermodynamic
      functions in terms of \f$ \mu, V \f$ and \f$ T \f$, so we have
      \f[
      s=s(\mu,T,V) \quad \mathrm{and} \quad n=n(\mu,T,V) \, .
      \f]
      We can then construct a function
      \f[
      s=s[\mu(n,T,V),T,V]
      \f]
      and then write the required derivative directly
      \f[
      \left(\frac{\partial s}{\partial T}\right)_{n,V} =
      \left(\frac{\partial s}{\partial \mu}\right)_{T,V}
      \left(\frac{\partial \mu}{\partial T}\right)_{n,V} +
      \left(\frac{\partial s}{\partial T}\right)_{\mu,V} \, .
      \f]
      Now we use the identity
      \f[
      \left(\frac{\partial \mu}{\partial T}\right)_{n,V} = -
      \left(\frac{\partial n}{\partial T}\right)_{\mu,V} 
      \left(\frac{\partial n}{\partial \mu}\right)_{T,V}^{-1} \, ,
      \f]
      and the Maxwell relation above to give
      \f[
      C_V = \frac{T}{n}
      \left[ 
      \left(\frac{\partial s}{\partial T}\right)_{\mu,V}
      -\left(\frac{\partial n}{\partial T}\right)_{\mu,V}^2
      \left(\frac{\partial n}{\partial \mu}\right)_{T,V}^{-1}
      \right]
      \f]
      which expresses the specific heat in terms of the three
      derivatives which are given.

      For, \f$ c_P \f$, defined as
      \f[
      c_P = \frac{T}{N} \left( \frac{\partial S}{\partial T} 
      \right)_{N,P}
      \f] 
      (which is also unitless) we can write functions
      \f[
      S=S(N,T,V) \qquad \mathrm{and} \qquad V=V(N,P,T)
      \f]
      which imply
      \f[
      \left( \frac{\partial S}{\partial T} \right)_{N,P} =
      \left( \frac{\partial S}{\partial T} \right)_{N,V} +
      \left( \frac{\partial S}{\partial V} \right)_{N,T}
      \left( \frac{\partial V}{\partial T} \right)_{N,P} \, .
      \f]
      Thus we require the derivatives
      \f[
      \left( \frac{\partial S}{\partial T} \right)_{N,V} ,
      \left( \frac{\partial S}{\partial V} \right)_{N,T} ,
      \qquad\mathrm{and}\qquad
      \left( \frac{\partial V}{\partial T} \right)_{N,P}
      \, .
      \f]

      To compute the new entropy derivatives, we can write
      \f[
      S=S(\mu(N,T,V),T,V)
      \f]
      to get
      \f[
      \left( \frac{\partial S}{\partial T} \right)_{N,V} =
      \left( \frac{\partial S}{\partial \mu} \right)_{T,V}
      \left( \frac{\partial \mu}{\partial T} \right)_{N,V} +
      \left( \frac{\partial S}{\partial T} \right)_{\mu,V} \, ,
      \f]
      and
      \f[
      \left( \frac{\partial S}{\partial V} \right)_{N,T} =
      \left( \frac{\partial S}{\partial \mu} \right)_{T,V}
      \left( \frac{\partial \mu}{\partial V} \right)_{N,T} +
      \left( \frac{\partial S}{\partial V} \right)_{\mu,T} \, .
      \f]
      These require the chemical potential derivatives which have
      associated Maxwell relations
      \f[
      \left( \frac{\partial \mu}{\partial T} \right)_{N,V} =
      -\left( \frac{\partial S}{\partial N} \right)_{T,V} 
      \qquad\mathrm{and}\qquad
      \left( \frac{\partial \mu}{\partial V} \right)_{N,T} =
      -\left( \frac{\partial P}{\partial N} \right)_{T,V} \, .
      \f]
      Finally, we can rewrite the derivatives on the right hand sides
      in terms of derivatives of functions of \f$ \mu, V \f$ and
      \f$ T \f$,
      \f[
      \left( \frac{\partial S}{\partial N} \right)_{T,V} =
      \left( \frac{\partial S}{\partial \mu} \right)_{T,V} 
      \left( \frac{\partial N}{\partial \mu} \right)_{T,V}^{-1} \, ,
      \f]
      and
      \f[
      \left( \frac{\partial P}{\partial N} \right)_{T,V} =
      \left( \frac{\partial P}{\partial \mu} \right)_{T,V} 
      \left( \frac{\partial N}{\partial \mu} \right)_{T,V}^{-1} \, .
      \f]

      The volume derivative,
      \f[
      \left( \frac{\partial V}{\partial T} \right)_{N,P} \, ,
      \f]
      is related to the coefficient of thermal expansion, sometimes 
      called \f$ \alpha \f$,
      \f[
      \alpha \equiv \frac{1}{V}
      \left( \frac{\partial V}{\partial T} \right)_{N,P} \, .
      \f]
      We can rewrite the derivative 
      \f[
      \left( \frac{\partial V}{\partial T} \right)_{N,P} =
      -\left( \frac{\partial P}{\partial T} \right)_{N,V} 
      \left( \frac{\partial P}{\partial V} \right)_{N,T}^{-1} \, .
      \f]
      The first term can be computed from the Maxwell relation
      \f[
      \left( \frac{\partial P}{\partial T} \right)_{N,V} = 
      \left( \frac{\partial S}{\partial V} \right)_{N,T} \, ,
      \f]
      where the entropy derivative was computed above. The second term
      (related to the inverse of the isothermal compressibility, \f$
      \kappa_T \equiv (-1/V) (\partial V/\partial P)_{T,N} \f$ can be
      computed from the function \f$ P = P[\mu(N,V,T),V,T] \f$
      \f[
      \left( \frac{\partial P}{\partial V} \right)_{N,T} = 
      \left( \frac{\partial P}{\partial \mu} \right)_{T,V} 
      \left( \frac{\partial \mu}{\partial V} \right)_{N,T} +
      \left( \frac{\partial P}{\partial V} \right)_{\mu,T} 
      \f]
      where the chemical potential derivative was computed above.

      The results above can be collected to give
      \f[
      \left( \frac{\partial S}{\partial T} \right)_{N,P} =
      \left( \frac{\partial S}{\partial T} \right)_{\mu,V} +
      \frac{S^2}{N^2}
      \left( \frac{\partial N}{\partial \mu} \right)_{T,V} -
      \frac{2 S}{N}
      \left( \frac{\partial N}{\partial T} \right)_{\mu,V} \, ,
      \f]
      which implies
      \f[
      c_P = 
      \frac{T}{n}
      \left( \frac{\partial s}{\partial T} \right)_{\mu,V} +
      \frac{s^2 T}{n^3}
      \left( \frac{\partial n}{\partial \mu} \right)_{T,V} -
      \frac{2 s T}{n^2}
      \left( \frac{\partial n}{\partial T} \right)_{\mu,V} \, ,
      \f]

      This derivation also gives the well-known relationship between
      the specific heats at constant volume and constant pressure,
      \f[
      c_P = c_V + \frac{T \alpha^2}{n \kappa_T} \, .
      \f]

      \hline

      \future Right now the part fields (n, ed, etc.) are
      initialized to zero but the derivatives aren't. This
      should be more consistent. 
  */
  class part_deriv : public part {
    
  public:
    
    /// Make a particle of mass \c mass and degeneracy \c dof.
  part_deriv(double mass=0.0, double dof=0.0) : part(mass,dof) {
    }

    /// Derivative of number density with respect to chemical potential
    double dndmu;
    
    /// Derivative of number density with respect to temperature
    double dndT;

    /// Derivative of entropy density with respect to temperature
    double dsdT;

    /// Derivative of number density with respect to the effective mass
    double dndm;

  };
  
  /** \brief A fermion with derivative information
   */
  class fermion_deriv : public part_deriv {
    
  public:
    
    /// Make a particle of mass \c mass and degeneracy \c dof.
  fermion_deriv(double mass=0.0, double dof=0.0) : part_deriv(mass,dof) {
    }
    
    /// Fermi momentum
    double kf;
    
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
      in Constantiou et al. 2014 which could be implemented here.
  */
  class fermion_deriv_thermo {

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
    */
    virtual bool calc_mu_deg(fermion_deriv &f, double temper,
			     double prec);
    
  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
