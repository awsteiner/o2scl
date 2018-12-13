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
#ifndef O2SCL_FERMION_DERIV_NR_H
#define O2SCL_FERMION_DERIV_NR_H

/** \file fermion_deriv_nr.h
    \brief File defining \ref o2scl::fermion_deriv_nr
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/root_cern.h>

#include <o2scl/part_deriv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a nonrelativistic fermion

      This does not include the rest mass energy in the chemical 
      potential or the rest mass energy density in the energy density
      to alleviate numerical precision problems at low densities

      This implements an equation of state for a nonrelativistic fermion
      using direct integration. After subtracting the rest mass from
      the chemical potentials, the distribution function is
      \f[
      \left\{1+\exp\left[\left(\frac{k^2}
      {2 m^{*}}-\nu\right)/T\right]\right\}^{-1}
      \f]
      where \f$ \nu \f$ is the effective chemical potential, \f$ m \f$ is
      the rest mass, and \f$ m^{*} \f$ is the effective mass.
      For later use, we define \f$ E^{*} = k^2/2/m^{*} \f$ .

      Uncertainties are given in \ref unc.

      \b Evaluation \b of \b the \b derivatives

      The relevant derivatives of the distribution function are
      \f[
      \frac{\partial f}{\partial T}=
      f(1-f)\frac{E^{*}-\nu}{T^2}
      \f]
      \f[
      \frac{\partial f}{\partial \nu}=
      f(1-f)\frac{1}{T}
      \f]
      \f[
      \frac{\partial f}{\partial k}=
      -f(1-f)\frac{k}{m^{*} T}
      \f]
      \f[
      \frac{\partial f}{\partial m^{*}}=
      f(1-f)\frac{k^2}{2 m^{*2} T}
      \f]

      We also need the derivative of the entropy integrand w.r.t. the 
      distribution function, which is quite simple
      \f[
      {\cal S}\equiv f \ln f +(1-f) \ln (1-f) \qquad
      \frac{\partial {\cal S}}{\partial f} = \ln 
      \left(\frac{f}{1-f}\right) = 
      \left(\frac{\nu-E^{*}}{T}\right)
      \f]
      where the entropy density is
      \f[
      s = - \frac{g}{2 \pi^2} \int_0^{\infty} {\cal S} k^2 d k
      \f]

      The derivatives can be integrated directly
      or they may be converted to integrals
      over the distribution function through an integration by parts
      \f[
      \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) g(k)\right|_{k=a}^{k=b}
      - \int_a^b g(k) \frac{d f(k)}{dk} dk 
      \f]
      using the distribution function for \f$ f(k) \f$ and 0 and \f$
      \infty \f$ as the limits, we have
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{d g(k)}{dk} f dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} g(k) f (1-f) \frac{k}{E^{*} T} dk 
      \f]
      as long as \f$ g(k) \f$ vanishes at \f$ k=0 \f$ .
      Rewriting,
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} h(k) f (1-f) dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} f \frac{T m^{*}}{k} 
      \left[ h^{\prime} - \frac{h}{k}\right] d k
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
      m^{*} f dk
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-\nu)}{T^2} 
      f (1-f) dk
      \f]
      Using \f$ h(k)=k^2(E^{*}-\nu)/T^2 \f$ we get
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
      \left[m^{*} \left(E^{*}-\nu\right) -k^2\right] d k
      \f]

      3) The derivative of the entropy wrt the chemical potential
      \f[
      \left(\frac{d s}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-\nu)}{T^2} dk
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
      \frac{(E^{*}-\nu)^2}{T^3} dk
      \f]
      Using \f$ h(k)=k^2 (E^{*}-\nu)^2/T^3 \f$ 
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      f \frac{m^{*}}{T^2} \left[\left( E^{*}-\nu \right)^2 
      +\frac{2 k^2}{m^{*}} \left(E^{*}-\nu\right)\right] d k
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      \frac{k^2}{2 m^{* 2} T} f (1-f) k^2 dk
      \f]
      Using \f$ h(k)=k^4/(2 m^{* 2} T) \f$ we get
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} f
      \frac{3 k^2}{2 m^{*}} d k 
      \f]

      <b>Conversion to unitless variables:</b>

      After integrating by parts 
      \f$ u = k^2/2/m^{*}/T \f$ and \f$ y=\mu/T \f$, so
      \f[
      k d k = m^{*} T d u 
      \f]
      or
      \f[
      d k = \frac{m^{*} T}{\sqrt{2 m^{*} T u}} d u =
      \sqrt{\frac{m^{*} T}{2 u}} d u 
      \f]
      
      1) The derivative of the density wrt the chemical potential
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g m^{* 3/2} \sqrt{T}}{2^{3/2} \pi^2} \int_0^{\infty} 
      u^{-1/2} f d u
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g m^{* 3/2} \sqrt{T}}
      {2^{3/2} \pi^2} \int_0^{\infty} f d u
      \left[ 3 u^{1/2} - y u^{-1/2}\right]
      \f]

      4) The derivative of the entropy wrt the temperature
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g m^{* 3/2} T^{1/2}}{2^{3/2} \pi^2} \int_0^{\infty} 
      f \left[ 5 u^{3/2} - 6 y u^{1/2} + y^2 u^{-1/2}\right] d u
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{3 g m{* 1/2} T^{3/2}}{2^{3/2} \pi^2} 
      \int_0^{\infty} u^{1/2} f d u
      \f]

      \todo Make a calc_density_zerot function
      \todo There is a lot of code duplication with fermion_nonrel
      which needs to be fixed.
  */
  class fermion_deriv_nr : public fermion_deriv_thermo {

  public:

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_deriv_nr();

    virtual ~fermion_deriv_nr();
    
    /** \brief The limit for the Fermi functions (default 20.0)
    
	fermion_deriv_nr will ignore corrections smaller than about
	 \f$ \exp(-\mathrm{f{l}imit}) \f$ . 
     */
    double flimit;
    
    /// Storage for the most recently calculated uncertainties 
    fermion_deriv unc;
    
    /** \brief Calculate properties as function of chemical potential
    */
    virtual int calc_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties as function of density
     */
    virtual int calc_density(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual int pair_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
     */
    virtual int pair_density(fermion_deriv &f, double temper);

    /// Calculate effective chemical potential from density
    virtual int nu_from_n(fermion_deriv &f, double temper);

    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_root(root<> &rp) {
      density_root=&rp;
      return;
    }

    /// The default solver for npen_density() and pair_density()
    root_cern<> def_density_root;

    /// Return string denoting type ("fermion_deriv_nr")
    virtual const char *type() { return "fermion_deriv_nr"; };

  protected:

#ifndef DOXYGEN_INTERNAL

    /// Solver to compute chemical potential from density
    root<> *density_root;
    
    /// Function to compute chemical potential from density
    double solve_fun(double x, double nog, double msT);

    /** \brief Function to compute chemical potential from density
	when antiparticles are included
     */
    double pair_fun(double x, fermion_deriv &f, double T);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
