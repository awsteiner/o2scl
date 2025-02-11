/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner

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
/** \file constants.h
    \brief File defining numerical constants
*/
#ifndef O2SCL_CONSTANTS_H
#define O2SCL_CONSTANTS_H

#include <cmath>
#include <iostream>

#include <o2scl/string_conv.h>

#define BOOST_DISABLE_ASSERTS
#include <boost/math/constants/constants.hpp>

/** \brief Constants
    
    \verbatim embed:rst
    Previous versions contained constants from [Luzum11]_ and [Mohr12]_ .
    \endverbatim

    CODATA 2022 values are from physics.nist.gov/constants. IAU 2015
    values are the nominal values from arXiv:1510.07674 and
    arXiv:1605.09788 .
*/
namespace o2scl_const {

  /// \name Unit prefixes
  //@{
  const double quetta=1e30;
  const double ronna=1e27;
  const double yotta=1e24;
  const double zetta=1e21;
  const double exa=1e18;
  const double peta=1e15;
  const double tera=1e12;
  const double giga=1e9;
  const double mega=1e6;
  const double kilo=1e3;
  const double milli=1e-3;
  const double micro=1e-6;
  const double nano=1e-9;
  const double pico=1e-12;
  const double femto=1e-15;
  const double atto=1e-18;
  const double zepto=1e-21;
  const double yocto=1e-24;
  const double ronto=1e-27;
  const double quecto=1e-30;
  //@}

  /// \name Unit systems
  //@{
  /// MKS units
  static const size_t o2scl_mks=1;
  /// CGS units
  static const size_t o2scl_cgs=2;
  //@}

  /// \name Mathematical constants
  //@{
  /// \f$ \pi \f$ 
  template<class fp_t> fp_t pi_f() {
    return boost::math::constants::pi<fp_t>();
  }
  /// \f$ \pi^2 \f$ 
  template<class fp_t> fp_t pi2_f() {
    return boost::math::constants::pi_sqr<fp_t>();
  }
  /// \f$ \sqrt{\pi} \f$ 
  template<class fp_t> fp_t root_pi_f() {
    return boost::math::constants::root_pi<fp_t>();
  }
  /// \f$ \zeta(2) \f$ 
  template<class fp_t> fp_t zeta2_f() {
    return boost::math::constants::zeta_two<fp_t>();
  }
  /// \f$ \zeta(3) \f$ 
  template<class fp_t> fp_t zeta3_f() {
    return boost::math::constants::zeta_three<fp_t>();
  }
  /// The Euler-Mascheroni constant
  template<class fp_t> fp_t euler_f() {
    return boost::math::constants::euler<fp_t>();
  }
  //@}

  /// \name Physical constants
  //@{
  /** \brief Fine structure constant (CODATA 2022 value)
   */
  template<class fp_t> fp_t fine_structure_f() {
    fp_t num=72973525643;
    fp_t den=10000000000000;
    fp_t ret=num/den;
    return ret;
  }
  
  /** \brief Avogadro's number (exact)
   */
  template<class fp_t> fp_t avogadro_f() {
    fp_t ret=602214076e15;
    return ret;
  }
  
  /** \brief Speed of light (exact) 
   */
  template<class fp_t> fp_t speed_of_light_f(size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      fp_t result=29979245800;
      return result;
    }
    fp_t result=299792458;
    return result;
  }
  
  /** \brief Planck constant (exact)
   */
  template<class fp_t> fp_t planck_f(size_t system=o2scl_mks) {
    fp_t numer=662607015;
    fp_t denom=100000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-27;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-34;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Reduced Planck constant (derived; exact)
   */
  template<class fp_t> fp_t hbar_f(size_t system=o2scl_mks) {
    return planck_f<fp_t>(system)/2/
      boost::math::constants::pi<fp_t>();
  }
  
  /** \brief Reduced Planck's constant times speed of light
      \f$ \hbar c \f$ (derived; exact)
   */
  template<class fp_t> fp_t hbarc_f(size_t system=o2scl_mks) {
    return hbar_f<fp_t>(system)*speed_of_light_f<fp_t>(system);
  }

  /** \brief Boltzmann's constant (exact) */
  template<class fp_t> fp_t boltzmann_f(size_t system=o2scl_mks) {
    fp_t numer=1380649;
    fp_t denom=1000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-16;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-23;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /** \brief Gravitational constant (CODATA 2022 value)
   */
  template<class fp_t> fp_t gravitational_constant_f
  (size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      fp_t numer=667430;
      fp_t denom=10000000000000;
      return numer/denom;
    }
    fp_t numer=667430;
    fp_t denom=10000000000000000;
    return numer/denom;
  }
  
  /// Elementary charge in C (exact)
  template<class fp_t> fp_t elem_charge_f() {
    fp_t numer=1602176634;
    fp_t denom=1000000000;
    fp_t frac=(numer/denom);
    fp_t base=10;
    fp_t exp=-19;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Electron volt (exact)
  template<class fp_t> fp_t electron_volt_f(size_t system=o2scl_mks) {
    fp_t numer=1602176634;
    fp_t denom=1000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-12;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-19;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /** \brief Reduced Planck constant times speed of light 
      in \f$ \mathrm{MeV}~\mathrm{fm} \f$ (derived; exact)
  */
  template<class fp_t> fp_t hc_mev_fm_f() {
    fp_t hbarc=hbar_f<fp_t>()*speed_of_light_f<fp_t>()/
      elem_charge_f<fp_t>()*1e9;
    return hbarc;
  }

  /** \brief Stefan-Boltzmann constant (\f$ \mathrm{kg} / \mathrm{K}^4 
      \mathrm{s}^3 \f$ in MKS and \f$ \mathrm{g} / \mathrm{K}^4 
      \mathrm{s}^3 \f$ in CGS; derived; exact)
  */
  template<class fp_t> fp_t stefan_boltz_cons_f(size_t system=o2scl_mks) {
    return pi2_f<fp_t>()*boltzmann_f<fp_t>(system)*
      boltzmann_f<fp_t>(system)*boltzmann_f<fp_t>(system)*
      boltzmann_f<fp_t>(system)/60/
      hbar_f<fp_t>(system)/hbar_f<fp_t>(system)/hbar_f<fp_t>(system)/
      speed_of_light_f<fp_t>(system)/speed_of_light_f<fp_t>(system);
  }
  
  /// \f$ \sin^2 \theta_W \f$ (unitless; PDG 2024 value)
  template<class fp_t> fp_t sin2_theta_weak_f() {
    fp_t numer=23129;
    fp_t denom=100000;
    fp_t ret=numer/denom;
    return ret;
  }
  
  /** \brief Fermi coupling constant in \f$ \mathrm{GeV}^{-2} \f$
      (PDG 2024 value)
  */
  template<class fp_t> fp_t gfermi_gev2_f() {
    fp_t numer=11663788;
    fp_t denom=1000000000000;
    fp_t ret=numer/denom;
    return ret;
  }

  /** \brief Fermi coupling constant in ( \f$ \mathrm{s}^4 /
      \mathrm{m}^4 \mathrm{kg}^2 \f$ in MKS and \f$ \mathrm{s}^4 /
      \mathrm{cm}^4 \mathrm{g}^2 \f$ in CGS; derived from \ref
      gfermi_gev2_f())
  */
  template<class fp_t> fp_t gfermi_f(size_t system=o2scl_mks) {
    return gfermi_gev2_f<fp_t>()/electron_volt_f<fp_t>(system)/
      electron_volt_f<fp_t>(system)/1e18;
  }
  //@}

  /// \name Astrophysical constants
  //@{
  /** \brief Astronomical unit (IAU 2009 value; now exact)
   */
  template<class fp_t> fp_t astronomical_unit_f
  (size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      return 14959787070000;
    }
    return 149597870700;
  }

  /** \brief Parsec (derived; exact)
   */
  template<class fp_t> fp_t parsec_f
  (size_t system=o2scl_mks) {
    return astronomical_unit_f<fp_t>(system)*648000/pi_f<fp_t>();
  }
  
  /** \brief Acccleration due to gravity ( \f$ \mathrm{m} /
      \mathrm{s}^2 \f$ in MKS and \f$ \mathrm{cm} / \mathrm{s}^2 \f$
      in CGS ) (CODATA 2018; now exact)
  */
  template<class fp_t> fp_t grav_accel_f
  (size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      fp_t numer=980665;
      fp_t denom=1000;
      return numer/denom;
    }
    fp_t numer=980665;
    fp_t denom=100000;
    return numer/denom;
  }

  /** \brief Sidereal year in s (PDG 2024 value)
  */
  template<class fp_t> fp_t sidereal_year_f(size_t system=o2scl_mks) {
    fp_t numer=315581498;
    fp_t denom=10;
    return numer/denom;
  }
  
  /** \brief Tropical year in s (PDG 2024 value)
  */
  template<class fp_t> fp_t tropical_year_f(size_t system=o2scl_mks) {
    fp_t numer=315569251;
    fp_t denom=10;
    return numer/denom;
  }
  
  /// Julian year in s (exact)
  template<class fp_t> fp_t julian_year_f(size_t system=o2scl_mks) {
    fp_t numer=36525;
    fp_t numer2=86400;
    fp_t denom=100;
    fp_t temp=numer/denom;
    return temp*numer2;
  }

  /// Light year in \f$ \mathrm{cm} \f$ (derived; exact)
  template<class fp_t> fp_t light_year_f(size_t system=o2scl_mks) {
    return julian_year_f<fp_t>(system)*speed_of_light_f<fp_t>(system);
  }
  //@}
  
  /// \name Particle masses
  //@{
  /// Neutron mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_neutron_f(size_t system=o2scl_mks) {
    fp_t numer=167492750056;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Proton mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_proton_f(size_t system=o2scl_mks) {
    fp_t numer=167262192595;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Proton mass in AMU (CODATA 2022 value)
  template<class fp_t> fp_t mass_proton_amu_f() {
    fp_t numer=10072764665789;
    fp_t denom=10000000000000;
    fp_t frac=(numer/denom);
    return frac;
  }

  /// Electron mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_electron_f(size_t system=o2scl_mks) {
    fp_t numer=91093837139;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-28;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-31;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Muon mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_muon_f(size_t system=o2scl_mks) {
    fp_t numer=1883531627;
    fp_t denom=1000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-25;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-28;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Tau mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_tau_f(size_t system=o2scl_mks) {
    fp_t numer=316754;
    fp_t denom=100000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  //@}

  /// \name Nuclear masses
  //@{
  /// Deuteron mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_deuteron_f(size_t system=o2scl_mks) {
    fp_t numer=33435837768;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Triton mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_triton_f(size_t system=o2scl_mks) {
    fp_t numer=50073567512;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Helion mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_helion_f(size_t system=o2scl_mks) {
    fp_t numer=50064127862;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Alpha mass (CODATA 2022 value)
  template<class fp_t> fp_t mass_alpha_f(size_t system=o2scl_mks) {
    fp_t numer=66446573450;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Unified atomic mass (CODATA 2018 value)
  template<class fp_t> fp_t unified_atomic_mass_f(size_t system=o2scl_mks) {
    fp_t numer=166053906892;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  //@}

  /// Bohr radius (CODATA 2022 value)
  template<class fp_t> fp_t bohr_radius_f(size_t system=o2scl_mks) {
    fp_t numer=529177210544;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-11;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /** \brief Thomson cross section (MKS in \f$ m^2 \f$ and 
      CGS in \f$ cm^2 \f$ ; CODATA 2022 value)
  */
  template<class fp_t> fp_t thomson_csec_f(size_t system=o2scl_mks) {
    fp_t numer=66524587051;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-25;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-29;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// \name Chemical constants
  //@{
  /** \brief Rydberg constant (\f$ \mathrm{kg} \mathrm{m}^2 /
      \mathrm{s}^2 \f$ in MKS and \f$ \mathrm{g} \mathrm{cm}^2 /
      \mathrm{s}^2 \f$ in CGS; CODATA 2022 value)

      This is referred to as the "Rydberg constant" by GSL and is
      equal to \f$ h c R_{\infty} \f$. The value \f$ R_{\infty} \f$ is
      the inverse Rydberg wavelength (also sometimes referred to as
      the Rydberg constant).
  */
  template<class fp_t> fp_t rydberg_f(size_t system=o2scl_mks) {
    fp_t numer=21798723611030;
    fp_t denom=10000000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-11;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-18;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /** \brief Molar gas constant, \f$ R\f$ , (\f$ \mathrm{kg}
      \mathrm{m}^2 / \mathrm{K} \mathrm{mol} \mathrm{s}^2 \f$ in MKS
      \f$ \mathrm{g} \mathrm{cm}^2 / \mathrm{K} \mathrm{mol}
      \mathrm{s}^2 \f$ in CGS; derived; exact)
  */
  template<class fp_t> fp_t molar_gas_f(size_t system=o2scl_mks) {
    return avogadro_f<fp_t>()*boltzmann_f<fp_t>(system);
  }

  /** \brief Molar volume of ideal gas at standard T and P (\f$
      \mathrm{m}/\mathrm{mol} \f$ in MKS and (\f$
      \mathrm{cm}/\mathrm{mol} \f$ in CGS; CODATA 2022 value)
  */
  template<class fp_t> fp_t std_gas_volume_f(size_t system=o2scl_mks) {
    fp_t numer=2271095464;
    fp_t denom=1000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=4;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-2;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  //@}
  
  /// \name Electromagnetic constants
  //@{
  /** \brief Electron magnetic moment (\f$ \mathrm{J}/\mathrm{T} \f$ 
      in MKS; CODATA 2022 value)

      \note This quantity is defined without a minus sign. 
  */
  template<class fp_t> fp_t electron_mag_mom_f(size_t system=o2scl_mks) {
    fp_t numer=92847646917;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-21;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Proton magnetic moment (\f$ \mathrm{J}/\mathrm{T} \f$ 
      in MKS; CODATA 2022 value)
  */
  template<class fp_t> fp_t proton_mag_mom_f(size_t system=o2scl_mks) {
    fp_t numer=141060679545;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-23;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-26;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Roentgen (\f$ \mathrm{A}~\mathrm{s}/\mathrm{kg} \f$ \f$)
   */
  template<class fp_t> fp_t roentgen_f(size_t system=o2scl_mks) {
    fp_t numer=258;
    fp_t denom=100;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-7;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-4;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Bohr magneton (\f$ \mathrm{J}/\mathrm{T} \f$ 
      in MKS; CODATA 2022 value)
  */
  template<class fp_t> fp_t bohr_magneton_f(size_t system=o2scl_mks) {
    fp_t numer=92740100657;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-21;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Nuclear magneton (\f$ \mathrm{J}/\mathrm{T} \f$ 
      in MKS; CODATA 2022 value)
  */
  template<class fp_t> fp_t nuclear_magneton_f(size_t system=o2scl_mks) {
    fp_t numer=50507837393;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=-24;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-27;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Faraday constant in \f$ \mathrm{A}~\mathrm{s} / \mathrm{mol} \f$
      (derived; exact)
  */
  template<class fp_t> fp_t faraday_f(size_t system=o2scl_mks) {
    fp_t base=o2scl_const::avogadro_f<double>()*
      o2scl_const::electron_volt_f<double>(system);
    if (system==o2scl_cgs) {
      return base/100000000;
    }
    return base;
  }

  /** \brief Vacuum permittivity in \f$ \mathrm{F}/\mathrm{m} \f$
      (CODATA 2022)
   */
  template<class fp_t> fp_t vacuum_permittivity_f(size_t system=o2scl_mks) {
    fp_t numer=88541878188;
    fp_t denom=10000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      O2SCL_ERR("Not supported.",o2scl::exc_efailed);
    }
    fp_t base=10;
    fp_t exp=-12;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /** \brief Vacuum permeability in \f$ \mathrm{kg}~
      \mathrm{m} / \mathrm{A}^2~\mathrm{s}^2 \f$
      (CODATA 2022)
   */
  template<class fp_t> fp_t vacuum_permeability_f(size_t system=o2scl_mks) {
    fp_t numer=125663706127;
    fp_t denom=100000000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      O2SCL_ERR("Not supported.",o2scl::exc_efailed);
    }
    fp_t base=10;
    fp_t exp=-6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Electron charge squared in Gaussian units (derived)

      In Gaussian Units:
      \f{eqnarray*}
      &\vec{\nabla} \cdot \vec{E} = 4 \pi \rho \, ,
      \quad
      \vec{E}=-\vec{\nabla} \Phi \, ,
      \quad
      \nabla^2 \Phi = - 4 \pi \rho \, ,
      &\\&
      F=\frac{q_1 q_2}{r^2} \, ,
      \quad
      W=\frac{1}{2} \int \rho V d^3 x
      =\frac{1}{8 \pi} \int | \vec{E} |^2 d^3 x \, ,
      \quad 
      \alpha=\frac{e^2}{\hbar c}=\frac{1}{137}&
      \f}
  */
  template<class fp_t> fp_t e2_gaussian_f() {
    return o2scl_const::hc_mev_fm_f<fp_t>()*
      o2scl_const::fine_structure_f<fp_t>();
  }

  /** \brief Electron charge sqaured in 
      Heaviside-Lorentz units where \f$\hbar=c=1\f$ (derived)

      In Heaviside-Lorentz units:
      \f{eqnarray*}
      &\vec{\nabla} \cdot \vec{E} = \rho \, ,
      \quad
      \vec{E}=-\vec{\nabla} \Phi \, ,
      \quad
      \nabla^2 \Phi = - \rho \, ,
      &\\&
      F=\frac{q_1 q_2}{4 \pi r^2} \, ,
      \quad
      W=\frac{1}{2} \int \rho V d^3 x
      =\frac{1}{2} \int | \vec{E} |^2 d^3 x \, ,
      \quad
      \alpha=\frac{e^2}{4 \pi}=\frac{1}{137}&
      \f}
  */      
  template<class fp_t> fp_t e2_hlorentz_f() {
    return o2scl_const::fine_structure_f<fp_t>()*4*pi_f<fp_t>();
  }

  /** \brief Electron charge squared in SI(MKS) units (derived)

      In MKS units:
      \f{eqnarray*}
      &\vec{\nabla} \cdot \vec{E} = \rho \, ,
      \quad
      \vec{E}=-\vec{\nabla} \Phi \, ,
      \quad
      \nabla^2 \Phi = - \rho \, ,
      &\\&
      F=\frac{1}{4 \pi \varepsilon_0}\frac{q_1 q_2}{r^2} \, ,
      \quad
      W=\frac{1}{2} \int \rho V d^3 x
      =\frac{\varepsilon_0}{2} \int | \vec{E} |^2 d^3 x \, ,
      \quad
      \alpha=\frac{e^2}{4 \pi \varepsilon_0 \hbar c}=\frac{1}{137}&
      \f}

      Note the conversion formulas
      \f[
      q_HL=\sqrt{4 \pi} q_G = \frac{1}{\sqrt{\varepsilon_0}} q_{SI}
      \f]
      as mentioned, e.g. in pg. 13 of D. Griffiths Intro to Elem. 
      Particles.
  */      
  template<class fp_t> fp_t elem_charge_squared_f() {
    return elem_charge_f<fp_t>()*elem_charge_f<fp_t>();
  }

  /** \brief 1 \f$\mathrm{Gauss}\f$ times the electron charge 
      in Gaussian units in \f$\mathrm{fm}^{-2}\f$
  */
  template<class fp_t> fp_t ec_gauss_fm2_f() {
    fp_t base=10;
    fp_t exp=-34;
    fp_t powt=pow(base,exp);
    return elem_charge_f<fp_t>()*powt/hbar_f<fp_t>(o2scl_mks);
  }

  /** \brief Conversion factor from \f$ \mathrm{Gauss}^2 \f$ to
      \f$\mathrm{fm}^{-4}\f$ in Gaussian units.

      This is useful, e.g. in converting magnetic field squared
      to an energy density.
  */
  template<class fp_t> fp_t gauss2_fm4_f() {
    return ec_gauss_fm2_f<fp_t>()*ec_gauss_fm2_f<fp_t>()/
      o2scl_const::fine_structure_f<fp_t>();
  }
  //@}
  
  /// \name Time measurements (s)
  //@{
  /// Minute
  template<class fp_t> fp_t minute_f(size_t system=o2scl_mks) {
    return 60;
  }

  /// Hour
  template<class fp_t> fp_t hour_f(size_t system=o2scl_mks) {
    return 3600;
  }

  /// Day
  template<class fp_t> fp_t day_f(size_t system=o2scl_mks) {
    return 24*3600;
  }

  /// Week
  template<class fp_t> fp_t week_f(size_t system=o2scl_mks) {
    return 168*3600;
  }
  //@}

  /// \name Length measurements (MKS in m and CGS in cm)
  //@{
  /// Inch
  template<class fp_t> fp_t inch_f(size_t system=o2scl_mks) {
    fp_t numer=254;
    if (system==o2scl_cgs) {
      fp_t denom=100;
      return numer/denom;
    }
    fp_t denom=10000;
    return numer/denom;
  }

  /// Foot
  template<class fp_t> fp_t foot_f(size_t system=o2scl_mks) {
    fp_t numer=3048;
    if (system==o2scl_cgs) {
      fp_t denom=100;
      return numer/denom;
    }
    fp_t denom=10000;
    return numer/denom;
  }

  /// Yard
  template<class fp_t> fp_t yard_f(size_t system=o2scl_mks) {
    fp_t numer=9144;
    if (system==o2scl_cgs) {
      fp_t denom=100;
      return numer/denom;
    }
    fp_t denom=10000;
    return numer/denom;
  }

  /// Mile
  template<class fp_t> fp_t mile_f(size_t system=o2scl_mks) {
    fp_t numer=1609344;
    if (system==o2scl_cgs) {
      fp_t denom=10;
      return numer/denom;
    }
    fp_t denom=1000;
    return numer/denom;
  }

  /// Nautical mile
  template<class fp_t> fp_t nautical_mile_f(size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      return 185200;
    }
    return 1852;
  }

  /// Fathom
  template<class fp_t> fp_t fathom_f(size_t system=o2scl_mks) {
    fp_t numer=18288;
    if (system==o2scl_cgs) {
      fp_t denom=100;
      return numer/denom;
    }
    fp_t denom=10000;
    return numer/denom;
  }

  /// Mil
  template<class fp_t> fp_t mil_f(size_t system=o2scl_mks) {
    fp_t numer=254;
    if (system==o2scl_cgs) {
      fp_t denom=100000;
      return numer/denom;
    }
    fp_t denom=10000000;
    return numer/denom;
  }

  /// Point
  template<class fp_t> fp_t point_f(size_t system=o2scl_mks) {
    fp_t numer=3175;
    if (system==o2scl_cgs) {
      fp_t denom=90000;
      return numer/denom;
    }
    fp_t denom=9000000;
    return numer/denom;
  }

  /// Texpoint
  template<class fp_t> fp_t texpoint_f(size_t system=o2scl_mks) {
    fp_t numer=351459803515;
    if (system==o2scl_cgs) {
      fp_t denom=10000000000000;
      return numer/denom;
    }
    fp_t denom=1000000000000000;
    return numer/denom;
  }

  /// Micron
  template<class fp_t> fp_t micron_f(size_t system=o2scl_mks) {
    fp_t numer=1;
    if (system==o2scl_cgs) {
      fp_t denom=10000;
      return numer/denom;
    }
    fp_t denom=1000000;
    return numer/denom;
  }

  /// Angstrom
  template<class fp_t> fp_t angstrom_f(size_t system=o2scl_mks) {
    fp_t numer=1;
    if (system==o2scl_cgs) {
      fp_t denom=100000000;
      return numer/denom;
    }
    fp_t denom=10000000000;
    return numer/denom;
  }
  //@}

  /// \name Particle masses from PDG 2024
  //@{
  /** \brief \f$ \Lambda \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
  */
  template<class fp_t> fp_t mass_lambda_MeV_f() {
    fp_t numer=1115683;
    fp_t denom=1000;
    return numer/denom;
  }

  /** \brief \f$ \Sigma^{-} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
   */
  template<class fp_t> fp_t mass_sigma_minus_MeV_f() {
    fp_t numer=1197449;
    fp_t denom=1000;
    return numer/denom;
  }

  /** \brief \f$ \Sigma^{0} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
   */
  template<class fp_t> fp_t mass_sigma_zero_MeV_f() {
    fp_t numer=1192642;
    fp_t denom=1000;
    return numer/denom;
  }

  /** \brief \f$ \Sigma^{+} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
   */
  template<class fp_t> fp_t mass_sigma_plus_MeV_f() {
    fp_t numer=118937;
    fp_t denom=100;
    return numer/denom;
  }

  /** \brief \f$ \Xi^{0} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
   */
  template<class fp_t> fp_t mass_cascade_zero_MeV_f() {
    fp_t numer=131486;
    fp_t denom=100;
    return numer/denom;
  }
  
  /** \brief \f$ \Xi^{-} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR FIT")
   */
  template<class fp_t> fp_t mass_cascade_minus_MeV_f() {
    fp_t numer=132171;
    fp_t denom=100;
    return numer/denom;
  }
  
  /** \brief Up quark mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR EVALUATION")
   */
  template<class fp_t> fp_t mass_up_MeV_f() {
    fp_t numer=216;
    fp_t denom=100;
    return numer/denom;
  }
  
  /** \brief Down quark mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR EVALUATION")
   */
  template<class fp_t> fp_t mass_down_MeV_f() {
    fp_t numer=470;
    fp_t denom=100;
    return numer/denom;
  }
  
  /** \brief Strange quark mass in \f$ \mathrm{MeV} \f$
      (used value labeled "OUR EVALUATION")
   */
  template<class fp_t> fp_t mass_strange_MeV_f() {
    fp_t numer=935;
    fp_t denom=10;
    return numer/denom;
  }
  //@}

  /// \name Solar system properties
  //@{
  /** \brief Solar mass times gravitational constant (\f$ \mathrm{m}^3
      / \mathrm{s}^2 \f$ in MKS and \f$ \mathrm{cm}^3 / \mathrm{s}^2
      \f$ in CGS; IAU 2015 value, see https://arxiv.org/abs/1510.07674)

      Note that this value differs slightly in Barycentric Coordinate
      Time and Barycentric Dynamical Time. This is the IAU's nominal
      value.
  */
  template<class fp_t> fp_t solar_mass_param_f(size_t system=o2scl_mks) {
    fp_t numer=13271244;
    fp_t denom=10000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=26;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=20;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Mass of the sun in g (derived)
  template<class fp_t> fp_t solar_mass_f(size_t system=o2scl_mks) {
    return solar_mass_param_f<fp_t>(system)/
      gravitational_constant_f<fp_t>(system);
  }
  
  /// Schwarzchild radius (MKS in m and CGS in cm) (derived)
  template<class fp_t> fp_t schwarzchild_radius_f
  (size_t system=o2scl_mks) {
    return 2*solar_mass_param_f<fp_t>(system)/
      speed_of_light_f<fp_t>(system)/speed_of_light_f<fp_t>(system);
  }
  
  /// Radius of the sun (MKS in m and CGS in cm) (IAU 2015 nominal value)
  template<class fp_t> fp_t solar_radius_f(size_t system=o2scl_mks) {
    fp_t numer=6957;
    fp_t denom=1000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=10;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=8;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Temperature of the sun's photosphere in K (IAU 2015 nominal value)
  template<class fp_t> fp_t solar_temperature_f(size_t system=o2scl_mks) {
    return 5772;
  }
  
  /// Luminosity of sun in erg/s (IAU 2015 nominal value)
  template<class fp_t> fp_t solar_luminosity_f(size_t system=o2scl_mks) {
    fp_t numer=3828;
    fp_t denom=1000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=40;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=33;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Solar total irradiance in W/m^2 (IAU 2015 nominal value)
  template<class fp_t> fp_t solar_irradiance_f() {
    fp_t numer=1361;
    return numer;
  }
  
  /** \brief Earth mass times gravitational constant (\f$ \mathrm{m}^3
      / \mathrm{s}^2 \f$ in MKS and \f$ \mathrm{cm}^3 / \mathrm{s}^2
      \f$ in CGS; IAU 2015 nominal values)
  */
  template<class fp_t> fp_t earth_mass_param_f(size_t system=o2scl_mks) {
    fp_t numer=3986004;
    fp_t denom=1000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=20;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=14;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Mass of the earth in g (derived)
  template<class fp_t> fp_t earth_mass_f(size_t system=o2scl_mks) {
    return earth_mass_param_f<fp_t>(system)/
      gravitational_constant_f<fp_t>(system);
  }
  
  /// Equatorial radius of earth (MKS in m and CGS in cm) (IAU 2015 value)
  template<class fp_t> fp_t earth_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=63781;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of earth (MKS in m and CGS in cm) (IAU 2015 value)
  template<class fp_t> fp_t earth_radius_pol_f(size_t system=o2scl_mks) {
    fp_t numer=63568;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /** \brief Jupter mass times gravitational constant (\f$
      \mathrm{m}^3 / \mathrm{s}^2 \f$ in MKS and \f$ \mathrm{cm}^3 /
      \mathrm{s}^2 \f$ in CGS; IAU 2015 nominal values)
  */
  template<class fp_t> fp_t jupiter_mass_param_f
  (size_t system=o2scl_mks) {
    fp_t numer=12668653;
    fp_t denom=10000000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=23;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=17;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Mass of jupiter in g (derived)
  template<class fp_t> fp_t jupiter_mass_f(size_t system=o2scl_mks) {
    return jupiter_mass_param_f<fp_t>(system)/
      gravitational_constant_f<fp_t>(system);
  }
  
  /// Equatorial radius of jupiter (MKS in m and CGS in cm) (IAU 2015 value)
  template<class fp_t> fp_t jupiter_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=71492;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of jupiter (MKS in m and CGS in cm) (IAU 2015 value)
  template<class fp_t> fp_t jupiter_radius_pol_f
  (size_t system=o2scl_mks) {
    fp_t numer=66854;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Mass of mercury in g
  template<class fp_t> fp_t mercury_mass_f(size_t system=o2scl_mks) {
    fp_t numer=33011;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=26;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Radius of mercury (MKS in m and CGS in cm)
  template<class fp_t> fp_t mercury_radius_f(size_t system=o2scl_mks) {
    fp_t numer=24397;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of venus in g
  template<class fp_t> fp_t venus_mass_f(size_t system=o2scl_mks) {
    fp_t numer=78675;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=27;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Radius of venus (MKS in m and CGS in cm)
  template<class fp_t> fp_t venus_radius_f(size_t system=o2scl_mks) {
    fp_t numer=60518;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of mars in g
  template<class fp_t> fp_t mars_mass_f(size_t system=o2scl_mks) {
    fp_t numer=64171;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=26;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Equatorial radius of mars (MKS in m and CGS in cm)
  template<class fp_t> fp_t mars_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=33962;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of mars (MKS in m and CGS in cm)
  template<class fp_t> fp_t mars_radius_pol_f(size_t system=o2scl_mks) {
    fp_t numer=33762;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of saturn in g
  template<class fp_t> fp_t saturn_mass_f(size_t system=o2scl_mks) {
    fp_t numer=56834;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=29;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Equatorial radius of saturn (MKS in m and CGS in cm) 
  template<class fp_t> fp_t saturn_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=60268;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of saturn (MKS in m and CGS in cm)
  template<class fp_t> fp_t saturn_radius_pol_f(size_t system=o2scl_mks) {
    fp_t numer=54364;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of uranus in g
  template<class fp_t> fp_t uranus_mass_f(size_t system=o2scl_mks) {
    fp_t numer=86810;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=28;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Equatorial radius of uranus (MKS in m and CGS in cm) 
  template<class fp_t> fp_t uranus_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=25559;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of uranus (MKS in m and CGS in cm)
  template<class fp_t> fp_t uranus_radius_pol_f(size_t system=o2scl_mks) {
    fp_t numer=24973;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of neptune in g
  template<class fp_t> fp_t neptune_mass_f(size_t system=o2scl_mks) {
    fp_t numer=102413;
    fp_t denom=100000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=29;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Equatorial radius of neptune (MKS in m and CGS in cm) 
  template<class fp_t> fp_t neptune_radius_eq_f
  (size_t system=o2scl_mks) {
    fp_t numer=24764;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Polar radius of neptune (MKS in m and CGS in cm)
  template<class fp_t> fp_t neptune_radius_pol_f(size_t system=o2scl_mks) {
    fp_t numer=24341;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=9;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=7;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }

  /// Mass of pluto in g
  template<class fp_t> fp_t pluto_mass_f(size_t system=o2scl_mks) {
    fp_t numer=1303;
    fp_t denom=1000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=25;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=-24;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  
  /// Radius of pluto (MKS in m and CGS in cm) 
  template<class fp_t> fp_t pluto_radius_f(size_t system=o2scl_mks) {
    fp_t numer=11883;
    fp_t denom=10000;
    fp_t frac=(numer/denom);
    if (system==o2scl_cgs) {
      fp_t base=10;
      fp_t exp=8;
      fp_t powt=pow(base,exp);
      fp_t result=frac*powt;
      return result;
    }
    fp_t base=10;
    fp_t exp=6;
    fp_t powt=pow(base,exp);
    fp_t result=frac*powt;
    return result;
  }
  //@}

  /// \name A few double-precision constants for convenience
  //@{
  /// \f$ \pi \f$ 
  const double pi=boost::math::constants::pi<double>();
  /// \f$ \pi^2 \f$ 
  const double pi2=boost::math::constants::pi_sqr<double>();
  /// \f$ \sqrt{\pi} \f$ 
  const double root_pi=boost::math::constants::root_pi<double>();
  /// \f$ \hbar c \f$ in MeV fm (exact)
  const double hc_mev_fm=hc_mev_fm_f<double>();
  /// \f$ \hbar c \f$ in MeV cm (exact)
  const double hc_mev_cm=hc_mev_fm*1.0e-13;
  //@}
  
}

#endif
