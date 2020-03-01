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
/** \file constants.h
    \brief File defining numerical constants
*/
#ifndef O2SCL_CONSTANTS_H
#define O2SCL_CONSTANTS_H

#include <cmath>

#include <boost/math/constants/constants.hpp>

/** \brief Constants
    
    CODATA 2014 values were from \ref Mohr16. CODATA 2018 values are
    from physics.nist.gov/constants. IAU 2015 values are the nominal
    values from arXiv:1510.07674 and arXiv:1605.09788 .
*/
namespace o2scl_const {

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

  /** \brief Fine structure constant (CODATA 2018 value)
   */
  const double fine_structure=7.2973525693e-3;
  /** \brief Avogadro's number (CODATA 2018 value; exact)
   */
  const double avogadro=6.02214076e23;

  /// \f$ \pi \f$ 
  const double pi=boost::math::constants::pi<double>();
  /// \f$ \pi^2 \f$ 
  const double pi2=boost::math::constants::pi_sqr<double>();
  /// \f$ \zeta(3/2) \f$
  const double zeta32=2.6123753486854883433;
  /// \f$ \zeta(2) \f$
  const double zeta2=boost::math::constants::zeta_two<double>();
  /// \f$ \zeta(5/2) \f$
  const double zeta52=1.3414872572509171798;
  /// \f$ \zeta(3) \f$
  const double zeta3=boost::math::constants::zeta_three<double>();
  /// \f$ \zeta(5) \f$
  const double zeta5=1.0369277551433699263;
  /// \f$ \zeta(7) \f$
  const double zeta7=1.0083492773819228268;
  
  /// \f$ \sin^2 \theta_W \f$ (PDG 2018 value)
  const double sin2_theta_weak=0.23122;
  /** \brief Fermi coupling constant in \f$ \mathrm{GeV}^{-2} \f$
      (CODATA 2018 value)
  */
  const double gfermi_gev2=1.1663787e-5;
  
  /// MKS units
  static const size_t o2scl_mks=1;
  /// CGS units
  static const size_t o2scl_cgs=2;
  
  /// Planck constant
  template<class fp_t> fp_t planck_f(size_t system=o2scl_mks) {
    fp_t numer=662607015;
    fp_t denom=100000000;
    if (system==o2scl_cgs) {
      fp_t result=(numer/denom)*1.0e-27;
      return result;
    }
    fp_t result=(numer/denom)*1.0e-34;
    return result;
  }

  /// Reduced Planck constant
  template<class fp_t> fp_t hbar_f(size_t system=o2scl_mks) {
    return planck_f<fp_t>(system)/2/
      boost::math::constants::pi<fp_t>();
  }

  /** \brief Speed of light */
  template<class fp_t> fp_t speed_of_light_f(size_t system=o2scl_mks) {
    if (system==o2scl_cgs) {
      fp_t result=29979245800;
      return result;
    }
    fp_t result=299792458;
    return result;
  }

  /** \brief Reduced Planck's constant times speed of light
      \f$ \hbar c \f$
   */
  template<class fp_t> fp_t hbarc_f(size_t system=o2scl_mks) {
    return hbar_f<fp_t>(system)*speed_of_light_f<fp_t>(system);
  }
  
  /// Elementary charge
  template<class fp_t> fp_t elem_charge_f() {
    fp_t numer=1602176634;
    fp_t denom=1000000000;
    fp_t result=(numer/denom)*1.0e-19;
    return result;
  }

  /** \brief Reduced Planck constant times speed of light 
      in \f$ \mathrm{MeV}~\mathrm{fm} \f$
  */
  template<class fp_t> fp_t hc_mev_fm_f() {
    fp_t hbarc=hbar_f<fp_t>()*speed_of_light_f<fp_t>()/
      elem_charge_f<fp_t>()*1.0e9;
    return hbarc;
  }
  
  /// \f$ \hbar c \f$ in MeV fm (exact)
  const double hc_mev_fm=hc_mev_fm_f<double>();

  /// \f$ \hbar c \f$ in MeV cm (exact)
  const double hc_mev_cm=hc_mev_fm*1.0e-13;
}

/** \brief Constants in CGS units 

    CODATA 2014 values were from \ref Mohr16. CODATA 2018 values are
    from physics.nist.gov/constants. IAU 2015 values are the nominal
    values from arXiv:1510.07674 and arXiv:1605.09788 .
*/
namespace o2scl_cgs {

  /// \name Fundamental constants
  //@{
  /// Speed of light in \f$ \mathrm{cm}/\mathrm{s} \f$ (exact)
  const double speed_of_light=2.99792458e10;
  /// Newtonian constant of gravitation in cm^3 / g s^2 (CODATA 2018 value)
  const double gravitational_constant=6.67430e-8;
  /// Planck constant in g cm^2 / s (CODATA 2018 value; exact)
  const double plancks_constant_h=6.62607015e-27;
  /// Planck constant divided by 2 pi in g cm^2 / s (derived)
  const double plancks_constant_hbar=o2scl_cgs::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Electron volt in g cm^2 / s^2 (CODATA 2018; exact)
  const double electron_volt=1.602176634e-12;
  /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2018; exact)
  const double boltzmann=1.380649e-16;
  /// Bohr radius in cm (CODATA 2018 value)
  const double bohr_radius=5.29177210903e-9;
  /// Stefan-Boltzmann constant in g / K^4 s^3 (CODATA 2018; derived; exact)
  const double stefan_boltzmann_constant=o2scl_const::pi*o2scl_const::pi*
    o2scl_cgs::boltzmann*o2scl_cgs::boltzmann*o2scl_cgs::boltzmann*
    o2scl_cgs::boltzmann/60.0/o2scl_cgs::plancks_constant_hbar/
    o2scl_cgs::plancks_constant_hbar/o2scl_cgs::plancks_constant_hbar/
    o2scl_cgs::speed_of_light/o2scl_cgs::speed_of_light;
  /// Thomson cross section in cm^2 (CODATA 2018 value)
  const double thomson_cross_section=6.6524587321e-25;
  /** \brief Fermi coupling constant in s^4 / cm^4 g^2
      (derived from CODATA 2018 value)
  */
  const double gfermi=o2scl_const::gfermi_gev2*1.0e-18/
    o2scl_cgs::electron_volt/o2scl_cgs::electron_volt;
  //@}

  /// \name Solar system properties
  //@{
  /** \brief Solar mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)

      Note that this value differs slightly in Barycentric Coordinate
      Time and Barycentric Dynamical Time. This is the IAU's nominal
      value.
  */
  const double solar_mass_parameter=1.3271244e26;
  /// Mass of the sun in g (derived)
  const double solar_mass=solar_mass_parameter/gravitational_constant;
  /// Radius of the sun in cm (IAU 2015 value)
  const double solar_radius=6.957e10;
  /// Temperature of the sun's photosphere in K (IAU 2015 value)
  const double solar_temperature=5772.0;
  /// Luminosity of sun in erg/s (IAU 2015 value)
  const double solar_luminosity=3.828e40;
  
  /** \brief Earth mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)
  */
  const double earth_mass_parameter=3.986004e20;
  /// Mass of the earth in g (derived)
  const double earth_mass=earth_mass_parameter/gravitational_constant;
  /// Equatorial radius of earth in cm (IAU 2015 value)
  const double earth_radius_equatorial=6.3781e8;
  /// Polar radius of earth in cm (IAU 2015 value)
  const double earth_radius_polar=6.3568e8;
  
  /** \brief Jupter mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)
  */
  const double jupiter_mass_parameter=1.2668653e23;
  /// Mass of jupiter in g (derived)
  const double jupiter_mass=jupiter_mass_parameter/gravitational_constant;
  /// Equatorial radius of jupiter in cm (IAU 2015 value)
  const double jupiter_radius_equatorial=7.1492e9;
  /// Polar radius of jupiter in cm (IAU 2015 value)
  const double jupiter_radius_polar=6.6854e9;
  
  /// Mass of mercury in g
  const double mercury_mass=3.3011e26;
  /// Radius of mercury in cm
  const double mercury_radius=2.4397e8;

  /// Mass of venus in g
  const double venus_mass=7.8675e27;
  /// Radius of venus in cm
  const double venus_radius=6.0518e8;

  /// Mass of mars in g
  const double mars_mass=6.4171e26;
  /// Equatorial radius of mars in cm
  const double mars_radius_equatorial=3.3962e8;
  /// Polar radius of mars in cm
  const double mars_radius_polar=3.3762e8;

  /// Mass of saturn in g
  const double saturn_mass=5.6834e29;
  /// Equatorial radius of saturn in cm 
  const double saturn_radius_equatorial=6.0268e9;
  /// Polar radius of saturn in cm
  const double saturn_radius_polar=5.4364e9;

  /// Mass of uranus in g
  const double uranus_mass=8.6810e28;
  /// Equatorial radius of uranus in cm 
  const double uranus_radius_equatorial=2.5559e9;
  /// Polar radius of uranus in cm
  const double uranus_radius_polar=2.4973e9;

  /// Mass of neptune in g
  const double neptune_mass=1.02413e29;
  /// Equatorial radius of neptune in cm 
  const double neptune_radius_equatorial=2.4764e9;
  /// Polar radius of neptune in cm
  const double neptune_radius_polar=2.4341e9;

  /// Mass of pluto in g
  const double pluto_mass=1.303e25;
  /// Radius of pluto in cm 
  const double pluto_radius=1.1883e8;
  //@}

  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in cm (IAU 2009 value; now exact)
  const double astronomical_unit=1.495978707e13;
  /// Parsec in \f$ \mathrm{cm} \f$ (derived; exact)
  const double parsec=o2scl_cgs::astronomical_unit*648000.0/o2scl_const::pi;
  /// Acccleration due to gravity in cm / s^2 (CODATA 2018; now exact)
  const double grav_accel=9.80665e2;
  /// Schwarzchild radius in cm (derived)
  const double schwarzchild_radius=2.0*o2scl_cgs::solar_mass_parameter/
    o2scl_cgs::speed_of_light/o2scl_cgs::speed_of_light;
  /** \brief Sidereal year in s 
      (from http://hpiers.obspm.fr/eop-pc/models/constants.html)
  */
  const double sidereal_year=365.256363004*8.64e4;
  /** \brief Tropical year in s 
      (from http://hpiers.obspm.fr/eop-pc/models/constants.html)
  */
  const double tropical_year=365.242190402*8.64e4;
  /// Julian year in s (exact)
  const double julian_year=365.25*8.64e4;
  /// Light year in \f$ \mathrm{cm} \f$ (derived; exact)
  const double light_year=o2scl_cgs::julian_year*o2scl_cgs::speed_of_light;
  //@}

  /// \name Particle masses
  //@{
  /// Electron mass in g (CODATA 2018 value)
  const double mass_electron=9.1093837015e-28;
  /// Muon mass in g (CODATA 2018 value)
  const double mass_muon=1.883531627e-25;
  /// Muon mass in g (CODATA 2018 value)
  const double mass_tau=3.16754e-24;
  /// Proton mass in g (CODATA 2018 value)
  const double mass_proton=1.67262192369e-24;
  /// Neutron mass in g (CODATA 2018 value)
  const double mass_neutron=1.67492749804e-24;
  //@}

  /// \name Nuclear masses
  //@{
  /// Deuteron mass in kg (CODATA 2018 value)
  const double mass_deuteron=3.3435837724e-24;
  /// Triton mass in kg (CODATA 2018 value)
  const double mass_triton=5.0073567446e-24;
  /// Helion mass in kg (CODATA 2018 value)
  const double mass_helion=5.0064127796e-24;
  /// Alpha particle mass in kg (CODATA 2018 value)
  const double mass_alpha=6.6446573357e-24;
  /// Atomic mass constant in g (CODATA 2018 value)
  const double unified_atomic_mass=1.6605390666e-24;
  //@}

  /// \name Chemical constants
  //@{
  /// Rydberg constant in g cm^2 / s^2 (CODATA 2018 value)
  const double rydberg=2.1798723611035e-11;
  /** \brief Molar gas constant, "R", in g cm^2 / K mol s^2 
      (CODATA 2018; exact; derived)
  */
  const double molar_gas=o2scl_const::avogadro*o2scl_cgs::boltzmann;
  /** \brief Molar volume of ideal gas at standard T and P in 
      cm^3 / mol (CODATA 2018 value)
  */
  const double standard_gas_volume=2.271095464e4;
  //@}

  /// \name Unit conversions
  //@{
  /// s
  const double minute=6e1;
  /// s
  const double hour=3.6e3;
  /// s
  const double day=8.64e4;
  /// s
  const double week=6.048e5;
  /// cm
  const double inch=2.54e0;
  /// cm
  const double foot=3.048e1;
  /// cm
  const double yard=9.144e1;
  /// cm
  const double mile=1.609344e5;
  /// cm
  const double nautical_mile=1.852e5;
  /// cm
  const double fathom=1.8288e2;
  /// cm
  const double mil=2.54e-3;
  /// cm
  const double point=3175.0/90000.0;
  /// cm
  const double texpoint=3.51459803515e-2;
  /// cm
  const double micron=1e-4;
  /// cm
  const double angstrom=1e-8;
  /// cm^2
  const double hectare=1e8;
  /// cm^2
  const double acre=4.04685642241e7;
  /// cm^2
  const double barn=1e-24;
  /// cm^3
  const double liter=1e3;
  /// cm^3
  const double us_gallon=3.78541178402e3;
  /// cm^3
  const double quart=9.46352946004e2;
  /// cm^3
  const double pint=4.73176473002e2;
  /// cm^3
  const double cup=2.36588236501e2;
  /// cm^3
  const double fluid_ounce=2.95735295626e1;
  /// cm^3
  const double tablespoon=1.47867647813e1;
  /// cm^3
  const double teaspoon=4.92892159375e0;
  /// cm^3
  const double canadian_gallon=4.54609e3;
  /// cm^3
  const double uk_gallon=4.546092e3;
  /// cm / s
  const double miles_per_hour=4.4704e1;
  /// cm / s
  const double kilometers_per_hour=250.0/9.0;
  /// cm / s
  const double knot=463.0/9.0;
  /// g
  const double pound_mass=4.5359237e2;
  /// g
  const double ounce_mass=2.8349523125e1;
  /// g
  const double ton=9.0718474e5;
  /// g
  const double metric_ton=1e6;
  /// g
  const double uk_ton=1.0160469088e6;
  /// g
  const double troy_ounce=3.1103475e1;
  /// g
  const double carat=2e-1;
  /// cm g / s^2
  const double gram_force=9.80665e2;
  /// cm g / s^2
  const double pound_force=4.44822161526e5;
  /// cm g / s^2
  const double kilopound_force=4.44822161526e8;
  /// cm g / s^2
  const double poundal=1.38255e4;
  /// g cm^2 / s^2
  const double calorie=4.1868e7;
  /// g cm^2 / s^2
  const double btu=1.05505585262e10;
  /// g cm^2 / s^2
  const double therm=1.05506e15;
  /// g cm^2 / s^3
  const double horsepower=7.457e9;
  /// g / cm s^2
  const double bar=1e6;
  /// g / cm s^2
  const double std_atmosphere=1.01325e6;
  /// g / cm s^2
  const double torr=1.33322368421e3;
  /// g / cm s^2
  const double meter_of_mercury=1.33322368421e6;
  /// g / cm s^2
  const double inch_of_mercury=3.38638815789e4;
  /// g / cm s^2
  const double inch_of_water=2.490889e3;
  /// g / cm s^2
  const double psi=6.89475729317e4;
  /// g / cm s
  const double poise=1e0;
  /// cm^2 / s
  const double stokes=1e0;
  /// cd / cm^2
  const double stilb=1e0;
  /// cd sr
  const double lumen=1e0;
  /// cd sr / cm^2
  const double lux=1e-4;
  /// cd sr / cm^2
  const double phot=1e0;
  /// cd sr / cm^2
  const double footcandle=1.076e-3;
  /// cd sr / cm^2
  const double lambert=1e0;
  /// cd sr / cm^2
  const double footlambert=1.07639104e-3;
  /// 1 / s
  const double curie=3.7e10;
  /// cm^2 / s^2
  const double rad=1e2;
  /// cm g / s^2
  const double newton=1e5;
  /// cm g / s^2
  const double dyne=1e0;
  /// g cm^2 / s^2
  const double joule=1e7;
  /// g cm^2 / s^2
  const double erg=1e0;
  //@}

  /// \name Electromagnetic constants
  //@{
  /// A s / g
  const double roentgen=2.58e-7;
  //@}
  
}

/** \brief Constants in CGSM units
    
    Where possible, constants here are defined in terms of the values
    in \ref o2scl_cgs, in order to make it easier to update these
    values. See also the documentation at \ref o2scl_cgs .
*/
namespace o2scl_cgsm {
  
  /// \name Fundamental constants
  //@{
  /// Speed of light in cm / s
  const double speed_of_light=o2scl_cgs::speed_of_light;
  /// Newtonian constant of gravitation in cm^3 / g s^2
  const double gravitational_constant=o2scl_cgs::gravitational_constant;
  /// Planck constant in g cm^2 / s
  const double plancks_constant_h=o2scl_cgs::plancks_constant_h;
  /// Planck constant divided by 2 pi in g cm^2 / s
  const double plancks_constant_hbar=o2scl_cgs::plancks_constant_hbar;
  /// Electron volt in g cm^2 / s^2
  const double electron_volt=o2scl_cgs::electron_volt;
  /// Boltzmann constant in g cm^2 / K s^2
  const double boltzmann=o2scl_cgs::boltzmann;
  /// Bohr radius in cm
  const double bohr_radius=o2scl_cgs::bohr_radius;
  /// Stefan-Boltzmann constant in g / K^4 s^3
  const double stefan_boltzmann_constant=
    o2scl_cgs::stefan_boltzmann_constant;
  /// Thomson cross section in cm^2
  const double thomson_cross_section=o2scl_cgs::thomson_cross_section;
  /// Fermi coupling constant in s^4 / cm^4 g^2
  const double gfermi=o2scl_cgs::gfermi;
  //@}

  /// \name Solar system properties
  //@{
  /** \brief Solar mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)

      Note that this value differs slightly in Barycentric Coordinate
      Time and Barycentric Dynamical Time. This is the IAU's nominal
      value.
  */
  const double solar_mass_parameter=o2scl_cgs::solar_mass_parameter;
  /// Mass of the sun in g (derived)
  const double solar_mass=o2scl_cgs::solar_mass;
  /// Radius of the sun in cm (IAU 2015 value)
  const double solar_radius=o2scl_cgs::solar_radius;
  /// Temperature of the sun's photosphere in K (IAU 2015 value)
  const double solar_temperature=o2scl_cgs::solar_temperature;
  /// Luminosity of sun in erg/s (IAU 2015 value)
  const double solar_luminosity=o2scl_cgs::solar_luminosity;
  
  /** \brief Earth mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)
  */
  const double earth_mass_parameter=o2scl_cgs::earth_mass_parameter;
  /// Mass of the earth in g (derived)
  const double earth_mass=o2scl_cgs::earth_mass;
  /// Equatorial radius of earth in cm (IAU 2015 value)
  const double earth_radius_equatorial=o2scl_cgs::earth_radius_equatorial;
  /// Polar radius of earth in cm (IAU 2015 value)
  const double earth_radius_polar=o2scl_cgs::earth_radius_polar;
  
  /** \brief Jupter mass times gravitational constant in cm^3 / s^2
      (IAU 2015 value)
  */
  const double jupiter_mass_parameter=o2scl_cgs::jupiter_mass_parameter;
  /// Mass of jupiter in g (derived)
  const double jupiter_mass=o2scl_cgs::jupiter_mass;
  /// Equatorial radius of jupiter in cm (IAU 2015 value)
  const double jupiter_radius_equatorial=o2scl_cgs::jupiter_radius_equatorial;
  /// Polar radius of jupiter in cm (IAU 2015 value)
  const double jupiter_radius_polar=o2scl_cgs::jupiter_radius_polar;
  
  /// Mass of mercury in g
  const double mercury_mass=o2scl_cgs::mercury_mass;
  /// Radius of mercury in cm
  const double mercury_radius=o2scl_cgs::mercury_radius;

  /// Mass of venus in g
  const double venus_mass=o2scl_cgs::venus_mass;
  /// Radius of venus in cm
  const double venus_radius=o2scl_cgs::venus_radius;

  /// Mass of mars in g
  const double mars_mass=o2scl_cgs::mars_mass;
  /// Equatorial radius of mars in cm
  const double mars_radius_equatorial=o2scl_cgs::mars_radius_equatorial;
  /// Polar radius of mars in cm
  const double mars_radius_polar=o2scl_cgs::mars_radius_polar;

  /// Mass of saturn in g
  const double saturn_mass=o2scl_cgs::saturn_mass;
  /// Equatorial radius of saturn in cm 
  const double saturn_radius_equatorial=o2scl_cgs::saturn_radius_equatorial;
  /// Polar radius of saturn in cm
  const double saturn_radius_polar=o2scl_cgs::saturn_radius_polar;

  /// Mass of uranus in g
  const double uranus_mass=o2scl_cgs::uranus_mass;
  /// Equatorial radius of uranus in cm 
  const double uranus_radius_equatorial=o2scl_cgs::uranus_radius_equatorial;
  /// Polar radius of uranus in cm
  const double uranus_radius_polar=o2scl_cgs::uranus_radius_polar;

  /// Mass of neptune in g
  const double neptune_mass=o2scl_cgs::neptune_mass;
  /// Equatorial radius of neptune in cm 
  const double neptune_radius_equatorial=o2scl_cgs::neptune_radius_equatorial;
  /// Polar radius of neptune in cm
  const double neptune_radius_polar=o2scl_cgs::neptune_radius_polar;

  /// Mass of pluto in g
  const double pluto_mass=o2scl_cgs::pluto_mass;
  /// Radius of pluto in cm 
  const double pluto_radius=o2scl_cgs::pluto_radius;
  //@}
  
  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in cm (IAU 2009 value; now exact)
  const double astronomical_unit=o2scl_cgs::astronomical_unit;
  /// Parsec in \f$ \mathrm{cm} \f$ (derived; exact)
  const double parsec=o2scl_cgs::parsec;
  /// Acccleration due to gravity in cm / s^2 (CODATA 2018; now exact)
  const double grav_accel=o2scl_cgs::grav_accel;
  /// Schwarzchild radius in cm (derived)
  const double schwarzchild_radius=o2scl_cgs::schwarzchild_radius;
  /** \brief Sidereal year in s 
      (from http://hpiers.obspm.fr/eop-pc/models/constants.html)
  */
  const double sidereal_year=o2scl_cgs::sidereal_year;
  /** \brief Tropical year in s 
      (from http://hpiers.obspm.fr/eop-pc/models/constants.html)
  */
  const double tropical_year=o2scl_cgs::tropical_year;
  /// Julian year in s (exact)
  const double julian_year=o2scl_cgs::julian_year;
  /// Light year in \f$ \mathrm{cm} \f$ (derived; exact)
  const double light_year=o2scl_cgs::light_year;
  //@}
  
  /// \name Particle masses
  //@{
  /// Electron mass in g
  const double mass_electron=o2scl_cgs::mass_electron;
  /// Muon mass in g
  const double mass_muon=o2scl_cgs::mass_muon;
  /// Muon mass in g
  const double mass_tau=o2scl_cgs::mass_tau;
  /// Proton mass in g
  const double mass_proton=o2scl_cgs::mass_proton;
  /// Neutron mass in g
  const double mass_neutron=o2scl_cgs::mass_neutron;
  //@}

  /// \name Nuclear masses
  //@{
  /// Deuteron mass in g
  const double mass_deuteron=o2scl_cgs::mass_deuteron;
  /// Triton mass in g
  const double mass_triton=o2scl_cgs::mass_triton;
  /// Helion mass in g
  const double mass_helion=o2scl_cgs::mass_helion;
  /// Alpha particle mass in g
  const double mass_alpha=o2scl_cgs::mass_alpha;
  /// Atomic mass constant in g
  const double unified_atomic_mass=o2scl_cgs::unified_atomic_mass;
  //@}

  /// \name Chemical constants
  //@{
  /// Rydberg constant in g cm^2 / s^2
  const double rydberg=o2scl_cgs::rydberg;
  /// Molar gas constant, "R", in g cm^2 / K mol s^2
  const double molar_gas=o2scl_cgs::molar_gas;
  /// Molar volume of ideal gas at standard T and P in cm^3 / mol
  const double standard_gas_volume=o2scl_cgs::standard_gas_volume;
  //@}

  /// \name Unit conversions
  //@{
  /// s
  const double minute=o2scl_cgs::minute;
  /// s
  const double hour=o2scl_cgs::hour;
  /// s
  const double day=o2scl_cgs::day;
  /// s
  const double week=o2scl_cgs::week;
  /// cm
  const double inch=o2scl_cgs::inch;
  /// cm
  const double foot=o2scl_cgs::foot;
  /// cm
  const double yard=o2scl_cgs::yard;
  /// cm
  const double mile=o2scl_cgs::mile;
  /// cm
  const double nautical_mile=o2scl_cgs::nautical_mile;
  /// cm
  const double fathom=o2scl_cgs::fathom;
  /// cm
  const double mil=o2scl_cgs::mil;
  /// cm
  const double point=o2scl_cgs::point;
  /// cm
  const double texpoint=o2scl_cgs::texpoint;
  /// cm
  const double micron=o2scl_cgs::micron;
  /// cm
  const double angstrom=o2scl_cgs::angstrom;
  /// cm^2
  const double hectare=o2scl_cgs::hectare;
  /// cm^2
  const double acre=o2scl_cgs::acre;
  /// cm^2
  const double barn=o2scl_cgs::barn;
  /// cm^3
  const double liter=o2scl_cgs::liter;
  /// cm^3
  const double us_gallon=o2scl_cgs::us_gallon;
  /// cm^3
  const double quart=o2scl_cgs::quart;
  /// cm^3
  const double pint=o2scl_cgs::pint;
  /// cm^3
  const double cup=o2scl_cgs::cup;
  /// cm^3
  const double fluid_ounce=o2scl_cgs::fluid_ounce;
  /// cm^3
  const double tablespoon=o2scl_cgs::tablespoon;
  /// cm^3
  const double teaspoon=o2scl_cgs::teaspoon;
  /// cm^3
  const double canadian_gallon=o2scl_cgs::canadian_gallon;
  /// cm^3
  const double uk_gallon=o2scl_cgs::uk_gallon;
  /// cm / s
  const double miles_per_hour=o2scl_cgs::miles_per_hour;
  /// cm / s
  const double kilometers_per_hour=o2scl_cgs::kilometers_per_hour;
  /// cm / s
  const double knot=o2scl_cgs::knot;
  /// g
  const double pound_mass=o2scl_cgs::pound_mass;
  /// g
  const double ounce_mass=o2scl_cgs::ounce_mass;
  /// g
  const double ton=o2scl_cgs::ton;
  /// g
  const double metric_ton=o2scl_cgs::metric_ton;
  /// g
  const double uk_ton=o2scl_cgs::uk_ton;
  /// g
  const double troy_ounce=o2scl_cgs::troy_ounce;
  /// g
  const double carat=o2scl_cgs::carat;
  /// cm g / s^2
  const double gram_force=o2scl_cgs::gram_force;
  /// cm g / s^2
  const double pound_force=o2scl_cgs::pound_force;
  /// cm g / s^2
  const double kilopound_force=o2scl_cgs::kilopound_force;
  /// cm g / s^2
  const double poundal=o2scl_cgs::poundal;
  /// g cm^2 / s^2
  const double calorie=o2scl_cgs::calorie;
  /// g cm^2 / s^2
  const double btu=o2scl_cgs::btu;
  /// g cm^2 / s^2
  const double therm=o2scl_cgs::therm;
  /// g cm^2 / s^3
  const double horsepower=o2scl_cgs::horsepower;
  /// g / cm s^2
  const double bar=o2scl_cgs::bar;
  /// g / cm s^2
  const double std_atmosphere=o2scl_cgs::std_atmosphere;
  /// g / cm s^2
  const double torr=o2scl_cgs::torr;
  /// g / cm s^2
  const double meter_of_mercury=o2scl_cgs::meter_of_mercury;
  /// g / cm s^2
  const double inch_of_mercury=o2scl_cgs::inch_of_mercury;
  /// g / cm s^2
  const double inch_of_water=o2scl_cgs::inch_of_water;
  /// g / cm s^2
  const double psi=o2scl_cgs::psi;
  /// g / cm s
  const double poise=o2scl_cgs::poise;
  /// cm^2 / s
  const double stokes=o2scl_cgs::stokes;
  /// cd / cm^2
  const double stilb=o2scl_cgs::stilb;
  /// cd sr
  const double lumen=o2scl_cgs::lumen;
  /// cd sr / cm^2
  const double lux=o2scl_cgs::lux;
  /// cd sr / cm^2
  const double phot=o2scl_cgs::phot;
  /// cd sr / cm^2
  const double footcandle=o2scl_cgs::footcandle;
  /// cd sr / cm^2
  const double lambert=o2scl_cgs::lambert;
  /// cd sr / cm^2
  const double footlambert=o2scl_cgs::footlambert;
  /// 1 / s
  const double curie=o2scl_cgs::curie;
  /// cm^2 / s^2
  const double rad=o2scl_cgs::rad;
  /// cm g / s^2
  const double newton=o2scl_cgs::newton;
  /// cm g / s^2
  const double dyne=o2scl_cgs::dyne;
  /// g cm^2 / s^2
  const double joule=o2scl_cgs::joule;
  /// g cm^2 / s^2
  const double erg=o2scl_cgs::erg;
  //@}

  /// \name Electromagnetic constants
  //@{
  /// Electron magnetic moment in abamp cm^2 (CODATA 2018 value)
  const double electron_magnetic_moment=9.2847647043e-21;
  /// Proton magnetic moment in abamp cm^2 (CODATA 2018 value)
  const double proton_magnetic_moment=1.41060679736e-23;
  /// Roentgen abamp s / g
  const double roentgen=o2scl_cgs::roentgen/10.0;
  /// Bohr magneton in abamp cm^2 (CODATA 2018 value)
  const double bohr_magneton=9.2740100783e-21;
  /// Nuclear magneton in abamp cm^2 (CODATA 2018 value)
  const double nuclear_magneton=5.0507837461e-24;
  /// Faraday constant in abamp s / mol (CODATA 2018 value; derived; exact)
  const double faraday=o2scl_const::avogadro*o2scl_cgs::electron_volt/1.0e8;
  /// Electron charge in abamp s (derived)
  const double electron_charge=electron_volt*1.0e-8;
  //@}
}

/** \brief Constants in MKS units

    Where possible, constants here are defined in terms of the values
    in \ref o2scl_cgs, in order to make it easier to update these
    values. See also the documentation at \ref o2scl_cgs .
*/
namespace o2scl_mks {

  /// \name Fundamental constants
  //@{
  /// Speed of light in m / s
  const double speed_of_light=o2scl_cgs::speed_of_light/1.0e2;
  /// Newtonian constant of gravitation in m^3 / kg s^2
  const double gravitational_constant=o2scl_cgs::gravitational_constant/1.0e3;
  /// Planck constant in kg m^2 / s
  const double plancks_constant_h=o2scl_cgs::plancks_constant_h/1.0e7;
  /// Planck constant divided by 2 pi in kg m^2 / s
  const double plancks_constant_hbar=o2scl_cgs::plancks_constant_hbar/1.0e7;
  /// Electron volt in kg m^2 / s^2
  const double electron_volt=o2scl_cgs::electron_volt/1.0e7;
  /// Boltzmann constant in kg m^2 / K s^2
  const double boltzmann=o2scl_cgs::boltzmann/1.0e7;
  /// Bohr radius in m
  const double bohr_radius=o2scl_cgs::bohr_radius/1.0e2;
  /// Stefan-Boltzmann constant in kg / K^4 s^3
  const double stefan_boltzmann_constant=
    o2scl_cgs::stefan_boltzmann_constant/1.0e3;
  /// Thomson cross section in m^2
  const double thomson_cross_section=o2scl_cgs::thomson_cross_section/1.0e4;
  /// Fermi coupling constant in s^4 / m^4 kg^2
  const double gfermi=o2scl_cgs::gfermi*1.0e14;
  //@{

  /// \name Solar system properties
  //@{
  /** \brief Solar mass times gravitational constant in km^3 / s^2
  */
  const double solar_mass_parameter=o2scl_cgs::solar_mass_parameter;
  /// Mass of the sun in kg
  const double solar_mass=o2scl_cgs::solar_mass/1.0e3;
  /// Radius of the sun in m
  const double solar_radius=o2scl_cgs::solar_radius/1.0e2;
  /// Temperature of the sun's photosphere in K 
  const double solar_temperature=o2scl_cgs::solar_temperature;
  /// Luminosity of sun in erg/s
  const double solar_luminosity=o2scl_cgs::solar_luminosity;
  
  /** \brief Earth mass times gravitational constant in m^3 / s^2
  */
  const double earth_mass_parameter=o2scl_cgs::earth_mass_parameter;
  /// Mass of the earth in kg
  const double earth_mass=o2scl_cgs::earth_mass/1.0e3;
  /// Equatorial radius of earth in m
  const double earth_radius_equatorial=o2scl_cgs::earth_radius_equatorial/1.0e2;
  /// Polar radius of earth in m
  const double earth_radius_polar=o2scl_cgs::earth_radius_polar/1.0e2;
  
  /** \brief Jupter mass times gravitational constant in m^3 / s^2
  */
  const double jupiter_mass_parameter=o2scl_cgs::jupiter_mass_parameter;
  /// Mass of jupiter in kg
  const double jupiter_mass=o2scl_cgs::jupiter_mass/1.0e3/1.0e2;
  /// Equatorial radius of jupiter in m
  const double jupiter_radius_equatorial=o2scl_cgs::jupiter_radius_equatorial;
  /// Polar radius of jupiter in m
  const double jupiter_radius_polar=o2scl_cgs::jupiter_radius_polar/1.0e2;
  
  /// Mass of mercury in kg
  const double mercury_mass=o2scl_cgs::mercury_mass/1.0e3;
  /// Radius of mercury in m
  const double mercury_radius=o2scl_cgs::mercury_radius/1.0e2;

  /// Mass of venus in kg
  const double venus_mass=o2scl_cgs::venus_mass/1.0e3;
  /// Radius of venus in m
  const double venus_radius=o2scl_cgs::venus_radius/1.0e2;

  /// Mass of mars in kg
  const double mars_mass=o2scl_cgs::mars_mass/1.0e3;
  /// Equatorial radius of mars in m
  const double mars_radius_equatorial=o2scl_cgs::mars_radius_equatorial/1.0e2;
  /// Polar radius of mars in m
  const double mars_radius_polar=o2scl_cgs::mars_radius_polar/1.0e2;

  /// Mass of saturn in kg
  const double saturn_mass=o2scl_cgs::saturn_mass/1.0e3;
  /// Equatorial radius of saturn in m 
  const double saturn_radius_equatorial=
    o2scl_cgs::saturn_radius_equatorial/1.0e2;
  /// Polar radius of saturn in m
  const double saturn_radius_polar=o2scl_cgs::saturn_radius_polar/1.0e2;

  /// Mass of uranus in kg
  const double uranus_mass=o2scl_cgs::uranus_mass/1.0e3;
  /// Equatorial radius of uranus in m 
  const double uranus_radius_equatorial=
    o2scl_cgs::uranus_radius_equatorial/1.0e2;
  /// Polar radius of uranus in m
  const double uranus_radius_polar=o2scl_cgs::uranus_radius_polar/1.0e2;

  /// Mass of neptune in kg
  const double neptune_mass=o2scl_cgs::neptune_mass/1.0e3;
  /// Equatorial radius of neptune in m 
  const double neptune_radius_equatorial=
    o2scl_cgs::neptune_radius_equatorial/1.0e2;
  /// Polar radius of neptune in m
  const double neptune_radius_polar=o2scl_cgs::neptune_radius_polar/1.0e2;

  /// Mass of pluto in kg
  const double pluto_mass=o2scl_cgs::pluto_mass/1.0e3;
  /// Radius of pluto in m 
  const double pluto_radius=o2scl_cgs::pluto_radius/1.0e2;
  //@}
  
  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in m
  const double astronomical_unit=o2scl_cgs::astronomical_unit/1.0e2;
  /// Parsec in \f$ \mathrm{m} \f$
  const double parsec=o2scl_cgs::parsec/1.0e2;
  /// Acccleration due to gravity in m / s^2
  const double grav_accel=o2scl_cgs::grav_accel/1.0e2;
  /// Schwarzchild radius in m
  const double schwarzchild_radius=o2scl_cgs::schwarzchild_radius/1.0e2;
  /** \brief Sidereal year in s 
  */
  const double sidereal_year=o2scl_cgs::sidereal_year;
  /** \brief Tropical year in s 
   */
  const double tropical_year=o2scl_cgs::tropical_year;
  /// Julian year in s
  const double julian_year=o2scl_cgs::julian_year;
  /// Light year in \f$ \mathrm{m} \f$
  const double light_year=o2scl_cgs::light_year/1.0e2;
  //@}
  
  /// \name Particle masses
  //@{
  /// Electron mass in kg
  const double mass_electron=o2scl_cgs::mass_electron/1.0e3;
  /// Muon mass in kg
  const double mass_muon=o2scl_cgs::mass_muon/1.0e3;
  /// Muon mass in kg
  const double mass_tau=o2scl_cgs::mass_tau/1.0e3;
  /// Proton mass in kg
  const double mass_proton=o2scl_cgs::mass_proton/1.0e3;
  /// Neutron mass in kg
  const double mass_neutron=o2scl_cgs::mass_neutron/1.0e3;
  //@}

  /// \name Nuclear masses
  //@{
  /// Deuteron mass in kg
  const double mass_deuteron=o2scl_cgs::mass_deuteron/1.0e3;
  /// Triton mass in kg
  const double mass_triton=o2scl_cgs::mass_triton/1.0e3;
  /// Helion mass in kg
  const double mass_helion=o2scl_cgs::mass_helion/1.0e3;
  /// Alpha particle mass in kg
  const double mass_alpha=o2scl_cgs::mass_alpha/1.0e3;
  /// Atomic mass constant in kg 
  const double unified_atomic_mass=o2scl_cgs::unified_atomic_mass/1.0e3;
  //@}

  /// \name Chemical constants
  //@{
  /// Rydberg constant in kg m^2 / s^2
  const double rydberg=o2scl_cgs::rydberg/1.0e7;
  /// kg m^2 / K mol s^2
  const double molar_gas=o2scl_cgs::molar_gas/1.0e7;
  /// m^3 / mol
  const double standard_gas_volume=o2scl_cgs::standard_gas_volume/1.0e6;
  //@}

  /// \name Unit conversions
  //@{
  /// s
  const double minute=o2scl_cgs::minute;
  /// s
  const double hour=o2scl_cgs::hour;
  /// s
  const double day=o2scl_cgs::day;
  /// s
  const double week=o2scl_cgs::week;
  /// m
  const double inch=o2scl_cgs::inch*1.0e-2;
  /// m
  const double foot=o2scl_cgs::foot*1.0e-2;
  /// m
  const double yard=o2scl_cgs::yard*1.0e-2;
  /// m
  const double mile=o2scl_cgs::mile*1.0e-2;
  /// m
  const double nautical_mile=o2scl_cgs::nautical_mile*1.0e-2;
  /// m
  const double fathom=o2scl_cgs::fathom*1.0e-2;
  /// m
  const double mil=o2scl_cgs::mil*1.0e-2;
  /// m
  const double point=o2scl_cgs::point*1.0e-2;
  /// m
  const double texpoint=o2scl_cgs::texpoint*1.0e-2;
  /// m
  const double micron=o2scl_cgs::micron*1.0e-2;
  /// m
  const double angstrom=o2scl_cgs::angstrom*1.0e-2;
  /// m^2
  const double hectare=o2scl_cgs::hectare*1.0e-4;
  /// m^2
  const double acre=o2scl_cgs::acre*1.0e-4;
  /// m^2
  const double barn=o2scl_cgs::barn*1.0e-4;
  /// m^3
  const double liter=o2scl_cgs::liter*1.0e-6;
  /// m^3
  const double us_gallon=o2scl_cgs::us_gallon*1.0e-6;
  /// m^3
  const double quart=o2scl_cgs::quart*1.0e-6;
  /// m^3
  const double pint=o2scl_cgs::pint*1.0e-6;
  /// m^3
  const double cup=o2scl_cgs::cup*1.0e-6;
  /// m^3
  const double fluid_ounce=o2scl_cgs::fluid_ounce*1.0e-6;
  /// m^3
  const double tablespoon=o2scl_cgs::tablespoon*1.0e-6;
  /// m^3
  const double teaspoon=o2scl_cgs::teaspoon*1.0e-6;
  /// m^3
  const double canadian_gallon=o2scl_cgs::canadian_gallon*1.0e-6;
  /// m^3
  const double uk_gallon=o2scl_cgs::uk_gallon*1.0e-6;
  /// m / s
  const double miles_per_hour=o2scl_cgs::miles_per_hour*1.0e-2;
  /// m / s
  const double kilometers_per_hour=o2scl_cgs::kilometers_per_hour*1.0e-2;
  /// m / s
  const double knot=o2scl_cgs::knot*1.0e-2;
  /// kg
  const double pound_mass=o2scl_cgs::pound_mass*1.0e-3;
  /// kg
  const double ounce_mass=o2scl_cgs::ounce_mass*1.0e-3;
  /// kg
  const double ton=o2scl_cgs::ton*1.0e-3;
  /// kg
  const double metric_ton=o2scl_cgs::metric_ton*1.0e-3;
  /// kg
  const double uk_ton=o2scl_cgs::uk_ton*1.0e-3;
  /// kg
  const double troy_ounce=o2scl_cgs::troy_ounce*1.0e-3;
  /// kg
  const double carat=o2scl_cgs::carat*1.0e-3;
  /// kg m / s^2
  const double gram_force=o2scl_cgs::gram_force*1.0e-5;
  /// kg m / s^2
  const double pound_force=o2scl_cgs::pound_force*1.0e-5;
  /// kg m / s^2
  const double kilopound_force=o2scl_cgs::kilopound_force*1.0e-5;
  /// kg m / s^2
  const double poundal=o2scl_cgs::poundal*1.0e-5;
  /// kg m^2 / s^2
  const double calorie=o2scl_cgs::calorie*1.0e-7;
  /// kg m^2 / s^2
  const double btu=o2scl_cgs::btu*1.0e-7;
  /// kg m^2 / s^2
  const double therm=o2scl_cgs::therm*1.0e-7;
  /// kg m^2 / s^3
  const double horsepower=o2scl_cgs::horsepower*1.0e-7;
  /// kg / m s^2
  const double bar=o2scl_cgs::bar*1.0e-1;
  /// kg / m s^2
  const double std_atmosphere=o2scl_cgs::std_atmosphere*1.0e-1;
  /// kg / m s^2
  const double torr=o2scl_cgs::torr*1.0e-1;
  /// kg / m s^2
  const double meter_of_mercury=o2scl_cgs::meter_of_mercury*1.0e-1;
  /// kg / m s^2
  const double inch_of_mercury=o2scl_cgs::inch_of_mercury*1.0e-1;
  /// kg / m s^2
  const double inch_of_water=o2scl_cgs::inch_of_water*1.0e-1;
  /// kg / m s^2
  const double psi=o2scl_cgs::psi*1.0e-1;
  /// kg m^-1 s^-1
  const double poise=o2scl_cgs::poise*1.0e-1;
  /// m^2 / s
  const double stokes=o2scl_cgs::stokes*1.0e-4;
  /// kg / A s^2
  const double gauss=1.0e-4;
  /// cd / m^2
  const double stilb=o2scl_cgs::stilb*1.0e4;
  /// cd sr
  const double lumen=o2scl_cgs::lumen;
  /// cd sr / m^2
  const double lux=o2scl_cgs::lux*1.0e4;
  /// cd sr / m^2
  const double phot=o2scl_cgs::phot*1.0e4;
  /// cd sr / m^2
  const double footcandle=o2scl_cgs::footcandle*1.0e4;
  /// cd sr / m^2
  const double lambert=o2scl_cgs::lambert*1.0e4;
  /// cd sr / m^2
  const double footlambert=o2scl_cgs::footlambert*1.0e4;
  /// 1 / s
  const double curie=o2scl_cgs::curie;
  /// m^2 / s^2
  const double rad=o2scl_cgs::rad*1.0e-4;
  /// kg m / s^2
  const double newton=o2scl_cgs::newton*1.0e-5;
  /// kg m / s^2
  const double dyne=o2scl_cgs::dyne*1.0e-5;
  /// kg m^2 / s^2
  const double joule=o2scl_cgs::joule*1.0e-7;
  /// kg m^2 / s^2
  const double erg=o2scl_cgs::erg*1.0e-7;
  //@}
  
  /// \name ELectromagnetic constants
  //@{
  /// A m^2
  const double electron_magnetic_moment=
    o2scl_cgsm::electron_magnetic_moment/1.0e3;
  /// A m^2
  const double proton_magnetic_moment=
    o2scl_cgsm::proton_magnetic_moment/1.0e3;
  /// A s / kg
  const double roentgen=o2scl_cgs::roentgen*1.0e3;
  /// Bohr magneton in A m^2
  const double bohr_magneton=o2scl_cgsm::bohr_magneton/1.0e3;
  /// A m^2
  const double nuclear_magneton=o2scl_cgsm::nuclear_magneton/1.0e3;
  /// A^2 s^4 / kg m^3 (derived)
  const double vacuum_permittivity=1.0/o2scl_mks::speed_of_light/
    o2scl_mks::speed_of_light/4.0e-7/o2scl_const::pi;
  /** \brief Vacuum permeability in kg m / A^2 s^2 
      (being redefined as of 5/20, this value is from Wikipedia)
  */
  const double vacuum_permeability=4.0*o2scl_const::pi*1.00000000082e-7;
  /// A s / mol
  const double faraday=o2scl_cgsm::faraday*10.0;
  /// A s (derived)
  const double electron_charge=o2scl_cgsm::electron_charge*1.0e1;
  //@}

}

// Other derived values to add to the namespace
namespace o2scl_const {

  /** \name Squared electron charge
   */
  //@{
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
  const double e2_gaussian=o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure;

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
  const double e2_hlorentz=o2scl_const::fine_structure*4.0*pi;

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
      as mentioned, e.g. in pg. 13 of D. Griffiths Intro to Elem. Particles.
  */      
  const double e2_mks=o2scl_mks::electron_charge;
  //@}

  /** \brief 1 \f$\mathrm{Gauss}\f$ times the electron charge 
      in Gaussian units in \f$\mathrm{fm}^{-2}\f$
  */
  const double ec_gauss_fm2=o2scl_mks::electron_charge*1.0e-34/
    o2scl_mks::plancks_constant_hbar;

  /** \brief Conversion factor from \f$ \mathrm{Gauss}^2 \f$ to
      \f$\mathrm{fm}^{-4}\f$ in Gaussian units.

      This is useful, e.g. in converting magnetic field squared
      to an energy density.
  */
  const double gauss2_fm4=ec_gauss_fm2*ec_gauss_fm2/
    o2scl_const::fine_structure;

  /// \name Particle masses from PDG (4/30/19)
  //@{
  /** \brief \f$ \Lambda \f$ hyperon mass in \f$ \mathrm{MeV} \f$
  */
  const double mass_lambda_MeV=1115.683;

  /** \brief \f$ \Sigma^{-} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_minus_MeV=1197.449;

  /** \brief \f$ \Sigma^{0} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_zero_MeV=1192.642;

  /** \brief \f$ \Sigma^{+} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_plus_MeV=1189.37;

  /** \brief \f$ \Xi^{0} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_cascade_zero_MeV=1314.86;
  
  /** \brief \f$ \Xi^{-} \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_cascade_minus_MeV=1321.71;
  
  /** \brief Up quark mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_up_MeV=2.2;
  
  /** \brief Down quark mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_down_MeV=4.7;
  
  /** \brief Strange quark mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_strange_MeV=95.0;
  //@}
}


#endif
