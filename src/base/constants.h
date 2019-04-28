/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

/** \brief Constants
    
    CODATA 2014 values are in \ref Mohr16.
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

  /** \brief Fine structure constant (CODATA 2014 value)
   */
  const double fine_structure=7.2973525664e-3;
  /** \brief Avogadro's number (CODATA 2014 value)
   */
  const double avogadro=6.022140857e23;

  /// \f$ \pi \f$ 
  const double pi=acos(-1.0);
  /// \f$ \pi^2 \f$ 
  const double pi2=pi*pi;
  /// \f$ \zeta(3/2) \f$
  const double zeta32=2.6123753486854883433;
  /// \f$ \zeta(2) \f$
  const double zeta2=1.6449340668482264365;
  /// \f$ \zeta(5/2) \f$
  const double zeta52=1.3414872572509171798;
  /// \f$ \zeta(3) \f$
  const double zeta3=1.2020569031595942854;
  /// \f$ \zeta(5) \f$
  const double zeta5=1.0369277551433699263;
  /// \f$ \zeta(7) \f$
  const double zeta7=1.0083492773819228268;
  
  /// \f$ \sin^2 \theta_W \f$ (PDG 2018 value)
  const double sin2_theta_weak=0.23122;
}
  
/** \brief Constants in CGS units 
    
    CODATA 2014 values are from \ref Mohr16 . The solar mass and solar
    mass parameter are from 2018 value at http://asa.usno.navy.mil/ .

*/
namespace o2scl_cgs {

  /// \name Fundamental constants
  //@{
  /// Speed of light in \f$ \mathrm{cm}/\mathrm{s} \f$ (exact)
  const double speed_of_light=2.99792458e10;
  /// Newtonian constant of gravitation in cm^3 / g s^2 (CODATA 2014 value)
  const double gravitational_constant=6.67408e-8;
  /// Planck constant in g cm^2 / s (CODATA 2014 value)
  const double plancks_constant_h=6.62607004e-27;
  /// Planck constant divided by 2 pi in g cm^2 / s (derived)
  const double plancks_constant_hbar=o2scl_cgs::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Electron volt in g cm^2 / s^2 (CODATA 2014 value)
  const double electron_volt=1.6021766208e-12;
  /// Bohr radius in cm (CODATA 2014 value)
  const double bohr_radius=5.2917721067e-9;
  /// Stefan-Boltzmann constant in g / K^4 s^3 (CODATA 2014 value)
  const double stefan_boltzmann_constant=5.670367e-5;
  /// Thomson cross section in cm^2 (CODATA 2014 value)
  const double thomson_cross_section=6.652457158e-25;
  /** \brief Fermi coupling constant in s^4 / cm^4 g^2, 
      defined as \f$ 1.1663787 \times 10^{-5}~\mathrm{GeV}^{-2} \f$
      (CODATA 2014 value)
  */
  const double Gfermi=1.1663787e-23/o2scl_cgs::electron_volt/
    o2scl_cgs::electron_volt;
  /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2014 value)
  const double boltzmann=1.38064852e-16;
  //@}

  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in cm (IAU 2009 value; now exact)
  const double astronomical_unit=1.495978707e13;
  /// Parsec in \f$ \mathrm{cm} \f$ (derived)
  const double parsec=o2scl_cgs::astronomical_unit*64800.0/o2scl_const::pi;
  /// Acccleration due to gravity in cm / s^2
  const double grav_accel=9.80665e2;
  /** \brief Solar mass times gravitational constant in cm^3 / s^2
      (from navy.mil)
  */
  const double solar_mass_parameter=1.32712441e26;
  /// Mass of the sun in g (from navy.mil)
  const double solar_mass=1.9884e33;
  /// Schwarzchild radius in cm (derived)
  const double schwarzchild_radius=2.0*o2scl_cgs::solar_mass_parameter/
    o2scl_cgs::speed_of_light/o2scl_cgs::speed_of_light;
  /// Sidereal year in s 
  const double sidereal_year=365.256363004*8.64e4;
  /// Tropical year in s 
  const double tropical_year=365.242190402*8.64e4;
  /// Julian year in s 
  const double julian_year=365.25*8.64e4;
  /// Light year in \f$ \mathrm{cm} \f$ (derived)
  const double light_year=o2scl_cgs::julian_year*o2scl_cgs::speed_of_light;
  //@}

  /// \name Particle masses
  //@{
  /// Electron mass in g (CODATA 2014 value)
  const double mass_electron=9.10938356e-28;
  /// Muon mass in g (CODATA 2014 value)
  const double mass_muon=1.883531594e-25;
  /// Proton mass in g (CODATA 2014 value)
  const double mass_proton=1.672621898e-24;
  /// Neutron mass in g (CODATA 2014 value)
  const double mass_neutron=1.674927471e-24;
  //@}

  /// \name Nuclear masses
  //@{
  /// Deuteron mass in kg (CODATA 2014 value)
  const double mass_deuteron=3.343583719e-24;
  /// Triton mass in kg (CODATA 2014 value)
  const double mass_triton=5.007356665e-24;
  /// Helion mass in kg (CODATA 2014 value)
  const double mass_helion=5.006412700e-24;
  /// Alpha particle mass in kg (CODATA 2014 value)
  const double mass_alpha=6.64465723e-24;
  /// Atomic mass constant in g (CODATA 2014 value)
  const double unified_atomic_mass=1.66053904e-24;
  //@}

  /// \name Chemical constants
  //@{
  /// Rydberg constant in g cm^2 / s^2 (CODATA 2014 value)
  const double rydberg=2.179872325e-11;
  /// Molar gas constant, "R", in g cm^2 / K mol s^2 (CODATA 2014 value)
  const double molar_gas=8.3144598e7;
  /** \brief Molar volume of ideal gas at standard T and P in 
      cm^3 / mol (CODATA 2014 value)
  */
  const double standard_gas_volume=2.2710947e4;
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
  /// Bohr radius in cm
  const double bohr_radius=o2scl_cgs::bohr_radius;
  /// Stefan-Boltzmann constant in g / K^4 s^3
  const double stefan_boltzmann_constant=
    o2scl_cgs::stefan_boltzmann_constant;
  /// Thomson cross section in cm^2
  const double thomson_cross_section=o2scl_cgs::thomson_cross_section;
  /// Fermi coupling constant in s^4 / cm^4 g^2
  const double Gfermi=o2scl_cgs::Gfermi;
  /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2010 value)
  const double boltzmann=o2scl_cgs::boltzmann;
  //@}

  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in cm
  const double astronomical_unit=o2scl_cgs::astronomical_unit;
  /// Parsec in \f$ \mathrm{cm} \f$
  const double parsec=o2scl_cgs::parsec;
  /// Acceleration due to gravity in cm / s^2
  const double grav_accel=o2scl_cgs::grav_accel;
  /// Solar mass times gravitational constant in cm^3 / s^2
  const double solar_mass_parameter=o2scl_cgs::solar_mass_parameter;
  /// Solar mass in g
  const double solar_mass=o2scl_cgs::solar_mass;
  /// Schwarzchild radius in cm
  const double schwarzchild_radius=o2scl_cgs::schwarzchild_radius;
  /// Sidereal year in s 
  const double sidereal_year=o2scl_cgs::sidereal_year;
  /// Tropical year in s 
  const double tropical_year=o2scl_cgs::tropical_year;
  /// Julian year in s 
  const double julian_year=o2scl_cgs::julian_year;
  /// Light year in \f$ \mathrm{cm} \f$
  const double light_year=o2scl_cgs::light_year;
  //@}

  /// \name Particle masses
  //@{
  /// Electron mass in g
  const double mass_electron=o2scl_cgs::mass_electron;
  /// Muon mass in g
  const double mass_muon=o2scl_cgs::mass_muon;
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
  /// Electron magnetic moment in abamp cm^2 (CODATA 2014 value)
  const double electron_magnetic_moment=9.284764620e-21;
  /// Proton magnetic moment in abamp cm^2 (CODATA 2014 value)
  const double proton_magnetic_moment=1.4106067873e-23;
  /// Roentgen abamp s / g
  const double roentgen=o2scl_cgs::roentgen/10.0;
  /// Bohr magneton in abamp cm^2 (CODATA 2014 value)
  const double bohr_magneton=9.274009994e-21;
  /// Nuclear magneton in abamp cm^2 (CODATA 2014 value)
  const double nuclear_magneton=5.050783699e-24;
  /// Faraday constant in abamp s / mol (CODATA 2014 value)
  const double faraday=9.648533289e3;
  /// Electron charge in abamp s (derived)
  const double electron_charge=electron_volt/10.0;
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
  /// Bohr radius in m
  const double bohr_radius=o2scl_cgs::bohr_radius/1.0e2;
  /// Thomson cross section in m^2
  const double thomson_cross_section=o2scl_cgs::thomson_cross_section/1.0e4;
  /// Fermi coupling constant in s^4 / m^4 kg^2 (FIXME)
  const double Gfermi=o2scl_cgs::Gfermi;
  /// Boltzmann constant in kg m^2 / K s^2
  const double boltzmann=o2scl_cgs::boltzmann/1.0e7;
  /// Stefan-Boltzmann constant in kg / K^4 s^3
  const double stefan_boltzmann_constant=
    o2scl_cgs::stefan_boltzmann_constant/1.0e3;
  //@{

  /// \name Astrophysical constants
  //@{
  /// Astronomical unit in m 
  const double astronomical_unit=o2scl_cgs::astronomical_unit/1.0e2;
  /// Parsec in m
  const double parsec=o2scl_cgs::parsec/1.0e2;
  /// Acceleration due to gravity in m / s^2
  const double grav_accel=o2scl_cgs::grav_accel/1.0e2;
  /// Solar mass times gravitational constant in m^3 / s^2
  const double solar_mass_parameter=o2scl_cgs::solar_mass_parameter/1.0e6;
  /// Mass of the sun in kg 
  const double solar_mass=o2scl_cgs::solar_mass/1.0e3;
  /// Schwarzchild radius in m
  const double schwarzchild_radius=o2scl_cgs::schwarzchild_radius/1.0e2;
  /// Sidereal year in s 
  const double sidereal_year=365.256363004*8.64e4;
  /// Tropical year in s 
  const double tropical_year=365.242190402*8.64e4;
  /// Julian year in s 
  const double julian_year=365.25*8.64e4;
  /// Light year in \f$ \mathrm{m} \f$
  const double light_year=o2scl_cgs::light_year/1.0e2;
  //@}

  /// \name Particle masses
  //@{
  /// Electron mass in kg
  const double mass_electron=o2scl_cgs::mass_electron/1.0e3;
  /// Muon mass in kg
  const double mass_muon=o2scl_cgs::mass_muon/1.0e3;
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
  const double minute=6e1;
  /// s
  const double hour=3.6e3;
  /// s
  const double day=8.64e4;
  /// s
  const double week=6.048e5;
  /// m
  const double inch=2.54e-2;
  /// m
  const double foot=3.048e-1;
  /// m
  const double yard=9.144e-1;
  /// m
  const double mile=1.609344e3;
  /// m
  const double nautical_mile=1.852e3;
  /// m
  const double fathom=1.8288e0;
  /// m
  const double mil=2.54e-5;
  /// m
  const double point=3.52777777778e-4;
  /// m
  const double texpoint=3.51459803515e-4;
  /// m
  const double micron=1e-6;
  /// m
  const double angstrom=1e-10;
  /// m^2
  const double hectare=1e4;
  /// m^2
  const double acre=4.04685642241e3;
  /// m^2
  const double barn=1e-28;
  /// m^3
  const double liter=1e-3;
  /// m^3
  const double us_gallon=3.78541178402e-3;
  /// m^3
  const double quart=9.46352946004e-4;
  /// m^3
  const double pint=4.73176473002e-4;
  /// m^3
  const double cup=2.36588236501e-4;
  /// m^3
  const double fluid_ounce=2.95735295626e-5;
  /// m^3
  const double tablespoon=1.47867647813e-5;
  /// m^3
  const double teaspoon=4.92892159375e-6;
  /// m^3
  const double canadian_gallon=4.54609e-3;
  /// m^3
  const double uk_gallon=4.546092e-3;
  /// m / s
  const double miles_per_hour=4.4704e-1;
  /// m / s
  const double kilometers_per_hour=2.77777777778e-1;
  /// m / s
  const double knot=5.14444444444e-1;
  /// kg
  const double pound_mass=4.5359237e-1;
  /// kg
  const double ounce_mass=2.8349523125e-2;
  /// kg
  const double ton=9.0718474e2;
  /// kg
  const double metric_ton=1e3;
  /// kg
  const double uk_ton=1.0160469088e3;
  /// kg
  const double troy_ounce=3.1103475e-2;
  /// kg
  const double carat=2e-4;
  /// kg m / s^2
  const double gram_force=9.80665e-3;
  /// kg m / s^2
  const double pound_force=4.44822161526e0;
  /// kg m / s^2
  const double kilopound_force=4.44822161526e3;
  /// kg m / s^2
  const double poundal=1.38255e-1;
  /// kg m^2 / s^2
  const double calorie=4.1868e0;
  /// kg m^2 / s^2
  const double btu=1.05505585262e3;
  /// kg m^2 / s^2
  const double therm=1.05506e8;
  /// kg m^2 / s^3
  const double horsepower=7.457e2;
  /// kg / m s^2
  const double bar=1e5;
  /// kg / m s^2
  const double std_atmosphere=1.01325e5;
  /// kg / m s^2
  const double torr=1.33322368421e2;
  /// kg / m s^2
  const double meter_of_mercury=1.33322368421e5;
  /// kg / m s^2
  const double inch_of_mercury=3.38638815789e3;
  /// kg / m s^2
  const double inch_of_water=2.490889e2;
  /// kg / m s^2
  const double psi=6.89475729317e3;
  /// kg m^-1 s^-1
  const double poise=1e-1;
  /// m^2 / s
  const double stokes=1e-4;
  /// kg / A s^2
  const double gauss=1e-4;
  /// cd / m^2
  const double stilb=1e4;
  /// cd sr
  const double lumen=1e0;
  /// cd sr / m^2
  const double lux=1e0;
  /// cd sr / m^2
  const double phot=1e4;
  /// cd sr / m^2
  const double footcandle=1.076e1;
  /// cd sr / m^2
  const double lambert=1e4;
  /// cd sr / m^2
  const double footlambert=1.07639104e1;
  /// 1 / s
  const double curie=3.7e10;
  /// m^2 / s^2
  const double rad=1e-2;
  /// kg m / s^2
  const double newton=1e0;
  /// kg m / s^2
  const double dyne=1e-5;
  /// kg m^2 / s^2
  const double joule=1e0;
  /// kg m^2 / s^2
  const double erg=1e-7;
  //@}
  
  /// \name ELectromagnetic constants
  //@{
  /// A m^2 (CODATA 2014 value)
  const double electron_magnetic_moment=9.284764620e-24;
  /// A m^2 (CODATA 2014 value)
  const double proton_magnetic_moment=1.4106067873e-26;
  /// A s / kg
  const double roentgen=2.58e-4;
  /// Bohr magneton in A m^2 (CODATA 2014 value)
  const double bohr_magneton=9.274009994e-24;
  /// A m^2 (CODATA 2014 value)
  const double nuclear_magneton=5.050783699e-27;
  /// A^2 s^4 / kg m^3 (derived)
  const double vacuum_permittivity=1.0/o2scl_mks::speed_of_light/
    o2scl_mks::speed_of_light/4.0e-7/o2scl_const::pi;
  /// kg m / A^2 s^2
  const double vacuum_permeability=1.25663706144e-6;
  /// A s / mol
  const double faraday=o2scl_cgsm::faraday*10.0;
  /// A s (derived)
  const double electron_charge=o2scl_mks::electron_volt;
  //@}

}

/** \brief Constants in MKSA units
    
    Where possible, constants here are defined in terms of the values
    in \ref o2scl_cgs, in order to make it easier to update these
    values. See also the documentation at \ref o2scl_cgs .
*/
namespace o2scl_mksa {
  /// m / s
  const double speed_of_light=2.99792458e8;
  /// Newtonian constant of gravitation in m^3 / kg s^2 (CODATA 2010 value)
  const double gravitational_constant=6.67384e-11;
  /// Planck constant in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_h=6.62606957e-34;
  /// Planck constant divided by 2 pi in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_hbar=o2scl_mksa::plancks_constant_h/
    2.0/o2scl_const::pi;
  
  /// m
  const double schwarzchild_radius=2.95325008e3;
  /// Astronomical unit in m (IAU 2009 value; now exact)
  const double astronomical_unit=1.49597870700e11;
  /// m
  const double light_year=9.46053620707e15;
  /// m
  const double parsec=3.08567758135e16;
  /// m / s^2
  const double grav_accel=9.80665e0;
  /// Electron volt in kg m^2 / s^2 (CODATA 2010 value)
  const double electron_volt=1.602176565e-19;
  /// Electron mass in kg (CODATA 2010 value)
  const double mass_electron=9.10938291e-31;
  /// Muon mass in kg (CODATA 2010 value)
  const double mass_muon=1.883531475e-28;
  /// Proton mass in kg (CODATA 2010 value)
  const double mass_proton=1.672621777e-27;
  /// Neutron mass in kg (CODATA 2010 value)
  const double mass_neutron=1.674927351e-27;

  /// Deuteron mass in kg (CODATA 2010 value)
  const double mass_deuteron=3.34358348e-27;
  /// Triton mass in kg (CODATA 2010 value)
  const double mass_triton=5.00735630e-27;
  /// Helion mass in kg (CODATA 2010 value)
  const double mass_helion=5.00641234e-27;
  /// Alpha particle mass in kg (CODATA 2010 value)
  const double mass_alpha=6.64465675e-27;

  /// Rydberg constant in kg m^2 / s^2 (CODATA 2010 value)
  const double rydberg=2.179872171e-18;
  /// Boltzmann constant in kg m^2 / K s^2 (CODATA 2010 value)
  const double boltzmann=1.3806488e-23;
  /// Bohr magneton in A m^2 (CODATA 2010 value)
  const double bohr_magneton=9.27400968e-24;
  /// A m^2
  const double nuclear_magneton=5.05078317e-27;
  /// A m^2
  const double electron_magnetic_moment=9.28476362e-24;
  /// A m^2
  const double proton_magnetic_moment=1.410606633e-26;
  /// kg m^2 / K mol s^2
  const double molar_gas=8.314472e0;
  /// m^3 / mol
  const double standard_gas_volume=2.2710981e-2;
  /// s
  const double minute=6e1;
  /// s
  const double hour=3.6e3;
  /// s
  const double day=8.64e4;
  /// s
  const double week=6.048e5;
  /// m
  const double inch=2.54e-2;
  /// m
  const double foot=3.048e-1;
  /// m
  const double yard=9.144e-1;
  /// m
  const double mile=1.609344e3;
  /// m
  const double nautical_mile=1.852e3;
  /// m
  const double fathom=1.8288e0;
  /// m
  const double mil=2.54e-5;
  /// m
  const double point=3.52777777778e-4;
  /// m
  const double texpoint=3.51459803515e-4;
  /// m
  const double micron=1e-6;
  /// m
  const double angstrom=1e-10;
  /// m^2
  const double hectare=1e4;
  /// m^2
  const double acre=4.04685642241e3;
  /// m^2
  const double barn=1e-28;
  /// m^3
  const double liter=1e-3;
  /// m^3
  const double us_gallon=3.78541178402e-3;
  /// m^3
  const double quart=9.46352946004e-4;
  /// m^3
  const double pint=4.73176473002e-4;
  /// m^3
  const double cup=2.36588236501e-4;
  /// m^3
  const double fluid_ounce=2.95735295626e-5;
  /// m^3
  const double tablespoon=1.47867647813e-5;
  /// m^3
  const double teaspoon=4.92892159375e-6;
  /// m^3
  const double canadian_gallon=4.54609e-3;
  /// m^3
  const double uk_gallon=4.546092e-3;
  /// m / s
  const double miles_per_hour=4.4704e-1;
  /// m / s
  const double kilometers_per_hour=2.77777777778e-1;
  /// m / s
  const double knot=5.14444444444e-1;
  /// kg
  const double pound_mass=4.5359237e-1;
  /// kg
  const double ounce_mass=2.8349523125e-2;
  /// kg
  const double ton=9.0718474e2;
  /// kg
  const double metric_ton=1e3;
  /// kg
  const double uk_ton=1.0160469088e3;
  /// kg
  const double troy_ounce=3.1103475e-2;
  /// kg
  const double carat=2e-4;
  /// Atomic mass constant in kg (CODATA 2010 value)
  const double unified_atomic_mass=1.660538921e-27;
  /// kg m / s^2
  const double gram_force=9.80665e-3;
  /// kg m / s^2
  const double pound_force=4.44822161526e0;
  /// kg m / s^2
  const double kilopound_force=4.44822161526e3;
  /// kg m / s^2
  const double poundal=1.38255e-1;
  /// kg m^2 / s^2
  const double calorie=4.1868e0;
  /// kg m^2 / s^2
  const double btu=1.05505585262e3;
  /// kg m^2 / s^2
  const double therm=1.05506e8;
  /// kg m^2 / s^3
  const double horsepower=7.457e2;
  /// kg / m s^2
  const double bar=1e5;
  /// kg / m s^2
  const double std_atmosphere=1.01325e5;
  /// kg / m s^2
  const double torr=1.33322368421e2;
  /// kg / m s^2
  const double meter_of_mercury=1.33322368421e5;
  /// kg / m s^2
  const double inch_of_mercury=3.38638815789e3;
  /// kg / m s^2
  const double inch_of_water=2.490889e2;
  /// kg / m s^2
  const double psi=6.89475729317e3;
  /// kg m^-1 s^-1
  const double poise=1e-1;
  /// m^2 / s
  const double stokes=1e-4;
  /// A s / mol
  const double faraday=9.64853429775e4;
  /// A s
  const double electron_charge=1.602176565e-19;
  /// kg / A s^2
  const double gauss=1e-4;
  /// cd / m^2
  const double stilb=1e4;
  /// cd sr
  const double lumen=1e0;
  /// cd sr / m^2
  const double lux=1e0;
  /// cd sr / m^2
  const double phot=1e4;
  /// cd sr / m^2
  const double footcandle=1.076e1;
  /// cd sr / m^2
  const double lambert=1e4;
  /// cd sr / m^2
  const double footlambert=1.07639104e1;
  /// 1 / s
  const double curie=3.7e10;
  /// A s / kg
  const double roentgen=2.58e-4;
  /// m^2 / s^2
  const double rad=1e-2;
  /// kg
  const double solar_mass=1.9884e30;
  /// m
  const double bohr_radius=5.291772083e-11;
  /// kg m / s^2
  const double newton=1e0;
  /// kg m / s^2
  const double dyne=1e-5;
  /// kg m^2 / s^2
  const double joule=1e0;
  /// kg m^2 / s^2
  const double erg=1e-7;
  /// kg / K^4 s^3
  const double stefan_boltzmann_constant=5.67039934436e-8;
  /// m^2
  const double thomson_cross_section=6.65245853542e-29;
  /// A^2 s^4 / kg m^3
  const double vacuum_permittivity=8.854187817e-12;
  /// kg m / A^2 s^2
  const double vacuum_permeability=1.25663706144e-6;
  /** \brief Fermi coupling constant in s^4 / cm^4 g^2, 
      defined as \f$ 1.166364 \times 10^{-5}~\mathrm{GeV}^{-2} \f$
      (CODATA 2010 value)
  */
  const double Gfermi=1.166364e-23/electron_volt/electron_volt;
}

// Other derived values to add to the namespace
namespace o2scl_const {

  /// \f$ \hbar c \f$ in MeV fm (derived)
  const double hc_mev_fm=o2scl_mks::plancks_constant_hbar*
    o2scl_mks::speed_of_light/o2scl_mks::electron_volt*1.0e9;

  /// \f$ \hbar c \f$ in MeV cm (derived)
  const double hc_mev_cm=hc_mev_fm*1.0e-13;

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

  /** \brief Electron charge squared in SI(MKSA) units (derived)

      In MKSA units:
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
  const double e2_mksa=o2scl_mksa::electron_charge;
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


  /** \brief \f$ \Lambda \f$ hyperon mass in \f$ \mathrm{MeV} \f$
      
      Value from PDG live (5/1/14).
  */
  const double mass_lambda_MeV=1115.683;

  /** \brief \f$ \Sigma \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_minus_MeV=1197.449;

  /** \brief \f$ \Sigma \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_zero_MeV=1192.642;

  /** \brief \f$ \Sigma \f$ hyperon mass in \f$ \mathrm{MeV} \f$
   */
  const double mass_sigma_plus_MeV=1189.37;
}


#endif
