/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
    
    CODATA 2010 values are in \ref Mohr12.
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

  /** \brief Fine structure constant (CODATA 2010 value)
  */
  const double fine_structure=7.2973525698e-3;
  /** \brief Avogadro's number (CODATA 2010 value)
  */
  const double avogadro=6.02214129e23;

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

  /// \f$ \sin^2 \theta_W \f$ (CODATA 2010 value)
  const double sin2_theta_weak=0.2224;
}
  
/** \brief Constants in CGS units 
    
    CODATA 2010 values are in \ref Mohr12. IAU 2009 values
    are from \ref Luzum11 . Solar mass from 
    http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf
*/
namespace o2scl_cgs {
  /// cm
  const double schwarzchild_radius=2.95325008e5;
  /// cm / s
  const double speed_of_light=2.99792458e10;
  /// Newtonian constant of gravitation in cm^3 / g s^2 (CODATA 2010 value)
  const double gravitational_constant=6.67384e-8;
  /// Planck constant in g cm^2 / s (CODATA 2010 value)
  const double plancks_constant_h=6.62606957e-27;
  /// Planck constant divided by 2 pi in g cm^2 / s (derived)
  const double plancks_constant_hbar=o2scl_cgs::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Astronomical unit in cm (IAU 2009 value)
  const double astronomical_unit=1.49597870700e13;
  /// cm
  const double light_year=9.46053620707e17;
  /// cm
  const double parsec=3.08567758135e18;
  /// cm / s^2
  const double grav_accel=9.80665e2;
  /// Electron volt in g cm^2 / s^2 (CODATA 2010 value)
  const double electron_volt=1.602176565e-12;
  /// Electron mass in g (CODATA 2010 value)
  const double mass_electron=9.10938291e-28;
  /// Muon mass in g (CODATA 2010 value)
  const double mass_muon=1.883531475e-25;
  /// Proton mass in g (CODATA 2010 value)
  const double mass_proton=1.672621777e-24;
  /// Neutron mass in g (CODATA 2010 value)
  const double mass_neutron=1.674927351e-24;

  /// Deuteron mass in kg (CODATA 2010 value)
  const double mass_deuteron=3.34358348e-24;
  /// Triton mass in kg (CODATA 2010 value)
  const double mass_triton=5.00735630e-24;
  /// Helion mass in kg (CODATA 2010 value)
  const double mass_helion=5.00641234e-24;
  /// Alpha particle mass in kg (CODATA 2010 value)
  const double mass_alpha=6.64465675e-24;

  /// Rydberg constant in g cm^2 / s^2 (CODATA 2010 value)
  const double rydberg=2.179872171e-11;
  /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2010 value)
  const double boltzmann=1.3806488e-16;
  /// g cm^2 / K mol s^2
  const double molar_gas=8.314472e7;
  /// cm^3 / mol
  const double standard_gas_volume=2.2710981e4;
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
  const double point=3.52777777778e-2;
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
  const double kilometers_per_hour=2.77777777778e1;
  /// cm / s
  const double knot=5.14444444444e1;
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
  /// Atomic mass constant in g (CODATA 2010 value)
  const double unified_atomic_mass=1.660538921e-24;
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
  /// A s / g
  const double roentgen=2.58e-7;
  /// cm^2 / s^2
  const double rad=1e2;
  /// g
  const double solar_mass=1.9884e33;
  /// cm
  const double bohr_radius=5.291772083e-9;
  /// cm g / s^2
  const double newton=1e5;
  /// cm g / s^2
  const double dyne=1e0;
  /// g cm^2 / s^2
  const double joule=1e7;
  /// g cm^2 / s^2
  const double erg=1e0;
  /// g / K^4 s^3
  const double stefan_boltzmann_constant=5.67039934436e-5;
  /// cm^2
  const double thomson_cross_section=6.65245853542e-25;
  /** \brief Fermi coupling constant in s^4 / cm^4 g^2, 
      defined as \f$ 1.166364 \times 10^{-5}~\mathrm{GeV}^{-2} \f$
      (CODATA 2010 value)
  */
  const double Gfermi=1.166364e-23/electron_volt/electron_volt;
}

/** \brief Constants in CGSM units
    
    CODATA 2010 values are in \ref Mohr12. IAU 2009 values
    are from \ref Luzum11 . Solar mass from 
    http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf
*/
namespace o2scl_cgsm {
  /// cm
  const double schwarzchild_radius=2.95325008e5;
  /// cm / s
  const double speed_of_light=2.99792458e10;
  /// Newtonian constant of gravitation in cm^3 / g s^2 (CODATA 2010 value)
  const double gravitational_constant=6.67384e-8;
  /// Planck constant in g cm^2 / s (CODATA 2010 value)
  const double plancks_constant_h=6.62606957e-27;
  /// Planck constant divided by 2 pi in g cm^2 / s (derived)
  const double plancks_constant_hbar=o2scl_cgsm::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Astronomical unit in cm (IAU 2009 value)
  const double astronomical_unit=1.49597870700e13;
  /// cm
  const double light_year=9.46053620707e17;
  /// cm
  const double parsec=3.08567758135e18;
  /// cm / s^2
  const double grav_accel=9.80665e2;
  /// Electron volt in g cm^2 / s^2 (CODATA 2010 value)
  const double electron_volt=1.602176565e-12;
  /// Electron mass in g (CODATA 2010 value)
  const double mass_electron=9.10938291e-28;
  /// Muon mass in g (CODATA 2010 value)
  const double mass_muon=1.883531475e-25;
  /// Proton mass in g (CODATA 2010 value)
  const double mass_proton=1.672621777e-24;
  /// Neutron mass in g (CODATA 2010 value)
  const double mass_neutron=1.674927351e-24;

  /// Deuteron mass in kg (CODATA 2010 value)
  const double mass_deuteron=3.34358348e-24;
  /// Triton mass in kg (CODATA 2010 value)
  const double mass_triton=5.00735630e-24;
  /// Helion mass in kg (CODATA 2010 value)
  const double mass_helion=5.00641234e-24;
  /// Alpha particle mass in kg (CODATA 2010 value)
  const double mass_alpha=6.64465675e-24;

  /// Rydberg constant in g cm^2 / s^2 (CODATA 2010 value)
  const double rydberg=2.179872171e-11;
  /// Boltzmann constant in g cm^2 / K s^2 (CODATA 2010 value)
  const double boltzmann=1.3806488e-16;
  /// abamp cm^2
  const double electron_magnetic_moment=9.28476362e-21;
  /// abamp cm^2
  const double proton_magnetic_moment=1.410606633e-23;
  /// g cm^2 / K mol s^2
  const double molar_gas=8.314472e7;
  /// cm^3 / mol
  const double standard_gas_volume=2.2710981e4;
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
  const double point=3.52777777778e-2;
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
  const double kilometers_per_hour=2.77777777778e1;
  /// cm / s
  const double knot=5.14444444444e1;
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
  /// Atomic mass constant in g (CODATA 2010 value)
  const double unified_atomic_mass=1.660538921e-24;
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
  /// abamp s / mol
  const double faraday=9.64853429775e3;
  /// abamp s
  const double electron_charge=1.602176565e-20;
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
  /// abamp s / g
  const double roentgen=2.58e-8;
  /// cm^2 / s^2
  const double rad=1e2;
  /// g
  const double solar_mass=1.9884e33;
  /// cm
  const double bohr_radius=5.291772083e-9;
  /// cm g / s^2
  const double newton=1e5;
  /// cm g / s^2
  const double dyne=1e0;
  /// g cm^2 / s^2
  const double joule=1e7;
  /// g cm^2 / s^2
  const double erg=1e0;
  /// g / K^4 s^3
  const double stefan_boltzmann_constant=5.67039934436e-5;
  /// cm^2
  const double thomson_cross_section=6.65245853542e-25;
  /// Bohr magneton in abamp cm^2 (CODATA 2010 value)
  const double bohr_magneton=9.27400968e-21;
  /// abamp cm^2
  const double nuclear_magneton=5.05078317e-24;
  /** \brief Fermi coupling constant in s^4 / cm^4 g^2, 
      defined as \f$ 1.166364 \times 10^{-5}~\mathrm{GeV}^{-2} \f$
      (CODATA 2010 value)
  */
  const double Gfermi=1.166364e-23/electron_volt/electron_volt;
}

/** \brief Constants in MKS units
    
    CODATA 2010 values are in \ref Mohr12. IAU 2009 values
    are from \ref Luzum11 . Solar mass from 
    http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf
*/
namespace o2scl_mks {
  /// m
  const double schwarzchild_radius=2.95325008e3;
  /// m / s
  const double speed_of_light=2.99792458e8;
  /// Newtonian constant of gravitation in m^3 / kg s^2 (CODATA 2010 value)
  const double gravitational_constant=6.67384e-11;
  /// Planck constant in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_h=6.62606957e-34;
  /// Planck constant divided by 2 pi in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_hbar=o2scl_mks::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Astronomical unit in m (IAU 2009 value)
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

/** \brief Constants in MKSA units
    
    CODATA 2010 values are in \ref Mohr12. IAU 2009 values
    are from \ref Luzum11 . Solar mass from 
    http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf
*/
namespace o2scl_mksa {
  /// m
  const double schwarzchild_radius=2.95325008e3;
  /// m / s
  const double speed_of_light=2.99792458e8;
  /// Newtonian constant of gravitation in m^3 / kg s^2 (CODATA 2010 value)
  const double gravitational_constant=6.67384e-11;
  /// Planck constant in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_h=6.62606957e-34;
  /// Planck constant divided by 2 pi in kg m^2 / s (CODATA 2010 value)
  const double plancks_constant_hbar=o2scl_mksa::plancks_constant_h/
    2.0/o2scl_const::pi;
  /// Astronomical unit in m (IAU 2009 value)
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
}


#endif
