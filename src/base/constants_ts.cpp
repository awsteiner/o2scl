/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2019-2023, Andrew W. Steiner
  
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
#include <gsl/gsl_const_cgs.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_mks.h>
#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <o2scl/string_conv.h>
#include <climits>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  t.test_rel(GSL_CONST_CGS_SPEED_OF_LIGHT,
	     o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_cgs),
             1.0e-7,"CGS speed_of_light");
  t.test_rel(GSL_CONST_MKS_SPEED_OF_LIGHT,
	     o2scl_const::speed_of_light_f<double>(o2scl_const::o2scl_mks),
             1.0e-7,"MKS speed_of_light");
  t.test_rel(GSL_CONST_CGS_GRAVITATIONAL_CONSTANT,
	     o2scl_const::gravitational_constant_f<double>
             (o2scl_const::o2scl_cgs),3.0e-4,
	     "CGS gravitational_constant");
  t.test_rel(GSL_CONST_MKS_GRAVITATIONAL_CONSTANT,
	     o2scl_const::gravitational_constant_f<double>
             (o2scl_const::o2scl_mks),3.0e-4,
	     "MKS gravitational_constant");
  t.test_rel(GSL_CONST_CGS_PLANCKS_CONSTANT_H,
             o2scl_const::planck_f<double>(o2scl_const::o2scl_cgs),3.0e-7,
             "CGS plancks_constant_h");
  t.test_rel(GSL_CONST_MKS_PLANCKS_CONSTANT_H,
             o2scl_const::planck_f<double>(o2scl_const::o2scl_mks),3.0e-7,
             "MKS plancks_constant_h");
  t.test_rel(GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR,
	     o2scl_const::hbar_f<double>(o2scl_const::o2scl_cgs),3.0e-7,
	     "CGS plancks_constant_hbar");
  t.test_rel(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR,
	     o2scl_const::hbar_f<double>(o2scl_const::o2scl_mks),3.0e-7,
	     "MKS plancks_constant_hbar");
  t.test_rel(GSL_CONST_CGS_ELECTRON_VOLT,
	     o2scl_const::electron_volt_f<double>(o2scl_const::o2scl_cgs),
             1.0e-7,"CGS electron_volt");
  t.test_rel(GSL_CONST_MKS_ELECTRON_VOLT,
	     o2scl_const::electron_volt_f<double>(o2scl_const::o2scl_mks),
             1.0e-7,"MKS electron_volt");
  t.test_rel(GSL_CONST_CGS_BOHR_RADIUS,
	     o2scl_const::bohr_radius_f<double>(o2scl_const::o2scl_cgs),
             1.0e-7,"CGS bohr_radius");
  t.test_rel(GSL_CONST_MKS_BOHR_RADIUS,
	     o2scl_const::bohr_radius_f<double>(o2scl_const::o2scl_mks),
             1.0e-7,"MKS bohr_radius");
  t.test_rel(GSL_CONST_CGS_STEFAN_BOLTZMANN_CONSTANT,
	     o2scl_const::stefan_boltz_cons_f<double>
             (o2scl_const::o2scl_cgs),1.0e-5,"CGS stefan_boltzmann");
  t.test_rel(GSL_CONST_MKS_STEFAN_BOLTZMANN_CONSTANT,
	     o2scl_const::stefan_boltz_cons_f<double>
             (o2scl_const::o2scl_mks),1.0e-5,"MKS stefan_boltzmann");
  t.test_rel(GSL_CONST_CGS_THOMSON_CROSS_SECTION,
	     o2scl_const::thomson_csec_f<double>(o2scl_const::o2scl_cgs),
             4.0e-7,"CGS thomson_cross_section");
  t.test_rel(GSL_CONST_MKS_THOMSON_CROSS_SECTION,
	     o2scl_const::thomson_csec_f<double>(o2scl_const::o2scl_mks),
             4.0e-7,"MKS thomson_cross_section");
  t.test_rel(GSL_CONST_CGS_BOLTZMANN,
	     o2scl_const::boltzmann_f<double>(o2scl_const::o2scl_cgs),
             4.0e-6,"CGS boltzmann");
  t.test_rel(GSL_CONST_MKS_BOLTZMANN,
	     o2scl_const::boltzmann_f<double>(o2scl_const::o2scl_mks),
             4.0e-6,"MKS boltzmann");
  t.test_rel(GSL_CONST_CGS_ASTRONOMICAL_UNIT,
	     o2scl_const::astronomical_unit_f<double>(o2scl_const::o2scl_cgs),
             1.0e-10,"CGS astronomical_unit");
  t.test_rel(GSL_CONST_MKS_ASTRONOMICAL_UNIT,
	     o2scl_const::astronomical_unit_f<double>(o2scl_const::o2scl_mks),
             1.0e-10,"MKS astronomical_unit");

  t.report();
  
  exit(-1);
  
  t.test_rel(GSL_CONST_CGS_PARSEC,
	     o2scl_cgs::parsec,1.0e-7,
	     "CGS parsec");
  t.test_rel(GSL_CONST_CGS_GRAV_ACCEL,
	     o2scl_cgs::grav_accel,1.0e-7,
	     "CGS grav_accel");
  t.test_rel(GSL_CONST_CGS_SOLAR_MASS,
	     o2scl_cgs::solar_mass,5.0e-4,
	     "CGS solar_mass");
  t.test_rel(GSL_CONST_CGS_LIGHT_YEAR,
	     o2scl_cgs::light_year,4.0e-5,
	     "CGS light_year");
  t.test_rel(GSL_CONST_CGS_MASS_ELECTRON,
	     o2scl_cgs::mass_electron,4.0e-7,
	     "CGS mass_electron");
  t.test_rel(GSL_CONST_CGS_MASS_MUON,
	     o2scl_cgs::mass_muon,4.0e-7,
	     "CGS mass_muon");
  t.test_rel(GSL_CONST_CGS_MASS_PROTON,
	     o2scl_cgs::mass_proton,4.0e-7,
	     "CGS mass_proton");
  t.test_rel(GSL_CONST_CGS_MASS_NEUTRON,
	     o2scl_cgs::mass_neutron,4.0e-7,
	     "CGS mass_neutron");
  t.test_rel(GSL_CONST_CGS_UNIFIED_ATOMIC_MASS,
	     o2scl_cgs::unified_atomic_mass,4.0e-7,
	     "CGS unified_atomic_mass");
  t.test_rel(GSL_CONST_CGS_RYDBERG,
	     o2scl_cgs::rydberg,4.0e-7,
	     "CGS rydberg");
  t.test_rel(GSL_CONST_CGS_MOLAR_GAS,
	     o2scl_cgs::molar_gas,4.0e-6,
	     "CGS molar_gas");
  t.test_rel(GSL_CONST_CGS_STANDARD_GAS_VOLUME,
	     o2scl_cgs::standard_gas_volume,4.0e-6,
	     "CGS standard_gas_volume");
  t.test_rel(GSL_CONST_CGS_MINUTE,
	     o2scl_cgs::minute,1.0e-7,
	     "CGS minute");
  t.test_rel(GSL_CONST_CGS_HOUR,
	     o2scl_cgs::hour,1.0e-7,
	     "CGS hour");
  t.test_rel(GSL_CONST_CGS_DAY,
	     o2scl_cgs::day,1.0e-7,
	     "CGS day");
  t.test_rel(GSL_CONST_CGS_WEEK,
	     o2scl_cgs::week,1.0e-7,
	     "CGS week");
  t.test_rel(GSL_CONST_CGS_INCH,
	     o2scl_cgs::inch,1.0e-7,
	     "CGS inch");
  t.test_rel(GSL_CONST_CGS_FOOT,
	     o2scl_cgs::foot,1.0e-7,
	     "CGS foot");
  t.test_rel(GSL_CONST_CGS_YARD,
	     o2scl_cgs::yard,1.0e-7,
	     "CGS yard");
  t.test_rel(GSL_CONST_CGS_MILE,
	     o2scl_cgs::mile,1.0e-7,
	     "CGS mile");
  t.test_rel(GSL_CONST_CGS_NAUTICAL_MILE,
	     o2scl_cgs::nautical_mile,1.0e-7,
	     "CGS nautical_mile");
  t.test_rel(GSL_CONST_CGS_FATHOM,
	     o2scl_cgs::fathom,1.0e-7,
	     "CGS fathom");
  t.test_rel(GSL_CONST_CGS_MIL,
	     o2scl_cgs::mil,1.0e-7,
	     "CGS mil");
  t.test_rel(GSL_CONST_CGS_POINT,
	     o2scl_cgs::point,1.0e-7,
	     "CGS point");
  t.test_rel(GSL_CONST_CGS_TEXPOINT,
	     o2scl_cgs::texpoint,1.0e-7,
	     "CGS texpoint");
  t.test_rel(GSL_CONST_CGS_MICRON,
	     o2scl_cgs::micron,1.0e-7,
	     "CGS micron");
  t.test_rel(GSL_CONST_CGS_ANGSTROM,
	     o2scl_cgs::angstrom,1.0e-7,
	     "CGS angstrom");
  t.test_rel(GSL_CONST_CGS_HECTARE,
	     o2scl_cgs::hectare,1.0e-7,
	     "CGS hectare");
  t.test_rel(GSL_CONST_CGS_ACRE,
	     o2scl_cgs::acre,1.0e-7,
	     "CGS acre");
  t.test_rel(GSL_CONST_CGS_BARN,
	     o2scl_cgs::barn,1.0e-7,
	     "CGS barn");
  t.test_rel(GSL_CONST_CGS_LITER,
	     o2scl_cgs::liter,1.0e-7,
	     "CGS liter");
  t.test_rel(GSL_CONST_CGS_US_GALLON,
	     o2scl_cgs::us_gallon,1.0e-7,
	     "CGS us_gallon");
  t.test_rel(GSL_CONST_CGS_QUART,
	     o2scl_cgs::quart,1.0e-7,
	     "CGS quart");
  t.test_rel(GSL_CONST_CGS_PINT,
	     o2scl_cgs::pint,1.0e-7,
	     "CGS pint");
  t.test_rel(GSL_CONST_CGS_CUP,
	     o2scl_cgs::cup,1.0e-7,
	     "CGS cup");
  t.test_rel(GSL_CONST_CGS_FLUID_OUNCE,
	     o2scl_cgs::fluid_ounce,1.0e-7,
	     "CGS fluid_ounce");
  t.test_rel(GSL_CONST_CGS_TABLESPOON,
	     o2scl_cgs::tablespoon,1.0e-7,
	     "CGS tablespoon");
  t.test_rel(GSL_CONST_CGS_TEASPOON,
	     o2scl_cgs::teaspoon,1.0e-7,
	     "CGS teaspoon");
  t.test_rel(GSL_CONST_CGS_CANADIAN_GALLON,
	     o2scl_cgs::canadian_gallon,1.0e-7,
	     "CGS canadian_gallon");
  t.test_rel(GSL_CONST_CGS_UK_GALLON,
	     o2scl_cgs::uk_gallon,1.0e-7,
	     "CGS uk_gallon");
  t.test_rel(GSL_CONST_CGS_MILES_PER_HOUR,
	     o2scl_cgs::miles_per_hour,1.0e-7,
	     "CGS miles_per_hour");
  t.test_rel(GSL_CONST_CGS_KILOMETERS_PER_HOUR,
	     o2scl_cgs::kilometers_per_hour,1.0e-7,
	     "CGS kilometers_per_hour");
  t.test_rel(GSL_CONST_CGS_KNOT,
	     o2scl_cgs::knot,1.0e-7,
	     "CGS knot");
  t.test_rel(GSL_CONST_CGS_POUND_MASS,
	     o2scl_cgs::pound_mass,1.0e-7,
	     "CGS pound_mass");
  t.test_rel(GSL_CONST_CGS_OUNCE_MASS,
	     o2scl_cgs::ounce_mass,1.0e-7,
	     "CGS ounce_mass");
  t.test_rel(GSL_CONST_CGS_TON,
	     o2scl_cgs::ton,1.0e-7,
	     "CGS ton");
  t.test_rel(GSL_CONST_CGS_METRIC_TON,
	     o2scl_cgs::metric_ton,1.0e-7,
	     "CGS metric_ton");
  t.test_rel(GSL_CONST_CGS_UK_TON,
	     o2scl_cgs::uk_ton,1.0e-7,
	     "CGS uk_ton");
  t.test_rel(GSL_CONST_CGS_TROY_OUNCE,
	     o2scl_cgs::troy_ounce,1.0e-7,
	     "CGS troy_ounce");
  t.test_rel(GSL_CONST_CGS_CARAT,
	     o2scl_cgs::carat,1.0e-7,
	     "CGS carat");
  t.test_rel(GSL_CONST_CGS_GRAM_FORCE,
	     o2scl_cgs::gram_force,1.0e-7,
	     "CGS gram_force");
  t.test_rel(GSL_CONST_CGS_POUND_FORCE,
	     o2scl_cgs::pound_force,1.0e-7,
	     "CGS pound_force");
  t.test_rel(GSL_CONST_CGS_KILOPOUND_FORCE,
	     o2scl_cgs::kilopound_force,1.0e-7,
	     "CGS kilopound_force");
  t.test_rel(GSL_CONST_CGS_POUNDAL,
	     o2scl_cgs::poundal,1.0e-7,
	     "CGS poundal");
  t.test_rel(GSL_CONST_CGS_CALORIE,
	     o2scl_cgs::calorie,1.0e-7,
	     "CGS calorie");
  t.test_rel(GSL_CONST_CGS_BTU,
	     o2scl_cgs::btu,1.0e-7,
	     "CGS btu");
  t.test_rel(GSL_CONST_CGS_THERM,
	     o2scl_cgs::therm,1.0e-7,
	     "CGS therm");
  t.test_rel(GSL_CONST_CGS_HORSEPOWER,
	     o2scl_cgs::horsepower,1.0e-7,
	     "CGS horsepower");
  t.test_rel(GSL_CONST_CGS_BAR,
	     o2scl_cgs::bar,1.0e-7,
	     "CGS bar");
  t.test_rel(GSL_CONST_CGS_STD_ATMOSPHERE,
	     o2scl_cgs::std_atmosphere,1.0e-7,
	     "CGS std_atmosphere");
  t.test_rel(GSL_CONST_CGS_TORR,
	     o2scl_cgs::torr,1.0e-7,
	     "CGS torr");
  t.test_rel(GSL_CONST_CGS_METER_OF_MERCURY,
	     o2scl_cgs::meter_of_mercury,1.0e-7,
	     "CGS meter_of_mercury");
  t.test_rel(GSL_CONST_CGS_INCH_OF_MERCURY,
	     o2scl_cgs::inch_of_mercury,1.0e-7,
	     "CGS inch_of_mercury");
  t.test_rel(GSL_CONST_CGS_INCH_OF_WATER,
	     o2scl_cgs::inch_of_water,1.0e-7,
	     "CGS inch_of_water");
  t.test_rel(GSL_CONST_CGS_PSI,
	     o2scl_cgs::psi,1.0e-7,
	     "CGS psi");
  t.test_rel(GSL_CONST_CGS_POISE,
	     o2scl_cgs::poise,1.0e-7,
	     "CGS poise");
  t.test_rel(GSL_CONST_CGS_STOKES,
	     o2scl_cgs::stokes,1.0e-7,
	     "CGS stokes");
  t.test_rel(GSL_CONST_CGS_STILB,
	     o2scl_cgs::stilb,1.0e-7,
	     "CGS stilb");
  t.test_rel(GSL_CONST_CGS_LUMEN,
	     o2scl_cgs::lumen,1.0e-7,
	     "CGS lumen");
  t.test_rel(GSL_CONST_CGS_LUX,
	     o2scl_cgs::lux,1.0e-7,
	     "CGS lux");
  t.test_rel(GSL_CONST_CGS_PHOT,
	     o2scl_cgs::phot,1.0e-7,
	     "CGS phot");
  t.test_rel(GSL_CONST_CGS_FOOTCANDLE,
	     o2scl_cgs::footcandle,1.0e-7,
	     "CGS footcandle");
  t.test_rel(GSL_CONST_CGS_LAMBERT,
	     o2scl_cgs::lambert,1.0e-7,
	     "CGS lambert");
  t.test_rel(GSL_CONST_CGS_FOOTLAMBERT,
	     o2scl_cgs::footlambert,1.0e-7,
	     "CGS footlambert");
  t.test_rel(GSL_CONST_CGS_CURIE,
	     o2scl_cgs::curie,1.0e-7,
	     "CGS curie");
  t.test_rel(GSL_CONST_CGS_RAD,
	     o2scl_cgs::rad,1.0e-7,
	     "CGS rad");
  t.test_rel(GSL_CONST_CGS_NEWTON,
	     o2scl_cgs::newton,1.0e-7,
	     "CGS newton");
  t.test_rel(GSL_CONST_CGS_DYNE,
	     o2scl_cgs::dyne,1.0e-7,
	     "CGS dyne");
  t.test_rel(GSL_CONST_CGS_JOULE,
	     o2scl_cgs::joule,1.0e-7,
	     "CGS joule");
  t.test_rel(GSL_CONST_CGS_ERG,
	     o2scl_cgs::erg,1.0e-7,
	     "CGS erg");
  t.test_rel(GSL_CONST_CGS_ROENTGEN,
	     o2scl_cgs::roentgen,1.0e-7,
	     "CGS roentgen");

  t.test_rel(GSL_CONST_CGSM_SPEED_OF_LIGHT,
	     o2scl_cgsm::speed_of_light,1.0e-7,
	     "CGSM speed_of_light");
  t.test_rel(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT,
	     o2scl_cgsm::gravitational_constant,3.0e-4,
	     "CGSM gravitational_constant");
  t.test_rel(GSL_CONST_CGSM_PLANCKS_CONSTANT_H,
	     o2scl_cgsm::plancks_constant_h,3.0e-7,
	     "CGSM plancks_constant_h");
  t.test_rel(GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR,
	     o2scl_cgsm::plancks_constant_hbar,3.0e-7,
	     "CGSM plancks_constant_hbar");
  t.test_rel(GSL_CONST_CGSM_ELECTRON_VOLT,
	     o2scl_cgsm::electron_volt,1.0e-7,
	     "CGSM electron_volt");
  t.test_rel(GSL_CONST_CGSM_BOHR_RADIUS,
	     o2scl_cgsm::bohr_radius,1.0e-7,
	     "CGSM bohr_radius");
  t.test_rel(GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT,
	     o2scl_cgsm::stefan_boltzmann_constant,1.0e-5,
	     "CGSM stefan_boltzmann");
  t.test_rel(GSL_CONST_CGSM_THOMSON_CROSS_SECTION,
	     o2scl_cgsm::thomson_cross_section,4.0e-7,
	     "CGSM thomson_cross_section");
  t.test_rel(GSL_CONST_CGSM_BOLTZMANN,
	     o2scl_cgsm::boltzmann,4.0e-6,
	     "CGSM boltzmann");
  t.test_rel(GSL_CONST_CGSM_ASTRONOMICAL_UNIT,
	     o2scl_cgsm::astronomical_unit,1.0e-7,
	     "CGSM astronomical_unit");
  t.test_rel(GSL_CONST_CGSM_PARSEC,
	     o2scl_cgsm::parsec,1.0e-7,
	     "CGSM parsec");
  t.test_rel(GSL_CONST_CGSM_GRAV_ACCEL,
	     o2scl_cgsm::grav_accel,1.0e-7,
	     "CGSM grav_accel");
  t.test_rel(GSL_CONST_CGSM_SOLAR_MASS,
	     o2scl_cgsm::solar_mass,5.0e-4,
	     "CGSM solar_mass");
  t.test_rel(GSL_CONST_CGSM_LIGHT_YEAR,
	     o2scl_cgsm::light_year,4.0e-5,
	     "CGSM light_year");
  t.test_rel(GSL_CONST_CGSM_MASS_ELECTRON,
	     o2scl_cgsm::mass_electron,4.0e-7,
	     "CGSM mass_electron");
  t.test_rel(GSL_CONST_CGSM_MASS_MUON,
	     o2scl_cgsm::mass_muon,4.0e-7,
	     "CGSM mass_muon");
  t.test_rel(GSL_CONST_CGSM_MASS_PROTON,
	     o2scl_cgsm::mass_proton,4.0e-7,
	     "CGSM mass_proton");
  t.test_rel(GSL_CONST_CGSM_MASS_NEUTRON,
	     o2scl_cgsm::mass_neutron,4.0e-7,
	     "CGSM mass_neutron");
  t.test_rel(GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS,
	     o2scl_cgsm::unified_atomic_mass,4.0e-7,
	     "CGSM unified_atomic_mass");
  t.test_rel(GSL_CONST_CGSM_RYDBERG,
	     o2scl_cgsm::rydberg,4.0e-7,
	     "CGSM rydberg");
  t.test_rel(GSL_CONST_CGSM_MOLAR_GAS,
	     o2scl_cgsm::molar_gas,4.0e-6,
	     "CGSM molar_gas");
  t.test_rel(GSL_CONST_CGSM_STANDARD_GAS_VOLUME,
	     o2scl_cgsm::standard_gas_volume,4.0e-6,
	     "CGSM standard_gas_volume");
  t.test_rel(GSL_CONST_CGSM_MINUTE,
	     o2scl_cgsm::minute,1.0e-7,
	     "CGSM minute");
  t.test_rel(GSL_CONST_CGSM_HOUR,
	     o2scl_cgsm::hour,1.0e-7,
	     "CGSM hour");
  t.test_rel(GSL_CONST_CGSM_DAY,
	     o2scl_cgsm::day,1.0e-7,
	     "CGSM day");
  t.test_rel(GSL_CONST_CGSM_WEEK,
	     o2scl_cgsm::week,1.0e-7,
	     "CGSM week");
  t.test_rel(GSL_CONST_CGSM_INCH,
	     o2scl_cgsm::inch,1.0e-7,
	     "CGSM inch");
  t.test_rel(GSL_CONST_CGSM_FOOT,
	     o2scl_cgsm::foot,1.0e-7,
	     "CGSM foot");
  t.test_rel(GSL_CONST_CGSM_YARD,
	     o2scl_cgsm::yard,1.0e-7,
	     "CGSM yard");
  t.test_rel(GSL_CONST_CGSM_MILE,
	     o2scl_cgsm::mile,1.0e-7,
	     "CGSM mile");
  t.test_rel(GSL_CONST_CGSM_NAUTICAL_MILE,
	     o2scl_cgsm::nautical_mile,1.0e-7,
	     "CGSM nautical_mile");
  t.test_rel(GSL_CONST_CGSM_FATHOM,
	     o2scl_cgsm::fathom,1.0e-7,
	     "CGSM fathom");
  t.test_rel(GSL_CONST_CGSM_MIL,
	     o2scl_cgsm::mil,1.0e-7,
	     "CGSM mil");
  t.test_rel(GSL_CONST_CGSM_POINT,
	     o2scl_cgsm::point,1.0e-7,
	     "CGSM point");
  t.test_rel(GSL_CONST_CGSM_TEXPOINT,
	     o2scl_cgsm::texpoint,1.0e-7,
	     "CGSM texpoint");
  t.test_rel(GSL_CONST_CGSM_MICRON,
	     o2scl_cgsm::micron,1.0e-7,
	     "CGSM micron");
  t.test_rel(GSL_CONST_CGSM_ANGSTROM,
	     o2scl_cgsm::angstrom,1.0e-7,
	     "CGSM angstrom");
  t.test_rel(GSL_CONST_CGSM_HECTARE,
	     o2scl_cgsm::hectare,1.0e-7,
	     "CGSM hectare");
  t.test_rel(GSL_CONST_CGSM_ACRE,
	     o2scl_cgsm::acre,1.0e-7,
	     "CGSM acre");
  t.test_rel(GSL_CONST_CGSM_BARN,
	     o2scl_cgsm::barn,1.0e-7,
	     "CGSM barn");
  t.test_rel(GSL_CONST_CGSM_LITER,
	     o2scl_cgsm::liter,1.0e-7,
	     "CGSM liter");
  t.test_rel(GSL_CONST_CGSM_US_GALLON,
	     o2scl_cgsm::us_gallon,1.0e-7,
	     "CGSM us_gallon");
  t.test_rel(GSL_CONST_CGSM_QUART,
	     o2scl_cgsm::quart,1.0e-7,
	     "CGSM quart");
  t.test_rel(GSL_CONST_CGSM_PINT,
	     o2scl_cgsm::pint,1.0e-7,
	     "CGSM pint");
  t.test_rel(GSL_CONST_CGSM_CUP,
	     o2scl_cgsm::cup,1.0e-7,
	     "CGSM cup");
  t.test_rel(GSL_CONST_CGSM_FLUID_OUNCE,
	     o2scl_cgsm::fluid_ounce,1.0e-7,
	     "CGSM fluid_ounce");
  t.test_rel(GSL_CONST_CGSM_TABLESPOON,
	     o2scl_cgsm::tablespoon,1.0e-7,
	     "CGSM tablespoon");
  t.test_rel(GSL_CONST_CGSM_TEASPOON,
	     o2scl_cgsm::teaspoon,1.0e-7,
	     "CGSM teaspoon");
  t.test_rel(GSL_CONST_CGSM_CANADIAN_GALLON,
	     o2scl_cgsm::canadian_gallon,1.0e-7,
	     "CGSM canadian_gallon");
  t.test_rel(GSL_CONST_CGSM_UK_GALLON,
	     o2scl_cgsm::uk_gallon,1.0e-7,
	     "CGSM uk_gallon");
  t.test_rel(GSL_CONST_CGSM_MILES_PER_HOUR,
	     o2scl_cgsm::miles_per_hour,1.0e-7,
	     "CGSM miles_per_hour");
  t.test_rel(GSL_CONST_CGSM_KILOMETERS_PER_HOUR,
	     o2scl_cgsm::kilometers_per_hour,1.0e-7,
	     "CGSM kilometers_per_hour");
  t.test_rel(GSL_CONST_CGSM_KNOT,
	     o2scl_cgsm::knot,1.0e-7,
	     "CGSM knot");
  t.test_rel(GSL_CONST_CGSM_POUND_MASS,
	     o2scl_cgsm::pound_mass,1.0e-7,
	     "CGSM pound_mass");
  t.test_rel(GSL_CONST_CGSM_OUNCE_MASS,
	     o2scl_cgsm::ounce_mass,1.0e-7,
	     "CGSM ounce_mass");
  t.test_rel(GSL_CONST_CGSM_TON,
	     o2scl_cgsm::ton,1.0e-7,
	     "CGSM ton");
  t.test_rel(GSL_CONST_CGSM_METRIC_TON,
	     o2scl_cgsm::metric_ton,1.0e-7,
	     "CGSM metric_ton");
  t.test_rel(GSL_CONST_CGSM_UK_TON,
	     o2scl_cgsm::uk_ton,1.0e-7,
	     "CGSM uk_ton");
  t.test_rel(GSL_CONST_CGSM_TROY_OUNCE,
	     o2scl_cgsm::troy_ounce,1.0e-7,
	     "CGSM troy_ounce");
  t.test_rel(GSL_CONST_CGSM_CARAT,
	     o2scl_cgsm::carat,1.0e-7,
	     "CGSM carat");
  t.test_rel(GSL_CONST_CGSM_GRAM_FORCE,
	     o2scl_cgsm::gram_force,1.0e-7,
	     "CGSM gram_force");
  t.test_rel(GSL_CONST_CGSM_POUND_FORCE,
	     o2scl_cgsm::pound_force,1.0e-7,
	     "CGSM pound_force");
  t.test_rel(GSL_CONST_CGSM_KILOPOUND_FORCE,
	     o2scl_cgsm::kilopound_force,1.0e-7,
	     "CGSM kilopound_force");
  t.test_rel(GSL_CONST_CGSM_POUNDAL,
	     o2scl_cgsm::poundal,1.0e-7,
	     "CGSM poundal");
  t.test_rel(GSL_CONST_CGSM_CALORIE,
	     o2scl_cgsm::calorie,1.0e-7,
	     "CGSM calorie");
  t.test_rel(GSL_CONST_CGSM_BTU,
	     o2scl_cgsm::btu,1.0e-7,
	     "CGSM btu");
  t.test_rel(GSL_CONST_CGSM_THERM,
	     o2scl_cgsm::therm,1.0e-7,
	     "CGSM therm");
  t.test_rel(GSL_CONST_CGSM_HORSEPOWER,
	     o2scl_cgsm::horsepower,1.0e-7,
	     "CGSM horsepower");
  t.test_rel(GSL_CONST_CGSM_BAR,
	     o2scl_cgsm::bar,1.0e-7,
	     "CGSM bar");
  t.test_rel(GSL_CONST_CGSM_STD_ATMOSPHERE,
	     o2scl_cgsm::std_atmosphere,1.0e-7,
	     "CGSM std_atmosphere");
  t.test_rel(GSL_CONST_CGSM_TORR,
	     o2scl_cgsm::torr,1.0e-7,
	     "CGSM torr");
  t.test_rel(GSL_CONST_CGSM_METER_OF_MERCURY,
	     o2scl_cgsm::meter_of_mercury,1.0e-7,
	     "CGSM meter_of_mercury");
  t.test_rel(GSL_CONST_CGSM_INCH_OF_MERCURY,
	     o2scl_cgsm::inch_of_mercury,1.0e-7,
	     "CGSM inch_of_mercury");
  t.test_rel(GSL_CONST_CGSM_INCH_OF_WATER,
	     o2scl_cgsm::inch_of_water,1.0e-7,
	     "CGSM inch_of_water");
  t.test_rel(GSL_CONST_CGSM_PSI,
	     o2scl_cgsm::psi,1.0e-7,
	     "CGSM psi");
  t.test_rel(GSL_CONST_CGSM_POISE,
	     o2scl_cgsm::poise,1.0e-7,
	     "CGSM poise");
  t.test_rel(GSL_CONST_CGSM_STOKES,
	     o2scl_cgsm::stokes,1.0e-7,
	     "CGSM stokes");
  t.test_rel(GSL_CONST_CGSM_STILB,
	     o2scl_cgsm::stilb,1.0e-7,
	     "CGSM stilb");
  t.test_rel(GSL_CONST_CGSM_LUMEN,
	     o2scl_cgsm::lumen,1.0e-7,
	     "CGSM lumen");
  t.test_rel(GSL_CONST_CGSM_LUX,
	     o2scl_cgsm::lux,1.0e-7,
	     "CGSM lux");
  t.test_rel(GSL_CONST_CGSM_PHOT,
	     o2scl_cgsm::phot,1.0e-7,
	     "CGSM phot");
  t.test_rel(GSL_CONST_CGSM_FOOTCANDLE,
	     o2scl_cgsm::footcandle,1.0e-7,
	     "CGSM footcandle");
  t.test_rel(GSL_CONST_CGSM_LAMBERT,
	     o2scl_cgsm::lambert,1.0e-7,
	     "CGSM lambert");
  t.test_rel(GSL_CONST_CGSM_FOOTLAMBERT,
	     o2scl_cgsm::footlambert,1.0e-7,
	     "CGSM footlambert");
  t.test_rel(GSL_CONST_CGSM_CURIE,
	     o2scl_cgsm::curie,1.0e-7,
	     "CGSM curie");
  t.test_rel(GSL_CONST_CGSM_RAD,
	     o2scl_cgsm::rad,1.0e-7,
	     "CGSM rad");
  t.test_rel(GSL_CONST_CGSM_NEWTON,
	     o2scl_cgsm::newton,1.0e-7,
	     "CGSM newton");
  t.test_rel(GSL_CONST_CGSM_DYNE,
	     o2scl_cgsm::dyne,1.0e-7,
	     "CGSM dyne");
  t.test_rel(GSL_CONST_CGSM_JOULE,
	     o2scl_cgsm::joule,1.0e-7,
	     "CGSM joule");
  t.test_rel(GSL_CONST_CGSM_ERG,
	     o2scl_cgsm::erg,1.0e-7,
	     "CGSM erg");
  t.test_rel(GSL_CONST_CGSM_ELECTRON_MAGNETIC_MOMENT,
	     o2scl_cgsm::electron_magnetic_moment,4.0e-7,
	     "CGSM electron_magnetic_moment");
  t.test_rel(GSL_CONST_CGSM_PROTON_MAGNETIC_MOMENT,
	     o2scl_cgsm::proton_magnetic_moment,4.0e-7,
	     "CGSM proton_magnetic_moment");
  t.test_rel(GSL_CONST_CGSM_ROENTGEN,
	     o2scl_cgsm::roentgen,1.0e-7,
	     "CGSM roentgen");
  t.test_rel(GSL_CONST_CGSM_BOHR_MAGNETON,
	     o2scl_cgsm::bohr_magneton,4.0e-7,
	     "CGSM bohr_magneton");
  t.test_rel(GSL_CONST_CGSM_NUCLEAR_MAGNETON,
	     o2scl_cgsm::nuclear_magneton,4.0e-7,
	     "CGSM nuclear_magneton");
  t.test_rel(GSL_CONST_CGSM_FARADAY,
	     o2scl_cgsm::faraday,3.0e-7,
	     "CGSM faraday");
  t.test_rel(GSL_CONST_CGSM_ELECTRON_CHARGE,
	     o2scl_cgsm::electron_charge,1.0e-7,
	     "CGSM electron_charge");

  t.test_rel(GSL_CONST_MKS_SPEED_OF_LIGHT,
	     o2scl_mks::speed_of_light,1.0e-7,
	     "MKS speed_of_light");
  t.test_rel(GSL_CONST_MKS_GRAVITATIONAL_CONSTANT,
	     o2scl_mks::gravitational_constant,3.0e-4,
	     "MKS gravitational_constant");
  t.test_rel(GSL_CONST_MKS_PLANCKS_CONSTANT_H,
	     o2scl_mks::plancks_constant_h,3.0e-7,
	     "MKS plancks_constant_h");
  t.test_rel(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR,
	     o2scl_mks::plancks_constant_hbar,3.0e-7,
	     "MKS plancks_constant_hbar");
  t.test_rel(GSL_CONST_MKS_ELECTRON_VOLT,
	     o2scl_mks::electron_volt,1.0e-7,
	     "MKS electron_volt");
  t.test_rel(GSL_CONST_MKS_BOHR_RADIUS,
	     o2scl_mks::bohr_radius,1.0e-7,
	     "MKS bohr_radius");
  t.test_rel(GSL_CONST_MKS_THOMSON_CROSS_SECTION,
	     o2scl_mks::thomson_cross_section,4.0e-7,
	     "MKS thomson_cross_section");
  t.test_rel(GSL_CONST_MKS_BOLTZMANN,
	     o2scl_mks::boltzmann,4.0e-6,
	     "MKS boltzmann");
  t.test_rel(GSL_CONST_MKS_STEFAN_BOLTZMANN_CONSTANT,
	     o2scl_mks::stefan_boltzmann_constant,1.0e-5,
	     "MKS stefan_boltzmann");
  t.test_rel(GSL_CONST_MKS_ASTRONOMICAL_UNIT,
	     o2scl_mks::astronomical_unit,1.0e-7,
	     "MKS astronomical_unit");
  t.test_rel(GSL_CONST_MKS_PARSEC,
	     o2scl_mks::parsec,1.0e-7,
	     "MKS parsec");
  t.test_rel(GSL_CONST_MKS_GRAV_ACCEL,
	     o2scl_mks::grav_accel,1.0e-7,
	     "MKS grav_accel");
  t.test_rel(GSL_CONST_MKS_SOLAR_MASS,
	     o2scl_mks::solar_mass,5.0e-4,
	     "MKS solar_mass");
  t.test_rel(GSL_CONST_MKS_LIGHT_YEAR,
	     o2scl_mks::light_year,4.0e-5,
	     "MKS light_year");
  t.test_rel(GSL_CONST_MKS_MASS_ELECTRON,
	     o2scl_mks::mass_electron,4.0e-7,
	     "MKS mass_electron");
  t.test_rel(GSL_CONST_MKS_MASS_MUON,
	     o2scl_mks::mass_muon,4.0e-7,
	     "MKS mass_muon");
  t.test_rel(GSL_CONST_MKS_MASS_PROTON,
	     o2scl_mks::mass_proton,4.0e-7,
	     "MKS mass_proton");
  t.test_rel(GSL_CONST_MKS_MASS_NEUTRON,
	     o2scl_mks::mass_neutron,4.0e-7,
	     "MKS mass_neutron");
  t.test_rel(GSL_CONST_MKS_UNIFIED_ATOMIC_MASS,
	     o2scl_mks::unified_atomic_mass,4.0e-7,
	     "MKS unified_atomic_mass");
  t.test_rel(GSL_CONST_MKS_RYDBERG,
	     o2scl_mks::rydberg,4.0e-7,
	     "MKS rydberg");
  t.test_rel(GSL_CONST_MKS_MOLAR_GAS,
	     o2scl_mks::molar_gas,4.0e-6,
	     "MKS molar_gas");
  t.test_rel(GSL_CONST_MKS_STANDARD_GAS_VOLUME,
	     o2scl_mks::standard_gas_volume,4.0e-6,
	     "MKS standard_gas_volume");
  t.test_rel(GSL_CONST_MKS_MINUTE,
	     o2scl_mks::minute,1.0e-7,
	     "MKS minute");
  t.test_rel(GSL_CONST_MKS_HOUR,
	     o2scl_mks::hour,1.0e-7,
	     "MKS hour");
  t.test_rel(GSL_CONST_MKS_DAY,
	     o2scl_mks::day,1.0e-7,
	     "MKS day");
  t.test_rel(GSL_CONST_MKS_WEEK,
	     o2scl_mks::week,1.0e-7,
	     "MKS week");
  t.test_rel(GSL_CONST_MKS_INCH,
	     o2scl_mks::inch,1.0e-7,
	     "MKS inch");
  t.test_rel(GSL_CONST_MKS_FOOT,
	     o2scl_mks::foot,1.0e-7,
	     "MKS foot");
  t.test_rel(GSL_CONST_MKS_YARD,
	     o2scl_mks::yard,1.0e-7,
	     "MKS yard");
  t.test_rel(GSL_CONST_MKS_MILE,
	     o2scl_mks::mile,1.0e-7,
	     "MKS mile");
  t.test_rel(GSL_CONST_MKS_NAUTICAL_MILE,
	     o2scl_mks::nautical_mile,1.0e-7,
	     "MKS nautical_mile");
  t.test_rel(GSL_CONST_MKS_FATHOM,
	     o2scl_mks::fathom,1.0e-7,
	     "MKS fathom");
  t.test_rel(GSL_CONST_MKS_MIL,
	     o2scl_mks::mil,1.0e-7,
	     "MKS mil");
  t.test_rel(GSL_CONST_MKS_POINT,
	     o2scl_mks::point,1.0e-7,
	     "MKS point");
  t.test_rel(GSL_CONST_MKS_TEXPOINT,
	     o2scl_mks::texpoint,1.0e-7,
	     "MKS texpoint");
  t.test_rel(GSL_CONST_MKS_MICRON,
	     o2scl_mks::micron,1.0e-7,
	     "MKS micron");
  t.test_rel(GSL_CONST_MKS_ANGSTROM,
	     o2scl_mks::angstrom,1.0e-7,
	     "MKS angstrom");
  t.test_rel(GSL_CONST_MKS_HECTARE,
	     o2scl_mks::hectare,1.0e-7,
	     "MKS hectare");
  t.test_rel(GSL_CONST_MKS_ACRE,
	     o2scl_mks::acre,1.0e-7,
	     "MKS acre");
  t.test_rel(GSL_CONST_MKS_BARN,
	     o2scl_mks::barn,1.0e-7,
	     "MKS barn");
  t.test_rel(GSL_CONST_MKS_LITER,
	     o2scl_mks::liter,1.0e-7,
	     "MKS liter");
  t.test_rel(GSL_CONST_MKS_US_GALLON,
	     o2scl_mks::us_gallon,1.0e-7,
	     "MKS us_gallon");
  t.test_rel(GSL_CONST_MKS_QUART,
	     o2scl_mks::quart,1.0e-7,
	     "MKS quart");
  t.test_rel(GSL_CONST_MKS_PINT,
	     o2scl_mks::pint,1.0e-7,
	     "MKS pint");
  t.test_rel(GSL_CONST_MKS_CUP,
	     o2scl_mks::cup,1.0e-7,
	     "MKS cup");
  t.test_rel(GSL_CONST_MKS_FLUID_OUNCE,
	     o2scl_mks::fluid_ounce,1.0e-7,
	     "MKS fluid_ounce");
  t.test_rel(GSL_CONST_MKS_TABLESPOON,
	     o2scl_mks::tablespoon,1.0e-7,
	     "MKS tablespoon");
  t.test_rel(GSL_CONST_MKS_TEASPOON,
	     o2scl_mks::teaspoon,1.0e-7,
	     "MKS teaspoon");
  t.test_rel(GSL_CONST_MKS_CANADIAN_GALLON,
	     o2scl_mks::canadian_gallon,1.0e-7,
	     "MKS canadian_gallon");
  t.test_rel(GSL_CONST_MKS_UK_GALLON,
	     o2scl_mks::uk_gallon,1.0e-7,
	     "MKS uk_gallon");
  t.test_rel(GSL_CONST_MKS_MILES_PER_HOUR,
	     o2scl_mks::miles_per_hour,1.0e-7,
	     "MKS miles_per_hour");
  t.test_rel(GSL_CONST_MKS_KILOMETERS_PER_HOUR,
	     o2scl_mks::kilometers_per_hour,1.0e-7,
	     "MKS kilometers_per_hour");
  t.test_rel(GSL_CONST_MKS_KNOT,
	     o2scl_mks::knot,1.0e-7,
	     "MKS knot");
  t.test_rel(GSL_CONST_MKS_POUND_MASS,
	     o2scl_mks::pound_mass,1.0e-7,
	     "MKS pound_mass");
  t.test_rel(GSL_CONST_MKS_OUNCE_MASS,
	     o2scl_mks::ounce_mass,1.0e-7,
	     "MKS ounce_mass");
  t.test_rel(GSL_CONST_MKS_TON,
	     o2scl_mks::ton,1.0e-7,
	     "MKS ton");
  t.test_rel(GSL_CONST_MKS_METRIC_TON,
	     o2scl_mks::metric_ton,1.0e-7,
	     "MKS metric_ton");
  t.test_rel(GSL_CONST_MKS_UK_TON,
	     o2scl_mks::uk_ton,1.0e-7,
	     "MKS uk_ton");
  t.test_rel(GSL_CONST_MKS_TROY_OUNCE,
	     o2scl_mks::troy_ounce,1.0e-7,
	     "MKS troy_ounce");
  t.test_rel(GSL_CONST_MKS_CARAT,
	     o2scl_mks::carat,1.0e-7,
	     "MKS carat");
  t.test_rel(GSL_CONST_MKS_GRAM_FORCE,
	     o2scl_mks::gram_force,1.0e-7,
	     "MKS gram_force");
  t.test_rel(GSL_CONST_MKS_POUND_FORCE,
	     o2scl_mks::pound_force,1.0e-7,
	     "MKS pound_force");
  t.test_rel(GSL_CONST_MKS_KILOPOUND_FORCE,
	     o2scl_mks::kilopound_force,1.0e-7,
	     "MKS kilopound_force");
  t.test_rel(GSL_CONST_MKS_POUNDAL,
	     o2scl_mks::poundal,1.0e-7,
	     "MKS poundal");
  t.test_rel(GSL_CONST_MKS_CALORIE,
	     o2scl_mks::calorie,1.0e-7,
	     "MKS calorie");
  t.test_rel(GSL_CONST_MKS_BTU,
	     o2scl_mks::btu,1.0e-7,
	     "MKS btu");
  t.test_rel(GSL_CONST_MKS_THERM,
	     o2scl_mks::therm,1.0e-7,
	     "MKS therm");
  t.test_rel(GSL_CONST_MKS_HORSEPOWER,
	     o2scl_mks::horsepower,1.0e-7,
	     "MKS horsepower");
  t.test_rel(GSL_CONST_MKS_BAR,
	     o2scl_mks::bar,1.0e-7,
	     "MKS bar");
  t.test_rel(GSL_CONST_MKS_STD_ATMOSPHERE,
	     o2scl_mks::std_atmosphere,1.0e-7,
	     "MKS std_atmosphere");
  t.test_rel(GSL_CONST_MKS_TORR,
	     o2scl_mks::torr,1.0e-7,
	     "MKS torr");
  t.test_rel(GSL_CONST_MKS_METER_OF_MERCURY,
	     o2scl_mks::meter_of_mercury,1.0e-7,
	     "MKS meter_of_mercury");
  t.test_rel(GSL_CONST_MKS_INCH_OF_MERCURY,
	     o2scl_mks::inch_of_mercury,1.0e-7,
	     "MKS inch_of_mercury");
  t.test_rel(GSL_CONST_MKS_INCH_OF_WATER,
	     o2scl_mks::inch_of_water,1.0e-7,
	     "MKS inch_of_water");
  t.test_rel(GSL_CONST_MKS_PSI,
	     o2scl_mks::psi,1.0e-7,
	     "MKS psi");
  t.test_rel(GSL_CONST_MKS_POISE,
	     o2scl_mks::poise,1.0e-7,
	     "MKS poise");
  t.test_rel(GSL_CONST_MKS_STOKES,
	     o2scl_mks::stokes,1.0e-7,
	     "MKS stokes");
  t.test_rel(GSL_CONST_MKS_GAUSS,
	     o2scl_mks::gauss,1.0e-7,
	     "MKS gauss");
  t.test_rel(GSL_CONST_MKS_STILB,
	     o2scl_mks::stilb,1.0e-7,
	     "MKS stilb");
  t.test_rel(GSL_CONST_MKS_LUMEN,
	     o2scl_mks::lumen,1.0e-7,
	     "MKS lumen");
  t.test_rel(GSL_CONST_MKS_LUX,
	     o2scl_mks::lux,1.0e-7,
	     "MKS lux");
  t.test_rel(GSL_CONST_MKS_PHOT,
	     o2scl_mks::phot,1.0e-7,
	     "MKS phot");
  t.test_rel(GSL_CONST_MKS_FOOTCANDLE,
	     o2scl_mks::footcandle,1.0e-7,
	     "MKS footcandle");
  t.test_rel(GSL_CONST_MKS_LAMBERT,
	     o2scl_mks::lambert,1.0e-7,
	     "MKS lambert");
  t.test_rel(GSL_CONST_MKS_FOOTLAMBERT,
	     o2scl_mks::footlambert,1.0e-7,
	     "MKS footlambert");
  t.test_rel(GSL_CONST_MKS_CURIE,
	     o2scl_mks::curie,1.0e-7,
	     "MKS curie");
  t.test_rel(GSL_CONST_MKS_RAD,
	     o2scl_mks::rad,1.0e-7,
	     "MKS rad");
  t.test_rel(GSL_CONST_MKS_NEWTON,
	     o2scl_mks::newton,1.0e-7,
	     "MKS newton");
  t.test_rel(GSL_CONST_MKS_DYNE,
	     o2scl_mks::dyne,1.0e-7,
	     "MKS dyne");
  t.test_rel(GSL_CONST_MKS_JOULE,
	     o2scl_mks::joule,1.0e-7,
	     "MKS joule");
  t.test_rel(GSL_CONST_MKS_ERG,
	     o2scl_mks::erg,1.0e-7,
	     "MKS erg");
  t.test_rel(GSL_CONST_MKS_ELECTRON_MAGNETIC_MOMENT,
	     o2scl_mks::electron_magnetic_moment,4.0e-7,
	     "MKS electron_magnetic_moment");
  t.test_rel(GSL_CONST_MKS_PROTON_MAGNETIC_MOMENT,
	     o2scl_mks::proton_magnetic_moment,4.0e-7,
	     "MKS proton_magnetic_moment");
  t.test_rel(GSL_CONST_MKS_ROENTGEN,
	     o2scl_mks::roentgen,1.0e-7,
	     "MKS roentgen");
  t.test_rel(GSL_CONST_MKS_BOHR_MAGNETON,
	     o2scl_mks::bohr_magneton,4.0e-7,
	     "MKS bohr_magneton");
  t.test_rel(GSL_CONST_MKS_NUCLEAR_MAGNETON,
	     o2scl_mks::nuclear_magneton,4.0e-7,
	     "MKS nuclear_magneton");
  t.test_rel(GSL_CONST_MKS_VACUUM_PERMITTIVITY,
	     o2scl_mks::vacuum_permittivity,1.0e-7,
	     "MKS vacuum_permittivity");
  t.test_rel(GSL_CONST_MKS_VACUUM_PERMEABILITY,
	     o2scl_mks::vacuum_permeability,1.0e-7,
	     "MKS vacuum_permeability");
  t.test_rel(GSL_CONST_MKS_FARADAY,
	     o2scl_mks::faraday,3.0e-7,
	     "MKS faraday");
  t.test_rel(GSL_CONST_MKS_ELECTRON_CHARGE,
	     o2scl_mks::electron_charge,1.0e-7,
	     "MKS electron_charge");

  t.report();

  return 0;
}
