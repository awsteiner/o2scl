/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

lib_settings_class o2scl::o2scl_settings;

lib_settings_class::lib_settings_class() {
  data_dir=((std::string)(O2SCL_DATA_DIR));
  cup=&def_cu;

  // Default conversions 

  def_cu.insert_cache("kg","1/fm",1.0e-15/o2scl_mks::plancks_constant_hbar*
		      o2scl_mks::speed_of_light);
  def_cu.insert_cache("kg","MeV",1.0e-15/o2scl_mks::plancks_constant_hbar*
		      o2scl_mks::speed_of_light*o2scl_const::hc_mev_fm);

  // For the TOV solver
  def_cu.insert_cache("g/cm^3","Msun/km^3",1.0e12/o2scl_mks::solar_mass);
  def_cu.insert_cache("erg/cm^3","Msun/km^3",1.0e12/o2scl_cgs::speed_of_light/
		      o2scl_cgs::speed_of_light/o2scl_mks::solar_mass);
  def_cu.insert_cache("dyne/cm^2","Msun/km^3",1.0e12/o2scl_cgs::speed_of_light/
		      o2scl_cgs::speed_of_light/o2scl_mks::solar_mass);
  def_cu.insert_cache("MeV/fm^3","Msun/km^3",
		      o2scl_cgs::electron_volt/o2scl_cgs::speed_of_light/
		      o2scl_cgs::speed_of_light/o2scl_mks::solar_mass*1.0e57);
  def_cu.insert_cache("1/fm^4","Msun/km^3",o2scl_const::hc_mev_fm*
		      o2scl_cgs::electron_volt/o2scl_cgs::speed_of_light/
		      o2scl_cgs::speed_of_light/o2scl_mks::solar_mass*1.0e57);

  def_cu.insert_cache("1/fm^3","1/cm^3",1.0e39);
  def_cu.insert_cache("1/fm^3","1/m^3",1.0e45);

}

lib_settings_class::~lib_settings_class() {
}

bool lib_settings_class::range_check() {
#if !O2SCL_NO_RANGE_CHECK
  return true;
#else
  return false;
#endif
}

std::string lib_settings_class::o2scl_version() {
  return PACKAGE_VERSION;
}

std::string lib_settings_class::o2scl_name() {
  return PACKAGE_NAME;
}

std::string lib_settings_class::o2scl_package() {
  return PACKAGE;
}

std::string lib_settings_class::o2scl_bugreport() {
  return PACKAGE_BUGREPORT;
}

std::string lib_settings_class::o2scl_string() {
  return PACKAGE_STRING;
}

std::string lib_settings_class::o2scl_tarname() {
  return PACKAGE_TARNAME;
}

std::string lib_settings_class::date_compiled() {
  return __DATE__;
}

std::string lib_settings_class::time_compiled() {
  return __TIME__;
}

bool lib_settings_class::eos_installed() {
#ifdef O2SCL_EOS
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::part_installed() {
#ifdef O2SCL_PART
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::hdf_support() {
#ifdef O2SCL_HDF
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::armadillo_support() {
#ifdef O2SCL_ARMA
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::eigen_support() {
#ifdef O2SCL_EIGEN
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::openmp_support() {
#ifdef O2SCL_OPENMP
  return true;
#else
  return false;
#endif
}

void lib_settings_class::config_h_report() {
#ifdef HAVE_ACOSH
  cout << "HAVE_ACOSH: " << HAVE_ACOSH << endl;
#else
  cout << "HAVE_ACOSH: <not defined>" << endl;
#endif
#ifdef HAVE_ASINH
  cout << "HAVE_ASINH: " << HAVE_ASINH << endl;
#else
  cout << "HAVE_ASINH: <not defined>" << endl;
#endif
#ifdef HAVE_ATANH
  cout << "HAVE_ATANH: " << HAVE_ATANH << endl;
#else
  cout << "HAVE_ATANH: <not defined>" << endl;
#endif
#ifdef HAVE_DLFCN_H
  cout << "HAVE_DLFCN_H: " << HAVE_DLFCN_H << endl;
#else
  cout << "HAVE_DLFCN_H: <not defined>" << endl;
#endif
#ifdef HAVE_FINITE
  cout << "HAVE_FINITE: " << HAVE_FINITE << endl;
#else
  cout << "HAVE_FINITE: <not defined>" << endl;
#endif
#ifdef HAVE_FLOOR
  cout << "HAVE_FLOOR: " << HAVE_FLOOR << endl;
#else
  cout << "HAVE_FLOOR: <not defined>" << endl;
#endif
#ifdef HAVE_INTTYPES_H
  cout << "HAVE_INTTYPES_H: " << HAVE_INTTYPES_H << endl;
#else
  cout << "HAVE_INTTYPES_H: <not defined>" << endl;
#endif
#ifdef HAVE_LIBCBLAS
  cout << "HAVE_LIBCBLAS: " << HAVE_LIBCBLAS << endl;
#else
  cout << "HAVE_LIBCBLAS: <not defined>" << endl;
#endif
#ifdef HAVE_LIBGSL
  cout << "HAVE_LIBGSL: " << HAVE_LIBGSL << endl;
#else
  cout << "HAVE_LIBGSL: <not defined>" << endl;
#endif
#ifdef HAVE_LIBGSLCBLAS
  cout << "HAVE_LIBGSLCBLAS: " << HAVE_LIBGSLCBLAS << endl;
#else
  cout << "HAVE_LIBGSLCBLAS: <not defined>" << endl;
#endif
#ifdef HAVE_LIBHDF5
  cout << "HAVE_LIBHDF5: " << HAVE_LIBHDF5 << endl;
#else
  cout << "HAVE_LIBHDF5: <not defined>" << endl;
#endif
#ifdef HAVE_LIBHDF5_HL
  cout << "HAVE_LIBHDF5_HL: " << HAVE_LIBHDF5_HL << endl;
#else
  cout << "HAVE_LIBHDF5_HL: <not defined>" << endl;
#endif
#ifdef HAVE_MALLOC
  cout << "HAVE_MALLOC: " << HAVE_MALLOC << endl;
#else
  cout << "HAVE_MALLOC: <not defined>" << endl;
#endif
#ifdef HAVE_MALLOC_H
  cout << "HAVE_MALLOC_H: " << HAVE_MALLOC_H << endl;
#else
  cout << "HAVE_MALLOC_H: <not defined>" << endl;
#endif
#ifdef HAVE_MEMORY_H
  cout << "HAVE_MEMORY_H: " << HAVE_MEMORY_H << endl;
#else
  cout << "HAVE_MEMORY_H: <not defined>" << endl;
#endif
#ifdef HAVE_POW
  cout << "HAVE_POW: " << HAVE_POW << endl;
#else
  cout << "HAVE_POW: <not defined>" << endl;
#endif
#ifdef HAVE_SELECT
  cout << "HAVE_SELECT: " << HAVE_SELECT << endl;
#else
  cout << "HAVE_SELECT: <not defined>" << endl;
#endif
#ifdef HAVE_SQRT
  cout << "HAVE_SQRT: " << HAVE_SQRT << endl;
#else
  cout << "HAVE_SQRT: <not defined>" << endl;
#endif
#ifdef HAVE_STDBOOL_H
  cout << "HAVE_STDBOOL_H: " << HAVE_STDBOOL_H << endl;
#else
  cout << "HAVE_STDBOOL_H: <not defined>" << endl;
#endif
#ifdef HAVE_STDINT_H
  cout << "HAVE_STDINT_H: " << HAVE_STDINT_H << endl;
#else
  cout << "HAVE_STDINT_H: <not defined>" << endl;
#endif
#ifdef HAVE_STDLIB_H
  cout << "HAVE_STDLIB_H: " << HAVE_STDLIB_H << endl;
#else
  cout << "HAVE_STDLIB_H: <not defined>" << endl;
#endif
#ifdef HAVE_STRCHR
  cout << "HAVE_STRCHR: " << HAVE_STRCHR << endl;
#else
  cout << "HAVE_STRCHR: <not defined>" << endl;
#endif
#ifdef HAVE_STRINGS_H
  cout << "HAVE_STRINGS_H: " << HAVE_STRINGS_H << endl;
#else
  cout << "HAVE_STRINGS_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_SELECT_H
  cout << "HAVE_SYS_SELECT_H: " << HAVE_SYS_SELECT_H << endl;
#else
  cout << "HAVE_SYS_SELECT_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_SOCKET_H
  cout << "HAVE_SYS_SOCKET_H: " << HAVE_SYS_SOCKET_H << endl;
#else
  cout << "HAVE_SYS_SOCKET_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_STAT_H
  cout << "HAVE_SYS_STAT_H: " << HAVE_SYS_STAT_H << endl;
#else
  cout << "HAVE_SYS_STAT_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_TYPES_H
  cout << "HAVE_SYS_TYPES_H: " << HAVE_SYS_TYPES_H << endl;
#else
  cout << "HAVE_SYS_TYPES_H: <not defined>" << endl;
#endif
#ifdef HAVE_UNISTD_H
  cout << "HAVE_UNISTD_H: " << HAVE_UNISTD_H << endl;
#else
  cout << "HAVE_UNISTD_H: <not defined>" << endl;
#endif
#ifdef HAVE__BOOL
  cout << "HAVE__BOOL: " << HAVE__BOOL << endl;
#else
  cout << "HAVE__BOOL: <not defined>" << endl;
#endif
  return;
}
