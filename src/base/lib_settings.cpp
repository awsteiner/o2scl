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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/lib_settings.h>
#include <o2scl/prev_commit.h>

#ifdef O2SCL_HDF
#ifdef O2SCL_PLAIN_HDF5_HEADER
#include <hdf5.h>
#else
#ifdef O2SCL_LINUX
#include <hdf5/serial/hdf5.h>
#else
#include <hdf5.h>
#endif
#endif
#endif

using namespace std;
using namespace o2scl;

lib_settings_class o2scl::o2scl_settings;

lib_settings_class::lib_settings_class() {
  data_dir=((std::string)(O2SCL_DATA_DIR));
  cup=&def_cu;
  fcp=&def_fc;
  
  // Default conversions are given here. Obviously GNU units is better
  // at handling these things, but it's nice to have some of the easy
  // conversions in by default rather than worrying about opening a
  // pipe, etc.
  def_cu.default_conversions();

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
#ifdef O2SCL_UBUNTU_PKG
  return "<package>";
#else
  return __DATE__;
#endif
}

std::string lib_settings_class::time_compiled() {
#ifdef O2SCL_UBUNTU_PKG
  return "<package>";
#else
  return __TIME__;
#endif
}

bool lib_settings_class::eos_installed() {
#ifdef O2SCL_EOS
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::hdf5_compression_support() {
#ifdef O2SCL_HDF5_COMP
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

void lib_settings_class::hdf5_header_version(unsigned &maj,
					     unsigned &min, unsigned &rel) {
#ifdef O2SCL_HDF
  maj=H5_VERS_MAJOR;
  min=H5_VERS_MINOR;
  rel=H5_VERS_RELEASE;
#else
  maj=0;
  min=0;
  rel=0;
#endif
  return;
}

void lib_settings_class::hdf5_lib_version(unsigned &maj,
					  unsigned &min, unsigned &rel) {
#ifdef O2SCL_HDF
  H5get_libversion(&maj,&min,&rel);
#else
  maj=0;
  min=0;
  rel=0;
#endif
  return;
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

bool lib_settings_class::fftw_support() {
#ifdef O2SCL_FFTW
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::python_support() {
#ifdef O2SCL_PYTHON
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

bool lib_settings_class::readline_support() {
#ifdef O2SCL_READLINE
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::ncurses_support() {
#ifdef O2SCL_NCURSES
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::gsl2_support() {
#ifdef O2SCL_GSL2
  return true;
#else
  return false;
#endif
}

std::string lib_settings_class::system_type() {
#ifdef O2SCL_LINUX
  return ((string)"Linux");
#else
#ifdef O2SCL_OSX
  return ((string)"OSX");
#else
  return ((string)"unknown");
#endif
#endif
}

void lib_settings_class::config_h_report() {
  cout << "Previous commit hash: " << O2SCL_PREV_COMMIT_HASH << endl;
  cout << "Previous commit date: " << O2SCL_PREV_COMMIT_DATE << endl;
  cout << "Commit branch: " << O2SCL_BRANCH << endl;
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
#ifdef HAVE_POPEN
  cout << "HAVE_POPEN: " << HAVE_POPEN << endl;
#else
  cout << "HAVE_POPEN: <not defined>" << endl;
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
