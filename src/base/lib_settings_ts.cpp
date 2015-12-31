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
#include <iostream>

#include <o2scl/lib_settings.h>
#include <o2scl/test_mgr.h>
#include <o2scl/err_hnd.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  
  cout << "O2scl version: " << o2scl_settings.o2scl_version() << endl;
  cout << "Range checking: " << o2scl_settings.range_check() << endl;
  cout << "EOS library: " << o2scl_settings.eos_installed() << endl;
  cout << "Particle library: " << o2scl_settings.part_installed() << endl;
  cout << "HDF support: " << o2scl_settings.hdf_support() << endl;
  cout << "Armadillo support: " << o2scl_settings.armadillo_support() << endl;
  cout << "Eigen support: " << o2scl_settings.eigen_support() << endl;
  cout << "OpenMP support: " << o2scl_settings.openmp_support() << endl;
  cout << "Data directory: " << o2scl_settings.get_data_dir() << endl;
  cout << "Name: " << o2scl_settings.o2scl_name() << endl;
  cout << "Package: " << o2scl_settings.o2scl_package() << endl;
  cout << "Bug-report: " << o2scl_settings.o2scl_bugreport() << endl;
  cout << "String: " << o2scl_settings.o2scl_string() << endl;
  cout << "Tarname: " << o2scl_settings.o2scl_tarname() << endl;
  cout << "Compile date: " << o2scl_settings.date_compiled() << endl;
  cout << "Compile time: " << o2scl_settings.time_compiled() << endl;
  cout << endl;
  cout << "config.h: " << endl;
  o2scl_settings.config_h_report();

  t.report();
  return 0;
}

