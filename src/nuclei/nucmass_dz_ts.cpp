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
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  nucmass_ame ame;
  o2scl_hdf::ame_load(ame,"12");
  nucmass_ame_exp amex;
  o2scl_hdf::ame_load(amex,"12");

  nucmass_dz_fit dmf;
  nucmass_dz_fit_33 dmf33;
  nucmass_dz_table dz("96");
  nucmass_dz_table dz2("95");

  // Compare the fits with numerical values from the original Fortran
  // code

  t.test_rel(dmf33.binding_energy(6,7),-96.74419,5.0e-7,"be1");
  t.test_rel(dmf33.binding_energy(7,7),-106.2967,5.0e-7,"be1");
  t.test_rel(dmf33.binding_energy(7,8),-114.9327,5.0e-7,"be1");
  t.test_rel(dmf33.binding_energy(10,11),-167.2058,5.0e-7,"be2");
  t.test_rel(dmf33.binding_energy(10,10),-160.6454,5.0e-7,"be3");
  t.test_rel(dmf33.binding_energy(82,126),-1637.013,5.0e-7,"be4");
  t.test_rel(dmf33.binding_energy(87,53),-527.2083,5.0e-7,"be5");
  t.test_rel(dmf33.binding_energy(37,53),-776.3599,5.0e-7,"be6");

  t.test_rel(dmf.binding_energy(11,10),-162.6734,1.0e-6,"be1");
  t.test_rel(dmf.mass_excess(11,10),-1.780304,3.0e-5,"me1");
  t.test_rel(dmf.binding_energy(10,11),-167.2239,1.0e-6,"be2");
  t.test_rel(dmf.mass_excess(10,11),-5.548401,3.0e-5,"me2");
  t.test_rel(dmf.binding_energy(10,10),-160.8553,1.0e-6,"be3");
  t.test_rel(dmf.mass_excess(10,10),-7.25116,3.0e-5,"me3");
  t.test_rel(dmf.binding_energy(82,126),-1637.278,1.0e-6,"be4");
  t.test_rel(dmf.mass_excess(82,126),-22.58337,3.0e-5,"me4");
  t.test_rel(dmf.binding_energy(87,53),-521.1669,1.0e-6,"be5");
  t.test_rel(dmf.mass_excess(87,53),540.7618,1.0e-4,"me5");
  t.test_rel(dmf.binding_energy(37,53),-776.4578,1.0e-6,"be6");
  t.test_rel(dmf.mass_excess(37,53),-78.98053,1.0e-4,"me6");

  // Fit masses to Audi et al.
  nucmass_fit mf;
  nucdist_set(mf.dist,ame);

  double res;
  mf.eval(dmf,res);
  cout << "DZ 10-parameter fit: " << res << endl;
  t.test_rel(res,1.154636,1.0e-5,"dmf10");

  mf.eval(dmf33,res);
  cout << "DZ 33-parameter fit: " << res << endl;
  t.test_rel(res,0.94796338,1.0e-5,"dmf33");

  nucdist_set(mf.dist,amex);

  mf.eval(dmf,res);
  cout << "DZ 10-parameter fit: " << res << endl;
  t.test_rel(res,0.5890547,1.0e-5,"dmf10ex");
  
  mf.eval(dmf33,res);
  cout << "DZ 33-parameter fit: " << res << endl;
  t.test_rel(res,0.3967063,1.0e-5,"dmf33ex");

  if (false) {
    for(size_t i=0;i<5;i++) {
      mf.def_mmin.ntrial*=30;
      mf.def_mmin.verbose=1;
      mf.fit(dmf33,res);
      cout << res << endl;
    }
  }

  // Compare tables to Audi et al.
  mf.eval(dz2,res);
  t.test_rel(res,0.3941622,1.0e-5,"95table");
  cout << "DZ 1995 table: " << res << endl;

  mf.eval(dz,res);
  t.test_rel(res,0.5913276,1.0e-5,"96table");
  cout << "DZ 1996 table: " << res << endl;
  
  // Test the binding energy of lead from tables
  t.test_rel(dz.binding_energy(82,126)/208.0,-7.867,1.0e-3,"dz be");
  t.test_rel(dz2.binding_energy(82,126)/208.0,-7.867,1.0e-3,"dz2 be");

  // Test the binding energy and mass excess from get_nucleus()
  double mass_neutron=o2scl_mks::mass_neutron*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_proton=o2scl_mks::mass_proton*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_electron=o2scl_mks::mass_electron*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_amu=o2scl_mks::unified_atomic_mass*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);

  nucleus n;
  dz.get_nucleus(82,126,n);
  t.test_rel(n.be*o2scl_const::hc_mev_fm/208.0,-7.867,4.0e-3,"ptr be");
  t.test_rel(n.m-126.0*mass_neutron-82.0*mass_proton,n.be,1.0e-12,"ptr be2");
  t.test_rel(n.m+82.0*mass_electron-208.0*mass_amu,n.mex,1.0e-11,
	     "ptr mex");
  
  dz2.get_nucleus(82,126,n);
  t.test_rel(n.be*o2scl_const::hc_mev_fm/208.0,-7.867,4.0e-3,"ptr be3");
  t.test_rel(n.m-126.0*mass_neutron-82.0*mass_proton,n.be,1.0e-12,"ptr be4");
  t.test_rel(n.m+82.0*mass_electron-208.0*mass_amu,n.mex,1.0e-11,
	     "ptr mex2");

  t.report();
  return 0;
}

