/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/nucmass_frdm.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  nucmass_frdm mo;
  nucmass_ame au;
  nucmass_semi_empirical sm;

  o2scl_hdf::ame_load(au,"12");

  // Show that the nucmass_frdm gives reasonable (but not great)
  // values for the binding energy of Lead 208
  cout << "AME2003 : ";
  cout << au.mass_excess(82,126) << " ";
  cout << au.binding_energy(82,126) << " ";
  cout << au.total_mass(82,126) << endl;
  cout << "Simple  : ";
  cout << sm.mass_excess(82,126) << " ";
  cout << sm.binding_energy(82,126) << " ";
  cout << sm.total_mass(82,126) << endl;
  cout << "FRDM    : ";
  cout << mo.mass_excess(82,126) << " ";
  cout << mo.binding_energy(82,126) << " ";
  cout << mo.total_mass(82,126) << endl;

  // What is this for?
  cout << au.binding_energy_d(32.33,295)/o2scl_const::hc_mev_fm << endl;
  cout << sm.binding_energy_d(32.33,295)/o2scl_const::hc_mev_fm << endl;
  cout << mo.binding_energy_d(32.33,295)/o2scl_const::hc_mev_fm << endl;
  cout << mo.drip_binding_energy_d(32.33,295,0,0,0)/
    o2scl_const::hc_mev_fm << endl;
  cout << mo.mass_excess_d(32.33,295)/o2scl_const::hc_mev_fm << endl;
  cout << mo.drip_mass_excess_d(32.33,295,0,0,0)/o2scl_const::hc_mev_fm 
       << endl;
  
  // Explicitly add the microscopic part of the binding energy for 
  // Lead 208 from the table in Moller et al. (1995) 
  // and compare with the theoretical value quoted in that table
  t.test_rel(mo.mass_excess(82,126)-12.84,-21.15,5.0e-4,"Lead 208");

  // Compare nucmass_frdm with the macroscopic parts from
  // nucmass_mnmsk table and show that they're almost the same
  nucmass_mnmsk mm;
  o2scl_hdf::mnmsk_load(mm);
  nucmass_mnmsk::entry mme;
  double comp=0.0;
  size_t nnuc=0;
  vector<nucleus> fd;
  nucdist_set(fd,mm);
  for(vector<nucleus>::iterator ndi=fd.begin();ndi!=fd.end();ndi++) {
    if (ndi->N>=8 && ndi->Z>=8) {
      mme=mm.get_ZN(ndi->Z,ndi->N);
      comp+=pow(mo.mass_excess(ndi->Z,ndi->N)-(mme.Mth-mme.Emic),2.0);
      nnuc++;
    }
  }
  comp=sqrt(comp/nnuc);
  t.test_abs(comp,0.0,5.0e-3,"Compare with Mth-Emic from full FRDM");
  
  cout << nnuc << " " << comp << endl;

  // Fit the macroscopic part to experiment
  nucmass_fit mf;
  nucdist_set(mf.dist,au);
  //mf.set_exp_mass(au);
  double qual;
  mf.eval(mo,qual);
  cout << "Before fit: " << qual << endl;
  t.test_rel(qual,3.1860,1.0e-2,"FRDM pre-fit.");
  mf.def_mmin.ntrial*=10;
  mf.fit(mo,qual);
  cout << "After fit: " << qual << endl;
  t.test_rel(qual,2.3836,1.0e-2,"FRDM post-fit.");

  t.report();

  return 0;
}
