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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_two_interp.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/cloud_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(3);

  nucmass_two_interp nti;
  nti.set_default();

  nucmass_ame ame;
  ame.load("20");
  nucmass_mnmsk mo12;
  o2scl_hdf::mnmsk_load(mo12,"msis16");
  vector<nucleus> mdist;

  nucdist_set(mdist,mo12);

  for(size_t i=0;i<mdist.size();i++) {
    int Z=mdist[i].Z;
    int N=mdist[i].N;
    if (ame.is_included(Z,N)) {
      cout << nti.binding_energy(Z,N)/(Z+N) << " "
           << ame.binding_energy(Z,N)/(Z+N) << " "
           << mo12.binding_energy(Z,N)/(Z+N) << endl;
    } else {
      cout << nti.binding_energy(Z,N)/(Z+N) << " "
           << 0.0 << " "
           << mo12.binding_energy(Z,N)/(Z+N) << endl;
    }
    char ch;
    cin >> ch;
  }
  
  t.report();
  return 0;
}

