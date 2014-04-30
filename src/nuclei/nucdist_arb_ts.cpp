/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <o2scl/nucdist_arb.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  /// Load Audi et al. masses
  nucmass_ame ame;
  o2scl_hdf::ame_load(ame,"");

  /// Select only nuclei with N=21 and N=22
  nucdist_arb ad;
  ad.set_dist(ame,"(N>20) & (N<23)");

  size_t cnt=0;
  for(nucdist::iterator ndi=ad.begin();ndi!=ad.end();ndi++) {
    cout << ndi->Z << " " << ndi->N << endl;
    cnt++;
  }
  cout << cnt << " total nuclei." << endl;

  t.test_gen(cnt==40,"nucdist_arb count");

  t.report();
  return 0;
}

