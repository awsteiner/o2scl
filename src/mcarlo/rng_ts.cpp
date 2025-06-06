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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <random>

#include <o2scl/rng.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  rng<> nr;
  rng<long double> nr_ld;

  for(size_t j=0;j<100;j++) {
    vector<double> arr={0,1,2,3,4,5,6,7,8,9};
    vector_shuffle<vector<double>,double>(nr,arr.size(),arr);
    vector_out(cout,arr,true);
  }

  t.report();
  return 0;
}
