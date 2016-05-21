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

#include <iostream>

#include <o2scl/rng_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  rng_gsl nr(10);

  double a1=nr.random();
  double a2=nr.random();

  rng_gsl nr2(10);

  // Test to make sure that given the same seed, we 
  // come up with the same numbers.
  t.test_rel(a1,nr2.random(),1.0e-14,"First random number.");
  t.test_rel(a2,nr2.random(),1.0e-14,"Second random number.");

  rng_gsl nr3=nr;

  t.test_rel(nr3.random(),nr.random(),1.0e-14,"Copy constructor.");

  cout << "Random integers [0,10): " << flush;
  for(int i=0;i<10;i++) {
    cout << nr2.random_int(10) << " " << flush;
  }
  cout << endl;

  t.report();
  return 0;
}
