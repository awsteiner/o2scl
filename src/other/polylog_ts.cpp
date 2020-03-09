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

#include <o2scl/test_mgr.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/polylog.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);
  
#ifdef O2SCL_LD_TYPES

#ifdef O2SCL_NEW_BOOST_INTEGRATION
  
  fermion_nr_integ_gsl f1;
  fermion_nr_integ_direct f2;

  t.test_rel(f1.calc_1o2(0.5),f2.calc_1o2(0.5),4.0e-16,"fd 1");
  t.test_rel(f1.calc_m1o2(0.5),f2.calc_m1o2(0.5),4.0e-16,"fd 2");
  t.test_rel(f1.calc_3o2(0.5),f2.calc_3o2(0.5),4.0e-16,"fd 3");

  t.test_rel(f2.calc_polylog(2,-0.5),-0.448414206923646,1.0e-14,"pl 1");
  t.test_rel(f2.calc_polylog(2,-2.0),-1.43674636688368,1.0e-14,"pl 2");
  t.test_rel(f2.calc_polylog(3,-0.5),-0.472597844658897,1.0e-14,"pl 3");
  t.test_rel(f2.calc_polylog(3,-2.0),-1.66828336396657,1.0e-14,"pl 4");
  
#endif
  
#endif
  
  t.report();
  return 0;
}
