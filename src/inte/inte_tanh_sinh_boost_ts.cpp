/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019, Andrew W. Steiner
  
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
#include <o2scl/funct.h>
#include <o2scl/inte_tanh_sinh_boost.h>

using namespace std;
using namespace o2scl;

// This function oscillates quite rapidly near x=0
double test_func_1(double x) {
  return -sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  inte_tanh_sinh_boost<funct,61> itsb;

  double ans, exact, err;

  cout.setf(ios::scientific);

  funct tf1=test_func_1;

  // Compare with the exact result
  itsb.integ_err(tf1,0.0,1.0,ans,err);
  exact=cos(100.0)-cos(1/1.01);
  std::cout << ans << " " << err << std::endl;
  t.test_rel(ans,exact,1.0e-8,"qag test");

  
  t.report();
  return 0;
}

