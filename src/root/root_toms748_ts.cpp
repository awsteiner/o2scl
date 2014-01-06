/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2014, Andrew W. Steiner
  
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
#include <o2scl/root_toms748.h>
#include <o2scl/funct.h>

double gfn(double x) {
  return sin(x-0.2);
}

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

#ifdef O2SCL_CPP11

  double a, b;
  funct11 f=gfn;

  root_toms748<> rt;
  a=1.0e-5;
  b=1.0;
  rt.solve_bkt(a,b,f);
  t.test_rel(a,0.2,1.0e-6,"1");

#endif

  t.report();
  return 0;
}

