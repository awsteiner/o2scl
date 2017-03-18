/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017, Andrew W. Steiner
  
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
#include <o2scl/interp.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // ---------------------------------------------------------------
  // Create test data

  ubvector x(5), y(5), rx(5), ry(5);
  double xa[5], ya[5], rxa[5], rya[5];
  for(size_t i=0;i<5;i++) {
    x[i]=((double)i);
    y[i]=((double)i*i);
    rx[4-i]=((double)i);
    ry[4-i]=((double)i*i);
    xa[i]=((double)i);
    ya[i]=((double)i*i);
    rxa[4-i]=((double)i);
    rya[4-i]=((double)i*i);
  }

  // ---------------------------------------------------------------
  //

  interp_krige<> gi(itp_linear);

  // ---------------------------------------------------------------
  // Test normal interpolation

  double x0=2.5;
  double y0;
  
  y0=gi.eval(x0,5,x,y);
  t.test_rel(y0,6.5,1.0e-5,"intp 1a.");
  y0=gi.deriv(x0,5,x,y);
  t.test_rel(y0,5.0,1.0e-5,"intp 1b.");
  y0=gi.deriv2(x0,5,x,y);
  t.test_rel(y0,0.0,1.0e-5,"intp 1c.");
  y0=gi.integ(x0,3.5,5,x,y);
  t.test_rel(y0,9.25,1.0e-5,"intp 1d.");

  t.report();

  return 0;
}
