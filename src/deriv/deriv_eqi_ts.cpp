/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

#include <o2scl/deriv_eqi.h>
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>

using namespace std;
using namespace o2scl;

double testfun(double x) {
  return sin(x);
}

int main(void) {
  deriv_eqi<funct> de;
  test_mgr t;
  t.set_output_level(2);
  int i;
  int vp=0;

  cout.setf(ios::scientific);
  cout.precision(12);
  
  funct_fptr tf(testfun);
  
  de.xoff=0.0;

  de.set_npoints(2);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(3);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(4);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(5);
  cout << de.deriv(0.5,tf) << endl;
  cout << cos(0.5) << endl;
  cout << endl;
  
  de.xoff=0.2*de.h;
  
  de.set_npoints(2);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(3);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(4);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(5);
  cout << de.deriv(0.5,tf) << endl;
  cout << cos(0.5+0.2*de.h) << endl;
  cout << endl;

  de.xoff=de.h;
  
  de.set_npoints(2);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(3);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(4);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(5);
  cout << de.deriv(0.5,tf) << endl;
  cout << cos(0.5+de.h) << endl;
  cout << endl;
  
  de.xoff=-de.h;
  
  de.set_npoints(2);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(3);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(4);
  cout << de.deriv(0.5,tf) << endl;
  de.set_npoints(5);
  cout << de.deriv(0.5,tf) << endl;
  cout << cos(0.5-de.h) << endl;
  cout << endl;

  int bign=30;
  typedef boost::numeric::ublas::vector<double> ubvector;
  ubvector x(bign), y(bign), dydx(bign), dydx2(bign);
  
  cout.precision(6);

  for(i=0;i<bign;i++) {
    x[i]=((double)i)/3.0;
    y[i]=exp(x[i]);
    cout << x[i] << " " << y[i] << " ";
    cout << de.deriv(0.5,tf) << endl;
  }
  
  de.h=1.0e-4;
  de.xoff=0.0;
  cout << "Check second derivatives:" << endl;
  de.set_npoints2(3);
  cout << de.deriv2(0.5,tf) << endl;
  de.set_npoints2(4);
  cout << de.deriv2(0.5,tf) << endl;
  de.set_npoints2(5);
  cout << de.deriv2(0.5,tf) << endl;
  cout << -sin(0.5) << endl;

  t.report();
  return 0;
}

