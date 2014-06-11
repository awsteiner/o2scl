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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/deriv_gsl.h>

using namespace std;
using namespace o2scl;

double testfun(double x) {
  return sin(x);
}

class tempc {
public:
  double operator()(double x) {
    return sin(x);
  }
};

int main(void) {
  int vp=0;
  double res;
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  cout.precision(10);
  
  deriv_gsl<funct11> de;
  
  // Demonstrate that one can instantiate the class with
  // a user-defined type as well
  deriv_gsl<tempc> de2;

  funct11 tf=testfun;

  de.h=1.0e-4;

  res=de.deriv(0.5,tf);
  cout << "First derivative: " << endl;
  cout << res << " " << de.get_err() 
       << " " << cos(0.5) << endl;
  t.test_rel(res,cos(0.5),1.0e-11,"simple derivative");

  cout << "Second derivative: " << endl;
  res=de.deriv2(0.5,tf);
  cout << res << " " 
       << de.get_err() << " " << -sin(0.5) << endl;
  t.test_rel(res,-sin(0.5),1.0e-6,"second derivative");

  cout << "Third derivative: " << endl;
  res=de.deriv3(0.5,tf);
  cout << res << " " 
  << de.get_err() << " " << -cos(0.5) << endl;
  t.test_rel(res,-cos(0.5),1.0e-2,"third derivative");
  cout << endl;

  de.verbose=1;
  res=de.deriv(0.5,tf);

  t.report();
  return 0;
}
