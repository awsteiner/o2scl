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

#include <o2scl/deriv_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double x) {
  return sin(x);
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  int vp=0;
  
  funct_fptr tf(testfun);
  deriv_cern<funct> de;
  
  cout.setf(ios::scientific);
  cout.precision(10);
  cout << "result exact" << endl;
  cout << de.deriv(0.5,tf) << " " << cos(0.5) << endl;
  t.test_rel(de.deriv(0.5,tf),cos(0.5),1.0e-6,"1st derivative.");
  cout << de.deriv2(0.5,tf) << " " 
       << -sin(0.5) << endl;
  t.test_rel(de.deriv2(0.5,tf),-sin(0.5),1.0e-4,"2nd derivative.");
  cout << de.deriv3(0.5,tf) << " " 
       << -cos(0.5) << endl;
  t.test_rel(de.deriv3(0.5,tf),-cos(0.5),1.0e-4,"3rd derivative.");
  cout << endl;

  de.verbose=1;
  de.deriv(0.5,tf);
  
  t.report();
  return 0;
}

