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
#include <o2scl/funct.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &pa);

double testfun(double tx, double &pa) {
  double a=pa;
  return (a*sin(tx)/(tx+0.01));
}

int main(void) {
  double a=3.0, calc, exact, diff;
  inte_gauss56_cern<funct11> cg;
  test_mgr t;
  t.set_output_level(2);

  funct11 tf=std::bind(testfun,std::placeholders::_1,a);

  cout.setf(ios::scientific);
  cout.precision(10);

  calc=cg.integ(tf,0.0,1.0);
  exact=a*(0.900729064796877177);
  t.test_rel(calc,exact,1.0e-2,"inte_gauss56_cern");
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << endl;
  
  t.report();
  return 0;
}

