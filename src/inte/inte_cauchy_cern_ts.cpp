/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

#include <o2scl/funct.h>
#include <o2scl/inte_cauchy_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx);
double testfun2(double tx);
double testfun3(double tx);

double testfun(double tx) {
  return exp(-tx)/(tx-2.0);
}

double testfun2(double tx) {
  return log(tx+1.0)/(tx-1.0);
}

double testfun3(double tx) {
  return log(tx+1.0)/pow((exp(tx/10.0)-0.1),2.0);
}

int main(void) {
  double a=3.0, calc, exact, diff;
  inte_cauchy_cern<funct> cc;
  test_mgr t;
  t.set_output_level(2);

  funct tf=testfun;
  funct tf2=testfun2;
  funct tf3=testfun3;

  cout.setf(ios::scientific);
  cout.precision(10);
  
  cc.s=2.0;
  calc=cc.integ(tf,1.0,3.0 );
  t.test_rel(calc,-0.286166,1.0e-5,"inte_cauchy_cern 1");
  
  cc.s=1.0;
  calc=cc.integ(tf2,0.0,2.0 );
  t.test_rel(calc,1.03065,1.0e-5,"inte_cauchy_cern 2");

  cc.s=1.0;
  calc=cc.integ(tf3,0.0,2.0 );
  t.test_rel(calc,1.21569,1.0e-5,"inte_cauchy_cern 3");
  
  t.report();
  return 0;
}

