/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2013-2023, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/min_brent_boost.h>
 
using namespace std;
using namespace o2scl;
 
double minfun(double x);
double minfun2(double x);

double minfun(double x) {
  return -exp(-(x-0.5)*(x-0.5));
}

// A more pathological function with a hidden sharp minimum
double minfun2(double x) {
  return pow(fabs(x),0.01)-1.0;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  double x, min;
  funct f=minfun;
  funct f2=minfun2;
  min_brent_boost<> mb;
  
  x=0.2;
  mb.min_bkt(x,-1.0,1.0,min,f);
  t.test_rel(x,0.5,1.0e-1,"val");
  t.test_rel(min,-1.0,3.0e-4,"min");
  
  mb.min_bkt(x,-1.0,1.0,min,f2);
  t.test_rel(x,0.0,1.0e-4,"val");
  t.test_rel(min,-1.0,1.0,"min");

  t.report();
  return 0;
}

