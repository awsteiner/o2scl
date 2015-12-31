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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
 
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/min_cern.h>
 
using namespace std;
using namespace o2scl;
 
double minfun(double x);

double minfun(double x) {
  return -exp(-(x-0.5)*(x-0.5));
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  double x, min=0.0;
  int vp=0;
  
  cout.setf(ios::scientific);
  
  min_cern<funct11> mb;
  funct11 mf=minfun;

  x=0.2;
  mb.min_bkt(x,-1.0,1.0,min,mf);
  t.test_rel(x,0.5,1.0e-5,"val");
  t.test_rel(min,-1.0,1.0e-5,"min");

  /// Demonstrate that an initial guess for x is not needed
  x=2.0;
  mb.min_bkt(x,-1.0,1.0,min,mf);
  t.test_rel(x,0.5,1.0e-5,"val");
  t.test_rel(min,-1.0,1.0e-5,"min");

  t.report();
  return 0;
}

