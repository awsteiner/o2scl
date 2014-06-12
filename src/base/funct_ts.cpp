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
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double func(double x, double &pa) {
  return pa+x;
}

class ac {
public:
  double mfunc(double x, double &pa) {
    return pa+x;
  }
};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  /*
    double a=2.0;
    funct_fptr_param<double> f1(func,a);
    double x=3.2, y;
    
    y=0.0;
    y=f1(x);
    t.test_rel(y,5.2,1.0e-6,"fptr");
    
    y=0.0;
    ac c1;
    funct_mfptr_param<ac,double> f3(&c1,&ac::mfunc,a);
    y=f3(x);
    t.test_rel(y,5.2,1.0e-6,"mfptr");
  */
  
  t.report();
  return 0;
}

