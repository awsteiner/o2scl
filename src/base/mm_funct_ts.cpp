/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/mm_funct.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  vector<string> list(2), funcs(2);
  list[0]="x1";
  list[1]="x2";
  funcs[0]="x1*x2-x1";
  funcs[1]="x1*x2-x2";
  mm_funct_strings f(2,funcs,list);
  ubvector x(2), y(2);
  x[0]=1.5;
  x[1]=2.1;
  f(2,x,y);
  t.test_rel(x[0]*x[1]-x[0],y[0],1.0e-12,"mm_funct_strings");
  t.test_rel(x[0]*x[1]-x[1],y[1],1.0e-12,"mm_funct_strings");

#ifdef O2SCL_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("./");
  {
    // We use the brackets to force the mm_funct_python
    // destructor to run before py_final()
    mm_funct_python<ubvector> fp("mm_funct_test","fun");
    int mfp_ret=fp(2,x,y);
    ubvector y2(2);
    y2[0]=x[0]*o2scl_const::pi;
    y2[1]=x[1]*o2scl_const::pi;
    t.test_rel_vec(2,y,y2,1.0e-12,"mm_funct_python");
  }
  o2scl_settings.py_final();
  
#endif
  
  t.report();
  return 0;
}

