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
#include <o2scl/multi_funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/funct_to_fp.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  vector<string> list(2);
  list[0]="x1";
  list[1]="x2";
  multi_funct_strings f("x1*x2-x1",2,list);
  ubvector x(2);
  x[0]=1.5;
  x[1]=2.1;
  t.test_rel(x[0]*x[1]-x[0],f(2,x),1.0e-12,"multi_funct_strings");

#ifdef O2SCL_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("../../data/o2scl/python/");
  {
    // We use the brackets to force the multi_funct_python
    // destructor to run before py_final()
    multi_funct_python<ubvector> fp("multi_funct_test","fun");
    t.test_rel(fp(2,x),o2scl_const::pi*3.6,1.0e-12,"multi_funct_python");
  }
  o2scl_settings.py_final();
  
#endif
  
  t.report();
  return 0;
}
