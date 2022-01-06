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
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  funct_string f("pi*r^2","r");
  f.set_parm("pi",o2scl_const::pi);
  for(double r=1.0;r<=2.0;r+=0.1) {
    t.test_rel(f(r),o2scl_const::pi*r*r,1.0e-12,"funct_string");
  }

#ifdef O2SCL_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("./");
  {
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python fp("funct_test","fun");
    t.test_rel(fp(2.0),o2scl_const::pi*2.0,1.0e-12,"funct_python");
  }
  o2scl_settings.py_final();
  
#endif
  
  t.report();
  return 0;
}

