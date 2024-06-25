/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/err_hnd.h>
#include <o2scl/constants.h>
#include <o2scl/funct_multip.h>
#include <o2scl/set_mpfr.h>
#include <o2scl/funct_to_fp.h>
#include <o2scl/set_python.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class fmc {
public:
  template<class fp_t> fp_t func(fp_t x) {
    return log(1+x);
  }
};

class fmc2 {
  
public:
  
  template<class fp_t> fp_t param_f(fp_t x) {
    fp_t one=1;
    fp_t ten=10;
    fp_t ret=one+one/ten;
    return ret;
  }
  
  template<class fp_t, class fp2_t> fp_t func(fp_t x, fp2_t a) {
    fp_t a2=static_cast<fp_t>(a);
    return log(a2+x);
  }
  
};

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  // The class funct_string is from funct_to_fp.h
  funct_string f("pi*r^2","r");
  f.set_parm("pi",o2scl_const::pi);
  for(double r=1.0;r<=2.0;r+=0.1) {
    t.test_rel(f(r),o2scl_const::pi*r*r,1.0e-12,"funct_string");
  }

#ifdef O2SCL_SET_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("../../data/o2scl/python",3);
  {
    cout << "1." << endl;
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python fp("funct_test","fun",2);
    t.test_rel(fp(2.0),o2scl_const::pi*2.0,1.0e-12,"funct_python");
  }
  {
    cout << "2." << endl;
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python_method fp("funct_test","c","fun2",2);
    t.test_rel(fp(2.0),o2scl_const::pi*4.0,1.0e-12,"funct_python 2");
  }
  {
    cout << "3." << endl;
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python fp("o2sclpy","cpp_test",2);
    t.test_rel(fp(2.0),o2scl_const::pi*2.0,1.0e-12,"funct_python 3");
  }

  o2scl_settings.py_final();
  
#endif
  
  t.report();
  return 0;
}

