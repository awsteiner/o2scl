/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class fmc {
public:
  template<class fp_t> fp_t func0(fp_t x) {
    return sin(x);
  }
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

#ifdef O2SCL_SET_MULTIP  
  {
    cout << "0." << endl;
    fmc f2;
    double val=0.0, err;
    funct_multip_cdf fm2;
    
    // No parameters
    fm2.verbose=2;
    fm2.eval_tol_err([f2](auto &&tx) mutable { return f2.func0(tx); },
                     1.0e-4,val,err);
    t.test_rel(val,sin(1.0e-4),1.0e-15,"funct_multip_cdf");
    cout << endl;
  }
  {
    cout << "1." << endl;
    fmc f2;
    double val=0.0, err;
    funct_multip_cdf fm2;
    
    // No parameters
    fm2.verbose=2;
    fm2.eval_tol_err([f2](auto &&tx) mutable { return f2.func(tx); },
                     1.0e-4,val,err);
    t.test_rel(val,log1p(1.0e-4),1.0e-15,"funct_multip_cdf");
    cout << endl;
  }
#ifdef O2SCL_SET_MPFR
  {
    cout << "2." << endl;
    fmc f2;
    double val, err;
    funct_multip_mpfr fm2;
    
    // No parameters
    fm2.eval_tol_err([f2](auto &&tx) mutable { return f2.func(tx); },
                     1.0e-4,val,err);
    t.test_rel(val,log1p(1.0e-4),1.0e-15,"funct_multip_mpfr");
  }
#endif

  {
    // A parameter with a fixed type
    fmc f2;
    double val, err;
    funct_multip fm2;
    fmc2 f3;
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<25>> one=1;
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<25>> ten=10;
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<25>> param=one+one/ten;
    fm2.eval_tol_err([f3,param](auto &&tx) mutable
    { return f3.func(tx,param); },1.0e-4,val,err);
    t.test_rel(val,log1p(0.1001),1.0e-15,"funct_multip 2");
    
    // A fully templated parameter defined by a function 
    fm2.eval_tol_err([f3](auto &&tx) mutable
    { return f3.func(tx,f3.param_f(tx)); },1.0e-4,val,err);
    t.test_rel(val,log1p(0.1001),1.0e-15,"funct_multip 3");
  }
#endif
  
  t.report();
  return 0;
}

