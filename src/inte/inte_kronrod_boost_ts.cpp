/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2022, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_kronrod_boost.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;

typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

template<class fp_t> fp_t test_func(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  return -sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
}

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  inte_kronrod_boost<61> imkb;
  
  {
    
    double ans, exact, err;
    
    funct tf=test_func<double>;
    
    // Compare with the exact result
    imkb.integ_err(tf,0.0,1.0,ans,err);
    exact=cos(100.0)-cos(1/1.01);
    std::cout << ans << " " << err << std::endl;
    t.test_rel(ans,exact,1.0e-8,"qag test");
  }

  {
    
    long double ans, exact, err;
    
    funct_ld tf=test_func<long double>;
    
    // Compare with the exact result
    imkb.tol_rel=1.0e-13;
    imkb.set_max_depth(15);
    imkb.integ_err(tf,0.0L,1.0L,ans,err);
    exact=cos(100.0L)-cos(1.0L/1.01L);
    std::cout << ans << " " << err << std::endl;
    t.test_rel<long double>(ans,exact,1.0e-15,"qag test");
  }

  typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
  {
    
    cpp_dec_float_50 ans, exact, err;
    
    funct_cdf50 tf=test_func<cpp_dec_float_50>;
    
    // Compare with the exact result
    imkb.tol_rel=1.0e-30;
    imkb.set_max_depth(25);
    cpp_dec_float_50 zero=0;
    cpp_dec_float_50 one=1;
    cpp_dec_float_50 hundred=100;
    imkb.integ_err(tf,zero,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<cpp_dec_float_50>(ans,exact,1.0e-30,"qag test");
  }

  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);

    funct_multip_string fms;
    fms.set_function("-sin(one/(x+1/100))/(x+1/100)^2","x");
    fms.verbose=2;

    imkb.verbose=2;
    imkb.integ_err_multip([](auto &&t) mutable { return test_func(t); },
                          a,b,val,err2,1.0e-8);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    imkb.integ_err_multip([](auto &&t) mutable { return test_func(t); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    imkb.integ_err_multip([fms](auto &&t) mutable { return fms(t); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip");
  }
  
  t.report();
  return 0;
}

