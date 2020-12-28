/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2021, Andrew W. Steiner
  
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
#ifdef O2SCL_NEW_BOOST_INTEGRATION
#include <o2scl/inte_kronrod_boost.h>
#endif

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;

// This function oscillates quite rapidly near x=0
double test_func(double x) {
  return -sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
}

#ifdef O2SCL_LD_TYPES

long double test_func_ld(long double x) {
  return -sin(1.0L/(x+0.01L))*pow(x+0.01L,-2.0L);
}

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

cpp_dec_float_50 test_func_cdf(cpp_dec_float_50 x) {
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 hundred=100;
  return -sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
}

#endif

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);

#ifdef O2SCL_NEW_BOOST_INTEGRATION
  
  {
    inte_kronrod_boost<funct,61> ikb;
    
    double ans, exact, err;
    
    funct tf=test_func;
    
    // Compare with the exact result
    ikb.integ_err(tf,0.0,1.0,ans,err);
    exact=cos(100.0)-cos(1/1.01);
    std::cout << ans << " " << err << std::endl;
    t.test_rel(ans,exact,1.0e-8,"qag test");
  }

#ifdef O2SCL_LD_TYPES

  {
    inte_kronrod_boost<funct_ld,61,long double> ikb;
    
    long double ans, exact, err;
    
    funct_ld tf=test_func_ld;
    
    // Compare with the exact result
    ikb.tol_rel=1.0e-14;
    ikb.set_max_depth(15);
    ikb.integ_err(tf,0.0L,1.0L,ans,err);
    exact=cos(100.0L)-cos(1.0L/1.01L);
    std::cout << ans << " " << err << std::endl;
    t.test_rel<long double>(ans,exact,1.0e-15,"qag test");
  }

  typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
  {
    inte_kronrod_boost<funct_cdf50,61,cpp_dec_float_50> ikb;
    
    cpp_dec_float_50 ans, exact, err;
    
    funct_cdf50 tf=test_func_cdf;
    
    // Compare with the exact result
    ikb.tol_rel=1.0e-30;
    ikb.set_max_depth(25);
    cpp_dec_float_50 one=1;
    cpp_dec_float_50 hundred=100;
    ikb.integ_err(tf,0.0L,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<cpp_dec_float_50>(ans,exact,1.0e-30,"qag test");
  }

#endif
#endif
  
  t.report();
  return 0;
}

