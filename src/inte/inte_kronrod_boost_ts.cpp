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
  t.set_output_level(1);

  {
    inte_kronrod_boost<funct,61> ikb;
    
    double ans, exact, err;
    
    funct tf=test_func<double>;
    
    // Compare with the exact result
    ikb.integ_err(tf,0.0,1.0,ans,err);
    exact=cos(100.0)-cos(1/1.01);
    std::cout << ans << " " << err << std::endl;
    t.test_rel(ans,exact,1.0e-8,"qag test");
  }

  {
    inte_kronrod_boost<funct_ld,61,long double> ikb;
    
    long double ans, exact, err;
    
    funct_ld tf=test_func<long double>;
    
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
    
    funct_cdf50 tf=test_func<cpp_dec_float_50>;
    
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

  funct f1=test_func<double>;
  funct_ld f2=test_func<long double>;
  funct_cdf25 f3=test_func<cpp_dec_float_25>;
  funct_cdf35 f4=test_func<cpp_dec_float_35>;
  funct_cdf50 f5=test_func<cpp_dec_float_50>;
  funct_cdf100 f6=test_func<cpp_dec_float_100>;
  funct_multip_wrapper fmw(f1,f2,f3,f4,f5,f6);
  funct_multip<> fm(fmw);

  /*
  inte_multip_kronrod_boost<funct_multip<>> imkb;
  double a=0.0, b=1.0, res, err, exact;
  exact=cos(100.0)-cos(1/1.01);
  imkb.integ_err(fm,a,b,res,err);
  cout << res << " " << exact << endl;

  long double a_ld=0.0, b_ld=1.0, res_ld, err_ld, exact_ld;
  exact_ld=cos(100.0L)-cos(1.0L/1.01L);
  imkb.integ_err(fm,a_ld,b_ld,res_ld,err_ld);
  cout << res_ld << " " << exact_ld << endl;
  */
  
  t.report();
  return 0;
}

