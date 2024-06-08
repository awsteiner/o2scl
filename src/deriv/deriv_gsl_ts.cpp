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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/deriv_gsl.h>

// AWS, 7/22/22, commenting this out because this is a relatively
// new function in boost and not available everywhere yet
#ifndef O2SCL_OLD_BOOST
#include <boost/math/differentiation/finite_difference.hpp>
#endif
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;

typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;

template<class fp_t> fp_t sin_fun(fp_t x) {
  return sin(x);
}

template<class fp_t> fp_t sqrt_fun(fp_t x) {
  return sqrt(x);
}

template<class fp_t> fp_t difficult_fun(fp_t x) {
  return exp(x)/(cos(x)*cos(x)*cos(x)+sin(x)*sin(x)*sin(x));
}

template<class fp_t> fp_t difficult_deriv(fp_t x) {
  fp_t den=(cos(x)*cos(x)*cos(x)+sin(x)*sin(x)*sin(x));
  return exp(x)*(2*cos(3*x)+3*sin(x)+sin(3*x))/2/den/den;
}

class tempc {
public:
  double operator()(double x) {
    return std::sin(x);
  }
};

int main(void) {
  double res;
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  cout.precision(10);
  
  deriv_gsl<funct> de;
  
  // Demonstrate that one can instantiate the class with
  // a user-defined type as well
  deriv_gsl<tempc> de2;

  funct tf=sin_fun<double>;
  funct tf2=sqrt_fun<double>;

  de.h=1.0e-4;

  res=de.deriv(0.5,tf);
  cout << "First derivative: " << endl;
  cout << res << " " << de.get_err() 
       << " " << std::cos(0.5) << endl;
  t.test_rel(res,std::cos(0.5),1.0e-11,"simple derivative");

  // Show how to use boost to compute the same derivative
  double err;
  // AWS, 7/22/22, commenting this out because this is a relatively
  // new function in boost and not available everywhere yet
#ifndef O2SCL_OLD_BOOST
  cout << boost::math::differentiation::finite_difference_derivative
    (tf,0.5,&err) << " " << err << endl;
#endif
  
  cout << "Second derivative: " << endl;
  res=de.deriv2(0.5,tf);
  cout << res << " " 
       << de.get_err() << " " << -sin(0.5) << endl;
  t.test_rel(res,-sin(0.5),1.0e-6,"second derivative");

  cout << "Third derivative: " << endl;
  res=de.deriv3(0.5,tf);
  cout << res << " " 
  << de.get_err() << " " << -std::cos(0.5) << endl;
  t.test_rel(res,-std::cos(0.5),1.0e-2,"third derivative");
  cout << endl;

  // Test a derivative where the step size initially
  // takes the function into non-finite values
  de.verbose=1;
  res=de.deriv(1.0e-6,tf2);
  cout << "First derivative: " << endl;
  cout << res << " " << de.get_err() 
       << " " << 1.0/2.0/std::sqrt(1.0e-6) << endl;
  t.test_rel(res,1.0/2.0/std::sqrt(1.0e-6),1.0e-7,"first, non-finite");
  cout << endl;

  // Try a long double derivative
  deriv_gsl<funct_ld,long double> de_ld;
  funct_ld tf_ld=sin_fun<long double>;
  long double ld_res;

  ld_res=de_ld.deriv(0.5L,tf_ld);
  cout << "First derivative: " << endl;
  cout << ld_res << " " << de.get_err() 
       << " " << std::cos(0.5L) << endl;
  t.test_rel(ld_res,std::cos(0.5L),4.0e-14L,"simple derivative long double");

  deriv_gsl<funct_cdf50,cpp_dec_float_50> de_cdf;
  funct_cdf50 tf_cdf=sin_fun<cpp_dec_float_50>;
  cpp_dec_float_50 cdf_res;
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 two=2;
  cpp_dec_float_50 half=one/two;

  cdf_res=de_cdf.deriv(half,tf_cdf);
  cout << "First derivative: " << endl;
  cout << cdf_res << " " << de.get_err() 
       << " " << cos(half) << endl;
  t.test_rel_boost<cpp_dec_float_50>(cdf_res,cos(half),1.0e-50,
				     "simple derivative cpp_dec_float_50");

  // A difficult function
  funct df=difficult_fun<double>;
  res=de.deriv(5.5,df);
  cout << "First derivative: " << endl;
  cout << res << " " << de.get_err() 
       << " " << difficult_deriv(5.5) << endl;
  t.test_rel(res,difficult_deriv(5.5),1.0e-8,"simple derivative");
  cout << endl;

  // Automatic multiprecision derivative
  double val, err2;
  deriv_multip_gsl dmg2;
  dmg2.verbose=1;
  dmg2.deriv_err([](auto &&tx) mutable { return difficult_fun(tx); },5.5,
                 val,err2);
  t.test_rel((long double)val,difficult_deriv<long double>(5.5L),1.0e-15L,
             "multip 1");
  t.test_abs(0.0,err2/val,1.0e-15,"multip 2");
  
  t.report();
  return 0;
}
