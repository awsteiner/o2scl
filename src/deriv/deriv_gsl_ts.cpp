/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/deriv_gsl.h>

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;

double testfun(double x) {
  return std::sin(x);
}

double testfun2(double x) {
  return std::sqrt(x);
}

long double testfun_ld(long double x) {
  return std::sin(x);
}

#ifdef O2SCL_LD_TYPES

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 x) {
  return sin(x);
}

#endif

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

  funct tf=testfun;
  funct tf2=testfun2;

  de.h=1.0e-4;

  res=de.deriv(0.5,tf);
  cout << "First derivative: " << endl;
  cout << res << " " << de.get_err() 
       << " " << std::cos(0.5) << endl;
  t.test_rel(res,std::cos(0.5),1.0e-11,"simple derivative");

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

#ifdef O2SCL_LD_TYPES
  
  // Try a long double derivative
  deriv_gsl<funct_ld,long double> de_ld;
  funct_ld tf_ld=testfun_ld;
  long double ld_res;

  ld_res=de_ld.deriv(0.5L,tf_ld);
  cout << "First derivative: " << endl;
  cout << ld_res << " " << de.get_err() 
       << " " << std::cos(0.5L) << endl;
  t.test_rel(ld_res,std::cos(0.5L),4.0e-14L,"simple derivative");
  
  deriv_gsl<funct_cdf50,cpp_dec_float_50> de_cdf;
  funct_cdf50 tf_cdf=testfun_cdf;
  cpp_dec_float_50 cdf_res;
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 two=2;
  cpp_dec_float_50 half=one/two;

  cdf_res=de_cdf.deriv(half,tf_cdf);
  cout << "First derivative: " << endl;
  cout << cdf_res << " " << de.get_err() 
       << " " << cos(half) << endl;
  t.test_rel_boost<cpp_dec_float_50>(cdf_res,cos(half),
				     4.0e-20,"simple derivative");
  
#endif
  
  t.report();
  return 0;
}
