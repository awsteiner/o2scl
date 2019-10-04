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
#include <o2scl/funct.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &pa);

double testfun(double tx, double &pa) {
  return (pa*sin(tx)/(tx+0.01));
}

long double testfun_ld(long double tx, long double &pa) {
  return (pa*sinl(tx)/(tx+0.01L));
}

#ifdef O2SCL_LD_TYPES
__float128 testfun_f128(__float128 tx, __float128 &pa) {
  return (pa*sinl(tx)/(tx+0.01));
}

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 tx, cpp_dec_float_50 &pa) {
  return (pa*sin(tx)/(tx+0.01));
}

/*

// based on
// https://gcc.gnu.org/onlinedocs/libquadmath/quadmath_005fsnprintf.html
std::string f128_to_s(const __float128 &x) {
  char buf[128];
  int width=46;
  int n=quadmath_snprintf(buf,128,"%Qa",x);
  if (n>=128 || n<=0) {
    cout << "Fail." << endl;
    exit(-1);
  }
  std::string s=buf;
  return s;
}

*/

typedef std::function<__float128(__float128)> funct_f128;
typedef std::function<cpp_dec_float_50(cpp_dec_float_50)> funct_cdf;

#endif

int main(void) {
  cout.setf(ios::scientific);
  cout.precision(10);

  test_mgr t;
  t.set_output_level(2);

  double a=3.0, calc, exact, diff;
  inte_gauss56_cern<funct> cg;

  funct tf=std::bind(testfun,std::placeholders::_1,a);

  calc=cg.integ(tf,0.0,1.0);
  exact=a*(0.900729064796877177);
  t.test_rel(calc,exact,1.0e-2,"inte_gauss56_cern");
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << endl;

  // Moving to long double here doesn't really improve the accuracy
  // for this particular function at the moment, but it verifies that
  // inte_gauss56_cern works with the long double type.
  
  inte_gauss56_cern<funct_ld,long double,
		    inte_gauss56_coeffs_long_double> cg_ld;
  long double a_ld=3.0L, calc_ld, exact_ld, diff_ld;

  funct_ld tf_ld=std::bind(testfun_ld,std::placeholders::_1,a_ld);
  
  calc_ld=cg_ld.integ(tf_ld,0.0L,1.0L);
  exact_ld=a_ld*(0.900729064796877177036268L);
  t.test_rel<long double>(calc_ld,exact_ld,1.0e-2L,"inte_gauss56_cern ls");
  diff_ld=fabsl(calc_ld-exact_ld);
  cout << calc_ld << " " << exact_ld << " " << diff_ld << endl;

  // ------------------------------------------------------

  inte_gauss56_cern<funct_f128,__float128,
		    inte_gauss56_coeffs_float128> cg_f128;
  __float128 a_f128=3.0L, calc_f128, exact_f128, diff_f128;

  funct_f128 tf_f128=std::bind(testfun_f128,std::placeholders::_1,a_f128);
  
  calc_f128=cg_f128.integ(tf_f128,0.0L,1.0L);
  exact_f128=a_f128*(0.900729064796877177036268L);
  //t.test_rel<__float128>(calc_f128,exact_f128,1.0e-2L,
  //"inte_gauss56_cern ls");
  diff_f128=fabsl(calc_f128-exact_f128);
  cout << ((long double)calc_f128) << " "
       << ((long double)exact_f128) << " "
       << ((long double)diff_f128) << endl;
  
  // ------------------------------------------------------

  inte_gauss56_cern<funct_cdf,cpp_dec_float_50,
		    inte_gauss56_coeffs_cpp_dec_float_50> cg_cdf;
  cpp_dec_float_50 a_cdf=3.0L, calc_cdf, exact_cdf, diff_cdf;

  funct_cdf tf_cdf=std::bind(testfun_cdf,std::placeholders::_1,a_cdf);
  
  calc_cdf=cg_cdf.integ(tf_cdf,0.0L,1.0L);
  exact_cdf=a_cdf*(0.900729064796877177036268L);
  //t.test_rel<cpp_dec_float_50>(calc_cdf,exact_cdf,1.0e-2L,
  //"inte_gauss56_cern ls");
  diff_cdf=boost::multiprecision::abs(calc_cdf-exact_cdf);
  cout << calc_cdf << " " << exact_cdf << " " << diff_cdf << endl;
  
  t.report();
  return 0;
}

