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
#include <o2scl/funct.h>
#include <o2scl/funct_multip.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

double testfun(double tx, double &pa) {
  return (pa*sin(tx)/(tx+0.01));
}

long double testfun_ld(long double tx, long double &pa) {
  return (pa*sinl(tx)/(tx+0.01L));
}

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 tx, cpp_dec_float_50 &pa) {
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 hundred=100;
  return (pa*sin(tx)/(tx+one/hundred));
}

int main(void) {
  cout.setf(ios::scientific);
  cout.precision(10);

  test_mgr t;
  t.set_output_level(2);

  double a=3.0, calc, exact, diff;
  inte_gauss56_cern<funct,double> cg;

  funct tf=std::bind(testfun,std::placeholders::_1,a);

  calc=cg.integ(tf,0.0,1.0);
  exact=a*(0.900729064796877177);
  t.test_rel(calc,exact,1.0e-2,"inte_gauss56_cern");
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << endl;

  // Moving to long double here doesn't really improve the accuracy
  // for this particular function, but it verifies that
  // inte_gauss56_cern compiles and executes with the long double
  // and boost::multiprecision types
  
  inte_gauss56_cern<funct_ld,long double> cg_ld;
  long double a_ld=3.0L, calc_ld, exact_ld, diff_ld;

  funct_ld tf_ld=std::bind(testfun_ld,std::placeholders::_1,a_ld);
  
  calc_ld=cg_ld.integ(tf_ld,0.0L,1.0L);
  exact_ld=a_ld*(0.900729064796877177036268L);
  t.test_rel<long double>(calc_ld,exact_ld,1.0e-2L,"inte_gauss56_cern ls");
  diff_ld=fabsl(calc_ld-exact_ld);
  cout << calc_ld << " " << exact_ld << " " << diff_ld << endl;
  
  inte_gauss56_cern<funct_cdf50,cpp_dec_float_50> cg_cdf;
  cpp_dec_float_50 a_cdf=3.0L, calc_cdf, exact_cdf, diff_cdf;
  
  funct_cdf50 tf_cdf=std::bind(testfun_cdf,std::placeholders::_1,a_cdf);
  
  calc_cdf=cg_cdf.integ(tf_cdf,0.0L,1.0L);
  exact_cdf=a_cdf*(0.900729064796877177036268L);
  t.test_rel_boost<cpp_dec_float_50>(calc_cdf,exact_cdf,1.0e-2L,
				     "inte_gauss56_cern ls");
  diff_cdf=fabs(calc_cdf-exact_cdf);
  cout << calc_cdf << " " << exact_cdf << " " << diff_cdf << endl;
  
  t.report();
  return 0;
}

