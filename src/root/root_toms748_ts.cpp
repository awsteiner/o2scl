/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2020, Andrew W. Steiner
  
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
#include <o2scl/root_toms748.h>
#include <o2scl/funct.h>

double gfn(double x) {
  return sin(x-0.2);
}

#ifdef O2SCL_LD_TYPES

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

long double gfn_ld(long double x) {
  return sin(x-0.2);
}

cpp_dec_float_50 gfn_cdf(cpp_dec_float_50 x) {
  return sin(x-0.2);
}

#endif

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  double a, b;
  funct f=gfn;

  root_toms748<> rt;
  a=1.0e-5;
  b=1.0;
  rt.solve_bkt(a,b,f);
  t.test_rel(a,0.2,1.0e-6,"1");

#ifdef O2SCL_LD_TYPES

  long double a_ld, b_ld;
  funct_ld f_ld=gfn_ld;

  root_toms748<funct_ld,long double> rt_ld;
  a_ld=1.0e-5;
  b_ld=1.0;
  rt_ld.solve_bkt(a_ld,b_ld,f_ld);
  t.test_rel<long double>(a_ld,0.2,1.0e-6,"2");

  cpp_dec_float_50 a_cdf, b_cdf;
  funct_cdf50 f_cdf=gfn_cdf;

  root_toms748<funct_cdf50,cpp_dec_float_50> rt_cdf;
  a_cdf=1.0e-5;
  b_cdf=1.0;
  rt_cdf.solve_bkt(a_cdf,b_cdf,f_cdf);
  t.test_rel_boost<cpp_dec_float_50>(a_cdf,0.2,1.0e-6,"1");
  
#endif
  
  t.report();
  return 0;
}

