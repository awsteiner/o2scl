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
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <o2scl/funct.h>
#include <o2scl/inte_gauss_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &a);

double testfun(double tx, double &a) {
  return (a*sin(tx)/(tx+0.01));
}

double testfun2(double tx) {
  return 4.0*std::sqrt(1.0-tx*tx);
}

long double testfun2_ld(long double tx) {
  return 4.0*sqrtl(1.0-tx*tx);
}

#ifdef O2SCL_LD_TYPES

__float128 testfun2_f128(__float128 tx) {
  return 4.0*sqrtl(1.0-tx*tx);
}

cpp_dec_float_50 testfun2_cdf(cpp_dec_float_50 tx) {
  return 4.0*sqrt(1.0-tx*tx);
}

typedef std::function<__float128(__float128)> funct_f128;
typedef std::function<cpp_dec_float_50(cpp_dec_float_50)> funct_cdf;

#endif

int main(void) {
  
  cout.setf(ios::scientific);
  cout.precision(10);
  
  test_mgr t;
  t.set_output_level(2);

  {
    inte_gauss_cern<funct> cg;
    double a=3.0, calc, exact, diff;
  
    funct tf=std::bind(testfun,std::placeholders::_1,a);

    calc=cg.integ(tf,0.0,1.0);
    exact=a*(0.900729064796877177);
    t.test_rel(calc,exact,1.0e-8,"inte_gauss_cern 1");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;

    funct tf2=testfun2;

    calc=cg.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<double>();
    t.test_rel(calc,exact,1.0e-8,"inte_gauss_cern 2");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;
  }

  {
    inte_gauss_cern<funct_ld,long double,
		    inte_gauss_coeffs_long_double> cg_ld;
    cg_ld.tol_rel=1.0e-20;
    cg_ld.tol_abs=1.0e-20;
    long double a=3.0, calc, exact, diff;

    funct_ld tf2=testfun2_ld;
    
    calc=cg_ld.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<long double>();
    t.test_rel<long double>(calc,exact,1.0e-16L,"inte_gauss_cern ld");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;
  }
  
  {
    inte_gauss_cern<funct_f128,__float128,
		    inte_gauss_coeffs_float128> cg_f128;
    cg_f128.tol_rel=1.0e-20;
    cg_f128.tol_abs=1.0e-20;
    __float128 a=3.0, calc, exact, diff;

    funct_f128 tf2=testfun2_f128;
    
    calc=cg_f128.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<__float128>();
    //t.test_rel<__float128>(calc,exact,1.0e-16L,"inte_gauss_cern f128");
    diff=fabsl(calc-exact);
    cout << ((long double)calc) << " "
	 << ((long double)exact) << " "
	 << ((long double)diff) << endl;
  }
  
  {
    inte_gauss_cern<funct_cdf,cpp_dec_float_50,
		    inte_gauss_coeffs_cpp_dec_float_50> cg_cdf;
    cg_cdf.tol_rel=1.0e-22;
    cg_cdf.tol_abs=1.0e-22;
    cpp_dec_float_50 a=3.0, calc, exact, diff;

    funct_cdf tf2=testfun2_cdf;
    
    calc=cg_cdf.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<cpp_dec_float_50>();
    //t.test_rel<cpp_dec_float_50>(calc,exact,1.0e-16L,"inte_gauss_cern cdf");
    diff=abs(calc-exact);
    cout << std::setprecision(std::numeric_limits<
			      cpp_dec_float_50>::digits10);
    cout << "calc:  " << calc << endl;
    cout << "exact: " << exact << endl;
    cout << "diff:  " << diff << endl;
  }
  
  t.report();
  return 0;
}

