/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/polylog.h>

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
cpp_dec_float_25;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
cpp_dec_float_35;

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<75> >
cpp_dec_float_75;

typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);
  
  gen_test_number<> gn;
  gen_test_number<long double> gn_ld;
  gen_test_number<cpp_dec_float_25> gn_cdf25;
  
  gn.set_radix(1.7);
  gn_ld.set_radix(1.7);
  gn_cdf25.set_radix(1.7);
  
  fermi_dirac_integ_bf<double,30,40,50,cpp_dec_float_25,
                       cpp_dec_float_35,cpp_dec_float_50> fdib;
  // double is typically 1.0e-15, so we keep an extra two digits 
  fdib.set_tol(1.0e-17);

  fermi_dirac_integ_bf<long double,30,40,50,cpp_dec_float_25,
                       cpp_dec_float_35,cpp_dec_float_50> fdib2;
  // long double is typically 1.0e-18, but we do two digits larger
  fdib2.set_tol(1.0e-20);
  
  fermi_dirac_integ_bf<cpp_dec_float_25,30,40,50,cpp_dec_float_35,
                       cpp_dec_float_50,cpp_dec_float_75> fdib3;
  // cpp_dec_float_25 is typically 1.0e-25, but we do two digits larger
  fdib3.set_tol(1.0e-27);

  for(size_t i=0;i<90;i++) {
    
    double x=gn.gen(), res, err;
    int method, method2, method3;
    int iret=fdib.calc_1o2_ret_full(x,res,err,method);
    
    long double x_ld=gn_ld.gen(), res_ld, err_ld;
    int iret2=fdib2.calc_1o2_ret_full(x_ld,res_ld,err_ld,method2);
    
    cpp_dec_float_25 x_cdf25=gn_cdf25.gen(), res_cdf25, err_cdf25;
    int iret3=fdib3.calc_1o2_ret_full(x_cdf25,res_cdf25,err_cdf25,method3);
    
    cout.width(4);
    cout << i << " ";
    cout.setf(ios::showpos);
    cout << x << " ";
    cout.unsetf(ios::showpos);
    cout << iret << " " << method << " "
         << dtos(res,0) << " " << dtos(err,0) << endl;
    t.test_gen(iret==0,"double prec");
    t.test_abs(err/res,0.0,1.0e-17,"double prec2");
    cout << "                   ";
    cout << iret2 << " " << method2 << " "
         << dtos(res_ld,0) << " " << dtos(err_ld,0) << endl;
    t.test_gen(iret2==0,"long double prec");
    t.test_abs<long double>(err_ld/res_ld,0.0,1.0e-20,"long double prec2");
    cout << "                   ";
    cout << iret3 << " " << method3 << " "
         << dtos(res_cdf25,0) << " " << dtos(err_cdf25,0) << endl;
    t.test_gen(iret3==0,"long double prec");
    t.test_abs_boost<cpp_dec_float_25>(err_cdf25/res_cdf25,
                                       0.0,1.0e-27,"long double prec2");
  }
  cout << endl;
  
  t.report();
  return 0;
}
