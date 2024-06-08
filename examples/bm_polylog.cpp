/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#include <o2scl/test_mgr.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/polylog.h>

#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>

using namespace std;
using namespace o2scl;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
cpp_dec_float_25;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
cpp_dec_float_35;

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::mpfr_float_50 mpfr_float_50;

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<75> >
cpp_dec_float_75;

typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
typedef boost::multiprecision::mpfr_float_100 mpfr_float_100;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);
  
  gen_test_number<> gn;
  gen_test_number<long double> gn_ld;
  gen_test_number<cpp_dec_float_25> gn_cdf25;
  
  gn.set_radix(1.4);
  gn_ld.set_radix(1.4);
  gn_cdf25.set_radix(1.4);
  
  fermi_dirac_integ_bf<double,30,40,50,cpp_dec_float_25,
                       cpp_dec_float_35,cpp_dec_float_50> fdib;
  // double is typically 1.0e-15, so we keep an extra two digits 
  fdib.set_tol(1.0e-17);

  fermi_dirac_integ_bf<long double,30,40,50,cpp_dec_float_25,
                       cpp_dec_float_35,cpp_dec_float_50> fdib2;
  // long double is typically 1.0e-18, but we do two digits larger
  fdib2.set_tol(1.0e-20);
  
  fermi_dirac_integ_bf<long double,30,40,50,cpp_dec_float_25,
                       cpp_dec_float_35,mpfr_float_50> fdib2b;
  // long double is typically 1.0e-18, but we do two digits larger
  fdib2.set_tol(1.0e-20);
  
  fermi_dirac_integ_bf<cpp_dec_float_25,30,40,50,cpp_dec_float_35,
                       cpp_dec_float_50,cpp_dec_float_75> fdib3;
  // cpp_dec_float_25 is typically 1.0e-25, but we do two digits larger
  fdib3.set_tol(1.0e-27);

  fermi_dirac_integ_bf<cpp_dec_float_25,30,40,50,cpp_dec_float_35,
                       mpfr_float_50,mpfr_float_100> fdib3b;
  // cpp_dec_float_25 is typically 1.0e-25, but we do two digits larger
  fdib3b.set_tol(1.0e-29);

  bessel_K_exp_integ_gsl bkeg;
  
  bessel_K_exp_integ_bf<double,30,40,50,cpp_dec_float_25,
                        cpp_dec_float_35,cpp_dec_float_50> bkeb;
  // double is typically 1.0e-15, so we keep an extra two digits 
  bkeb.set_tol(1.0e-17);

  bessel_K_exp_integ_bf<long double,30,40,50,cpp_dec_float_25,
                        cpp_dec_float_35,cpp_dec_float_50> bkeb2;
  // long double is typically 1.0e-18, but we do two digits larger
  bkeb2.set_tol(1.0e-20);
  
  bessel_K_exp_integ_bf<cpp_dec_float_25,30,40,50,cpp_dec_float_35,
                        cpp_dec_float_50,cpp_dec_float_75> bkeb3;
  // cpp_dec_float_25 is typically 1.0e-25, but we do two digits larger
  bkeb3.set_tol(1.0e-27);

  if (false) {

    double x;
    long double x_ld;
    cpp_dec_float_25 x_cdf25;
    
    for(size_t i=1;i<100;i++) {
    
      if (i==1) {
        x=0;
        x_ld=0;
        x_cdf25=0;
      }

      if (i!=1) {
        x=-x;
        x_ld=-x_ld;
        x_cdf25=-x_cdf25;
      }
      
      if (i%2==0) {
        x+=6;
        x_ld+=6;
        x_cdf25+=6;
      }
      
      double res, err;
      int method, method2, method3, method2b, method3b;
      int iret=fdib.calc_1o2_ret_full(x,res,err,method);
    
      long double res_ld, err_ld, res_ldb, err_ldb;
      int iret2=fdib2.calc_1o2_ret_full(x_ld,res_ld,err_ld,method2);
      int iret2b=fdib2b.calc_1o2_ret_full(x_ld,res_ldb,err_ldb,method2b);
    
      cpp_dec_float_25 res_cdf25, err_cdf25;
      cpp_dec_float_25 res_cdf25b, err_cdf25b;
      int iret3=fdib3.calc_1o2_ret_full(x_cdf25,res_cdf25,err_cdf25,method3);
      int iret3b=fdib3b.calc_1o2_ret_full(x_cdf25,res_cdf25b,err_cdf25b,
                                          method3b);
    
      cout.width(4);
      cout << i << " ";
      cout.setf(ios::showpos);
      
      cout << x << " ";
      cout.unsetf(ios::showpos);
      cout << iret << " " << method << " "
           << dtos(res,0) << " " << dtos(err,0) << endl;
      t.test_gen(iret==0,"double ret");
      t.test_abs(err/res,0.0,1.0e-15,"double prec");
      
      cout << "                   ";
      cout << iret2 << " " << method2 << " "
           << dtos(res_ld,0) << " " << dtos(err_ld,0) << endl;
      t.test_gen(iret2==0,"long double ret");
      t.test_abs<long double>(err_ld/res_ld,0.0,1.0e-18,"long double prec");
      
      cout << "                   ";
      cout << iret2b << " " << method2b << " "
           << dtos(res_ldb,0) << " " << dtos(err_ldb,0) << endl;
      t.test_gen(iret2==0,"long double ret 2b");
      t.test_abs<long double>(err_ldb/res_ldb,0.0,1.0e-18,
                              "long double prec 2b");
      
      cout << "                   ";
      cout << iret3 << " " << method3 << " "
           << dtos(res_cdf25,0) << " " << dtos(err_cdf25,0) << endl;
      t.test_gen(iret3==0,"cdf25 ret");
      t.test_abs_boost<cpp_dec_float_25>(err_cdf25/res_cdf25,
                                         0.0,1.0e-25,"cdf25 prec");
      cout << "                   ";
      cout << iret3b << " " << method3b << " "
           << dtos(res_cdf25b,0) << " " << dtos(err_cdf25b,0) << endl;
      t.test_gen(iret3b==0,"cdf25 ret 3b");
      t.test_abs_boost<cpp_dec_float_25>(err_cdf25b/res_cdf25b,
                                         0.0,1.0e-25,"cdf25 prec 3b");
    }
    cout << endl;
  }
  
  gn.reset();
  gn_ld.reset();
  gn_cdf25.reset();

  for(size_t i=0;i<90;i++) {
    
    double x=gn.gen(), res, err;
    int method, method2, method3;
    if (x>0.0) {
      cout << x << " " << bkeg.K1exp(x) << endl;
      int iret=bkeb.K1exp_ret_full(x,res,err,method);
      
      cout.width(4);
      cout << i << " ";
      cout.setf(ios::showpos);
      cout << x << " ";
      cout.unsetf(ios::showpos);
      cout << iret << " " << method << " "
           << dtos(res,0) << " " << dtos(err,0) << endl;
      t.test_gen(iret==0,"double ret");
      t.test_abs(err/res,0.0,1.0e-17,"double prec");
      
      if (true) {
        
        long double x_ld=gn_ld.gen(), res_ld, err_ld;
        int iret2=bkeb2.K1exp_ret_full(x_ld,res_ld,err_ld,method2);
        
        cout << "                   ";
        cout << iret2 << " " << method2 << " "
             << dtos(res_ld,0) << " " << dtos(err_ld,0) << endl;
        t.test_gen(iret2==0,"long double ret");
        t.test_abs<long double>(err_ld/res_ld,0.0,1.0e-20,"long double prec");
      }

      if (false) {
        
        cpp_dec_float_25 x_cdf25=gn_cdf25.gen(), res_cdf25, err_cdf25;
        int iret3=bkeb3.K1exp_ret_full(x_cdf25,res_cdf25,err_cdf25,method3);
        
        cout << "                   ";
        cout << iret3 << " " << method3 << " "
             << dtos(res_cdf25,0) << " " << dtos(err_cdf25,0) << endl;
        t.test_gen(iret3==0,"cdf25 ret");
        t.test_abs_boost<cpp_dec_float_25>(err_cdf25/res_cdf25,
                                           0.0,1.0e-27,"cdf25 prec");
      }
    }
  }
  cout << endl;
  
  t.report();
  return 0;
}
