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
#include <o2scl/test_mgr.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/polylog.h>
#include <o2scl/set_mpfr.h>

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

typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;

#ifdef O2SCL_SET_MPFR
typedef boost::multiprecision::mpfr_float_50 mpfr_float_50;

typedef boost::multiprecision::mpfr_float_100 mpfr_float_100;
#endif

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);
  
  // Compare polylog values with hard-coded values
  polylog<> p;
  polylog_multip<double,long double> p2;
  t.test_rel(p.calc(2,-0.5),-0.448414206923646,4.0e-15,"pl 1");
  t.test_rel(p2.calc(2,-0.5),-0.448414206923646,4.0e-15,"pl 1b");
  t.test_rel(p.calc(2,-2.0),-1.43674636688368,4.0e-15,"pl 2");
  t.test_rel(p.calc(3,-0.5),-0.472597844658897,4.0e-15,"pl 3");
  t.test_rel(p.calc(3,-2.0),-1.66828336396657,4.0e-15,"pl 4");
  t.test_rel(p.calc(2,0.5),0.5822405264650125,4.0e-15,"pl 5");
  t.test_rel(p.calc(3,0.5),0.5372131936080402,4.0e-15,"pl 6");

  bessel_K_exp_integ_gsl be_gsl;
  bessel_K_exp_integ_boost<double,double> be_boost;
  bessel_K_exp_integ_boost<double,long double> be_boost2;
  bessel_K_exp_integ_direct<> be_d_ld;

  cout.precision(15);
  cout << be_gsl.K1exp(3.0e2) << endl;
  cout << be_boost.K1exp(3.0e2) << endl;
  cout << be_gsl.K1exp(7.0e2) << endl;
  cout << be_boost.K1exp(7.0e2) << endl;

  // Demonstrate we can use long double precision to get around
  // difficulties with large arguments
  cout << be_gsl.K1exp(1.0e3) << endl;
  cout << be_boost.K1exp(1.0e3) << endl;
  cout << be_boost2.K1exp(1.0e3) << endl;
  cout.precision(6);
  
  // Compare bessel_K_exp GSL and direct
  t.test_rel(be_gsl.K1exp(2.0),be_d_ld.K1exp(2.0),1.0e-15,"be_d_ld 1");
  t.test_rel(be_gsl.K2exp(2.0),be_d_ld.K2exp(2.0),1.0e-15,"be_d_ld 2");
  t.test_rel(be_gsl.K3exp(2.0),be_d_ld.K3exp(2.0),1.0e-15,"be_d_ld 3");

  gen_test_number<> gn;
  gen_test_number<long double> gn_ld;
  gen_test_number<cpp_dec_float_35> gn_cdf35;
  gen_test_number<cpp_dec_float_50> gn_cdf50;

  // Test Fermi-Dirac integrals
  
  if (true) {
    
    fermi_dirac_integ_gsl fdig;
    
    fermi_dirac_multip fdib2;
    fdib2.set_tol(1.0e-28);

    fermi_dirac_integ_direct<cpp_dec_float_35,funct_cdf50,20,
                             cpp_dec_float_50> fd_35_50;
    fd_35_50.set_tol(1.0e-28);
    
    for(size_t i=0;i<60;i+=2) {
      
      double x=gn.gen();
      cpp_dec_float_35 x_cdf35=gn_cdf35.gen();
      
      double y1=fdig.calc_3(x);
      fdib2.set_tol(1.0e-28);
      cpp_dec_float_35 y5=fdib2.calc_3(x_cdf35);
      cout.width(4);
      cout << i << " ";
      cout.setf(ios::showpos);
      cout << x << " ";
      cout.unsetf(ios::showpos);
      
      t.test_gen(abs(y1-y5)/y1<1.0e-15,"fd");
      
      // fd_35_50 starts to fail once i=49
      if (i<49) {
        cpp_dec_float_35 y4=fd_35_50.calc_3(x_cdf35);
        cout << y1 << " " << y4 << " " << y5 << " "
             << abs(y4-y5)/y4 << endl;
        t.test_gen(abs(y4-y5)/y4<1.0e-28,"fd 2");
      } else {
        cout << y1 << " " << 0.0 << " " << y5 << endl;
      }
      
    }
  }

  bessel_K_exp_integ_direct<long double,funct_cdf35,20,
			    cpp_dec_float_35> be_ld_35;
  be_ld_35.set_tol(1.0e-21);
  bessel_K_exp_integ_direct<cpp_dec_float_35,funct_cdf50,20,
                            cpp_dec_float_50> be_35_50;
  be_35_50.set_tol(1.0e-37);
  
  if (true) {

    gn.reset();
    gn_ld.reset();
    gn_cdf35.reset();
    
    // This section compares the GSL class bessel_K_exp_integ_gsl with
    // o2scl versions of with various accuracies and floating point
    // types. However, the o2scl version with the default types object
    // named 'be_d_ld', currently fails for large enough arguments.
    
    for(size_t i=0;i<4;i++) {
      
      double x=gn.gen();
      long double x_ld=gn_ld.gen();
      cpp_dec_float_35 x_cdf35=gn_cdf35.gen();

      if (x>0) {
      
        double y1=be_gsl.K1exp(x);
        double y2=be_d_ld.K1exp(x);
        long double y3=be_ld_35.K1exp(x_ld);
        cpp_dec_float_35 y4=be_35_50.K1exp(x_cdf35);
        cout.width(4);
        cout << i << " ";
        cout.setf(ios::showpos);
        cout << x << " ";
        cout.unsetf(ios::showpos);
        cout << dtos(y1,0) << endl;
        cout << "                 " << dtos(y2,0) << endl;
        cout << "                 " << dtos(y3,0) << endl;
        cout << "                 " << dtos(y4,0) << endl;
        double d21=y2-y1;
        long double d32=y3-y2;
        long double d31=y3-y1;
        cpp_dec_float_35 d43=y4-y3;
        cout << "                 " << abs(d21) << " " << abs(d31) << " "
             << abs(d32) << " " << abs(d43) << endl;
        t.test_gen(abs(d43)*1.0e2<abs(d32),
                   "bessel_K1_exp long double accuracy");
        t.test_gen(abs(d32)<=abs(d31),"bessel_K1_exp o2scl better than gsl");
      
        y1=be_gsl.K3exp(x);
        y2=be_d_ld.K3exp(x);
        y3=be_ld_35.K3exp(x_ld);
        y4=be_35_50.K3exp(x_cdf35);
        cout.width(4);
        cout << i << " ";
        cout.setf(ios::showpos);
        cout << x << " ";
        cout.unsetf(ios::showpos);
        cout << dtos(y1,0) << endl;
        cout << "                 " << dtos(y2,0) << endl;
        cout << "                 " << dtos(y3,0) << endl;
        cout << "                 " << dtos(y4,0) << endl;
        d21=y2-y1;
        d32=y3-y2;
        d31=y3-y1;
        d43=y4-y3;
        cout << "                 " << abs(d21) << " " << abs(d31) << " "
             << abs(d32) << " " << abs(d43) << endl;
        t.test_gen(abs(d43)*1.0e2<abs(d32),
                   "bessel_K3_exp long double accuracy");
        t.test_gen(abs(d32)<=abs(d31),
                   "bessel_K3_exp o2scl better than gsl");

      }
      
    }
  }
  
  t.report();
  
  return 0;
}
