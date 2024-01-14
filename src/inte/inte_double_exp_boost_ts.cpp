/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2019-2024, Andrew W. Steiner
  
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
#include <o2scl/funct.h>
#include <o2scl/funct_multip.h>
#include <o2scl/inte_double_exp_boost.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

template<class fp_t> fp_t test_func(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  return -sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
}

template<class fp_t> fp_t test_func_i(fp_t x) {
  return exp(-x*x);
}

template<class fp_t> fp_t exact_func() {
  fp_t one=1;
  fp_t hundred=100;
  fp_t exact=cos(hundred)-cos(one/(one+one/hundred));
  return exact;
}

template<class fp_t> fp_t test_func_2(fp_t x) {
  return exp(-x*x);
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  if (true) {
    
    inte_double_exp_boost<funct> ideb;

    double ans, exact, err;

    funct tf1=test_func<double>;
    funct tf2=test_func_2<double>;

    // Finite integral, moderately difficult integrand
    cout << "tanh-sinh, finite interval:" << endl;
    ideb.integ_err(tf1,0.0,1.0,ans,err);
    exact=exact_func<double>();
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    t.test_rel(ans,exact,1.0e-8,"tanh_sinh test");
    cout << endl;

    // Semi-infinite domains
    cout << "tanh-sinh, infinite upper limit:" << endl;
    exact=sqrt(acos(-1))/2.0;
    ideb.integ_err(tf2,0.0,std::numeric_limits<double>::infinity(),ans,err);
    t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 2");
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    cout << endl;

    cout << "tanh-sinh, infinite lower limit:" << endl;
    ideb.integ_err(tf2,-std::numeric_limits<double>::infinity(),0.0,ans,err);
    t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 3");
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    cout << endl;

    // Infinite domain
    cout << "tanh-sinh, infinite lower and upper limit:" << endl;
    exact=sqrt(acos(-1));
    ideb.integ_err(tf2,-std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   ans,err);
    t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 4");
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    cout << endl;

    // exp_sinh integration, native
    cout << "exp-sinh, infinite upper limit:" << endl;
    ideb.integ_iu_err(tf2,0.0,ans,err);
    exact=sqrt(acos(-1.0))/2.0;
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    t.test_rel(ans,exact,1.0e-8,"exp_sinh test 1");
    cout << endl;

    cout << "exp-sinh, infinite upper limit, explicit limits:" << endl;
    ideb.integ_err(tf2,0.0,std::numeric_limits<double>::infinity(),ans,err);
    exact=sqrt(acos(-1.0))/2.0;
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    t.test_rel(ans,exact,1.0e-8,"exp_sinh test 2");
    cout << endl;
  
    cout << "exp-sinh, infinite lower limit:" << endl;
    ideb.integ_err(tf2,-std::numeric_limits<double>::infinity(),0.0,ans,err);
    exact=sqrt(acos(-1.0))/2.0;
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    t.test_rel(ans,exact,1.0e-8,"exp_sinh test 3");
    cout << endl;
  
    cout << "sinh-sinh, infinite interval:" << endl;
    ideb.integ_i_err(tf2,ans,err);
    exact=sqrt(acos(-1.0));
    std::cout << ans << " " << err << " " << exact << " "
              << ideb.L1norm << std::endl;
    t.test_rel(ans,exact,1.0e-8,"sinh_sinh test");
    cout << endl;

#ifndef O2SCL_NO_BOOST_MULTIPRECISION

    funct_ld tf1_ld=test_func<long double>;
    funct_ld tf2_ld=test_func_2<long double>;
    funct_cdf50 tf1_cdf=test_func<cpp_dec_float_50>;
    funct_cdf50 tf2_cdf=test_func_2<cpp_dec_float_50>;
    funct_cdf35 tf1_cdf35=test_func<cpp_dec_float_35>;
    funct_cdf35 tf2_cdf35=test_func_2<cpp_dec_float_35>;

    inte_double_exp_boost<funct_ld,long double> ideb_ld;
    inte_double_exp_boost<funct_cdf35,cpp_dec_float_35> ideb_cdf35;
    inte_double_exp_boost<funct_cdf50,cpp_dec_float_50> ideb_cdf;

    // Finite integral, moderately difficult integrand
    cout << "tanh-sinh, finite interval, long double:" << endl;
    long double one_ld=1;
    long double zero_ld=0;
    long double hundred_ld=100;
    long double exact_ld=exact_func<long double>();
    long double ans_ld, err_ld;
    ideb_ld.integ_err(tf1_ld,zero_ld,one_ld,ans_ld,err_ld);
    t.test_rel<long double>(ans_ld,exact_ld,1.0e-16,"tanh_sinh test ld");
    std::cout << ans_ld << " " << err_ld << " " << exact_ld << " "
              << ideb_ld.L1norm << std::endl;
    cout << endl;

    // Finite integral, moderately difficult integrand
    cout << "tanh-sinh, finite interval, cpp_dec_float_35:" << endl;
    cpp_dec_float_35 zero_cdf35=0;
    cpp_dec_float_35 one_cdf35=1;
    cpp_dec_float_35 hundred_cdf35=100;
    cpp_dec_float_35 exact_cdf35=exact_func<cpp_dec_float_35>();
    cpp_dec_float_35 ans_cdf35, err_cdf35;
    ideb_cdf35.tol_rel=1.0e-30;
    ideb_cdf35.integ_err(tf1_cdf35,zero_cdf35,one_cdf35,ans_cdf35,err_cdf35);
    std::cout << ans_cdf35 << " " << err_cdf35 << std::endl;
    std::cout << dtos(exact_cdf35,0) << std::endl;
    std::cout << ideb_cdf35.L1norm << std::endl;
    t.test_rel_boost<cpp_dec_float_35>(ans_cdf35,exact_cdf35,1.0e-29,
                                       "tanh_sinh test cdf35");
    cout << endl;

    // Finite integral, moderately difficult integrand
    cout << "tanh-sinh, finite interval, cpp_dec_float_50:" << endl;
    cpp_dec_float_50 zero_cdf=0;
    cpp_dec_float_50 one_cdf=1;
    cpp_dec_float_50 hundred_cdf=100;
    cpp_dec_float_50 exact_cdf=exact_func<cpp_dec_float_50>();
    cpp_dec_float_50 ans_cdf, err_cdf;
    ideb_cdf.tol_rel=1.0e-40;
    ideb_cdf.integ_err(tf1_cdf,zero_cdf,one_cdf,ans_cdf,err_cdf);
    std::cout << ans_cdf << " " << err_cdf << std::endl;
    std::cout << dtos(exact_cdf,0) << std::endl;
    std::cout << ideb_cdf.L1norm << std::endl;
    t.test_rel_boost<cpp_dec_float_50>(ans_cdf,exact_cdf,1.0e-39,
                                       "tanh_sinh test cdf");
    cout << endl;

  }

  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);
    inte_double_exp_boost<> ideb;
    ideb.verbose=2;
    ideb.integ_err_multip([](auto &&tx) mutable { return test_func(tx); },
                           a,b,val,err2,1.0e-8);
    ideb.integ_err_multip([](auto &&tx) mutable { return test_func(tx); },
                           a,b,val,err2);
    t.test_rel(val,exact,1.0e-15,"multip");

    // Make sure infinite integrals work
    ideb.integ_iu_err_multip([](auto &&tx) mutable { return test_func_i(tx); },
                              a,val,err2);
    cout << val << endl;
    t.test_rel(val,root_pi/2.0,1.0e-10,"iu multip");
    cout << endl;
    
    ideb.integ_il_err_multip([](auto &&tx) mutable { return test_func_i(tx); },
                              a,val,err2);
    cout << val << endl;
    t.test_rel(val,root_pi/2.0,1.0e-10,"il multip");
    cout << endl;
    
    ideb.integ_i_err_multip([](auto &&tx) mutable { return test_func_i(tx); },
                             val,err2);
    cout << val << endl;
    t.test_rel(val,root_pi,1.0e-10,"i multip");
    cout << endl;
    
  }

  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);
    inte_double_exp_boost<cpp_dec_float_25,
                                 cpp_dec_float_35,cpp_dec_float_50,
                                 cpp_dec_float_100> ideb;
    ideb.verbose=0;
    ideb.integ_err_multip([](auto &&tx) mutable { return test_func(tx); },
                           a,b,val,err2);
    t.test_rel(val,exact,1.0e-15,"multip cpp_dec_float");
  }

#ifdef O2SCL_OSX
  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);
    inte_double_exp_boost<mpfr_25,mpfr_35,mpfr_50,mpfr_100> ideb;
                                 
    ideb.verbose=0;
    ideb.integ_err_multip([](auto &&tx) mutable { return test_func(tx); },
                           a,b,val,err2);
    t.test_rel(val,exact,1.0e-15,"multip mpfr");
  }
#endif

#endif
  
  t.report();
  return 0;
}

