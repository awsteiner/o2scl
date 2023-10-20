/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2019-2023, Andrew W. Steiner
  
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
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/funct_to_fp.h>
#include <o2scl/set_mpfr.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

template<class fp_t> fp_t test_func(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  fp_t res=-sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
  return res;
}

template<class fp_t> fp_t test_func_i(fp_t x) {
  return exp(-x*x);
}

template<class fp_t, class fp2_t> fp_t test_func_param
(fp_t x, fp2_t a) {
  fp_t one=1;
  fp_t a2=static_cast<fp_t>(a);
  return -sin(one/(x+a2))/(x+a2)/(x+a2);
}

template<class fp_t> fp_t param_f(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  fp_t ret=one/hundred;
  return ret;
}

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  inte_kronrod_boost<61> imkb;
#ifdef O2SCL_SET_MPFR
  inte_kronrod_boost<61,mpfr_25,mpfr_35,mpfr_50,mpfr_100> imkb_mpfr;
#endif
  
  {
    // Integrate test_func over [0,1] and compare to the exact result
    // at double precision
    
    double ans, exact, err;
    
    funct tf=test_func<double>;
    
    imkb.integ_err(tf,0.0,1.0,ans,err);
    exact=cos(100.0)-cos(1/1.01);
    std::cout << ans << " " << err << std::endl;
    t.test_rel(ans,exact,1.0e-8,"imkb double");
    cout << endl;
  }

  {
    // Integrate test_func over [0,1] and compare to the exact result
    // at long double precision
    
    long double ans, exact, err;
    
    funct_ld tf=test_func<long double>;
    
    imkb.tol_rel=1.0e-13;
    imkb.set_max_depth(15);
    imkb.integ_err(tf,0.0L,1.0L,ans,err);
    exact=cos(100.0L)-cos(1.0L/1.01L);
    std::cout << ans << " " << err << std::endl;
    t.test_rel<long double>(ans,exact,1.0e-15,"imkb long double");
    cout << endl;
  }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  {
    // Integrate test_func over [0,1] and compare to the exact result
    // at 50-digit precision
    
    cpp_dec_float_50 ans, exact, err;
    
    funct_cdf50 tf=test_func<cpp_dec_float_50>;
    
    // Compare with the exact result
    imkb.tol_rel=1.0e-30;
    imkb.set_max_depth(25);
    cpp_dec_float_50 zero=0;
    cpp_dec_float_50 one=1;
    cpp_dec_float_50 hundred=100;
    imkb.integ_err(tf,zero,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<cpp_dec_float_50>(ans,exact,1.0e-45,"imkb 50-digit");
    cout << endl;
  }

#ifdef O2SCL_SET_MPFR
  {
    // Integrate test_func over [0,1] and compare to the exact result
    // at 50-digit precision
    
    mpfr_50 ans, exact, err;
    
    funct_mpfr50 tf=test_func<mpfr_50>;
    
    // Compare with the exact result
    imkb.tol_rel=1.0e-30;
    imkb.set_max_depth(25);
    mpfr_50 zero=0;
    mpfr_50 one=1;
    mpfr_50 hundred=100;
    imkb.integ_err(tf,zero,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<mpfr_50>(ans,exact,1.0e-45,"imkb_mpfr 50-digit");
    cout << endl;

    // Return tol_rel to its default value
    imkb.tol_rel=0.0;
  }
  
  {
    // Integrate test_func over [0,1] and compare to the exact result
    // at 50-digit precision
    
    mpfr_50 ans, exact, err;
    
    funct_mpfr50 tf=test_func<mpfr_50>;
    
    // Compare with the exact result
    imkb.tol_rel=1.0e-30;
    imkb.set_max_depth(25);
    mpfr_50 zero=0;
    mpfr_50 one=1;
    mpfr_50 hundred=100;
    imkb.integ_err(tf,zero,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<mpfr_50>(ans,exact,1.0e-45,"imkb_mpfr 50-digit");
    cout << endl;

    // Return tol_rel to its default value
    imkb.tol_rel=0.0;
  }
#endif
  
  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);

    // Multiprecision integration
    
    imkb.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2,1.0e-8);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-8,"multip 1");
    cout << endl;

#ifdef O2SCL_OSX
    // AWS, 10/17/23, this doesn't work on the docker images, possibly
    // because boost was installed without quadmath or mpfr, but I'm not
    // sure, so I'm just commenting them out for now
    imkb.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip 2");
    cout << endl;
#endif

#ifdef O2SCL_SET_MPFR
    // AWS, 10/17/23, this doesn't work on the docker images, possibly
    // because boost was installed without quadmath or mpfr, but I'm not
    // sure, so I'm just commenting them out for now
    imkb_mpfr.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2,1.0e-8);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-8,"multip 1 mpfr");
    cout << endl;
    
#ifdef O2SCL_OSX
    imkb_mpfr.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip 2 mpfr");
    cout << endl;
#endif

    // AWS, 10/17/23, this doesn't work on the docker images, possibly
    // because boost was installed without quadmath or mpfr, but I'm not
    // sure, so I'm just commenting them out for now
    imkb_mpfr.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2,1.0e-8);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-8,"multip 1 mpfr");
    cout << endl;
    
#ifdef O2SCL_OSX
    imkb_mpfr.integ_err_multip([](auto &&tb) mutable { return test_func(tb); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip 2 mpfr");
    cout << endl;
#endif
#endif

    // Multiprecision integration with infinite limits
    
    imkb.integ_iu_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                             a,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi/2.0,1.0e-15,"multip 3");
    cout << endl;

    imkb.integ_il_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                             a,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi/2.0,1.0e-15,"multip 4");
    cout << endl;

    imkb.integ_i_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                            val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi,1.0e-15,"multip 5");
    cout << endl;

    // AWS 3/16/23: These tests may work (the integrator fails without
    // calling the error handler as requested), but I'm commenting
    // them out because they take a very long time to run.
    if (false) {
      // Try a function which is difficult to integrate
      imkb.err_nonconv=false;
      int ret=imkb.integ_il_err_multip([](auto &&tb) mutable
      { return test_func(tb); },a,val,err2);
      cout << "ret: " << ret << endl;
      t.test_gen(ret!=0,"fail il");
      cout << endl;

      ret=imkb.integ_i_err_multip([](auto &&tb) mutable
      { return test_func(tb); },a,val,err2);
      cout << "ret: " << ret << endl;
      t.test_gen(ret!=0,"fail i");
      cout << endl;
      imkb.err_nonconv=true;
    }
    
    // Multiprecision integration with a template function which has
    // a parameter with a fixed type

    cpp_dec_float_25 one=1;
    cpp_dec_float_25 hundred=100;
    cpp_dec_float_25 param=one/hundred;

    imkb.tol_rel=0.0;
    imkb.integ_err_multip([param](auto &&tb) mutable
    { return test_func_param(tb,param); },a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip param 1");
    cout << endl;
      
    // Multiprecision integration with a template function which has
    // a template parameter
      
    imkb.integ_err_multip([param](auto &&tb) mutable
    { return test_func_param(tb,param_f(tb)); },a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip param 2");
    cout << endl;
    
#ifdef O2SCL_OSX
    
    // Multiprecision integration with a funct_multip_string object
      
    funct_multip_string fms;
    fms.set_function("-sin(1/(x+1/100))/(x+1/100)^2","x");
    funct_multip_string *fmsp=&fms;
      
    imkb.integ_err_multip([fmsp](auto &&tb) mutable { return (*fmsp)(tb); },
                          a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip string");
    cout << endl;
    
#ifdef O2SCL_SET_MPFR
    imkb_mpfr.integ_err_multip([fmsp](auto &&tb) mutable {
      return (*fmsp)(tb); },a,b,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,exact,1.0e-15,"multip string mpfr");
    cout << endl;
#endif
    
#endif  
      
  }

#endif
  
  t.report();
  return 0;
}

