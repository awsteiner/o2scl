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

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/funct.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/test_mgr.h>

#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

template<class fp_t>
fp_t testfun(fp_t tx, fp_t &a) {
  return -cos(1/(tx+a))/(a+tx)/(a+tx);
}

template<class fp_t> fp_t test_func_i(fp_t x) {
  return exp(-x*x);
}

template<class fp_t>
fp_t testfun2(fp_t tx) {
  fp_t a=1;
  fp_t a2=100;
  fp_t a3=a/a2;
  return -cos(1/(tx+a3))/(a3+tx)/(a3+tx);
}

template<class fp_t> fp_t sin_recip(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  return sin(1/(-x+one/hundred))*pow(-x+one/hundred,-2);
}

template<class fp_t> fp_t mp_test_func(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  return -sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
}

typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

#ifdef O2SCL_MPFR
typedef boost::multiprecision::mpfr_float_50 mpfr_float_50;
#endif

template<class func_t=funct, class fp_t=double,
	 class def_inte_t, size_t nsub>
void test_iac(test_mgr &t, func_t &f, fp_t acc,
	      std::string comment, fp_t &diff, bool output_sub=false) {
  
  cout << comment << ":\n  ";
  fp_t zero=0;
  fp_t one=1;
  fp_t ten=10;
  fp_t hundred=100;
  fp_t a=one/hundred, calc, ei;
  fp_t exact=sin(one/(one+a))-sin(one/a);
  
  //inte_adapt_cern<func_t,def_inte_t,nsub,fp_t> iac;
  //iac.tol_rel=acc;
  //iac.tol_abs=acc;

  inte_adapt_cern imac;
  inte_subdiv<fp_t> is(nsub);
  imac.set_nsub(nsub);
  imac.tol_rel=(double)acc;
  imac.tol_abs=(double)acc;
  
  //iac.integ_err(f,0.0,one,calc,ei);
  imac.integ_err_is(f,zero,one,calc,ei,is);
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << " "
       << ei << endl;
  cout << "  subdivisions: ";
  cout << is.get_nsubdivisions() << endl;
  cout << endl;

  if (output_sub) {
    size_t n=is.get_nsubdivisions();
    typedef boost::numeric::ublas::vector<fp_t> ubvector;
    ubvector xlo(n), xhi(n), val(n), err(n);
    is.get_subdivisions(xlo,xhi,val,err);
    cout << "  xlo              xhi               ";
    cout << "val              err              " << endl;;
    for(size_t i=0;i<n;i+=10) {
      cout << "  " << xlo[i] << " " << xhi[i] << " ";
      cout.setf(ios::showpos);
      cout << val[i] << " ";
      cout.unsetf(ios::showpos);
      cout << err[i] << endl;
    }
    cout << endl;
  }

  return;
}

int main(void) {
  
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  cout.precision(10);

  {
    cout << "inte_adapt_cern, double, testfun:\n  ";
    
    double a=0.01, diff;
    funct tf=std::bind(testfun<double>,std::placeholders::_1,a);
    test_iac<funct,double,
	     inte_gauss56_cern<funct,double>,100>
      (t,tf,1.0e-8,"iac, double, testfun",diff,true);
    t.test_abs<double>(diff,0.0,1.0e-7,"inte_adapt_cern");

    cout << "inte_adapt_cern, long double, testfun:\n  ";
    
    long double a_ld=0.01L, diff_ld;
    funct_ld tf_ld=std::bind(testfun<long double>,
                             std::placeholders::_1,a_ld);
    test_iac<funct_ld,long double,
	     inte_gauss56_cern<funct_ld,long double>,1000>
      (t,tf_ld,1.0e-15,"iac, long double, testfun",diff_ld);
    t.test_abs<long double>(diff_ld,0.0,1.0e-14,"inte_adapt_cern_ld");

    cout << "inte_adapt_cern, cpp_dec_float_50, testfun:\n  ";
    
    cpp_dec_float_50 one=1.0, diff_cdf;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 a_cdf=one/hundred;
    funct_cdf50 tf_cdf=std::bind(testfun<cpp_dec_float_50>,
                                 std::placeholders::_1,a_cdf);
    test_iac<funct_cdf50,cpp_dec_float_50,
	     inte_gauss56_cern<funct_cdf50,cpp_dec_float_50>,10000>
      (t,tf_cdf,1.0e-30,"iac, cpp_dec_float_50, testfun",diff_cdf);
    t.test_abs_boost<cpp_dec_float_50>(diff_cdf,0.0,1.0e-29,
				       "inte_adapt_cern_cdf");
    
#ifdef O2SCL_MPFR
    
    cout << "inte_adapt_cern, mpfr_float_50, testfun:\n  ";
    
    mpfr_float_50 one_mpfr=1.0, diff_mpfr;
    mpfr_float_50 hundred_mpfr=100.0;
    mpfr_float_50 a_mpfr=one_mpfr/hundred_mpfr;
    funct_mpfr50 tf_mpfr=std::bind(testfun<mpfr_float_50>,
                                 std::placeholders::_1,a_mpfr);
    test_iac<funct_mpfr50,mpfr_float_50,
	     inte_gauss56_cern<funct_mpfr50,mpfr_float_50>,10000>
      (t,tf_mpfr,1.0e-30,"iac, mpfr_float_50, testfun",diff_mpfr);
    t.test_abs_boost<mpfr_float_50>(diff_mpfr,0.0,1.0e-29,
                                    "inte_adapt_cern_mpfr");

#endif
    
  }

  if (true) {

    cout << "inte_transform with inte_adapt_cern, double precision, "
         << "sin_recip:\n  " << endl;
    
    double calc, ei, diff;

    //inte_transform<funct,inte_adapt_cern<>,double> it_iac;
    inte_adapt_cern imac;
    double exact=1.0-cos(100.0/101.0);
    funct tf2=std::bind(sin_recip<double>,std::placeholders::_1);
    //it_iac.integ_il_err(tf2,-1.0,calc,ei);
    imac.integ_il_err(tf2,-1.0,calc,ei);
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << " " << ei << endl;
    t.test_rel<double>(calc,exact,1.0e-12,"it_iac double");
    cout << endl;
  
    cout << "inte_transform with inte_adapt_cern, long double precision, "
         << "sin_recip:\n  " << endl;

    /*
    inte_transform<funct_ld,inte_adapt_cern
		   <funct_ld,inte_gauss56_cern
		    <funct_ld,long double>,100,
		    long double>,long double> it_iac_ld;
    */
    //inte_transform<funct_ld,inte_adapt_cern,long double> it_imac_ld;

    long double exact_ld=1.0L-cos(100.0L/101.0L);
    funct_ld tf2_ld=std::bind(sin_recip<long double>,
                              std::placeholders::_1);
    long double calc_ld, ei_ld;
    //it_iac_ld.def_inte.tol_rel=1.0e-15;
    //it_iac_ld.def_inte.tol_abs=1.0e-15;
    imac.tol_rel=1.0e-15;
    imac.tol_abs=1.0e-15;
    //it_imac_ld.def_inte.tol_abs=1.0e-15;
    //it_iac_ld.integ_il_err(tf2_ld,-1.0L,calc_ld,ei_ld);
    imac.integ_il_err(tf2_ld,-1.0L,calc_ld,ei_ld);
    long double diff_ld=fabs(calc_ld-exact_ld);
    cout << calc_ld << " " << exact_ld << " "
	 << diff_ld << " " << ei_ld << endl;
    t.test_rel<double>(calc_ld,exact_ld,1.0e-15,"it_iac double");
    cout << endl;
    
    cout << "inte_transform with inte_adapt_cern, cpp_dec_float_50 precision, "
         << "sin_recip:\n  " << endl;

    /*
    inte_transform<funct_cdf50,inte_adapt_cern
		   <funct_cdf50,inte_gauss56_cern
		    <funct_cdf50,cpp_dec_float_50>,100,
		    cpp_dec_float_50>,cpp_dec_float_50> it_iac_cdf;
    */
    
    cpp_dec_float_50 one=1.0;
    cpp_dec_float_50 mone=-1.0;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 exact_cdf=one-cos(hundred/(hundred+one));
    funct_cdf50 tf2_cdf=std::bind(sin_recip<cpp_dec_float_50>,
                                  std::placeholders::_1);
    cpp_dec_float_50 calc_cdf, ei_cdf;
    //it_iac_cdf.def_inte.tol_rel=1.0e-30;
    //it_iac_cdf.def_inte.tol_abs=1.0e-30;
    imac.tol_rel=1.0e-30;
    imac.tol_abs=1.0e-30;
    //it_iac_cdf.integ_il_err(tf2_cdf,-one,calc_cdf,ei_cdf);
    imac.integ_il_err(tf2_cdf,mone,calc_cdf,ei_cdf);
    cpp_dec_float_50 diff_cdf=fabs(calc_cdf-exact_cdf);
    cout << calc_cdf << " " << exact_cdf << " "
	 << diff_cdf << " " << ei_cdf << endl;
    t.test_rel_boost<cpp_dec_float_50>(calc_cdf,exact_cdf,1.0e-29,
				       "it_iac cpp_dec_float_50");
    cout << endl;

#ifdef O2SCL_MPFR
    
    cout << "inte_transform with inte_adapt_cern, mpfr_float_50 precision, "
         << "sin_recip:\n  " << endl;

    /*
    inte_transform<funct_mpfr50,inte_adapt_cern
		   <funct_mpfr50,inte_gauss56_cern
		    <funct_mpfr50,mpfr_float_50>,100,
		    mpfr_float_50>,mpfr_float_50> it_iac_mp50;
    */
    
    mpfr_float_50 one_mp50=1.0;
    mpfr_float_50 mone_mp50=-1.0;
    mpfr_float_50 hundred_mp50=100.0;
    mpfr_float_50 exact_mp50=one_mp50-cos(hundred_mp50/(hundred_mp50+one_mp50));
    funct_mpfr50 tf2_mp50=std::bind(sin_recip<mpfr_float_50>,
                                  std::placeholders::_1);
    mpfr_float_50 calc_mp50, ei_mp50;
    //it_iac_mp50.def_inte.tol_rel=1.0e-30;
    //it_iac_mp50.def_inte.tol_abs=1.0e-30;
    imac.tol_rel=1.0e-30;
    imac.tol_abs=1.0e-30;
    //it_iac_mp50.integ_il_err(tf2_mp50,-one_mp50,calc_mp50,ei_mp50);
    imac.integ_il_err(tf2_mp50,mone_mp50,calc_mp50,ei_mp50);
    mpfr_float_50 diff_mp50=fabs(calc_mp50-exact_mp50);
    cout << calc_mp50 << " " << exact_mp50 << " "
	 << diff_mp50 << " " << ei_mp50 << endl;
    t.test_rel_boost<mpfr_float_50>(calc_mp50,exact_mp50,1.0e-29,
				       "it_iac mpfr_float_50");
    cout << endl;

#endif
    
  }

  {

    inte_adapt_cern imac;
    
    double ans, exact, err;
    
    funct tf=mp_test_func<double>;
    
    // Compare with the exact result
    imac.integ_err(tf,0.0,1.0,ans,err);
    exact=cos(100.0)-cos(1/1.01);
    std::cout << ans << " " << err << std::endl;
    t.test_rel(ans,exact,1.0e-8,"imac test");
  }

  {
    inte_adapt_cern imac;
    
    long double ans, exact, err;
    
    funct_ld tf=mp_test_func<long double>;
    
    // Compare with the exact result
    imac.tol_rel=1.0e-15;
    imac.set_nsub(1000);
    imac.integ_err(tf,0.0L,1.0L,ans,err);
    exact=cos(100.0L)-cos(1.0L/1.01L);
    std::cout << ans << " " << err << std::endl;
    t.test_rel<long double>(ans,exact,1.0e-15,"imac test");
  }

  {
    inte_adapt_cern imac;
    
    cpp_dec_float_50 ans, exact, err;
    
    funct_cdf50 tf=mp_test_func<cpp_dec_float_50>;
    
    // Compare with the exact result
    imac.set_nsub(20000);
    imac.tol_rel=1.0e-31;
    cpp_dec_float_50 zero=0;
    cpp_dec_float_50 one=1;
    cpp_dec_float_50 hundred=100;
    imac.integ_err(tf,zero,one,ans,err);
    exact=cos(hundred)-cos(hundred/(hundred+one));
    std::cout << ans << " " << err << std::endl;
    t.test_rel_boost<cpp_dec_float_50>(ans,exact,1.0e-30,"imac test");
  }

  {
    // Multiprecision integration with infinite limits
    
    inte_adapt_cern imac;
    
    double val, err2, a=0, b=1;
    
    imac.integ_iu_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                             a,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi/2.0,1.0e-15,"multip 3");
    cout << endl;

    imac.integ_il_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                             a,val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi/2.0,1.0e-15,"multip 4");
    cout << endl;

    imac.integ_i_err_multip([](auto &&tb) mutable { return test_func_i(tb); },
                            val,err2);
    cout << dtos(val,0) << " " << dtos(err2,0) << endl;
    t.test_rel(val,root_pi,1.0e-15,"multip 5");
    cout << endl;

    
  }
  t.report();
  return 0;
}

