/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

double testfun(double tx, double &a) {
  return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

long double testfun_ld(long double tx, long double &a) {
  return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

double sin_recip(double x) {
  return sin(1.0/(-x+0.01))*pow(-x+0.01,-2.0);
}

long double sin_recip_ld(long double x) {
  return sin(1.0L/(-x+0.01L))*pow(-x+0.01L,-2.0L);
}

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 tx, cpp_dec_float_50 &a) {
  cpp_dec_float_50 one=1.0;
  return -cos(one/(tx+a))/(a+tx)/(a+tx);
}

cpp_dec_float_50 sin_recip_cdf(cpp_dec_float_50 x) {
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 hundred=100;
  return sin(one/(-x+one/hundred))*pow(-x+one/hundred,-one-one);
}

#ifdef O2SCL_MPFR
typedef boost::multiprecision::mpfr_float_50 mpfr_float_50;

mpfr_float_50 testfun_mpfr(mpfr_float_50 tx, mpfr_float_50 &a) {
  mpfr_float_50 one=1.0;
  return -cos(one/(tx+a))/(a+tx)/(a+tx);
}

mpfr_float_50 sin_recip_mpfr(mpfr_float_50 x) {
  mpfr_float_50 one=1;
  mpfr_float_50 hundred=100;
  return sin(one/(-x+one/hundred))*pow(-x+one/hundred,-one-one);
}

#endif

template<class func_t=funct, class fp_t=double,
	 class def_inte_t, size_t nsub>
void test_iac(test_mgr &t, func_t &f, fp_t acc,
	      std::string comment, fp_t &diff, bool output_sub=false) {

  cout << comment << ":\n  ";
  fp_t one=1;
  fp_t ten=10;
  fp_t hundred=100;
  fp_t a=one/hundred, calc, ei;
  fp_t exact=sin(one/(one+a))-sin(one/a);
  
  inte_adapt_cern<func_t,def_inte_t,nsub,fp_t> iac;
  iac.tol_rel=acc;
  iac.tol_abs=acc;
  iac.integ_err(f,0.0,one,calc,ei);
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << " "
       << ei << endl;
  cout << "  subdivisions: ";
  cout << iac.get_nsubdivisions() << endl;
  cout << endl;

  if (output_sub) {
    size_t n=iac.get_nsubdivisions();
    typedef boost::numeric::ublas::vector<fp_t> ubvector;
    ubvector xlo(n), xhi(n), val(n), err(n);
    iac.get_subdivisions(xlo,xhi,val,err);
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
    funct tf=std::bind(testfun,std::placeholders::_1,a);
    test_iac<funct,double,
	     inte_gauss56_cern<funct,double,
			       inte_gauss56_coeffs_double>,100>
      (t,tf,1.0e-8,"iac, double, testfun",diff,true);
    t.test_abs<double>(diff,0.0,1.0e-7,"inte_adapt_cern");

    cout << "inte_adapt_cern, long double, testfun:\n  ";
    
    long double a_ld=0.01L, diff_ld;
    funct_ld tf_ld=std::bind(testfun_ld,std::placeholders::_1,a_ld);
    test_iac<funct_ld,long double,
	     inte_gauss56_cern<funct_ld,long double,
			       inte_gauss56_coeffs_long_double>,1000>
      (t,tf_ld,1.0e-15,"iac, long double, testfun",diff_ld);
    t.test_abs<long double>(diff_ld,0.0,1.0e-14,"inte_adapt_cern_ld");

    cout << "inte_adapt_cern, cpp_dec_float_50, testfun:\n  ";
    
    cpp_dec_float_50 one=1.0, diff_cdf;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 a_cdf=one/hundred;
    funct_cdf50 tf_cdf=std::bind(testfun_cdf,std::placeholders::_1,a_cdf);
    test_iac<funct_cdf50,cpp_dec_float_50,
	     inte_gauss56_cern<funct_cdf50,cpp_dec_float_50,
			       inte_gauss56_coeffs_float_50<cpp_dec_float_50>
                               >,10000>
      (t,tf_cdf,1.0e-30,"iac, cpp_dec_float_50, testfun",diff_cdf);
    t.test_abs_boost<cpp_dec_float_50>(diff_cdf,0.0,1.0e-29,
				       "inte_adapt_cern_cdf");

#ifdef O2SCL_MPFR
    
    cout << "inte_adapt_cern, mpfr_float_50, testfun:\n  ";
    
    mpfr_float_50 one_mpfr=1.0, diff_mpfr;
    mpfr_float_50 hundred_mpfr=100.0;
    mpfr_float_50 a_mpfr=one_mpfr/hundred_mpfr;
    funct_mp50 tf_mpfr=std::bind(testfun_mpfr,std::placeholders::_1,a_mpfr);
    test_iac<funct_mp50,mpfr_float_50,
	     inte_gauss56_cern<funct_mp50,mpfr_float_50,
			       inte_gauss56_coeffs_float_50<mpfr_float_50>
                               >,10000>
      (t,tf_mpfr,1.0e-30,"iac, mpfr_float_50, testfun",diff_mpfr);
    t.test_abs_boost<mpfr_float_50>(diff_mpfr,0.0,1.0e-29,
                                    "inte_adapt_cern_mpfr");

#endif
    
  }

  if (true) {

    cout << "inte_transform with inte_adapt_cern, double precision, "
         << "sin_recip:\n  " << endl;
    
    double calc, ei, diff;

    inte_transform<funct,inte_adapt_cern<>,double> it_iac;
    double exact=1.0-cos(100.0/101.0);
    funct tf2=std::bind(sin_recip,std::placeholders::_1);
    it_iac.integ_il_err(tf2,-1.0,calc,ei);
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << " " << ei << endl;
    t.test_rel<double>(calc,exact,1.0e-12,"it_iac double");
    cout << endl;
  
    cout << "inte_transform with inte_adapt_cern, long double precision, "
         << "sin_recip:\n  " << endl;
    
    inte_transform<funct_ld,inte_adapt_cern
		   <funct_ld,inte_gauss56_cern
		    <funct_ld,long double,
		     inte_gauss56_coeffs_long_double>,100,
		    long double>,long double> it_iac_ld;

    long double exact_ld=1.0L-cos(100.0L/101.0L);
    funct_ld tf2_ld=std::bind(sin_recip_ld,std::placeholders::_1);
    long double calc_ld, ei_ld;
    it_iac_ld.def_inte.tol_rel=1.0e-15;
    it_iac_ld.def_inte.tol_abs=1.0e-15;
    it_iac_ld.integ_il_err(tf2_ld,-1.0,calc_ld,ei_ld);
    long double diff_ld=fabs(calc_ld-exact_ld);
    cout << calc_ld << " " << exact_ld << " "
	 << diff_ld << " " << ei_ld << endl;
    t.test_rel<double>(calc_ld,exact_ld,1.0e-15,"it_iac double");
    cout << endl;
    
    cout << "inte_transform with inte_adapt_cern, cpp_dec_float_50 precision, "
         << "sin_recip:\n  " << endl;
    
    inte_transform<funct_cdf50,inte_adapt_cern
		   <funct_cdf50,inte_gauss56_cern
		    <funct_cdf50,cpp_dec_float_50,
		     inte_gauss56_coeffs_float_50<cpp_dec_float_50>>,100,
		    cpp_dec_float_50>,cpp_dec_float_50> it_iac_cdf;
    
    cpp_dec_float_50 one=1.0;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 exact_cdf=one-cos(hundred/(hundred+one));
    funct_cdf50 tf2_cdf=std::bind(sin_recip_cdf,std::placeholders::_1);
    cpp_dec_float_50 calc_cdf, ei_cdf;
    it_iac_cdf.def_inte.tol_rel=1.0e-30;
    it_iac_cdf.def_inte.tol_abs=1.0e-30;
    it_iac_cdf.integ_il_err(tf2_cdf,-one,calc_cdf,ei_cdf);
    cpp_dec_float_50 diff_cdf=fabs(calc_cdf-exact_cdf);
    cout << calc_cdf << " " << exact_cdf << " "
	 << diff_cdf << " " << ei_cdf << endl;
    t.test_rel_boost<cpp_dec_float_50>(calc_cdf,exact_cdf,1.0e-29,
				       "it_iac cpp_dec_float_50");
    cout << endl;

#ifdef O2SCL_MPFR
    
    cout << "inte_transform with inte_adapt_cern, mpfr_float_50 precision, "
         << "sin_recip:\n  " << endl;
    
    inte_transform<funct_mp50,inte_adapt_cern
		   <funct_mp50,inte_gauss56_cern
		    <funct_mp50,mpfr_float_50,
		     inte_gauss56_coeffs_float_50<mpfr_float_50>>,100,
		    mpfr_float_50>,mpfr_float_50> it_iac_mp50;
    
    mpfr_float_50 one_mp50=1.0;
    mpfr_float_50 hundred_mp50=100.0;
    mpfr_float_50 exact_mp50=one_mp50-cos(hundred_mp50/(hundred_mp50+one_mp50));
    funct_mp50 tf2_mp50=std::bind(sin_recip_mpfr,std::placeholders::_1);
    mpfr_float_50 calc_mp50, ei_mp50;
    it_iac_mp50.def_inte.tol_rel=1.0e-30;
    it_iac_mp50.def_inte.tol_abs=1.0e-30;
    it_iac_mp50.integ_il_err(tf2_mp50,-one_mp50,calc_mp50,ei_mp50);
    mpfr_float_50 diff_mp50=fabs(calc_mp50-exact_mp50);
    cout << calc_mp50 << " " << exact_mp50 << " "
	 << diff_mp50 << " " << ei_mp50 << endl;
    t.test_rel_boost<mpfr_float_50>(calc_mp50,exact_mp50,1.0e-29,
				       "it_iac mpfr_float_50");
    cout << endl;

#endif
    
  }
  
  t.report();
  return 0;
}

