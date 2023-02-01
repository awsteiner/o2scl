/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/funct_multip.h>
#include <o2scl/inte_double_exp_boost.h>

using namespace std;
using namespace o2scl;

typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

template<class fp_t> fp_t test_func(fp_t x) {
  fp_t one=1;
  fp_t hundred=100;
  return -sin(one/(x+one/hundred))/(x+one/hundred)/(x+one/hundred);
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

  inte_tanh_sinh_boost<funct,61> itsb;
  inte_exp_sinh_boost<funct,61> iesb;
  inte_sinh_sinh_boost<funct,61> issb;

  double ans, exact, err;

  funct tf1=test_func<double>;
  funct tf2=test_func_2<double>;

  // Finite integral, moderately difficult integrand
  cout << "tanh-sinh, finite interval:" << endl;
  itsb.integ_err(tf1,0.0,1.0,ans,err);
  exact=exact_func<double>();
  std::cout << ans << " " << err << " " << exact << " "
	    << itsb.L1norm << std::endl;
  t.test_rel(ans,exact,1.0e-8,"tanh_sinh test");
  cout << endl;

  // Semi-infinite domains
  cout << "tanh-sinh, infinite upper limit:" << endl;
  exact=sqrt(acos(-1))/2.0;
  itsb.integ_err(tf2,0.0,std::numeric_limits<double>::infinity(),ans,err);
  t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 2");
  std::cout << ans << " " << err << " " << exact << " "
	    << itsb.L1norm << std::endl;
  cout << endl;

  cout << "tanh-sinh, infinite lower limit:" << endl;
  itsb.integ_err(tf2,-std::numeric_limits<double>::infinity(),0.0,ans,err);
  t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 3");
  std::cout << ans << " " << err << " " << exact << " "
	    << itsb.L1norm << std::endl;
  cout << endl;

  // Infinite domain
  cout << "tanh-sinh, infinite lower and upper limit:" << endl;
  exact=sqrt(acos(-1));
  itsb.integ_err(tf2,-std::numeric_limits<double>::infinity(),
		 std::numeric_limits<double>::infinity(),
		 ans,err);
  t.test_rel(ans,exact,1.0e-8,"tanh_sinh test 4");
  std::cout << ans << " " << err << " " << exact << " "
	    << itsb.L1norm << std::endl;
  cout << endl;

  // exp_sinh integration, native
  cout << "exp-sinh, infinite upper limit:" << endl;
  iesb.integ_iu_err(tf2,0.0,ans,err);
  exact=sqrt(acos(-1.0))/2.0;
  std::cout << ans << " " << err << " " << exact << " "
	    << iesb.L1norm << std::endl;
  t.test_rel(ans,exact,1.0e-8,"exp_sinh test 1");
  cout << endl;

  cout << "exp-sinh, infinite upper limit, explicit limits:" << endl;
  iesb.integ_err(tf2,0.0,std::numeric_limits<double>::infinity(),ans,err);
  exact=sqrt(acos(-1.0))/2.0;
  std::cout << ans << " " << err << " " << exact << " "
	    << iesb.L1norm << std::endl;
  t.test_rel(ans,exact,1.0e-8,"exp_sinh test 2");
  cout << endl;
  
  cout << "exp-sinh, infinite lower limit:" << endl;
  iesb.integ_err(tf2,-std::numeric_limits<double>::infinity(),0.0,ans,err);
  exact=sqrt(acos(-1.0))/2.0;
  std::cout << ans << " " << err << " " << exact << " "
	    << iesb.L1norm << std::endl;
  t.test_rel(ans,exact,1.0e-8,"exp_sinh test 3");
  cout << endl;
  
  cout << "sinh-sinh, infinite interval:" << endl;
  issb.integ_i_err(tf2,ans,err);
  exact=sqrt(acos(-1.0));
  std::cout << ans << " " << err << " " << exact << " "
	    << issb.L1norm << std::endl;
  t.test_rel(ans,exact,1.0e-8,"sinh_sinh test");
  cout << endl;

  funct_ld tf1_ld=test_func<long double>;
  funct_ld tf2_ld=test_func_2<long double>;
  funct_cdf50 tf1_cdf=test_func<cpp_dec_float_50>;
  funct_cdf50 tf2_cdf=test_func_2<cpp_dec_float_50>;
  funct_cdf35 tf1_cdf35=test_func<cpp_dec_float_35>;
  funct_cdf35 tf2_cdf35=test_func_2<cpp_dec_float_35>;

  inte_tanh_sinh_boost<funct_ld,61,long double> itsb_ld;
  inte_tanh_sinh_boost<funct_cdf35,61,cpp_dec_float_35> itsb_cdf35;
  inte_tanh_sinh_boost<funct_cdf50,61,cpp_dec_float_50> itsb_cdf;

  // Finite integral, moderately difficult integrand
  cout << "tanh-sinh, finite interval, long double:" << endl;
  long double one_ld=1;
  long double hundred_ld=100;
  long double exact_ld=exact_func<long double>();
  long double ans_ld, err_ld;
  itsb_ld.integ_err(tf1_ld,0.0,1.0,ans_ld,err_ld);
  t.test_rel<long double>(ans_ld,exact_ld,1.0e-16,"tanh_sinh test ld");
  std::cout << ans_ld << " " << err_ld << " " << exact_ld << " "
	    << itsb_ld.L1norm << std::endl;
  cout << endl;

  // Finite integral, moderately difficult integrand
  cout << "tanh-sinh, finite interval, cpp_dec_float_35:" << endl;
  cpp_dec_float_35 one_cdf35=1;
  cpp_dec_float_35 hundred_cdf35=100;
  cpp_dec_float_35 exact_cdf35=exact_func<cpp_dec_float_35>();
  cpp_dec_float_35 ans_cdf35, err_cdf35;
  itsb_cdf35.tol_rel=1.0e-30;
  itsb_cdf35.integ_err(tf1_cdf35,0.0,1.0,ans_cdf35,err_cdf35);
  std::cout << ans_cdf35 << " " << err_cdf35 << std::endl;
  std::cout << dtos(exact_cdf35,0) << std::endl;
  std::cout << itsb_cdf35.L1norm << std::endl;
  t.test_rel_boost<cpp_dec_float_35>(ans_cdf35,exact_cdf35,1.0e-29,
				     "tanh_sinh test cdf35");
  cout << endl;

  // Finite integral, moderately difficult integrand
  cout << "tanh-sinh, finite interval, cpp_dec_float_50:" << endl;
  cpp_dec_float_50 one_cdf=1;
  cpp_dec_float_50 hundred_cdf=100;
  cpp_dec_float_50 exact_cdf=exact_func<cpp_dec_float_50>();
  cpp_dec_float_50 ans_cdf, err_cdf;
  itsb_cdf.tol_rel=1.0e-40;
  itsb_cdf.integ_err(tf1_cdf,0.0,1.0,ans_cdf,err_cdf);
  std::cout << ans_cdf << " " << err_cdf << std::endl;
  std::cout << dtos(exact_cdf,0) << std::endl;
  std::cout << itsb_cdf.L1norm << std::endl;
  t.test_rel_boost<cpp_dec_float_50>(ans_cdf,exact_cdf,1.0e-39,
				     "tanh_sinh test cdf");
  cout << endl;

#ifdef O2SCL_OSX
  {
    double val, err2, a=0, b=1;
    double exact=cos(100.0)-cos(1/1.01);
    inte_multip_tanh_sinh_boost<> imtsb;
    imtsb.verbose=2;
    imtsb.integ_err([](auto &&t) mutable { return test_func(t); },
                    a,b,val,err2,1.0e-8);
    imtsb.integ_err([](auto &&t) mutable { return test_func(t); },
                    a,b,val,err2);
    t.test_rel(val,exact,1.0e-15,"multip");
  }
#endif  
  
  t.report();
  return 0;
}

