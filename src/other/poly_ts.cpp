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
#include <string>
#include <ctime>
#include <o2scl/test_mgr.h>
#include <o2scl/misc.h>
#include <o2scl/poly.h>

using namespace std;
using namespace o2scl;

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::cpp_bin_float_50 cpp_bin_float_50;
typedef boost::multiprecision::cpp_complex_50 cpp_complex_50;

test_mgr tst;
int wid=21;

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_real_coeff
(quadratic_real_coeff<fp_t,cx_t> *po, 
 string str, fp_t alpha, fp_t e1, fp_t e2, 
 fp_t e3, fp_t e4, bool check_ret) {

  fp_t s1=0, s2=0, m1=0, m2=0;
  size_t wrong_ret=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count+=po->test_coeffs_zero_disc(alpha,s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  count+=po->test_real_roots(s1,s2,m1,m2,wrong_ret);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quadratic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quadratic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quadratic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quadratic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time;
  if (check_ret) {
    cout << " " << wrong_ret << endl;
  } else {
    cout << endl;
  }

  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_real_coeff_boost
(quadratic_real_coeff<fp_t,cx_t> *po, 
 string str, fp_t alpha, fp_t e1, fp_t e2, 
 fp_t e3, fp_t e4, bool check_ret) {

  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0, wrong_ret=0;
  
  count+=po->test_coeffs_zero_disc(alpha,s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  count+=po->test_real_roots(s1,s2,m1,m2,wrong_ret);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quadratic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quadratic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quadratic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quadratic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time;
  if (check_ret) {
    cout << " " << wrong_ret << endl;
  } else {
    cout << endl;
  }

  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_complex(quadratic_complex<fp_t,cx_t> *po,
                            string str, fp_t e1, fp_t e2, fp_t e3, fp_t e4) {
                            
  fp_t s1=0, s2=0, m1=0, m2=0;
  
  clock_t lt1=clock();
  
  size_t count=0;
  
  count+=po->test_complex_coeffs(s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  tst.test_abs<fp_t>(s1,0.0,e1,"quadratic_complex s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quadratic_complex s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quadratic_complex m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quadratic_complex m2");
  cout.width(wid);
  cout << str.c_str();
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_complex_boost(quadratic_complex<fp_t,cx_t> *po,
                                  string str, 
                                  double e1, double e2, double e3, double e4) {
  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count+=po->test_complex_coeffs(s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quadratic_complex s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quadratic_complex s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quadratic_complex m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quadratic_complex m2");
  cout.width(wid);
  cout << str.c_str();
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_real_coeff(cubic_real_coeff<fp_t,cx_t> *po,
                           string str, fp_t alpha, fp_t e1, fp_t e2,
                           fp_t e3, fp_t e4, bool check_ret) {
			   
  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0, wrong_ret=0;
  
  count=po->test_real_coeffs(alpha,s1,s2,m1,m2,wrong_ret);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"cubic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"cubic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"cubic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"cubic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time;
  if (check_ret) {
    cout << " " << wrong_ret << endl;
  } else {
    cout << endl;
  }
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_real_coeff_boost(cubic_real_coeff<fp_t,cx_t> *po,
                                 string str, fp_t alpha, fp_t e1,
                                 fp_t e2, fp_t e3, 
                                 fp_t e4, bool check_ret) {

  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0, wrong_ret=0;
  
  count=po->test_real_coeffs(alpha,s1,s2,m1,m2,wrong_ret);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"cubic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"cubic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"cubic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"cubic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time;
  if (check_ret) {
    cout << " " << wrong_ret << endl;
  } else {
    cout << endl;
  }
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_complex(cubic_complex<fp_t,cx_t> *po,
                        string str, fp_t e1, 
			fp_t e2, fp_t e3, fp_t e4) {

  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_complex_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"cubic_complex s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"cubic_complex s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"cubic_complex m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"cubic_complex m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_complex_boost(cubic_complex<fp_t,cx_t> *po,
                              string str, fp_t e1, 
                              fp_t e2, fp_t e3, fp_t e4) {

  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_complex_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;

  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"cubic_complex s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"cubic_complex s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"cubic_complex m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"cubic_complex m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double>
void test_quartic_real(quartic_real<fp_t> *po,
                       string str, fp_t alpha, 
		       fp_t e1, fp_t e2, fp_t e3, fp_t e4) {

  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_real_roots(alpha,s1,m1);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quartic_real s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quartic_real s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quartic_real m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quartic_real m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double>
void test_quartic_real_boost(quartic_real<fp_t> *po,
                             string str, fp_t alpha, 
                             fp_t e1, fp_t e2, fp_t e3, fp_t e4) {

  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_real_roots(alpha,s1,m1);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quartic_real s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quartic_real s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quartic_real m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quartic_real m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_real_coeff(quartic_real_coeff<fp_t,cx_t> *po,
                             string str, fp_t e1, fp_t e2, fp_t e3,
                             fp_t e4) {
  
  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_real_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_real_coeff_boost(quartic_real_coeff<fp_t,cx_t> *po,
                                   string str, fp_t e1, fp_t e2, fp_t e3,
                                   fp_t e4) {
  
  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_real_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_complex(quartic_complex<fp_t,cx_t> *po,
                          string str,
			  fp_t e1, fp_t e2, fp_t e3, fp_t e4) {
  
  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_complex_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_complex_boost(quartic_complex<fp_t,cx_t> *po,
                                string str,
                                fp_t e1, fp_t e2, fp_t e3, fp_t e4) {

  fp_t s1=0,s2=0,m1=0,m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count=po->test_complex_coeffs(s1,s2,m1,m2);
  
  clock_t lt2=clock();
  
  s1/=count;
  s2/=count;
  
  double time=((double)(lt2-lt1))/CLOCKS_PER_SEC;

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << time << endl;
  return;
}

int main(void) {
  tst.set_output_level(1);
  
  cout.setf(ios::left | ios::scientific);
  cout.precision(4);

  // Generic polynomial solvers
  poly_real_coeff_gsl<> p3;
  
  // Quadratic solvers
  quadratic_real_coeff_gsl t3;
  quadratic_real_coeff_gsl2<> t1;
  quadratic_complex_std<> t2;
  
  // Cubic solvers
  cubic_real_coeff_cern<> c1;
  cubic_real_coeff_gsl c4;
  cubic_real_coeff_gsl2<> c2;
  cubic_complex_std<> c3;
  cubic_real_coeff_multip c5;
  //c5.verbose=1;
  
  // Quartic solvers
  quartic_real_coeff_cern<> q1;
  quartic_real_std<> q4;
  quartic_complex_std<> q5;

  quadratic_real_coeff_gsl2<long double,std::complex<long double> > t1_ld;
  quadratic_real_coeff_gsl2<cpp_bin_float_50,cpp_complex_50> t1_cdf50;

  quadratic_complex_std<long double,std::complex<long double> > t2_ld;
  quadratic_complex_std<cpp_bin_float_50,cpp_complex_50> t2_cdf50;
  
  cubic_real_coeff_cern<long double,std::complex<long double> > c1_ld;
  cubic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> c1_cdf50;
  c1_ld.eps=1.0e-7;
  c1_ld.delta=1.0e-18;
  c1_cdf50.eps=1.0e-20;
  c1_cdf50.delta=1.0e-40;
  
  cubic_real_coeff_gsl2<long double,std::complex<long double> > c2_ld;
  cubic_real_coeff_gsl2<cpp_bin_float_50,cpp_complex_50> c2_cdf50;

  cubic_complex_std<long double,std::complex<long double> > c3_ld;
  cubic_complex_std<cpp_bin_float_50,cpp_complex_50> c3_cdf50;
  
  quartic_real_coeff_cern<long double,std::complex<long double> > q1_ld;
  quartic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> q1_cdf50;
  q1_ld.cub_obj.eps=1.0e-7;
  q1_ld.cub_obj.delta=1.0e-18;
  q1_cdf50.cub_obj.eps=1.0e-20;
  q1_cdf50.cub_obj.delta=1.0e-40;

  quartic_real_std<long double> q4_ld;
  quartic_real_std<cpp_bin_float_50> q4_cdf50;

  quartic_complex_std<long double,std::complex<long double> > q5_ld;
  quartic_complex_std<cpp_bin_float_50,cpp_complex_50> q5_cdf50;

  cout << "Quadratics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11,false);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11,true);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11,true);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11,false);
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0,
                                         1.0e-16,1.0e-16,1.0e-13,1.0e-14,true);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0,
                                         1.0e-16,1.0e-16,1.0e-13,1.0e-14,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl2_50",1.0,
     1.0e-47,1.0e-48,1.0e-45,1.0e-45,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",1.0,
     1.0e-47,1.0e-48,1.0e-45,1.0e-45,true);
  cout << endl;
    
  cout << "Quadratics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11,false);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11,true);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,5.0e-11,true);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11,false);
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0e-5,
                                         1.0e-17,1.0e-17,1.0e-13,1.0e-14,true);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0e-5,
                                         1.0e-17,1.0e-17,1.0e-13,1.0e-14,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl2_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-48,1.0e-48,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-48,1.0e-48,true);
  cout << endl;

  cout << "Quadratic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_complex(&t2,"quad_c_std",
			 1.0e-13,1.0e-9,1.0e-11,1.0e-7);
  test_quadratic_complex<long double>(&t2_ld,"quad_c_std_ld",
                                      1.0e-17,1.0e-13,1.0e-14,1.0e-10);
  test_quadratic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",
     1.0e-48,1.0e-44,1.0e-45,1.0e-42);
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(&c1,"cubic_rc_cern",1.0,
			1.0e0,1.0e6,1.0e1,1.0e7,true);
  test_cubic_real_coeff(&c4,"cubic_rc_gsl",1.0,
                        1.0e-1,1.0e-2,1.0e1,8.0e0,false);
  test_cubic_real_coeff(&c2,"cubic_rc_gsl2",1.0,
			1.0e-1,1.0e-2,1.0e1,8.0e0,true);
  test_cubic_real_coeff(&c3,"cubic_c_std",1.0,
			1.0e-1,1.0e-1,1.0e1,1.0e1,true);
  test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0,
			1.0e-1,4.0e-2,1.0e1,1.0e1,false);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c1_ld,"cubic_rc_cern_ld",
     1.0,5.0e-11,5.0e-11,5.0e-7,1.0e-10,true);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c2_ld,"cubic_rc_gsl2_ld",
     1.0,1.0e-12,1.0e-12,1.0e-9,1.0e-10,true);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0,1.0e-1,1.0e-1,1.0e1,1.0e1,true);
  test_cubic_real_coeff(&c5,"cubic_rc_mp",1.0,
                        1.0e-1,1.0e-2,1.0e1,8.0e0,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c1_cdf50,"cubic_rc_cern_50",
     1.0,1.0e-40,1.0e-40,1.0e-38,1.0e-37,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c2_cdf50,"cubic_rc_gsl2_50",
     1.0,1.0e-10,1.0e-10,1.0e-10,1.0e-10,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0,1.0e-1,1.0e-1,1.0e1,1.0e1,true);
  cout << endl;

  cout << "Cubics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(&c1,"cubic_rc_cern",1.0e-3,
			1.0e-5,1.0e-3,1.0e-3,1.0e-2,true);
  test_cubic_real_coeff(&c4,"cubic_rc_gsl",1.0e-3,
			1.0e-4,5.0e-4,1.0e-3,1.0e-3,false);
  test_cubic_real_coeff(&c2,"cubic_rc_gsl2",1.0e-3,
			1.0e-4,5.0e-4,1.0e-3,1.0e-3,true);
  test_cubic_real_coeff(&c3,"cubic_c_std",1.0e-3,
			1.0e-1,1.0e+2,4.0e0,1.0e+5,true);
  test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0e-3,
			1.0e-4,1.0e-3,1.0e-2,1.0e-2,false);
  //p3.check_refine=true;
  //test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0e-3,
  //1.0e-4,1.0e-3,1.0e-2,1.0e-2);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c1_ld,"cubic_rc_cern_ld",
     1.0e-3,1.0e-7,1.0e-5,1.0e-6,1.0e-4,true);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c2_ld,"cubic_rc_gsl2_ld",
     1.0e-3,1.0e-7,1.0e-5,1.0e-6,1.0e-4,true);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0e-3,1.0e1,1.0e-1,1.0e1,1.0e1,true);
  test_cubic_real_coeff(&c5,"cubic_rc_mp",1.0e-3,
                        1.0e-1,1.0e-2,1.0e1,8.0e0,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c1_cdf50,"cubic_rc_cern_50",
     1.0e-3,1.0e-38,1.0e-36,1.0e-37,1.0e-35,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c2_cdf50,"cubic_rc_gsl2_50",
     1.0e-3,1.0e-10,1.0e-10,1.0e-10,1.0e-10,true);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0e-3,1.0e-1,1.0e-1,1.0e1,1.0e1,true);
  cout << endl;
  
  cout << "Cubic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_complex(&c3,"cubic_c_std",
                     1.0e1,1.0e1,1.0e1,1.0e1);
  test_cubic_complex<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0e1,1.0e1,1.0e1,1.0e1);
  test_cubic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0e1,1.0e1,1.0e1,1.0e1);
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n";
  cout << " leading coefficient 1:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(&q1,"quartic_rc_cern",1.0,
		    5.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real(&q4,"quartic_r_std",1.0,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(&q5,"quartic_c_std",1.0,
		    1.0e-1,1.0,1.0e2,1.0);
  test_quartic_real(&p3,"poly_rc_gsl",1.0,
		    1.0e-14,1.0,1.0e-10,1.0);
  test_quartic_real<long double>
    (&q1_ld,"quartic_rc_cern_ld",1.0,
     6.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real<long double>
    (&q4_ld,"quartic_r_std_ld",1.0,
     5.0e-10,1.0,2.0e-6,1.0);
  test_quartic_real<long double>
    (&q5_ld,"quartic_c_std_ld",1.0,
     5.0e-12,1.0,2.0e-6,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q1_cdf50,"quartic_rc_cern_50",1.0,
     5.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q4_cdf50,"quartic_r_std_50",1.0,
     5.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q5_cdf50,"quartic_c_std_50",1.0,
     5.0e-12,1.0,2.0e-8,1.0);
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n"
       << " leading coefficient 1, coefficients of odd powers small:" 
       << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(&q1,"quartic_rc_cern",1.0e-5,
		    1.0e-2,1.0,1.0e2,1.0);
  test_quartic_real(&q4,"quartic_r_std",1.0e-5,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(&q5,"quartic_c_std",1.0e-5,
		    1.0e-1,1.0,5.0,1.0);
  test_quartic_real(&p3,"poly_rc_gsl",1.0e-5,
		    1.0e-14,1.0,1.0e-10,1.0);
  test_quartic_real<long double>
    (&q1_ld,"quartic_rc_cern_ld",1.0e-5,
     1.0e-11,1.0,2.0e-6,1.0);
  test_quartic_real<long double>
    (&q4_ld,"quartic_r_std_ld",1.0e-5,
     5.0e-1,1.0,1.0e2,1.0);
  test_quartic_real<long double>
    (&q5_ld,"quartic_c_std_ld",1.0e-5,
     5.0e-7,1.0,1.0,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q1_cdf50,"quartic_rc_cern_50",1.0e-5,
     5.0e-12,1.0,2.0e-6,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q4_cdf50,"quartic_r_std_50",1.0e-5,
     1.0e0,1.0,1.0e2,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q5_cdf50,"quartic_c_std_50",1.0e-5,
     1.0e-7,1.0,1.0,1.0);
  cout << endl;

  cout << "Quartics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real_coeff(&q1,"quartic_rc_cern",
			  1.0e0,1.0e4,1.0e2,1.0e5);
  test_quartic_real_coeff(&q5,"quartic_c_std",
			  1.0e1,1.0e1,1.0e3,1.0e4);
  test_quartic_real_coeff(&p3,"poly_rc_gsl",
			  1.0e-13,1.0e-6,1.0e-10,1.0e-4);
  test_quartic_real_coeff<long double,std::complex<long double> >
    (&q1_ld,"quartic_rc_cern_ld",
     1.0e-1,1.0e-1,1.0e1,1.0e2);
  test_quartic_real_coeff<long double,std::complex<long double> >
    (&q5_ld,"quartic_c_std_ld",
     1.0e-1,1.0e-1,1.0e1,1.0e2);
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&q1_cdf50,"quartic_rc_cern_50",
     1.0e-1,1.0e1,1.0e2,1.0e2);
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&q5_cdf50,"quartic_c_std_50",
     1.0e-1,1.0e1,1.0e2,1.0e2);
  cout << endl;
  
  cout << "Quartics with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_complex(&q5,"quartic_c_std",
		       1.0e-2,1.0e-2,1.0e1,1.0e2);
  test_quartic_complex<long double,std::complex<long double> >
    (&q5_ld,"quartic_c_std_ld",
     1.0e-2,1.0e-2,1.0e1,1.0e2);
  test_quartic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&q5_cdf50,"quartic_c_std_50",
     1.0e-2,1.0e-2,1.0e1,1.0e2);
  cout << endl;
  
  tst.report();

  return 0;
}


