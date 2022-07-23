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
#ifndef O2SCL_OLD_BOOST
#include <boost/multiprecision/cpp_complex.hpp>
#endif

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::cpp_bin_float_50 cpp_bin_float_50;
#ifndef O2SCL_OLD_BOOST
typedef boost::multiprecision::cpp_complex_50 cpp_complex_50;
#endif

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

template<class fp_t=double, class cx_t=std::complex<fp_t>,
         class poly_t=cubic_real_coeff<fp_t,cx_t>>
size_t test_cubic_real_coeffs2(poly_t *po,
                               fp_t alpha, fp_t &s1, fp_t &s2,
                               fp_t &m1, fp_t &m2, size_t &wrong_ret,
                               size_t n=16) {
  
  size_t count=0;
  
  gen_test_number<> ga, gb, gc, gd;
  
  for(size_t j1=0;j1<n;j1++) {
    fp_t ca=ga.gen()*alpha;
    gb.reset();
    for(size_t j2=0;j2<n;j2++) {
      fp_t cb=gb.gen();
      gc.reset();
      for(size_t j3=0;j3<n;j3++) {
        fp_t cc=gc.gen()*alpha;
        gd.reset();
        for(size_t j4=0;j4<n;j4++) {
          fp_t cd=gd.gen();
          
          if (fabs(ca)>0.0) {
            
            fp_t cr1;
            cx_t cr2, cr3;
            int ret=po->solve_rc(ca,cb,cc,cd,cr1,cr2,cr3);
            fp_t disc=po->disc3_r(ca,cb,cc,cd);
            if ((disc>=0.0 && ret==1) ||
                (disc<0.0 && ret==3) ||
                (ret!=1 && ret!=3)) {
              wrong_ret++;
            }
            
            cx_t cbp=-(cr1+cr2+cr3)*ca;
            cx_t ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
            cx_t cdp=-(cr1*cr2*cr3)*ca;
            
            fp_t czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
            cx_t czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
            cx_t czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
            
            fp_t q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
                         fabs(cc-ccp.real())*fabs(cc-ccp.real())+
                         fabs(cd-cdp.real())*fabs(cd-cdp.real()));
            fp_t q2=sqrt(fabs(czo1)*fabs(czo1)+abs(czo2)*abs(czo2)+
                         abs(czo3)*abs(czo3));
            
            s1+=q1;
            if (q1>m1) m1=q1;
            s2+=q2;
            if (q2>m2) m2=q2;
            
            count++;
          }
        }
      }
    }
  }
  
  return count;
}

template<class fp_t=double, class cx_t=std::complex<fp_t>,
         class poly_t=cubic_real_coeff<fp_t,cx_t>>
void test_cubic_real_coeff(poly_t *po,
                           string str, fp_t alpha, fp_t e1, fp_t e2,
                           fp_t e3, fp_t e4, size_t max_ret) {
			   
  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0, wrong_ret=0;
  
  count=test_cubic_real_coeffs2(po,alpha,s1,s2,m1,m2,wrong_ret);
  
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
  if (max_ret>0) {
    cout << " " << wrong_ret << endl;
    tst.test_gen(wrong_ret<max_ret,"cubic_real_coeff int");
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
  
  count=po->test_cubic_real_coeffs(alpha,s1,s2,m1,m2,wrong_ret);
  
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
  
  count=po->test_quartic_real_coeffs(s1,s2,m1,m2);
  
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
  
  count=po->test_quartic_real_coeffs(s1,s2,m1,m2);
  
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
#ifndef O2SCL_OLD_BOOST
  quadratic_real_coeff_gsl2<cpp_bin_float_50,cpp_complex_50> t1_cdf50;
#endif
  
  quadratic_complex_std<long double,std::complex<long double> > t2_ld;
#ifndef O2SCL_OLD_BOOST
  quadratic_complex_std<cpp_bin_float_50,cpp_complex_50> t2_cdf50;
#endif
  
  cubic_real_coeff_cern<long double,std::complex<long double> > c1_ld;
  c1_ld.eps=1.0e-7;
  c1_ld.delta=1.0e-18;
#ifndef O2SCL_OLD_BOOST
  cubic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> c1_cdf50;
  c1_cdf50.eps=1.0e-20;
  c1_cdf50.delta=1.0e-40;
#endif
  
  cubic_real_coeff_gsl2<long double,std::complex<long double> > c2_ld;
#ifndef O2SCL_OLD_BOOST
  cubic_real_coeff_gsl2<cpp_bin_float_50,cpp_complex_50> c2_cdf50;
#endif

  cubic_complex_std<long double,std::complex<long double> > c3_ld;
#ifndef O2SCL_OLD_BOOST
  cubic_complex_std<cpp_bin_float_50,cpp_complex_50> c3_cdf50;
#endif
  
  quartic_real_coeff_cern<long double,std::complex<long double> > q1_ld;
  q1_ld.cub_obj.eps=1.0e-7;
  q1_ld.cub_obj.delta=1.0e-18;
#ifndef O2SCL_OLD_BOOST
  quartic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> q1_cdf50;
  q1_cdf50.cub_obj.eps=1.0e-20;
  q1_cdf50.cub_obj.delta=1.0e-40;
#endif

  quartic_real_std<long double> q4_ld;
#ifndef O2SCL_OLD_BOOST
  quartic_real_std<cpp_bin_float_50> q4_cdf50;
#endif

  quartic_complex_std<long double,std::complex<long double> > q5_ld;
#ifndef O2SCL_OLD_BOOST
  quartic_complex_std<cpp_bin_float_50,cpp_complex_50> q5_cdf50;
#endif

  cout << "Quadratics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0,
                            1.0e-14,1.0e-14,1.0e-12,1.0e-12,false);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0,
                            1.0e-14,1.0e-14,1.0e-12,1.0e-12,true);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0,
                            1.0e-14,1.0e-14,1.0e-11,1.0e-11,true);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0,
                            1.0e-14,1.0e-14,1.0e-12,1.0e-12,false);
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0,
                                         1.0e-18,1.0e-18,1.0e-16,1.0e-16,true);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0,
                                         1.0e-17,1.0e-17,1.0e-14,1.0e-14,true);
#ifndef O2SCL_OLD_BOOST
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl2_50",1.0,
     1.0e-49,1.0e-49,1.0e-47,1.0e-47,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",1.0,
     1.0e-48,1.0e-48,1.0e-46,1.0e-46,true);
#endif
  cout << endl;
    
  cout << "Quadratics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-13,1.0e-12,false);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-13,1.0e-12,true);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-13,5.0e-12,true);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-13,1.0e-12,false);
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0e-5,
                                         1.0e-18,1.0e-17,1.0e-16,1.0e-16,true);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0e-5,
                                         1.0e-18,1.0e-17,1.0e-16,1.0e-16,true);
#ifndef O2SCL_OLD_BOOST
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl2_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-47,1.0e-47,true);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-47,1.0e-47,true);
#endif
  cout << endl;

  cout << "Quadratic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_complex(&t2,"quad_c_std",
			 1.0e-14,1.0e-14,1.0e-13,1.0e-3);
  test_quadratic_complex<long double>(&t2_ld,"quad_c_std_ld",
                                      1.0e-17,1.0e-17,1.0e-17,1.0e-16);
#ifndef O2SCL_OLD_BOOST
  test_quadratic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",
     1.0e-49,1.0e-49,1.0e-48,1.0e-48);
#endif
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(&c1,"cubic_rc_cern",1.0,
			1.0e-12,1.0e-12,1.0e-7,1.0e-9,11);
  test_cubic_real_coeff(&c4,"cubic_rc_gsl",1.0,
                        1.0e-14,1.0e-14,1.0e-13,1.0e-12,0);
  test_cubic_real_coeff(&c2,"cubic_rc_gsl2",1.0,
			1.0e-14,1.0e-14,1.0e-13,1.0e-12,9);
  test_cubic_real_coeff(&c5,"cubic_rc_mp",1.0,
                        1.0e-14,1.0e-14,1.0e-14,1.0e-12,9);
  test_cubic_real_coeff(&c3,"cubic_c_std",1.0,
			1.0e-10,1.0e-9,1.0e-6,1.0e-6,12300);
  test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0,
			1.0e-13,1.0e-13,1.0e-13,1.0e-12,0);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c1_ld,"cubic_rc_cern_ld",
     1.0,1.0e-9,1.0e-15,1.0e-5,1.0e-12,9);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c2_ld,"cubic_rc_gsl2_ld",
     1.0,1.0e-17,1.0e-17,1.0e-16,1.0e-15,43);
#ifndef O2SCL_FAST_TEST
  test_cubic_real_coeff<long double,std::complex<long double>,
                        cubic_real_coeff_multip>
    (&c5,"cubic_rc_mp_ld",
     1.0,1.0e-19,1.0e-17,1.0e-16,1.0e-15,9);
#endif
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0,1.0e-13,1.0e-13,1.0e-9,1.0e-9,11910);
#ifndef O2SCL_OLD_BOOST
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c1_cdf50,"cubic_rc_cern_50",
     1.0,1.0e-47,1.0e-46,1.0e-43,1.0e-43,9);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c2_cdf50,"cubic_rc_gsl2_50",
     1.0,1.0e-48,1.0e-48,1.0e-47,1.0e-46,19);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0,1.0e-44,1.0e-44,1.0e-40,1.0e-40,12516);
#endif
  cout << endl;

  cout << "Cubics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(&c1,"cubic_rc_cern",1.0e-3,
			1.0e-9,1.0e-7,1.0e-7,1.0e-6,1);
  test_cubic_real_coeff(&c4,"cubic_rc_gsl",1.0e-3,
			1.0e-9,1.0e-8,1.0e-7,1.0e-6,0);
  test_cubic_real_coeff(&c2,"cubic_rc_gsl2",1.0e-3,
			1.0e-9,1.0e-8,1.0e-7,1.0e-6,33);
  test_cubic_real_coeff(&c5,"cubic_rc_mp",1.0e-3,
			1.0e-14,1.0e-8,1.0e-14,1.0e-6,1);
  test_cubic_real_coeff(&c3,"cubic_c_std",1.0e-3,
			1.0e-8,1.0e-8,1.0e-5,1.0e-5,30600);
  test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0e-3,
			1.0e-13,1.0e-8,1.0e-12,1.0e-6,0);
  //p3.check_refine=true;
  //test_cubic_real_coeff(&p3,"poly_rc_gsl",1.0e-3,
  //1.0e-4,1.0e-3,1.0e-2,1.0e-2);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c1_ld,"cubic_rc_cern_ld",
     1.0e-3,1.0e-12,1.0e-11,1.0e-10,1.0e-9,1);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c2_ld,"cubic_rc_gsl2_ld",
     1.0e-3,1.0e-12,1.0e-11,1.0e-10,1.0e-9,154);
#ifndef O2SCL_FAST_TEST
  test_cubic_real_coeff<long double,std::complex<long double>,
                        cubic_real_coeff_multip>
    (&c5,"cubic_rc_mp_ld",
     1.0e-3,1.0e-18,1.0e-11,1.0e-17,1.0e-10,1);
#endif
  test_cubic_real_coeff<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0e-3,1.0e-11,1.0e-11,1.0e-8,1.0e-8,30400);
#ifndef O2SCL_OLD_BOOST
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c1_cdf50,"cubic_rc_cern_50",
     1.0e-3,1.0e-43,1.0e-42,1.0e-41,1.0e-40,1);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c2_cdf50,"cubic_rc_gsl2_50",
     1.0e-3,1.0e-43,1.0e-42,1.0e-41,1.0e-41,122);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0e-3,1.0e-43,1.0e-42,1.0e-39,1.0e-39,30593);
#endif
  cout << endl;
  
  cout << "Cubic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_complex(&c3,"cubic_c_std",
                     1.0e-13,1.0e-13,1.0e-12,1.0e-11);
  test_cubic_complex<long double,std::complex<long double> >
    (&c3_ld,"cubic_c_std_ld",
     1.0e-17,1.0e-16,1.0e-15,1.0e-15);
#ifndef O2SCL_OLD_BOOST
  test_cubic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&c3_cdf50,"cubic_c_std_50",
     1.0e-48,1.0e-47,1.0e-46,1.0e-46);
#endif
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n";
  cout << " leading coefficient 1:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(&q1,"quartic_rc_cern",1.0,
		    1.0e-13,1.0,1.0e-11,1.0);
  test_quartic_real(&q4,"quartic_r_std",1.0,
		    1.0e-13,1.0,1.0e-11,1.0);
  test_quartic_real(&q5,"quartic_c_std",1.0,
		    1.0e-14,1.0,1.0e-12,1.0);
  test_quartic_real(&p3,"poly_rc_gsl",1.0,
		    1.0e-13,1.0,1.0e-11,1.0);
  test_quartic_real<long double>
    (&q1_ld,"quartic_rc_cern_ld",1.0,
     1.0e-17,1.0,1.0e-15,1.0);
  test_quartic_real<long double>
    (&q4_ld,"quartic_r_std_ld",1.0,
     1.0e-9,1.0,1.0e-6,1.0);
  test_quartic_real<long double>
    (&q5_ld,"quartic_c_std_ld",1.0,
     1.0e-17,1.0,1.0e-15,1.0);
#ifndef O2SCL_OLD_BOOST
  test_quartic_real_boost<cpp_bin_float_50>
    (&q1_cdf50,"quartic_rc_cern_50",1.0,
     1.0e-48,1.0,1.0e-46,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q4_cdf50,"quartic_r_std_50",1.0,
     1.0e-25,1.0,1.0e-22,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q5_cdf50,"quartic_c_std_50",1.0,
     1.0e-49,1.0,1.0e-47,1.0);
#endif
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n"
       << " leading coefficient 1, coefficients of odd powers small:" 
       << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(&q1,"quartic_rc_cern",1.0e-5,
		    1.0e-1,1.0,1.0e2,1.0);
  test_quartic_real(&q4,"quartic_r_std",1.0e-5,
		    1.0,1.0,1.0e3,1.0);
  test_quartic_real(&q5,"quartic_c_std",1.0e-5,
		    1.0e-13,1.0,1.0e-12,1.0);
  test_quartic_real(&p3,"poly_rc_gsl",1.0e-5,
		    1.0e-13,1.0,1.0e-12,1.0);
  test_quartic_real<long double>
    (&q1_ld,"quartic_rc_cern_ld",1.0e-5,
     1.0e-10,1.0,1.0e-8,1.0);
  test_quartic_real<long double>
    (&q4_ld,"quartic_r_std_ld",1.0e-5,
     1.0,1.0,1.0e3,1.0);
  test_quartic_real<long double>
    (&q5_ld,"quartic_c_std_ld",1.0e-5,
     1.0e-17,1.0,1.0e-16,1.0);
#ifndef O2SCL_OLD_BOOST
  test_quartic_real_boost<cpp_bin_float_50>
    (&q1_cdf50,"quartic_rc_cern_50",1.0e-5,
     1.0e-39,1.0,1.0e-37,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q4_cdf50,"quartic_r_std_50",1.0e-5,
     1.0e1,1.0,1.0e4,1.0);
  test_quartic_real_boost<cpp_bin_float_50>
    (&q5_cdf50,"quartic_c_std_50",1.0e-5,
     1.0e-48,1.0,1.0e-47,1.0);
#endif
  cout << endl;

  cout << "Quartics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real_coeff(&q1,"quartic_rc_cern",
			  1.0e-1,1.0e-1,1.0e2,1.0e3);
  test_quartic_real_coeff(&q5,"quartic_c_std",
			  1.0e-3,1.0e-3,1.0e1,1.0e1);
  test_quartic_real_coeff(&p3,"poly_rc_gsl",
			  1.0e-13,1.0e-13,1.0e-12,1.0e-11);
  test_quartic_real_coeff<long double,std::complex<long double> >
    (&q1_ld,"quartic_rc_cern_ld",
     1.0e-1,1.0e-2,1.0e2,1.0e3);
  test_quartic_real_coeff<long double,std::complex<long double> >
    (&q5_ld,"quartic_c_std_ld",
     1.0e-3,1.0e-2,1.0e1,1.0e2);
#ifndef O2SCL_OLD_BOOST
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&q1_cdf50,"quartic_rc_cern_50",
     1.0e-1,1.0e-1,1.0e2,1.0e3);
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&q5_cdf50,"quartic_c_std_50",
     1.0e-3,1.0e-2,1.0e2,1.0e2);
#endif
  cout << endl;
  
  cout << "Quartics with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_complex(&q5,"quartic_c_std",
		       1.0e-13,1.0e-13,1.0e-9,1.0e-9);
  test_quartic_complex<long double,std::complex<long double> >
    (&q5_ld,"quartic_c_std_ld",
     1.0e-17,1.0e-16,1.0e-13,1.0e-12);
#ifndef O2SCL_OLD_BOOST
  test_quartic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (&q5_cdf50,"quartic_c_std_50",
     1.0e-48,1.0e-48,1.0e-44,1.0e-44);
#endif
  cout << endl;
  
  tst.report();

  return 0;
}


