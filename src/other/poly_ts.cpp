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
#include <string>
#include <ctime>
#include <o2scl/test_mgr.h>
#include <o2scl/misc.h>
#include <o2scl/poly.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::cpp_bin_float_50 cpp_bin_float_50;
typedef boost::multiprecision::cpp_complex_50 cpp_complex_50;
#endif

test_mgr tst;
int wid=21;

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_real_coeff
(quadratic_real_coeff<fp_t,cx_t> *po, 
 string str, fp_t alpha, fp_t e1, fp_t e2, 
 fp_t e3, fp_t e4) {

  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count+=po->test_coeffs_zero_disc(alpha,s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  count+=po->test_real_roots(s1,s2,m1,m2);
  
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
       << m2 << " " << time << endl;

  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_real_coeff_boost
(quadratic_real_coeff<fp_t,cx_t> *po, 
 string str, fp_t alpha, fp_t e1, fp_t e2, 
 fp_t e3, fp_t e4) {

  fp_t s1=0, s2=0, m1=0, m2=0;

  clock_t lt1=clock();
  
  size_t count=0;
  
  count+=po->test_coeffs_zero_disc(alpha,s1,s2,m1,m2);
  count+=po->test_complex_roots(s1,s2,m1,m2);
  count+=po->test_real_roots(s1,s2,m1,m2);
  
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
       << m2 << " " << time << endl;

  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quadratic_complex(size_t ne, quadratic_complex<fp_t,cx_t> *po,
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
void test_quadratic_complex_boost(size_t ne, quadratic_complex<fp_t,cx_t> *po,
                                  string str, 
                                  double e1, double e2, double e3, double e4) {
  fp_t s1,s2,m1,m2;

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
void test_cubic_real_coeff_base(size_t ne, cubic_real_coeff<fp_t,cx_t> *po,
                                string str, 
                                fp_t alpha, fp_t &s1, fp_t &s2,
                                fp_t &m1, fp_t &m2, clock_t &lt1,
                                clock_t &lt2) {
  
  cx_t cr2,cr3,czo2,czo3,cap,cbp,ccp,cdp;
  fp_t ca,cb,cc,cd,cr1,czo1;
  cx_t i(0.0,1.0);
  size_t j;
  fp_t q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  
  gen_test_number ga, gb, gc, gd;
  for(int j1=0;j1<16;j1++) {
    ca=ga.gen()*alpha;
    gb.reset();
    for(int j2=0;j2<16;j2++) {
      cb=gb.gen();
      gc.reset();
      for(int j3=0;j3<16;j3++) {
	cc=gc.gen()*alpha;
        gd.reset();
	for(int j4=0;j4<16;j4++) {
	  cd=gd.gen();
	  
	  if (fabs(ca)>0.0) {
	    po->solve_rc(ca,cb,cc,cd,cr1,cr2,cr3);
	    
	    cbp=-(cr1+cr2+cr3)*ca;
	    ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
	    cdp=-(cr1*cr2*cr3)*ca;
	    
	    czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
	    czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
	    czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
	    q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
		    fabs(cc-ccp.real())*fabs(cc-ccp.real())+
		    fabs(cd-cdp.real())*fabs(cd-cdp.real()));
	    q2=sqrt(fabs(czo1)*fabs(czo1)+abs(czo2)*abs(czo2)+
		    abs(czo3)*abs(czo3));
	    s1+=q1;
	    if (q1>m1) m1=q1;
	    s2+=q2;
	    if (q2>m2) m2=q2;

            if (q1>1.0e-2 || q2>1.0e-2) {
              cout << "Problem." << endl;
              cout << ca << " " << cb << " " << cc << " " << cd << endl;
              cout << cr1 << " " << cr2 << " " << cr3 << endl;
              exit(-1);
            }
            
	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((fp_t)ne);
  s2/=((fp_t)ne);

  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_real_coeff(size_t ne, cubic_real_coeff<fp_t,cx_t> *po,
                           string str, 
			   fp_t alpha, fp_t e1, fp_t e2, fp_t e3, 
			   fp_t e4) {

  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;

  test_cubic_real_coeff_base(ne,po,str,alpha,s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"cubic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"cubic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"cubic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"cubic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_real_coeff_boost(size_t ne, cubic_real_coeff<fp_t,cx_t> *po,
                                 string str, 
                                 fp_t alpha, fp_t e1, fp_t e2, fp_t e3, 
                                 fp_t e4) {

  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;

  test_cubic_real_coeff_base(ne,po,str,alpha,s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"cubic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"cubic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"cubic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"cubic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_complex_base(size_t ne, cubic_complex<fp_t,cx_t> *po,
                             string str, fp_t &s1, fp_t &s2,
                             fp_t &m1, fp_t &m2, clock_t &lt1,
                             clock_t &lt2) {
  
  cx_t ca,cb,cc,cd,cr1,cr2,cr3,czo1,czo2,czo3,cap,cbp,ccp,cdp;
  cx_t i(0.0,1.0);
  size_t j;
  fp_t q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  gen_test_number ga, gb, gc, gd;
  fp_t rca, rcb, rcc, rcd;
  for(int it=0;it<2;it++) {
    for(int j1=0;j1<9;j1++) {
      rca=ga.gen();
      gb.reset();
      for(int j2=0;j2<9;j2++) {
	rcb=gb.gen();
        gc.reset();
	for(int j3=0;j3<9;j3++) {
	  rcc=gc.gen();
          gd.reset();
	  for(int j4=0;j4<9;j4++) {
	    rcd=gd.gen();

	    if (it==0) {
	      ca=rca+i;
	      cb=rcb+i;
	      cc=rcc+i;
	      cd=rcd+i;
	    } else {
              fp_t one=1.0;
	      ca=one+i*rca;
	      cb=one+i*rcb;
	      cc=one+i*rcc;
	      cd=one+i*rcd;
	    }

            if (fabs(ca.real())>0.0 || fabs(ca.imag())>0.0) {
              po->solve_c(ca,cb,cc,cd,cr1,cr2,cr3);
              
              cbp=-(cr1+cr2+cr3)*ca;
              ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
              cdp=-(cr1*cr2*cr3)*ca;
              
              czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
              czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
              czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
              q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
                      abs(cc-ccp)*abs(cc-ccp)+
                      abs(cd-cdp)*abs(cd-cdp));
              q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
                      abs(czo3)*abs(czo3));

              s1+=q1;
              if (q1>m1) m1=q1;
              s2+=q2;
              if (q2>m2) m2=q2;

            }

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((fp_t)ne);
  s2/=((fp_t)ne);
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_complex(size_t ne, cubic_complex<fp_t,cx_t> *po,
                        string str, fp_t e1, 
			fp_t e2, fp_t e3, fp_t e4) {
  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;
  
  test_cubic_complex_base(ne,po,str,s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"cubic_complex s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"cubic_complex s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"cubic_complex m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"cubic_complex m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_cubic_complex_boost(size_t ne, cubic_complex<fp_t,cx_t> *po,
                              string str, fp_t e1, 
                              fp_t e2, fp_t e3, fp_t e4) {
  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;
  
  test_cubic_complex_base(ne,po,str,s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"cubic_complex s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"cubic_complex s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"cubic_complex m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"cubic_complex m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

void test_quartic_real(size_t ne, quartic_real<double> *po,
                       string str, double alpha, 
		       double e1, double e2, double e3, double e4) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  double cr1,cr2,cr3,cr4;
  double r1,r2,r3,r4;
  double ca,cb,cc,cd,ce,zo1,zo2,zo3,zo4;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();

  gen_test_number ga, gb, gc, gd;
  for(int j1=0;j1<9;j1++) {
    r1=ga.gen();
    gb.reset();
    for(int j2=0;j2<9;j2++) {
      r2=-r1+alpha*gb.gen();
      gc.reset();
      for(int j3=0;j3<9;j3++) {
	r3=gc.gen();
        gd.reset();
	for(int j4=0;j4<9;j4++) {
	  r4=-r3+alpha*gd.gen();
	  
	  ca=1.0;
	  cb=-(r1+r2+r3+r4);
	  cc=(r1*r2+r1*r3+r2*r3+r1*r4+r2*r4+r3*r4);
	  cd=-(r1*r2*r3+r1*r2*r4+r1*r3*r4+r2*r3*r4);
	  ce=r1*r2*r3*r4;
	  
	  po->solve_r(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	  
	  zo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	  zo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	  zo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	  zo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	  q1=sqrt(zo1*zo1+zo2*zo2+zo3*zo3+zo4*zo4);
	  q2=0.0;
	  s1+=q1;
	  if (q1>m1) m1=q1;
	  s2+=q2;
	  if (q2>m2) m2=q2;
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs(s1,0.0,e1,"quartic_real s1");
  tst.test_abs(s2,0.0,e2,"quartic_real s2");
  tst.test_abs(m1,0.0,e3,"quartic_real m1");
  tst.test_abs(m2,0.0,e4,"quartic_real m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_real_coeff_base(size_t ne, quartic_real_coeff<fp_t,cx_t> *po,
                                  string str, fp_t &s1, fp_t &s2, fp_t &m1,
                                  fp_t &m2, clock_t &lt1, clock_t &lt2) {

  cx_t cr1,cr2,cr3,cr4,czo1,czo2,czo3,czo4;  
  cx_t i(0.0,1.0),cap,cbp,ccp,cdp,cep;
  fp_t ca,cb,cc,cd,ce;
  size_t j;
  fp_t q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();

  gen_test_number ga, gb, gc, gd, ge;
  for(int j1=0;j1<9;j1++) {
    ca=ga.gen();
    gb.reset();
    for(int j2=0;j2<9;j2++) {
      cb=gb.gen();
      gc.reset();
      for(int j3=0;j3<9;j3++) {
	cc=gc.gen();
        gd.reset();
	for(int j4=0;j4<9;j4++) {
	  cd=gd.gen();
          ge.reset();
	  for(int j5=0;j5<9;j5++) {
	    ce=ge.gen();
	    
	    if (fabs(ca)>0.0) {
	      po->solve_rc(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	      
	      cbp=-(cr1+cr2+cr3+cr4)*ca;
	      ccp=(cr1*cr2+cr1*cr3+cr2*cr3+cr1*cr4+cr2*cr4+cr3*cr4)*ca;
	      cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+cr1*cr3*cr4+cr2*cr3*cr4)*ca;
	      cep=cr1*cr2*cr3*cr4*ca;
	      
	      czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	      czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	      czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	      czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	      q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		      abs(cc-ccp)*abs(cc-ccp)+
		      abs(cd-cdp)*abs(cd-cdp)+
		      abs(ce-cep)*abs(ce-cep));
	      q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		      abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
	      s1+=q1;
	      if (q1>m1) m1=q1;
	      s2+=q2;
	      if (q2>m2) m2=q2;
	    }

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((fp_t)ne);
  s2/=((fp_t)ne);
  
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_real_coeff(size_t ne, quartic_real_coeff<fp_t,cx_t> *po,
                             string str, fp_t e1, fp_t e2, fp_t e3,
                             fp_t e4) {
  
  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;
  
  test_quartic_real_coeff_base<fp_t,cx_t>(ne,po,str,
                                          s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_real_coeff_boost(size_t ne, quartic_real_coeff<fp_t,cx_t> *po,
                                   string str, fp_t e1, fp_t e2, fp_t e3,
                                   fp_t e4) {
  
  fp_t s1=0,s2=0,m1=0,m2=0;
  clock_t lt1, lt2;
  
  test_quartic_real_coeff_base<fp_t,cx_t>(ne,po,str,
                                          s1,s2,m1,m2,lt1,lt2);
  
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_complex_base(size_t ne, quartic_complex<fp_t,cx_t> *po, 
                               string str, fp_t &s1, fp_t &s2, fp_t &m1,
                               fp_t &m2, clock_t &lt1, clock_t &lt2) {
  
  cx_t ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4,czo1,czo2,czo3,czo4;  
  cx_t i(0.0,1.0),cap,cbp,ccp,cdp,cep;
  size_t j;
  fp_t q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  cx_t one(1,0);

  gen_test_number ga, gb, gc, gd, ge;
  fp_t rca, rcb, rcc, rcd, rce;
  for(int it=0;it<2;it++) {
    for(int j1=0;j1<9;j1++) {
      rca=ga.gen();
      gb.reset();
      for(int j2=0;j2<9;j2++) {
	rcb=gb.gen();
        gc.reset();
	for(int j3=0;j3<9;j3++) {
	  rcc=gc.gen();
          gd.reset();
	  for(int j4=0;j4<9;j4++) {
	    rcd=gd.gen();
            ge.reset();
	    for(int j5=0;j5<9;j5++) {
	      rce=ge.gen();

	      if (it==0) {
		ca=rca+i;
		cb=rcb+i;
		cc=rcc+i;
		cd=rcd+i;
		ce=rce+i;
	      } else {
		ca=one+i*rca;
		cb=one+i*rcb;
		cc=one+i*rcc;
		cd=one+i*rcd;
		ce=one+i*rce;
	      }

	      po->solve_c(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	    
	      cbp=-(cr1+cr2+cr3+cr4)*ca;
	      ccp=(cr1*cr2+cr1*cr3+cr2*cr3+cr1*cr4+cr2*cr4+cr3*cr4)*ca;
	      cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+cr1*cr3*cr4+cr2*cr3*cr4)*ca;
	      cep=cr1*cr2*cr3*cr4*ca;
	    
	      czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	      czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	      czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	      czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	      q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		      abs(cc-ccp)*abs(cc-ccp)+
		      abs(cd-cdp)*abs(cd-cdp)+
		      abs(ce-cep)*abs(ce-cep));
	      q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		      abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
              //std::cout << q1 << " " << q2 << std::endl;
	      if (!o2isfinite(q1) || !o2isfinite(q2)
                  || q1>1.0e-10 || q2>1.0e-10) {
                std::cout << "ca,cb,cc,cd,ce: "
                          << ca << " " << cb << " " << cc << " "
                          << cd << " " << ce << std::endl;
                std::cout << "cr1,cr2,cr3,cr4: "
                          << cr1 << " " << cr2 << " "
                          << cr3 << " " << cr4 << std::endl;
                exit(-1);
		O2SCL_ERR("Failure in test_quartic_complex().",
			  exc_esanity);
	      }
	      s1+=q1;
	      if (q1>m1) m1=q1;
	      s2+=q2;
	      if (q2>m2) m2=q2;
	    }

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((fp_t)ne);
  s2/=((fp_t)ne);
  
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_complex(size_t ne, quartic_complex<fp_t,cx_t> *po,
                          string str,
			  fp_t e1, fp_t e2, fp_t e3, fp_t e4) {
  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;

  test_quartic_complex_base(ne,po,str,s1,s2,m1,m1,lt1,lt2);

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

template<class fp_t=double, class cx_t=std::complex<fp_t> >
void test_quartic_complex_boost(size_t ne, quartic_complex<fp_t,cx_t> *po,
                                string str,
                                fp_t e1, fp_t e2, fp_t e3, fp_t e4) {
  fp_t s1,s2,m1,m2;
  clock_t lt1, lt2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;

  test_quartic_complex_base(ne,po,str,s1,s2,m1,m1,lt1,lt2);

  cout.width(wid);
  cout << str.c_str();
  tst.test_abs_boost<fp_t>(s1,0.0,e1,"quartic_real_coeff s1");
  tst.test_abs_boost<fp_t>(s2,0.0,e2,"quartic_real_coeff s2");
  tst.test_abs_boost<fp_t>(m1,0.0,e3,"quartic_real_coeff m1");
  tst.test_abs_boost<fp_t>(m2,0.0,e4,"quartic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
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
  
  // Quartic solvers
  quartic_real_coeff_cern<> q1;
  quartic_real_gsl q2;
  quartic_real_gsl2 q3;
  quartic_real_std<> q4;
  quartic_complex_std<> q5;

#ifdef O2SCL_LD_TYPES
  
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
  
  cubic_complex_std<long double,std::complex<long double> > c3_ld;
  cubic_complex_std<cpp_bin_float_50,cpp_complex_50> c3_cdf50;
  
  quartic_real_coeff_cern<long double,std::complex<long double> > q1_ld;
  quartic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> q1_cdf50;

  quartic_real_std<long double> q4_ld;
  quartic_real_std<cpp_bin_float_50> q4_cdf50;

  quartic_complex_std<long double,std::complex<long double> > q5_ld;
  quartic_complex_std<cpp_bin_float_50,cpp_complex_50> q5_cdf50;
  
#endif

  // I think this number is no longer used, except to 
  // give an overall scale for the timings
  size_t ne=10000;

  cout << "Quadratics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0,
                            1.0e-13,1.0e-13,1.0e-10,1.0e-11);
#ifdef O2SCL_LD_TYPES
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0,
                                         1.0e-16,1.0e-16,1.0e-13,1.0e-14);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0,
                                         1.0e-16,1.0e-16,1.0e-13,1.0e-14);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl_50",1.0,
     1.0e-47,1.0e-48,1.0e-45,1.0e-45);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_c_std_50",1.0,
     1.0e-47,1.0e-48,1.0e-45,1.0e-45);
#endif
  cout << endl;
    
  cout << "Quadratics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(&t3,"quad_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11);
  test_quadratic_real_coeff(&t1,"quad_rc_gsl2",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11);
  test_quadratic_real_coeff(&t2,"quad_c_std",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,5.0e-11);
  test_quadratic_real_coeff(&p3,"poly_rc_gsl",1.0e-5,
                            1.0e-14,1.0e-14,1.0e-10,1.0e-11);
#ifdef O2SCL_LD_TYPES
  test_quadratic_real_coeff<long double>(&t1_ld,"quad_rc_gsl2_ld",1.0e-5,
                                         1.0e-17,1.0e-17,1.0e-13,1.0e-14);
  test_quadratic_real_coeff<long double>(&t2_ld,"quad_c_std_ld",1.0e-5,
                                         1.0e-17,1.0e-17,1.0e-13,1.0e-14);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t1_cdf50,"quad_rc_gsl_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-48,1.0e-48);
  test_quadratic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (&t2_cdf50,"quad_rc_std_50",1.0e-5,
     1.0e-49,1.0e-49,1.0e-48,1.0e-48);
#endif
  cout << endl;

  cout << "Quadratic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_complex(ne,&t2,"quad_c_std",
			 1.0e-13,1.0e-9,1.0e-11,1.0e-7);
#ifdef O2SCL_LD_TYPES
  test_quadratic_complex<long double>(ne,&t2_ld,"quad_c_std_ld",
                                      1.0e-17,1.0e-13,1.0e-14,1.0e-10);
  test_quadratic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&t2_cdf50,"quad_c_std_50",
     1.0e-48,1.0e-44,1.0e-45,1.0e-42);
#endif
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(ne,&c1,"cubic_rc_cern",1.0,
			1.0e0,1.0e6,1.0e1,1.0e7);
  test_cubic_real_coeff(ne,&c4,"cubic_rc_gsl",1.0,
			1.0e-1,1.0e-2,1.0e1,8.0e0);
  test_cubic_real_coeff(ne,&c2,"cubic_rc_gsl2",1.0,
			1.0e-1,1.0e-2,1.0e1,8.0e0);
  test_cubic_real_coeff(ne,&c3,"cubic_c_std",1.0,
			1.0e-1,1.0e-1,1.0e1,1.0e1);
  test_cubic_real_coeff(ne,&p3,"poly_rc_gsl",1.0,
			1.0e-1,4.0e-2,1.0e1,1.0e1);
#ifdef O2SCL_LD_TYPES
  test_cubic_real_coeff<long double,std::complex<long double> >
    (ne,&c1_ld,"cubic_rc_cern_ld",
     1.0,1.0e-12,1.0e-12,1.0e-9,1.0e-10);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (ne,&c3_ld,"cubic_c_std_ld",
     1.0,1.0e-1,1.0e-1,1.0e1,1.0e1);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&c1_cdf50,"cubic_rc_cern_50",
     1.0,1.0e-40,1.0e-40,1.0e-38,1.0e-37);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&c3_cdf50,"cubic_c_std_50",
     1.0,1.0e-1,1.0e-1,1.0e1,1.0e1);
#endif
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(ne,&c1,"cubic_rc_cern",1.0e-3,
			1.0e-5,1.0e-3,1.0e-3,1.0e-2);
  test_cubic_real_coeff(ne,&c4,"cubic_rc_gsl",1.0e-3,
			1.0e-4,5.0e-4,1.0e-3,1.0e-3);
  test_cubic_real_coeff(ne,&c2,"cubic_rc_gsl2",1.0e-3,
			1.0e-4,5.0e-4,1.0e-3,1.0e-3);
  test_cubic_real_coeff(ne,&c3,"cubic_c_std",1.0e-3,
			1.0e-1,1.0e+2,4.0e0,1.0e+5);
  test_cubic_real_coeff(ne,&p3,"poly_rc_gsl",1.0e-3,
			1.0e-4,1.0e-3,1.0e-2,1.0e-2);
#ifdef O2SCL_LD_TYPES
  test_cubic_real_coeff<long double,std::complex<long double> >
    (ne,&c1_ld,"cubic_rc_cern_ld",
     1.0e-3,1.0e-7,1.0e-5,1.0e-6,1.0e-4);
  test_cubic_real_coeff<long double,std::complex<long double> >
    (ne,&c3_ld,"cubic_c_std_ld",
     1.0e-3,1.0e1,1.0e-1,1.0e1,1.0e1);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&c1_cdf50,"cubic_rc_cern_50",
     1.0e-3,1.0e-38,1.0e-36,1.0e-37,1.0e-35);
  test_cubic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&c3_cdf50,"cubic_c_std_50",
     1.0e-3,1.0e-1,1.0e-1,1.0e1,1.0e1);
#endif
  cout << endl;
  
  cout << "Cubic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_complex(ne,&c3,"cubic_c_std",
                     1.0e1,1.0e1,1.0e1,1.0e1);
#ifdef O2SCL_LD_TYPES
  test_cubic_complex<long double,std::complex<long double> >
    (ne,&c3_ld,"cubic_c_std_ld",
     1.0e1,1.0e1,1.0e1,1.0e1);
  test_cubic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&c3_cdf50,"cubic_c_std_50",
     1.0e1,1.0e1,1.0e1,1.0e1);
#endif
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n";
  cout << " leading coefficient 1:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(ne,&q1,"quartic_rc_cern",1.0,
		    5.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real(ne,&q2,"quartic_real_gsl",1.0,
		    5.0e-2,1.0,1.0e2,1.0);
  test_quartic_real(ne,&q3,"quartic_real_gsl2",1.0,
		    1.0e-2,1.0,1.0e2,1.0);
  test_quartic_real(ne,&q4,"quartic_real_std",1.0,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(ne,&q5,"quartic_c_std",1.0,
		    1.0e-1,1.0,1.0e2,1.0);
  test_quartic_real(ne,&p3,"poly_rc_gsl",1.0,
		    1.0e-14,1.0,1.0e-12,1.0);
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n"
       << " leading coefficient 1, coefficients of odd powers small:" 
       << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(ne,&q1,"quartic_rc_cern",1.0e-5,
		    1.0e-2,1.0,1.0e1,1.0);
  test_quartic_real(ne,&q2,"quartic_real_gsl",1.0e-5,
		    1.0e-3,1.0,1.0e1,1.0);
  test_quartic_real(ne,&q3,"quartic_real_gsl2",1.0e-5,
		    1.0e-4,1.0,1.0e-3,1.0);
  test_quartic_real(ne,&q4,"quartic_real_std",3.0e-4,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(ne,&q5,"quartic_c_std",1.0e-5,
		    1.0e-1,1.0,5.0,1.0);
  test_quartic_real(ne,&p3,"poly_rc_gsl",1.0e-5,
		    1.0e-14,1.0,1.0e-13,1.0);
  cout << endl;

  cout << "Quartics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real_coeff(ne,&q1,"cern_real_coeff",
			  1.0e0,1.0e4,1.0e2,1.0e5);
  test_quartic_real_coeff(ne,&q5,"quartic_c_std",
			  1.0e1,1.0e1,1.0e3,1.0e4);
  test_quartic_real_coeff(ne,&p3,"poly_rc_gsl",
			  1.0e-13,1.0e-6,1.0e-10,1.0e-4);
#ifdef O2SCL_LD_TYPES
  test_quartic_real_coeff<long double,std::complex<long double> >
    (ne,&q1_ld,"quartic_rc_cern_ld",
     1.0e-1,1.0e-1,1.0e1,1.0e2);
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&q1_cdf50,"quartic_rc_cern_50",
     1.0e-1,1.0e1,1.0e2,1.0e2);
#endif
  cout << endl;
  
  cout << "Quartics with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_complex(ne,&q5,"quartic_c_std",
		       1.0e-2,1.0e-2,1.0e1,1.0e2);
#ifdef O2SCL_LD_TYPES
  test_quartic_complex<long double,std::complex<long double> >
    (ne,&q5_ld,"quartic_c_std_ld",
     1.0e-2,1.0e-2,1.0e1,1.0e2);
  test_quartic_complex_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&q5_cdf50,"quartic_c_std_50",
     1.0e-2,1.0e-2,1.0e1,1.0e2);
#endif
  cout << endl;
  
  tst.report();

  return 0;
}


