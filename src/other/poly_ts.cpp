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

#define SWAP(a,b) do { double tmp=b ; b=a ; a=tmp ; } while(0)

// Test patch suggested by Lorenzo Moneta

int
gsl_poly_complex_solve_cubic_lm(double a, double b, double c, 
                                gsl_complex *z0, gsl_complex *z1, 
                                gsl_complex *z2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

//  double CR2 = 729 * r * r;
//  double CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0)
    {
      GSL_REAL (*z0) = -a / 3;
      GSL_IMAG (*z0) = 0;
      GSL_REAL (*z1) = -a / 3;
      GSL_IMAG (*z1) = 0;
      GSL_REAL (*z2) = -a / 3;
      GSL_IMAG (*z2) = 0;
      return 3;
    }
  else if (R2 == Q3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         will be considered to be a pair of complex roots z = x +/-
         epsilon i close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          GSL_REAL (*z0) = -2 * sqrtQ - a / 3;
          GSL_IMAG (*z0) = 0;
          GSL_REAL (*z1) = sqrtQ - a / 3;
          GSL_IMAG (*z1) = 0;
          GSL_REAL (*z2) = sqrtQ - a / 3;
          GSL_IMAG (*z2) = 0;
        }
      else
        {
          GSL_REAL (*z0) = -sqrtQ - a / 3;
          GSL_IMAG (*z0) = 0;
          GSL_REAL (*z1) = -sqrtQ - a / 3;
          GSL_IMAG (*z1) = 0;
          GSL_REAL (*z2) = 2 * sqrtQ - a / 3;
          GSL_IMAG (*z2) = 0;
        }
      return 3;
    }
  else if (R2 < Q3)  /* equivalent to R2 < Q3 */
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double ctheta = R / sqrtQ3;
      double theta = 0; 
      if ( ctheta <= -1.0) 
         theta = M_PI;
      else if ( ctheta < 1.0) 
         theta = acos (R / sqrtQ3);
      
      double norm = -2 * sqrtQ;
      double r0 = norm * cos (theta / 3) - a / 3;
      double r1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      double r2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

      /* Sort r0, r1, r2 into increasing order */

      if (r0 > r1)
        SWAP (r0, r1);

      if (r1 > r2)
        {
          SWAP (r1, r2);

          if (r0 > r1)
            SWAP (r0, r1);
        }

      GSL_REAL (*z0) = r0;
      GSL_IMAG (*z0) = 0;

      GSL_REAL (*z1) = r1;
      GSL_IMAG (*z1) = 0;

      GSL_REAL (*z2) = r2;
      GSL_IMAG (*z2) = 0;

      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0 / 3.0);
      double B = Q / A;

      if (A + B < 0)
        {
          GSL_REAL (*z0) = A + B - a / 3;
          GSL_IMAG (*z0) = 0;

          GSL_REAL (*z1) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z1) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z2) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z2) = (sqrt (3.0) / 2.0) * fabs(A - B);
        }
      else
        {
          GSL_REAL (*z0) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z0) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z1) = -0.5 * (A + B) - a / 3;
          GSL_IMAG (*z1) = (sqrt (3.0) / 2.0) * fabs(A - B);

          GSL_REAL (*z2) = A + B - a / 3;
          GSL_IMAG (*z2) = 0;
        }
      
      return 3;
    }
}


void test_quadratic_real_coeff(size_t ne, quadratic_real_coeff<> *po, 
			       string str, double alpha, double e1, 
			       double e2, double e3, double e4,
			       double &time_taken);
void test_quadratic_complex(size_t ne, quadratic_complex<> *po, string str, 
			    double e1, double e2, double e3, double e4);
void test_cubic_real_coeff(size_t ne, cubic_real_coeff<double> *po, string str, 
			   double alpha, double e1, double e2, double e3, 
			   double e4);
void test_cubic_complex(size_t ne, cubic_complex<> *po, string str, double e1, 
			double e2, double e3, double e4);
void test_quartic_real(size_t ne, quartic_real<double> *po, string str,
                       double alpha, 
		       double e1, double e2, double e3, double e4);
void test_quartic_complex(size_t ne, quartic_complex<> *po, string str, 
			  double e1, double e2, double e3, double e4);

void test_quadratic_real_coeff(size_t ne, quadratic_real_coeff<> *po, 
			       string str, double alpha, double e1, double e2, 
			       double e3, double e4, double &time_taken) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> cr2,cr3,czo2,czo3,cap,cbp,ccp,cdp,cr1,czo1;
  double ca,cb,cc;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();

  gen_test_number<40> ga, gb, gc;
  for(int j1=0;j1<40;j1++) {
    ca=ga.gen();
    for(int j2=0;j2<40;j2++) {
      cb=gb.gen()*alpha;
      for(int j3=0;j3<40;j3++) {
	cc=gc.gen();

	// Ensure there is a solution
	if (fabs(ca)>0.0) {

	  po->solve_rc(ca,cb,cc,cr1,cr2);
	  
	  cbp=-(cr1+cr2)*ca;
	  ccp=(cr1*cr2)*ca;
	  
	  czo1=(ca*cr1+cb)*cr1+cc;
	  czo2=(ca*cr2+cb)*cr2+cc;
	  q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
		  fabs(cc-ccp.real())*fabs(cc-ccp.real()));
	  q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
	  
	  s1+=q1;
	  if (q1>m1) m1=q1;
	  s2+=q2;
	  if (q2>m2) m2=q2;

	  if (!std::isfinite(q1) || !std::isfinite(q2) || 
	      fabs(q1)>1.0e-10 || fabs(q2)>1.0e-10) {
	    O2SCL_ERR("Failure in test_quadratic_real_coeff().",
		      exc_esanity);
	  }
	  
	}

      }
    }
  }

  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs(s1,0.0,e1,"quadratic_real_coeff s1");
  tst.test_abs(s2,0.0,e2,"quadratic_real_coeff s2");
  tst.test_abs(m1,0.0,e3,"quadratic_real_coeff m1");
  tst.test_abs(m2,0.0,e4,"quadratic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " ";
  time_taken=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  cout << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;

  return;
}

void test_quadratic_complex(size_t ne, quadratic_complex<> *po, string str, 
			    double e1, double e2, double e3, double e4) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,cr1,cr2,cr3,czo1,czo2,czo3,cap,cbp,ccp,cdp;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  
  gen_test_number<9> ga, gb, gc, gd, ge, gf;
  double rca, rcb, rcc, rcd, rce, rcf;
  for(int j1=0;j1<9;j1++) {
    rca=ga.gen();
    for(int j2=0;j2<9;j2++) {
      rcb=gb.gen();
      for(int j3=0;j3<9;j3++) {
	rcc=gc.gen();
	for(int j4=0;j4<9;j4++) {
	  rcd=gd.gen();
	  for(int j5=0;j5<9;j5++) {
	    rce=ge.gen();
	    for(int j6=0;j6<9;j6++) {
	      rcf=gf.gen();
	      
	      if (fabs(rca)>0.0 || fabs(rcb)>0.0) {
		ca=rca+i*rcb;
		cb=rcc+i*rcd;
		cc=rce+i*rcf;
    
		po->solve_c(ca,cb,cc,cr1,cr2);
	      
		cbp=-(cr1+cr2)*ca;
		ccp=(cr1*cr2)*ca;
	      
		czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
		czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
		q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
			abs(cc-ccp)*abs(cc-ccp));
		q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
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
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  tst.test_abs(s1,0.0,e1,"quadratic_complex s1");
  tst.test_abs(s2,0.0,e2,"quadratic_complex s2");
  tst.test_abs(m1,0.0,e3,"quadratic_complex m1");
  tst.test_abs(m2,0.0,e4,"quadratic_complex m2");
  cout.width(wid);
  cout << str.c_str();
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

void test_cubic_real_coeff(size_t ne, cubic_real_coeff<double> *po, string str, 
			   double alpha, double e1, double e2, double e3, 
			   double e4) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> cr2,cr3,czo2,czo3,cap,cbp,ccp,cdp;
  double ca,cb,cc,cd,cr1,czo1;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  
  gen_test_number<16> ga, gb, gc, gd;
  for(int j1=0;j1<16;j1++) {
    ca=ga.gen()*alpha;
    for(int j2=0;j2<16;j2++) {
      cb=gb.gen();
      for(int j3=0;j3<16;j3++) {
	cc=gc.gen()*alpha;
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
	    
	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs(s1,0.0,e1,"cubic_real_coeff s1");
  tst.test_abs(s2,0.0,e2,"cubic_real_coeff s2");
  tst.test_abs(m1,0.0,e3,"cubic_real_coeff m1");
  tst.test_abs(m2,0.0,e4,"cubic_real_coeff m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

void compare_gsl_cubic(size_t ne, string str, 
		       double alpha, int sw, cubic_real_coeff_gsl &gcrc) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> czo2,czo3,cap,cbp,ccp,cdp,czo1;
  double ca,cb,cc,cd;
  complex<double> i(0.0,1.0);
  gsl_complex z0,z1,z2;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();

  size_t nfails=0;
  const int ntest=16;

  gen_test_number<ntest> ga, gb, gc, gd;
  for(int j2=0;j2<ntest;j2++) {
    cb=gb.gen();
    for(int j3=0;j3<ntest;j3++) {
      cc=gc.gen()*alpha;
      for(int j4=0;j4<ntest;j4++) {
	cd=gd.gen();
	
	if (sw==0) {
	  gsl_poly_complex_solve_cubic(cb,cc,cd,&z0,&z1,&z2);
	} else if (sw==1) {
	  gcrc.gsl_poly_complex_solve_cubic2(cb,cc,cd,&z0,&z1,&z2);
	} else {
	  gsl_poly_complex_solve_cubic_lm(cb,cc,cd,&z0,&z1,&z2);
	}
	if (!std::isfinite(GSL_REAL(z0)) ||
	    !std::isfinite(GSL_IMAG(z0)) ||
	    !std::isfinite(GSL_REAL(z1)) ||
	    !std::isfinite(GSL_IMAG(z1)) ||
	    !std::isfinite(GSL_REAL(z2)) ||
	    !std::isfinite(GSL_IMAG(z2))) {
	  nfails++;
	}
	std::complex<double> cr1(GSL_REAL(z0),GSL_IMAG(z0));
	std::complex<double> cr2(GSL_REAL(z2),GSL_IMAG(z1));
	std::complex<double> cr3(GSL_REAL(z2),GSL_IMAG(z1));
	
	cbp=-(cr1+cr2+cr3);
	ccp=(cr1*cr2+cr1*cr3+cr2*cr3);
	cdp=-(cr1*cr2*cr3);
	
	czo1=((cr1+cb)*cr1+cc)*cr1+cd;
	czo2=((cr2+cb)*cr2+cc)*cr2+cd;
	czo3=((cr3+cb)*cr3+cc)*cr3+cd;
	q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
		fabs(cc-ccp.real())*fabs(cc-ccp.real())+
		fabs(cd-cdp.real())*fabs(cd-cdp.real()));
	q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		abs(czo3)*abs(czo3));

	if (std::isfinite(q1)) {
	  s1+=q1;
	  if (q1>m1) m1=q1;
	} 
	if (std::isfinite(q2)) {
	  s2+=q2;
	  if (q2>m2) m2=q2;
	}
	
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid-4);
  cout << str.c_str();
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << " "
       << nfails << endl;
  return;
}

void test_cubic_complex(size_t ne, cubic_complex<> *po, string str, double e1, 
			double e2, double e3, double e4) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,cr1,cr2,cr3,czo1,czo2,czo3,cap,cbp,ccp,cdp;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  gen_test_number<9> ga, gb, gc, gd;
  double rca, rcb, rcc, rcd;
  for(int it=0;it<2;it++) {
    for(int j1=0;j1<9;j1++) {
      rca=ga.gen();
      for(int j2=0;j2<9;j2++) {
	rcb=gb.gen();
	for(int j3=0;j3<9;j3++) {
	  rcc=gc.gen();
	  for(int j4=0;j4<9;j4++) {
	    rcd=gd.gen();

	    if (it==0) {
	      ca=rca+i;
	      cb=rcb+i;
	      cc=rcc+i;
	      cd=rcd+i;
	    } else {
	      ca=1.0+i*rca;
	      cb=1.0+i*rcb;
	      cc=1.0+i*rcc;
	      cd=1.0+i*rcd;
	    }
    
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
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs(s1,0.0,e1,"cubic_complex s1");
  tst.test_abs(s2,0.0,e2,"cubic_complex s2");
  tst.test_abs(m1,0.0,e3,"cubic_complex m1");
  tst.test_abs(m2,0.0,e4,"cubic_complex m2");
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

  gen_test_number<9> ga, gb, gc, gd;
  for(int j1=0;j1<9;j1++) {
    r1=ga.gen();
    for(int j2=0;j2<9;j2++) {
      r2=-r1+alpha*gb.gen();
      for(int j3=0;j3<9;j3++) {
	r3=gc.gen();
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

  gen_test_number<9> ga, gb, gc, gd, ge;
  for(int j1=0;j1<9;j1++) {
    ca=ga.gen();
    for(int j2=0;j2<9;j2++) {
      cb=gb.gen();
      for(int j3=0;j3<9;j3++) {
	cc=gc.gen();
	for(int j4=0;j4<9;j4++) {
	  cd=gd.gen();
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
  
  fp_t s1,s2,m1,m2;
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

void test_quartic_complex(size_t ne, quartic_complex<> *po, string str, 
			  double e1, double e2, double e3, double e4) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4,czo1,czo2,czo3,czo4;  
  complex<double> i(0.0,1.0),cap,cbp,ccp,cdp,cep;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();

  gen_test_number<9> ga, gb, gc, gd, ge;
  double rca, rcb, rcc, rcd, rce;
  for(int it=0;it<2;it++) {
    for(int j1=0;j1<9;j1++) {
      rca=ga.gen();
      for(int j2=0;j2<9;j2++) {
	rcb=gb.gen();
	for(int j3=0;j3<9;j3++) {
	  rcc=gc.gen();
	  for(int j4=0;j4<9;j4++) {
	    rcd=gd.gen();
	    for(int j5=0;j5<9;j5++) {
	      rce=ge.gen();

	      if (it==0) {
		ca=rca+i;
		cb=rcb+i;
		cc=rcc+i;
		cd=rcd+i;
		ce=rce+i;
	      } else {
		ca=1.0+i*rca;
		cb=1.0+i*rcb;
		cc=1.0+i*rcc;
		cd=1.0+i*rcd;
		ce=1.0+i*rce;
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
	      if (!std::isfinite(q1)) {
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
  s1/=((double)ne);
  s2/=((double)ne);
  cout.width(wid);
  cout << str.c_str();
  tst.test_abs(s1,0.0,e1,"quartic_complex s1");
  tst.test_abs(s2,0.0,e2,"quartic_complex s2");
  tst.test_abs(m1,0.0,e3,"quartic_complex m1");
  tst.test_abs(m2,0.0,e4,"quartic_complex m2");
  cout << ": " << s1 << " " << s2 << " " << m1 << " " 
       << m2 << " " << ((double)(lt2-lt1))/CLOCKS_PER_SEC << endl;
  return;
}

int main(void) {
  tst.set_output_level(1);
  
  cout.setf(ios::left | ios::scientific);
  cout.precision(4);

  // Generic polynomial solvers
  poly_real_coeff_gsl p3;
  
  // quadratic solvers
  quadratic_real_coeff_gsl<> t1;
  quadratic_complex_std<> t2;
  
  // cubic solvers
  cubic_real_coeff_cern<> c1;
  cubic_real_coeff_gsl c2;
  cubic_complex_std<> c3;
  
  // quartic solvers
  quartic_real_coeff_cern<> q1;
  quartic_real_gsl q2;
  quartic_real_gsl2 q3;
  quartic_real_simple q4;
  quartic_complex_simple<> q5;

#ifdef O2SCL_LD_TYPES
  quartic_real_coeff_gsl<cpp_bin_float_50,cpp_complex_50> t1_cdf50;
  quartic_real_coeff_cern<cpp_bin_float_50,cpp_complex_50> q1_cdf50;
#endif
  
  // I think this number is no longer used, except to 
  // give an overall scale for the timings
  size_t ne=10000;

  double tt;
  
  cout << "Quadratics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(ne,&t1,"gsl_quad_real_coeff",1.0,
			    1.0e-15,1.0e-13,1.0e-15,1.0e-9,tt);
  test_quadratic_real_coeff(ne,&t2,"quadratic_complex_std",1.0,
			    1.0e-13,1.0e-12,1.0e-9,1.0e-9,tt);
  test_quadratic_real_coeff(ne,&p3,"poly_real_coeff_gsl",1.0,
			    1.0e-13,1.0e-12,1.0e-9,1.0e-9,tt);
  cout << endl;

  cout << "Quadratics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_real_coeff(ne,&t1,"gsl_quad_real_coeff",1.0e-5,
			    1.0e-15,1.0e-15,1.0e-15,1.0e-15,tt);
  test_quadratic_real_coeff(ne,&t2,"quadratic_complex_std",1.0e-5,
			    1.0e-15,1.0e-15,1.0e-15,5.0e-15,tt);
  test_quadratic_real_coeff(ne,&p3,"poly_real_coeff_gsl",1.0e-5,
			    1.0e-15,1.0e-15,1.0e-15,1.0e-15,tt);
  cout << endl;

  cout << "Quadratic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quadratic_complex(ne,&t2,"quadratic_complex_std",
			 1.0e-12,1.0e-9,1.0e-12,1.0e-8);
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(ne,&c1,"cern_real_coeff",1.0,
			1.0e0,1.0e6,1.0e1,1.0e7);
  test_cubic_real_coeff(ne,&c2,"cubic_real_coeff_gsl",1.0,
			1.0e-1,1.0e-2,1.0e1,8.0e0);
  test_cubic_real_coeff(ne,&c3,"cubic_complex_std",1.0,
			1.0e-1,1.0e-2,1.0e1,1.0e1);
  test_cubic_real_coeff(ne,&p3,"poly_real_coeff_gsl",1.0,
			1.0e-1,4.0e-2,1.0e1,1.0e1);
  cout << endl;
  
  cout << "Cubics with real coefficients and complex roots -\n"
       << " coefficients of odd powers small:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_real_coeff(ne,&c1,"cern_real_coeff",1.0e-3,
			1.0e-5,1.0e-3,1.0e-3,1.0e-2);
  test_cubic_real_coeff(ne,&c2,"cubic_real_coeff_gsl",1.0e-3,
			1.0e-4,5.0e-4,1.0e-3,1.0e-3);
  test_cubic_real_coeff(ne,&c3,"cubic_complex_std",1.0e-3,
			1.0e-4,1.0e+2,1.0e-2,1.0e+5);
  test_cubic_real_coeff(ne,&p3,"poly_real_coeff_gsl",1.0e-3,
			1.0e-5,1.0e-3,1.0e-2,1.0e-2);
  cout << endl;
  
  cout << "Cubic with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_cubic_complex(ne,&c3,"cubic_complex_std",
		     1.0e1,1.0e1,1.0e1,1.0e1);
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n";
  cout << " leading coefficient 1:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(ne,&q1,"cern_real_coeff",1.0,
		    5.0e-12,1.0,2.0e-8,1.0);
  test_quartic_real(ne,&q2,"quartic_real_gsl",1.0,
		    5.0e-2,1.0,1.0e2,1.0);
  test_quartic_real(ne,&q3,"quartic_real_gsl2",1.0,
		    1.0e-2,1.0,1.0e1,1.0);
  test_quartic_real(ne,&q4,"quartic_real_simple",1.0,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(ne,&q5,"simple_quartic_cmplx",1.0,
		    1.0e-1,1.0,1.0e2,1.0);
  test_quartic_real(ne,&p3,"poly_real_coeff_gsl",1.0,
		    1.0e-14,1.0,1.0e-13,1.0);
  cout << endl;
  
  cout << "Quartics with real coefficients and real roots -\n"
       << " leading coefficient 1, coefficients of odd powers small:" 
       << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real(ne,&q1,"cern_real_coeff",1.0e-5,
		    1.0e-3,1.0,1.0e1,1.0);
  test_quartic_real(ne,&q2,"quartic_real_gsl",1.0e-5,
		    1.0e-3,1.0,1.0e1,1.0);
  test_quartic_real(ne,&q3,"quartic_real_gsl2",1.0e-5,
		    1.0e-5,1.0,1.0e-3,1.0);
  test_quartic_real(ne,&q4,"quartic_real_simple",3.0e-4,
		    2.0e5,1.0,3.0e7,1.0);
  test_quartic_real(ne,&q5,"simple_quartic_cmplx",1.0e-5,
		    1.0e-1,1.0,5.0,1.0);
  test_quartic_real(ne,&p3,"poly_real_coeff_gsl",1.0e-5,
		    1.0e-15,1.0,1.0e-14,1.0);
  cout << endl;

  cout << "Quartics with real coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_real_coeff(ne,&q1,"cern_real_coeff",
			  1.0e0,1.0e4,1.0e2,1.0e5);
  test_quartic_real_coeff(ne,&q5,"simple_quartic_cmplx",
			  1.0e1,1.0e1,1.0e3,1.0e4);
  test_quartic_real_coeff(ne,&p3,"poly_real_coeff_gsl",
			  1.0e-13,1.0e-6,1.0e-10,1.0e-4);
#ifdef O2SCL_LD_TYPES
  test_quartic_real_coeff_boost<cpp_bin_float_50,cpp_complex_50>
    (ne,&q1_cdf50,"cern_real_coeff_50",
     1.0e0,1.0e4,1.0e2,1.0e5);
#endif
  cout << endl;
  
  cout << "Quartics with complex coefficients and complex roots:" << endl;
  cout << "type                   Avg 1      Avg 2      Max 1"
       << "      Max 2      time" << endl;
  test_quartic_complex(ne,&q5,"simple_quartic_cmplx",
		       1.0e-2,1.0e-2,1.0e1,1.0e2);
  cout << endl;
  
  cout << "Compare gsl cubic functions:" << endl;
  cout << "On some systems, the GSL version gives solutions which are not" 
       << " finite:" << endl;
  compare_gsl_cubic(ne,"GSL",1.0e-9,0,c2);
  compare_gsl_cubic(ne,"Revised",1.0e-9,1,c2);
  compare_gsl_cubic(ne,"LM",1.0e-9,2,c2);

  tst.report();

  return 0;
}


