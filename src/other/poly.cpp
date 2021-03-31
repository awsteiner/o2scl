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
/* poly/solve_cubic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/poly.h>
// For is_finite()
#include <o2scl/misc.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int quadratic_real_coeff_gsl::solve_rc
(double a2, double b2, double c2,
 complex<double> &x1, complex<double> &x2) { 
  gsl_complex r1, r2;

  gsl_poly_complex_solve_quadratic(a2,b2,c2,&r1,&r2);
  
  x1=complex<double>(GSL_REAL(r1),GSL_IMAG(r1));
  x2=complex<double>(GSL_REAL(r2),GSL_IMAG(r2));

  return success;
}

inline double max(double a, double b) {
  if (a>b) return a;
  return b;
}

inline void swap(double &a, double &b) {
  double tmp=b;
  b=a;
  a=tmp;
  return;
}

poly_real_coeff_gsl::poly_real_coeff_gsl() {
  w2=gsl_poly_complex_workspace_alloc(3);
  w3=gsl_poly_complex_workspace_alloc(4);
  w4=gsl_poly_complex_workspace_alloc(5);
  gen_size=0;
}

poly_real_coeff_gsl::~poly_real_coeff_gsl() {
  gsl_poly_complex_workspace_free(w2);
  gsl_poly_complex_workspace_free(w3);
  gsl_poly_complex_workspace_free(w4);
  if (gen_size>0) {
    gsl_poly_complex_workspace_free(wgen);
  }
}

int poly_real_coeff_gsl::solve_rc(const double a2, const double b2, 
				  const double c2, complex<double> &r1, 
				  complex<double> &r2) {
  
  if (a2==0.0) {
    O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
  }

  double a[3]={c2,b2,a2};
  double z[4];
  complex<double> i(0.0,1.0);

  gsl_poly_complex_solve(a,3,w2,z);
  
  r1=z[0]+i*z[1];
  r2=z[2]+i*z[3];

  return success;
}

int poly_real_coeff_gsl::solve_rc(const double a3, const double b3, 
				  const double c3, const double d3, 
				  double &r1, 
				  complex<double> &r2, complex<double> &r3) {
  
  if (a3==0.0) {
    O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
  }

  double a[4]={d3,c3,b3,a3};  
  double z[6],s1,s2,s3;
  complex<double> i(0.0,1.0);

  gsl_poly_complex_solve(a,4,w3,z);
  
  s1=fabs(z[1]/z[0]);
  s2=fabs(z[3]/z[2]);
  s3=fabs(z[5]/z[4]);
  if (s1<s2 && s1<s3) {
    r1=z[0];
    r2=z[2]+i*z[3];
    r3=z[4]+i*z[5];
  } else if (s2<s1 && s2<s3) {
    r1=z[2];
    r2=z[0]+i*z[1];
    r3=z[4]+i*z[5];
  } else {
    r1=z[4];
    r2=z[0]+i*z[1];
    r3=z[2]+i*z[3];
  }

  return success;
}

int poly_real_coeff_gsl::solve_rc
(const double a4, const double b4, const double c4, const double d4, 
 const double e4, complex<double> &r1, complex<double> &r2, 
 complex<double> &r3, complex<double> &r4) {
  
  if (a4==0.0) {
    O2SCL_ERR2("Leading coefficient zero in ",
		  "poly_real_coeff_gsl::solve_rc().",
		  exc_einval);
  }

  double a[5]={e4,d4,c4,b4,a4};
  double z[8];
  complex<double> i(0.0,1.0);

  gsl_poly_complex_solve(a,5,w4,z);
  
  r1=z[0]+i*z[1];
  r2=z[2]+i*z[3];
  r3=z[4]+i*z[5];
  r4=z[6]+i*z[7];

  return success;
}

int poly_real_coeff_gsl::solve_rc_arr(int n, const double co[], 
				      std::complex<double> ro[]) {
  int j;
  typedef boost::numeric::ublas::vector<double> ubvector;
  ubvector a(n+1), z(2*n);
  complex<double> i(0.0,1.0);
  
  for(j=0;j<n+1;j++) {
    a[j]=co[n-j];
  }  
  if (gen_size!=n) {
    if (gen_size>0) {
      gsl_poly_complex_workspace_free(wgen);
    }
    wgen=gsl_poly_complex_workspace_alloc(n+1);
    gen_size=n;
  }
  if (a[n]==0.0) {
    O2SCL_ERR2("Leading coefficient zero in ",
	       "poly_real_coeff_gsl::solve_rc().",
	       exc_einval);
  }
  gsl_poly_complex_solve(&a[0],n+1,wgen,&z[0]);
  
  for(j=0;j<n;j++) {
    ro[j]=z[2*j]+i*z[2*j+1];
  }

  return success;
}

int cubic_real_coeff_gsl::gsl_poly_complex_solve_cubic2
(double a, double b, double c, gsl_complex *z0, gsl_complex *z1, 
 gsl_complex *z2) {
  
  double q=(a*a-3*b);
  double r=(2*a*a*a-9*a*b+27*c);
  
  double Q=q/9;
  double R=r/54;

  double Q3=Q*Q*Q;
  double R2=R*R;

  double CR2=729*r*r;
  double CQ3=2916*q*q*q;

  if (R == 0 && Q == 0) {
    GSL_REAL(*z0)=-a/3;
    GSL_IMAG(*z0)=0;
    GSL_REAL(*z1)=-a/3;
    GSL_IMAG(*z1)=0;
    GSL_REAL(*z2)=-a/3;
    GSL_IMAG(*z2)=0;
    return 3;
  } else if (CR2 == CQ3)  {

    /* this test is actually R2 == Q3, written in a form suitable
       for exact computation with integers */
    
    /* Due to finite precision some double roots may be missed, and
       will be considered to be a pair of complex roots z=x +/-
       epsilon i close to the real axis. */
    
    double sqrtQ=sqrt(Q);

    if (R > 0) {
      GSL_REAL(*z0)=-2*sqrtQ-a/3;
      GSL_IMAG(*z0)=0;
      GSL_REAL(*z1)=sqrtQ-a/3;
      GSL_IMAG(*z1)=0;
      GSL_REAL(*z2)=sqrtQ-a/3;
      GSL_IMAG(*z2)=0;
    } else {
      GSL_REAL(*z0)=-sqrtQ-a/3;
      GSL_IMAG(*z0)=0;
      GSL_REAL(*z1)=-sqrtQ-a/3;
      GSL_IMAG(*z1)=0;
      GSL_REAL(*z2)=2*sqrtQ-a/3;
      GSL_IMAG(*z2)=0;
    }
    return 3;

  } else if (CR2 < CQ3)  {

    /* equivalent to R2 < Q3 */

    double sqrtQ=sqrt(Q);
    double sqrtQ3=sqrtQ*sqrtQ*sqrtQ;
    double theta=acos(R/sqrtQ3);
    
    // Modified from the original GSL routine
    // Sometimes R/sqrtQ3 is slightly larger than one
    // when the coefficients are arranged just right
    // so theta becomes not finite.

    if (R/sqrtQ3>=1.0) theta=0.0;
    if (R/sqrtQ3<=-1.0) theta=o2scl_const::pi;

    double norm=-2*sqrtQ;
    double r0=norm*cos(theta/3)-a/3;
    double r1=norm*cos((theta+2.0*M_PI)/3)-a/3;
    double r2=norm*cos((theta-2.0*M_PI)/3)-a/3;

    /* Sort r0, r1, r2 into increasing order */

    if (r0 > r1)
      swap(r0,r1);
    if (r1 > r2) {
      swap(r1,r2);
      if (r0 > r1)
	swap(r0,r1);
    }
    
    GSL_REAL(*z0)=r0;
    GSL_IMAG(*z0)=0;
    
    GSL_REAL(*z1)=r1;
    GSL_IMAG(*z1)=0;
    
    GSL_REAL(*z2)=r2;
    GSL_IMAG(*z2)=0;
    
    return 3;

  } else {

    double sgnR=(R>=0?1:-1);
    
    double A;
    if (fabs(R)+sqrt(R2-Q3)<0.0) A=0.0;
    else A=-sgnR*pow(fabs(R)+sqrt(R2-Q3),1.0/3.0);
    
    // Modification from original GSL behavior: Just in case R2=Q3,
    // finite precision can cause the argument of the sqrt() function
    // to be negative. We correct for this here.
    if (R2<=Q3) {
      A=-sgnR*pow(fabs(R),1.0/3.0);
    }

    double B=Q/A;

    if (A+B < 0) {
      GSL_REAL(*z0)=A+B-a/3;
      GSL_IMAG(*z0)=0;
      
      GSL_REAL(*z1)=-0.5*(A+B)-a/3;
      GSL_IMAG(*z1)=-(sqrt(3.0)/2.0)*fabs(A-B);
      
      GSL_REAL(*z2)=-0.5*(A+B)-a/3;
      GSL_IMAG(*z2)=(sqrt(3.0)/2.0)*fabs(A-B);
    } else {
      GSL_REAL(*z0)=-0.5*(A+B)-a/3;
      GSL_IMAG(*z0)=-(sqrt(3.0)/2.0)*fabs(A-B);
      
      GSL_REAL(*z1)=-0.5*(A+B)-a/3;
      GSL_IMAG(*z1)=(sqrt(3.0)/2.0)*fabs(A-B);
      
      GSL_REAL(*z2)=A+B-a/3;
      GSL_IMAG(*z2)=0;
    }   
    return 3;
  }
}

int cubic_real_coeff_gsl::solve_rc
(const double a3, const double b3, const double c3, 
 const double d3, double &x1, std::complex<double> &x2, 
 std::complex<double> &x3) {
  double s1, s2, s3;
  gsl_complex r1, r2, r3;
  complex<double> i(0.0,1.0);

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_real_coeff_gsl::solve_rc().",
       exc_einval);
  }

  gsl_poly_complex_solve_cubic2(b3/a3,c3/a3,d3/a3,&r1,&r2,&r3);
  s1=fabs(GSL_IMAG(r1)/GSL_REAL(r1));
  s2=fabs(GSL_IMAG(r2)/GSL_REAL(r2));
  s3=fabs(GSL_IMAG(r3)/GSL_REAL(r3));
  if (s1<s2 && s1<s3) {
    x1=GSL_REAL(r1);
    x2=GSL_REAL(r2)+i*GSL_IMAG(r2);
    x3=GSL_REAL(r3)+i*GSL_IMAG(r3);
  } else if (s2<s1 && s2<s3) {
    x1=GSL_REAL(r2);
    x2=GSL_REAL(r1)+i*GSL_IMAG(r1);
    x3=GSL_REAL(r3)+i*GSL_IMAG(r3);
  } else {
    x1=GSL_REAL(r3);
    x2=GSL_REAL(r1)+i*GSL_IMAG(r1);
    x3=GSL_REAL(r2)+i*GSL_IMAG(r2);
  }

  return success;
}

int quartic_real_simple::solve_r
(const double a4, const double b4, const double c4, const double d4, 
 const double e4, double &x1, double &x2, double &x3, double &x4) {

  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_simple::solve_r().",
       exc_einval);
  }

  double a34, a24, a14, a04, a23, a13, a03;
  double u1, u2, u3;
  double u4, t1;
  double b2a, b2b, c2a, c2b;
  
  a34=b4/a4;
  a24=c4/a4;
  a14=d4/a4;
  a04=e4/a4;
  
  //---------------------------------------
  // Solve the resolvent cubic:
  
  a23=-a24;
  a13=(a14*a34-4.0*a04);
  a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);
  
  gsl_poly_solve_cubic(a23,a13,a03,&u1,&u2,&u3);
  
  u4=u2;
  
  //---------------------------------------
  // Now construct the two quadratics:
  
  t1=u4+a34*a34/4.0-a24;
  t1=sqrt(t1);
  
  b2a=-t1+a34/2.0;
  b2b=t1+a34/2.0;
  t1=u4*u4/4.0;
  
  // A temporary hack that fixes problems when numerical errors
  // make t1 slightly smaller that a04.
  if (fabs((u4*u4/4.0-a04)/a04)<cube_root_tol) {
    t1=0.0;
  } else {
    t1=sqrt(t1-a04);
  }
  
  c2a=u4/2.0-t1;
  c2b=u4/2.0+t1;
  
  if (fabs((b2a*c2b+c2a*b2b-d4)/d4)>1.0e-4) {
    t1=u4+a34*a34/4.0-a24;
    t1=-sqrt(t1);
    
    b2a=-t1+a34/2.0;
    b2b=t1+a34/2.0;
    
    t1=u4*u4/4.0;
    if (fabs((u4*u4/4.0-a04)/a04)<1.0e-6) {
      t1=0.0;
    } else {
      t1=sqrt(t1-a04);
    }
    c2a=u4/2.0-t1;
    c2b=u4/2.0+t1;
  }
  
  //---------------------------------------
  // The solutions to the two quadratics:
  
  gsl_poly_solve_quadratic(1.0,b2a,c2a,&x1,&x2);
  gsl_poly_solve_quadratic(1.0,b2b,c2b,&x3,&x4);

  return success;
}

int quartic_real_gsl::solve_r
(const double a4, const double b4, const double c4, 
 const double d4, const double e4, double &x1, double &x2, 
 double &x3, double &x4) {
  double a34, a24, a14, a04, a23, a13, a03;
  gsl_complex u1, u2, u3;
  gsl_complex u4, t1, t2, imagi;
  gsl_complex b2a, b2b, c2a, c2b;
  
  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_gsl::solve_r().",
       exc_einval);
  }

  GSL_REAL(imagi)=0.0;
  GSL_IMAG(imagi)=1.0;
  
  a34=b4/a4;
  a24=c4/a4;
  a14=d4/a4;
  a04=e4/a4;
  
  //---------------------------------------
  // Solve the resolvent cubic:
  
  a23=-a24;
  a13=(a14*a34-4.0*a04);
  a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);
  
  gsl_poly_complex_solve_cubic(a23,a13,a03,&u1,&u2,&u3);
  
  GSL_SET_COMPLEX(&u4,GSL_REAL(u2),GSL_IMAG(u2));
  
  //---------------------------------------
  // Now construct the two quadratics:
  
  t1=gsl_complex_add_real(u4,a34*a34/4.0-a24);
  t1=gsl_complex_sqrt(t1);
  
  b2a=gsl_complex_add_real(gsl_complex_negative(t1),a34/2.0);
  b2b=gsl_complex_add_real(t1,a34/2.0);
  
  t1=gsl_complex_mul_real(gsl_complex_mul(u4,u4),0.25);
  t1=gsl_complex_sqrt(gsl_complex_add_real(t1,-a04));
  c2a=gsl_complex_sub(gsl_complex_div_real(u4,2.0),t1);
  c2b=gsl_complex_add(gsl_complex_div_real(u4,2.0),t1);
  
  if (fabs((GSL_REAL(b2a)*GSL_REAL(c2b)+
	    GSL_REAL(c2a)*GSL_REAL(b2b)-d4)/d4)>cube_root_tol) {
    t1=gsl_complex_add_real(u4,a34*a34/4.0-a24);
    t1=gsl_complex_sqrt(t1);
    t2=gsl_complex_exp(gsl_complex_mul_real(imagi,pi));
    t1=gsl_complex_mul(t1,t2);
    
    b2a=gsl_complex_add_real(gsl_complex_negative(t1),a34/2.0);
    b2b=gsl_complex_add_real(t1,a34/2.0);
    
    t1=gsl_complex_mul_real(gsl_complex_mul(u4,u4),0.25);
    t1=gsl_complex_sqrt(gsl_complex_add_real(t1,-a04));
    c2a=gsl_complex_sub(gsl_complex_div_real(u4,2.0),t1);
    c2b=gsl_complex_add(gsl_complex_div_real(u4,2.0),t1);
    
  }
  
  //---------------------------------------
  // The solutions to the two quadratics:
    
  gsl_poly_solve_quadratic(1.0,GSL_REAL(b2a),GSL_REAL(c2a),&x1,&x2);
  gsl_poly_solve_quadratic(1.0,GSL_REAL(b2b),GSL_REAL(c2b),&x3,&x4);

  return success;
}

int quartic_real_gsl2::solve_r
(const double a4, const double b4, const double c4, 
 const double d4, const double e4, double &x1, double &x2, 
 double &x3, double &x4) {
 double q3, r3;
  double a34, a24, a14, a04, a23, a13, a03;
  gsl_complex imagi, t1, t2;
  gsl_complex b2a, c2a, b2b, c2b, dis;
  gsl_complex s1, s2, z1;
    
  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_gsl2::solve_r().",
       exc_einval);
  }

  GSL_REAL(imagi)=0.0;
  GSL_IMAG(imagi)=1.0;
    
  a34=b4/a4;
  a24=c4/a4;
  a14=d4/a4;
  a04=e4/a4;
    
  //---------------------------------------
  // Solve the resolvent cubic:
  
  a23=-a24;
  a13=(a14*a34-4.0*a04);
  a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);
    
  q3=a13/3.0-a23*a23/9.0;
  r3=(a13*a23-3.0*a03)/6.0-a23*a23*a23/27.0;
    
  dis=gsl_complex_sqrt_real(q3*q3*q3+r3*r3);
    
  s1=gsl_complex_add_real(dis,r3);
  s1=gsl_complex_pow_real(s1,1.0/3.0);
  s2=gsl_complex_negative(dis);
  s2=gsl_complex_add_real(s2,r3);
  s2=gsl_complex_pow_real(s2,1.0/3.0);
    
  z1=gsl_complex_add(s1,s2);
  z1=gsl_complex_sub_real(z1,a23/3.0);
    
  //---------------------------------------
  // Now construct the two quadratics:
    
  t1=gsl_complex_add_real(z1,a34*a34/4.0-a24);
  t1=gsl_complex_sqrt(t1);
  b2a=gsl_complex_add_real(gsl_complex_negative(t1),a34/2.0);
  b2b=gsl_complex_add_real(t1,a34/2.0);
    
  t1=gsl_complex_mul_real(gsl_complex_mul(z1,z1),0.25);
  t1=gsl_complex_sqrt(gsl_complex_add_real(t1,-a04));
  c2a=gsl_complex_sub(gsl_complex_div_real(z1,2.0),t1);
  c2b=gsl_complex_add(gsl_complex_div_real(z1,2.0),t1);
    
    
  //---------------------------------------
  // Check to see that we have made the correct choice for
  // the square roots
    
  t1=gsl_complex_mul(b2a,c2b);
  t1=gsl_complex_add(t1,gsl_complex_mul(b2b,c2a));
    
  if (fabs((GSL_REAL(t1)-a14)/a14)>cube_root_tol) {
      
    t1=gsl_complex_add_real(z1,a34*a34/4.0-a24);
    t1=gsl_complex_sqrt(t1);
    b2a=gsl_complex_add_real(t1,a34/2.0);
    b2b=gsl_complex_add_real(gsl_complex_negative(t1),a34/2.0);
      
    t1=gsl_complex_mul_real(gsl_complex_mul(z1,z1),0.25);
    t1=gsl_complex_sqrt(gsl_complex_add_real(t1,-a04));
    c2a=gsl_complex_sub(gsl_complex_div_real(z1,2.0),t1);
    c2b=gsl_complex_add(gsl_complex_div_real(z1,2.0),t1);
      
  }
    
  //---------------------------------------
  // Solve the quadratics
    
  t1=gsl_complex_mul(b2a,b2a);
  t1=gsl_complex_sub(t1,gsl_complex_mul_real(c2a,4.0));
  t1=gsl_complex_div_real(gsl_complex_sqrt(t1),2.0);
  t2=gsl_complex_div_real(gsl_complex_negative(b2a),2.0);
  x1=GSL_REAL(gsl_complex_add(t1,t2));
  x2=GSL_REAL(gsl_complex_sub(t2,t1));
    
  t1=gsl_complex_mul(b2b,b2b);
  t1=gsl_complex_sub(t1,gsl_complex_mul_real(c2b,4.0));
  t1=gsl_complex_div_real(gsl_complex_sqrt(t1),2.0);
  t2=gsl_complex_div_real(gsl_complex_negative(b2b),2.0);
  x3=GSL_REAL(gsl_complex_add(t1,t2));
  x4=GSL_REAL(gsl_complex_sub(t2,t1));
    
  return success;
}
