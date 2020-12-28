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

double quartic_real::disc_r(const double a, const double b, const double c, 
			    const double d, const double e) {
  double a2=a*a;
  double b2=b*b;
  double c2=c*c;
  double d2=d*d;
  double e2=e*e;
  
  double a3=a2*a;
  double b3=b2*b;
  double c3=c2*c;
  double d3=d2*d;
  double e3=e2*e;
  
  double b4=b2*b2;
  double c4=c2*c2;
  double d4=d2*d2;
  
  return 256.0*a3*e3-192.0*a2*b*d*e2-128.0*a2*c2*e2+144.0*a2*c*d2*e-
    27.0*a2*d4+144.0*a*b2*c*e2-6.0*a*b2*d2*e-80.0*a*b*c2*d*e+
    18.0*a*b*c*d3+16.0*a*c4*e-4.0*a*c3*d2-27.0*b4*e2+18.0*b3*c*d*e-
    4.0*b3*d3-4.0*b2*c3*e+b2*c2*d2;
}

int quadratic_complex_std::solve_c
(complex<double> a2, complex<double> b2, complex<double> c2,
 complex<double> &x1, complex<double> &x2) { 
  
  if (a2==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quadratic_complex_std::solve_c().",
       exc_einval);
  }
  x1=(-b2+sqrt(b2*b2-4.0*a2*c2))/2.0/a2;
  x2=(-b2-sqrt(b2*b2-4.0*a2*c2))/2.0/a2;
  return success;
}

int quadratic_real_coeff_gsl::solve_rc
(double a2, double b2, double c2,
 complex<double> &x1, complex<double> &x2) { 
  gsl_complex r1, r2;

  gsl_poly_complex_solve_quadratic(a2,b2,c2,&r1,&r2);
  
  x1=complex<double>(GSL_REAL(r1),GSL_IMAG(r1));
  x2=complex<double>(GSL_REAL(r2),GSL_IMAG(r2));

  return success;
}

inline double sign(double a, double b) {
  if (b>=0.0) return fabs(a);
  return -fabs(a);
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

int quartic_real_coeff_cern::solve_rc
(const double a4, const double b4, const double c4, const double d4, 
 const double e4, std::complex<double> &x1, std::complex<double> &x2, 
 std::complex<double> &x3, std::complex<double> &x4) {

  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_coeff_cern::solve_rc().",
       exc_einval);
  }

  int mt;
  double dc;
  complex<double> x[4];
  complex<double> i(0.0,1.0);

  rrteq4(b4/a4,c4/a4,d4/a4,e4/a4,x,dc,mt);
  x1=x[0];
  x2=x[1];
  x3=x[2];
  x4=x[3];

  return success;
}

// There are a couple differences with the original routine.
// The arrays z[] and u[] are now zero-indexed.
int quartic_real_coeff_cern::rrteq4(double a, double b, double c, double d, 
				    complex<double> z[], double &dc, 
				    int &mt) {
  complex<double> i(0.0,1.0), z0[5];
  complex<double> w1(0.0,0.0), w2(0.0,0.0), w3;
  double r4=1.0/4.0, r12=1.0/12.0;
  double q2=1.0/2.0, q4=1.0/4.0, q8=1.0/8.0;
  double q1=3.0/8.0, q3=3.0/16.0;
  double u[3], v[4], v1, v2;
  int j, k1=0, k2=0;
  
  // degenerate cases
  if (b==0 && c==0) {
    if (d==0) {
      mt=1;
      z[0]=-a;
      z[1]=0;
      z[2]=0;
      z[3]=0;
      dc=0;
      return success;
    } else if (a==0) {
      if (d>0) {
	mt=2;
	z[0]=sqrt(i*sqrt(d));
	z[1]=-z[0];
	z[3]=sqrt(-z[0]*z[0]);
	z[2]=-z[3];
      } else {
	mt=3;
	z[0]=sqrt(sqrt(-d));
	z[1]=-z[0];
	z[2]=sqrt(-z[0]*z[0]);
	z[3]=-z[2];
      }
      dc=-r12*d*r12*d*r12*d;
      return success;
    }
  }

  // Solve the resolvant cubic
  double aa=a*a;
  double pp=b-q1*aa;
  double qq=c-q2*a*(b-q4*aa);
  double rr=d-q4*(a*c-q4*aa*(b-q3*aa));
  double rc=q2*pp;
  double sc=q4*(q4*pp*pp-rr);
  double tc=-(q8*qq*q8*qq);
  
  cub_obj.rrteq3(rc,sc,tc,u,dc);
  
  double q=qq;
  double h=r4*a;
  if (dc==0) u[2]=u[1];
  if (dc<=0) {
    mt=2;
    v[1]=fabs(u[0]);
    v[2]=fabs(u[1]);
    v[3]=fabs(u[2]);
    v1=max(max(v[1],v[2]),v[3]);
    if (v1==v[1]) {
      k1=0;
      v2=max(v[2],v[3]);
    } else if (v1==v[2]) {
      k1=1;
      v2=max(v[1],v[3]);
    } else {
      k1=2;
      v2=max(v[1],v[2]);
    }
    if (v2==v[1]) {
      k2=0;
    } else if (v2==v[2]) {
      k2=1;
    } else {
      k2=2;
    }
    w1=sqrt(((complex<double>)(u[k1])));
    w2=sqrt(((complex<double>)(u[k2])));
  } else {
    mt=3;
    w1=sqrt(u[1]+i*u[2]);
    w2=sqrt(u[1]-i*u[2]);
  }
  w3=0;
  if (w1*w2!=0.0) w3=-q/(8.0*w1*w2);
  z0[1]=w1+w2+w3-h;
  z0[2]=-w1-w2+w3-h;
  z0[3]=-w1+w2-w3-h;
  z0[4]=w1-w2-w3-h;
  if (mt==2) {
    if (u[k1]>=0 && u[k2]>=0) {
      mt=1;
      for(j=1;j<=4;j++) {
	z[j-1]=z0[j].real();
      }
    } else if (u[k1]>=0 && u[k2]<0) {
      z[0]=z0[1];
      z[1]=z0[4];
      z[2]=z0[3];
      z[3]=z0[2];
    } else if (u[k1]<0 && u[k2]>=0) {
      z[0]=z0[1];
      z[1]=z0[3];
      z[2]=z0[4];
      z[3]=z0[2];
    } else if (u[k1]<0 && u[k2]<0) {
      z[0]=z0[1];
      z[1]=z0[2];
      z[2]=z0[4];
      z[3]=z0[3];
    }
  } else if (mt==3) {
    for(j=1;j<=2;j++) {
      z[j-1]=z0[j].real();
    }
    z[2]=z0[4];
    z[3]=z0[3];
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

int cubic_real_coeff_cern::solve_rc
(const double a3, const double b3, const double c3, const double d3, 
 double &r1, complex<double> &r2, complex<double> &r3) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_real_coeff_cern::solve_rc().",
       exc_einval);
  }

  double x[3],d;
  complex<double> i(0.0,1.0);

  rrteq3(b3/a3,c3/a3,d3/a3,x,d);
  if (d>0.0) {
    r1=x[0];
    r2=x[1]+i*x[2];
    r3=x[1]-i*x[2];
  } else {
    r1=x[0];
    r2=x[1];
    r3=x[2];
  }

  return success;
}

int cubic_real_coeff_cern::rrteq3(double r, double s, double t, 
				  double x[], double &d) {

  double delta2=delta;
  double r1=2.0/27.0, r2=0.5, r3=1.0/3.0;
  double w3=sqrt(3.0), r4=w3/2.0;
  double q1=2.0/27.0, q2=0.5, q3=1.0/3.0;
  double y[3];
  complex<double> z[3], i(0.0,1.0);
  double h2, h3;
  int j,k;
  
  if (s==0.0 && t==0.0) {
    x[0]=-r;
    x[1]=0.0;
    x[2]=0.0;
    d=0;
    return success;
  }
  double p=s-r3*r*r;
  double q=(r1*r*r-r3*s)*r+t;
  d=r2*r2*q*q+r3*p*r3*p*r3*p;
  if (fabs(d)<=eps) {
    double pp=s-q3*r*r;
    double qq=(q1*r*r-q3*s)*r+t;
    d=q2*q2*qq*qq+q3*pp*q3*pp*q3*pp;
    p=pp;
    q=qq;
  }
  double h=r3*r;
  double h1=r2*q;
  double u,v,d_new;

  // The discriminant in 'd' has units of [x]^6 so it is very
  // sensitive to the absolute magnitude of the roots. We attempt to
  // fix this by using the ratio instead of the sum.
  if (improve_scale) {
    double da=r2*r2*q*q;
    double db=r3*p*r3*p*r3*p;
    if (db==0.0) {
      delta2=0.0;
      d_new=da;
    } else if (db>0.0) {
      d_new=da/db+1.0;
    } else {
      d_new=-da/db-1.0;
    }
  } else {
    d_new=d;
  }

  if (d_new>delta2) {
    h2=sqrt(d);
    double u0=-h1+h2;
    double v0=-h1-h2;
    if (fabs(u0)==0.0) u=sign(0.0,u0);
    else u=sign(pow(fabs(u0),r3),u0);
    if (fabs(v0)==0.0) v=sign(0.0,v0);
    else v=sign(pow(fabs(v0),r3),v0);
    x[0]=u+v-h;
    x[1]=-r2*(u+v)-h;
    x[2]=r4*fabs(u-v);
    if (fabs(u0)<=eps || fabs(v0)<=eps) {
      y[0]=x[0];
      for(k=0;k<=1;k++) {
	y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3.0*y[k]+2.0*r)*y[k]+s);
      }
      x[0]=y[2];
     z[0]=x[1]+i*x[2];
      for(k=0;k<=1;k++) {
	z[k+1]=z[k]-(((z[k]+r)*z[k]+s)*z[k]+t)/((3.0*z[k]+2.0*r)*z[k]+s);
      }
      x[1]=z[2].real();
      x[2]=z[2].imag();
    }

  } else if (fabs(d_new)<=delta2) {

    d=0.0;
    if (fabs(h1)==0.0) u=sign(0.0,-h1);
    else u=sign(pow(fabs(h1),r3),-h1);
    x[0]=u+u-h;
    x[1]=-u-h;
    x[2]=x[1];
    if (fabs(h1)<=eps) {
      y[0]=x[0];
      for(k=0;k<=1;k++) {
	h1=(3.0*y[k]+2.0*r)*y[k]+s;
	if (fabs(h1)>delta2) {
	  y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/h1;
	} else {
	  x[0]=-r3*r;
	  x[1]=x[0];
	  x[2]=x[0];
	  return success;
	}
      }
      x[0]=y[2];
      x[1]=-r2*(r+x[0]);
      x[2]=x[1];
    }

  } else {

    h3=fabs(r3*p);
    h3=sqrt(h3*h3*h3);
    h2=r3*acos(-h1/h3);
    if (h3==0.0) h1=0.0;
    else h1=pow(h3,r3);
    u=h1*cos(h2);
    v=w3*h1*sin(h2);
    x[0]=u+u-h;
    x[1]=-u-v-h;
    x[2]=-u+v-h;
    if (h3<=eps || x[0]<=eps || x[1]<=eps || x[2]<=eps) {
      for(j=0;j<3;j++) {
	y[0]=x[j];
	for(k=0;k<=1;k++) {
	  y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/((3.0*y[k]+2.0*r)*y[k]+s);
	}
	x[j]=y[2];
      }
    }
  }

  return success;
}

int cubic_complex_std::solve_c
(const std::complex<double> a3, const std::complex<double> b3, 
 const std::complex<double> c3, const std::complex<double> d3, 
 std::complex<double> &x1, std::complex<double> &x2, 
 std::complex<double> &x3) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_std::complex::solve_c().",
       exc_einval);
  }

  complex<double> p3, q3, mo;
  complex<double> alpha, beta, cbrta, cbrtb;
  complex<double> e2, e4, cacb;
  double test, re_p3;

  if (a3==0.0) {
    quadratic_complex_std qsc;
    qsc.solve_c(b3,c3,d3,x1,x2);
    x3=0.0;
    return success;
  }

  mo=-1.0;
  mo=sqrt(mo);

  p3=(3.0*a3*c3-b3*b3)/9.0/a3/a3;
  q3=(2.0*b3*b3*b3-9.0*a3*b3*c3+27.0*a3*a3*d3)/27.0/a3/a3/a3;

  alpha=(-q3+sqrt(q3*q3+4.0*p3*p3*p3))/2.0;
  beta=(q3+sqrt(q3*q3+4.0*p3*p3*p3))/2.0;
  
  if (alpha.real()==0.0) cbrta=0.0;
  else cbrta=pow(alpha,1.0/3.0);
  if (beta.real()==0.0) cbrtb=0.0;
  else cbrtb=pow(beta,1.0/3.0);

  // It seems that if the real part of alpha is < 0 and the imaginary
  // part is zero, then cbrta is NaN (esp. w/Cygwin). We fix this
  // here:
  if (!std::isfinite(cbrta.real())) {
    cbrta=pow(-alpha,1.0/3.0)*exp(mo*pi/3.0);
  }
  if (!std::isfinite(cbrtb.real())) {
    cbrtb=pow(-beta,1.0/3.0)*exp(mo*pi/3.0);
  }

  e2=exp(mo*2.0*pi/3.0);
  e4=exp(mo*4.0*pi/3.0);

  // This next section is nessary to ensure that the code
  // selects the correct cube roots of alpha and beta.
  // I changed this because I wanted code that had no chance
  // of accidentally falling into an infinite loop.
  re_p3=p3.real();
  cacb=cbrta*cbrtb;
  test=fabs((cacb.real()-re_p3)/re_p3);
  if (fabs(((cacb*e2).real()-re_p3)/re_p3)<test) {
    cbrta*=e2;
    cacb=cbrta*cbrtb;
    test=fabs((cacb.real()-re_p3)/re_p3);
    if (fabs(((cacb*e2).real()-re_p3)/re_p3)<test) {
      cbrta*=e2;
    }
  } else {
    if (fabs(((cacb*e2*e2).real()-re_p3)/re_p3)<test) {
      cbrta*=e2*e2;
    }
  }

  x1=cbrta-cbrtb-b3/3.0/a3;
  x2=cbrta*e2-cbrtb*e4-b3/3.0/a3;
  x3=cbrta*e4-cbrtb*e2-b3/3.0/a3;

  return success;
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

int quartic_complex_simple::solve_c
(const complex<double> a4, const complex<double> b4, const complex<double> c4,
 const complex<double> d4, const complex<double> e4, complex<double> &x1,
 complex<double> &x2, complex<double> &x3, complex<double> &x4) {
  complex<double> p4, q4, r4;
  complex<double> a3, b3, c3, d3;
  complex<double> b2a, c2a, b2b, c2b;
  complex<double> u4, u41, u42;
  
  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_complex_simple::solve_c().",
       exc_einval);
  }

  p4=(8.0*a4*c4-3.0*b4*b4)/8.0/a4/a4;
  q4=(b4*b4*b4-4.0*a4*b4*c4+8.0*a4*a4*d4)/8.0/(a4*a4*a4);
  r4=(16.0*a4*b4*b4*c4+256.0*a4*a4*a4*e4-3.0*b4*b4*b4*b4-64.0*a4*a4*b4*d4)/
    256.0/(a4*a4*a4*a4);

  //---------------------------------------
  // Solve the resolvent cubic:

  a3=1.0;
  b3=-p4;
  c3=-4.0*r4;
  d3=4.0*p4*r4-q4*q4;

  cub_obj.solve_c(a3,b3,c3,d3,u4,u41,u42);

  //---------------------------------------

  // What to do when u4==p4?
  // Temporary hack:
  if (u4==p4) {
    b2a=0.0;
    b2b=0.0;
    c2a=u4/2.0;
    c2b=u4/2.0;
  } else {
    b2a=sqrt(u4-p4);
    b2b=-sqrt(u4-p4);
    c2a=-sqrt(u4-p4)*q4/2.0/(u4-p4)+u4/2.0;
    c2b=sqrt(u4-p4)*q4/2.0/(u4-p4)+u4/2.0;
  }

  x1=(-b2a+sqrt(b2a*b2a-4.0*c2a))/2.0-b4/4.0/a4;
  x2=(-b2a-sqrt(b2a*b2a-4.0*c2a))/2.0-b4/4.0/a4;
  x3=(-b2b+sqrt(b2b*b2b-4.0*c2b))/2.0-b4/4.0/a4;
  x4=(-b2b-sqrt(b2b*b2b-4.0*c2b))/2.0-b4/4.0/a4;
  
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

int quadratic_real_coeff::solve_r
(const double a3, const double b3, const double c3, 
 double &x1, double &x2) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quadratic_real_coeff::solve_r().",
       exc_einval);
  }

  complex<double> r1,r2,r3;
  int ret=solve_rc(a3,b3,c3,r1,r2);
  x1=r1.real();
  x2=r2.real();
  return ret;
}

int quadratic_complex::solve_r
(const double a3, const double b3, const double c3, 
 double &x1, double &x2) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quadratic_complex::solve_r().",
       exc_einval);
  }

  complex<double> r1,r2;
  int ret=solve_c(a3,b3,c3,r1,r2);
  x1=r1.real();
  x2=r2.real();
  return ret;
}

int quadratic_complex::solve_rc
(const double a3, const double b3, const double c3, 
 complex<double> &x1, complex<double> &x2) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quadratic_complex::solve_rc().",
       exc_einval);
  }

  int ret=solve_c(a3,b3,c3,x1,x2);
  return ret;
}

int cubic_real_coeff::solve_r
(const double a3, const double b3, const double c3, 
 const double d3, double &x1, double &x2, double &x3) {
  complex<double> r2,r3;

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_real_coeff::solve_r().",
       exc_einval);
  }

  int ret=solve_rc(a3,b3,c3,d3,x1,r2,r3);
  x2=r2.real();
  x3=r3.real();
  return ret;
}

int cubic_complex::solve_r
(const double a3, const double b3, const double c3, 
 const double d3, double &x1, double &x2, double &x3) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_complex::solve_r().",
       exc_einval);
  }

  complex<double> r1,r2,r3;
  int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
  x1=r1.real();
  x2=r2.real();
  x3=r3.real();
  return ret;
}

int cubic_complex::solve_rc
(const double a3, const double b3, const double c3, 
 const double d3, double &x1, complex<double> &x2, 
 complex<double> &x3) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in cubic_complex::solve_rc().",
       exc_einval);
  }

  complex<double> r1,r2,r3;
  int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
  double s1,s2,s3;
  s1=fabs(r1.imag()/r1.real());
  s2=fabs(r2.imag()/r2.real());
  s3=fabs(r3.imag()/r3.real());
  if (s1<s2 && s1<s3) {
    x1=r1.real();
    x2=r2;
    x3=r3;
  } else if (s2<s1 && s2<s3) {
    x1=r2.real();
    x2=r1;
    x3=r3;
  } else {
    x1=r3.real();
    x2=r1;
    x3=r2;
  }
  return ret;
}

int quartic_real_coeff::solve_r
(const double a3, const double b3, const double c3, 
 const double d3, const double e3, double &x1, double &x2, double &x3,
 double &x4) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_coeff::solve_r().",
       exc_einval);
  }

  complex<double> r1,r2,r3,r4;
  int ret=solve_rc(a3,b3,c3,d3,e3,r1,r2,r3,r4);
  x1=r1.real();
  x2=r2.real();
  x3=r3.real();
  x4=r4.real();
  return ret;
}

int quartic_complex::solve_rc
(const double a4, const double b4, const double c4, 
 const double d4, const double e4, 
 std::complex<double> &x1, std::complex<double> &x2, 
 std::complex<double> &x3, std::complex<double> &x4) {

  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in std::complex<double> &x4) {().",
       exc_einval);
  }

  return solve_c(a4,b4,c4,d4,e4,x1,x2,x3,x4);
}

int quartic_complex::solve_r
(const double a3, const double b3, const double c3, 
 const double d3, const double e3, double &x1, double &x2, double &x3,
 double &x4) {

  if (a3==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_complex::solve_r().",
       exc_einval);
  }

  complex<double> r1,r2,r3,r4;
  int ret=solve_c(a3,b3,c3,d3,e3,r1,r2,r3,r4);
  x1=r1.real();
  x2=r2.real();
  x3=r3.real();
  x4=r4.real();
  return ret;
}

