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

inline double max(double a, double b) {
  if (a>b) return a;
  return b;
}

int quadratic_real_gsl::solve_r(const double a2, const double b2,
                                const double c2, 
                                double &x1, double &x2) {
  return gsl_poly_solve_quadratic(a2,b2,c2,&x1,&x2);
}

int quadratic_real_coeff_gsl::solve_rc(const double a2, const double b2,
                                       const double c2, 
                                       std::complex<double> &x1,
                                       std::complex<double> &x2) {
  gsl_complex z0, z1;
  int ret=gsl_poly_complex_solve_quadratic(a2,b2,c2,&z0,&z1);
  x1.real(GSL_REAL(z0));
  x1.imag(GSL_IMAG(z0));
  x2.real(GSL_REAL(z1));
  x2.imag(GSL_IMAG(z1));
  return ret;
}

int cubic_real_coeff_gsl::solve_rc(const double a3, const double b3,
                                   const double c3, const double d3,
                                   double &x1,
                                   std::complex<double> &x2,
                                   std::complex<double> &x3) {
  
  gsl_complex z0, z1, z2;
  int ret=gsl_poly_complex_solve_cubic(b3/a3,c3/a3,d3/a3,&z0,&z1,&z2);
  if (GSL_IMAG(z0)==0.0) {
    x1=GSL_REAL(z0);
    x2.real(GSL_REAL(z1));
    x2.imag(GSL_IMAG(z1));
    x3.real(GSL_REAL(z2));
    x3.imag(GSL_IMAG(z2));
  } else if (GSL_IMAG(z1)==0.0) {
    x1=GSL_REAL(z1);
    x2.real(GSL_REAL(z0));
    x2.imag(GSL_IMAG(z0));
    x3.real(GSL_REAL(z2));
    x3.imag(GSL_IMAG(z2));
  } else if (GSL_IMAG(z2)==0.0) {
    x1=GSL_REAL(z2);
    x2.real(GSL_REAL(z0));
    x2.imag(GSL_IMAG(z0));
    x3.real(GSL_REAL(z1));
    x3.imag(GSL_IMAG(z1));
  } else {
    O2SCL_ERR("GSL returned three complex roots.",o2scl::exc_einval);
  }
  return ret;
}

int quartic_real_std::solve_r
(const double a4, const double b4, const double c4, const double d4, 
 const double e4, double &x1, double &x2, double &x3, double &x4) {

  if (a4==0.0) {
    O2SCL_ERR
      ("Leading coefficient zero in quartic_real_std::solve_r().",
       exc_einval);
  }

  double a34=b4/a4;
  double a24=c4/a4;
  double a14=d4/a4;
  double a04=e4/a4;
  
  //---------------------------------------
  // Solve the resolvent cubic:
  
  double a23=-a24;
  double a13=(a14*a34-4.0*a04);
  double a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);

  double u1, u2, u3;
  cub2.solve_r(1,a23,a13,a03,&u1,&u2,&u3);
  //gsl_poly_solve_cubic(a23,a13,a03,&u1,&u2,&u3);
  
  double u4=u2;
  
  //---------------------------------------
  // Now construct the two quadratics:
  
  double t1=u4+a34*a34/4.0-a24;
  t1=sqrt(t1);
  
  double b2a=-t1+a34/2.0;
  double b2b=t1+a34/2.0;
  t1=u4*u4/4.0;
  
  // A temporary hack that fixes problems when numerical errors
  // make t1 slightly smaller that a04.
  /*
  if (fabs((u4*u4/4.0-a04)/a04)<cube_root_tol) {
    t1=0.0;
  } else {
    t1=sqrt(t1-a04);
  }
  */
  if (t1>a04) {
    t1=sqrt(t1-a04);
  } else {
    t1=0;
  }
  
  double c2a=u4/2.0-t1;
  double c2b=u4/2.0+t1;
  
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

  quad2.solve_r(1,b2a,c2a,&x1,&x2);
  quad2.solve_r(1,b2b,c2b,&x3,&x4);
  
  //gsl_poly_solve_quadratic(1.0,b2a,c2a,&x1,&x2);
  //gsl_poly_solve_quadratic(1.0,b2b,c2b,&x3,&x4);

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
