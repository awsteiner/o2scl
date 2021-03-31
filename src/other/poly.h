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
#ifndef O2SCL_POLY_H
#define O2SCL_POLY_H

/** \file poly.h
    \brief Classes for solving polynomials

    \warning 
    One must be careful about using pow() in functions using
    complex<double> since pow(((complex<double>)0.0),3.0) returns
    (nan,nan). Instead, we should use pow(((complex<double>)0.0),3)
    which takes an integer for the second argument. The sqrt()
    function, always succeeds i.e. sqrt(((complex<double>)0.0))=0.0

    \future The quartics are tested only for a4=1, which should
    probably be generalized.
*/

#include <iostream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Solve a quadratic polynomial with real coefficients and 
      real roots [abstract base]
  */
  class quadratic_real {
    
  public:

    virtual ~quadratic_real() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const double a2, const double b2, const double c2, 
			double &x1, double &x2)=0;

    /** \brief Compute the quadratic discriminant, \f$ b^2-4ac \f$
     */
    virtual double disc_r(double a2, double b2, double c2) {
      return b2*b2-4.0*a2*c2;
    }
    
    /// Return a string denoting the type ("quadratic_real")
    const char *type() { return "quadratic_real"; }
  };

  /** \brief Solve a quadratic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  class quadratic_real_coeff : public quadratic_real {

  public:

    virtual ~quadratic_real_coeff() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const double a2, const double b2, const double c2, 
			double &x1, double &x2);

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const double a2, const double b2, const double c2, 
			 std::complex<double> &x1, 
			 std::complex<double> &x2)=0;

    /// Return a string denoting the type ("quadratic_real_coeff")
    const char *type() { return "quadratic_real_coeff"; }
  };

  /** \brief Solve a quadratic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  class quadratic_complex : public quadratic_real_coeff {
  public:

    virtual ~quadratic_complex() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const double a2, const double b2, const double c2, 
			double &x1, double &x2);

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const double a2, const double b2, const double c2, 
			 std::complex<double> &x1, std::complex<double> &x2);

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const std::complex<double> a2, 
			const std::complex<double> b2, 
			const std::complex<double> c2, 
			std::complex<double> &x1, std::complex<double> &x2)=0;

    /// Return a string denoting the type ("quadratic_complex")
    const char *type() { return "quadratic_complex"; }
  };  

  /** \brief Solve a cubic polynomial with real coefficients and real roots
      [abstract base]
  */
  class cubic_real {
  public:

    virtual ~cubic_real() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const double a3, const double b3, const double c3, 
			const double d3, double &x1, double &x2, 
			double &x3)=0;

    /** \brief Compute the cubic discriminant, 
	\f$ b^2 c^2 - 4 a c^3 - 4 b^3 d - 27 a^2 d^2 + 18 a b c d \f$
     */
    virtual double disc_r(const double a3, const double b3, const double c3, 
			  const double d3) {
      return b3*b3*c3*c3-4.0*a3*c3*c3*c3-4.0*b3*b3*b3*d3-
	27.0*a3*a3*d3*d3+18.0*a3*b3*c3*d3;
    }

    /// Return a string denoting the type ("cubic_real")
    const char *type() { return "cubic_real"; }
  };

  /** \brief Solve a cubic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  class cubic_real_coeff : public cubic_real {

  public:

    virtual ~cubic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const double a3, const double b3, const double c3, 
			const double d3, double &x1, double &x2, double &x3);

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, std::complex<double> &x2,
			 std::complex<double> &x3)=0;

    /// Return a string denoting the type ("cubic_real_coeff")
    const char *type() { return "cubic_real_coeff"; }
  };

  /** \brief Solve a cubic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  class cubic_complex : public cubic_real_coeff {

  public:

    virtual ~cubic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const double a3, const double b3, const double c3, 
			const double d3, double &x1, double &x2, double &x3);

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, std::complex<double> &x2,
			 std::complex<double> &x3) {
      if (a3==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cubic_complex::solve_rc().",
           exc_einval);
      }
      
      std::complex<double> r1,r2,r3;
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
    
    /** \brief Solves the complex polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the three complex solutions \f$ x=x_1 \f$ , 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_c(const std::complex<double> a3, 
			const std::complex<double> b3, 
			const std::complex<double> c3, 
			const std::complex<double> d3, 
			std::complex<double> &x1, std::complex<double> &x2, 
			std::complex<double> &x3)=0;

    /// Return a string denoting the type ("cubic_complex")
    const char *type() { return "cubic_complex"; }
  };

  /** \brief Solve a quartic polynomial with real coefficients and 
      real roots [abstract base]
  */
  class quartic_real {

  public:

    virtual ~quartic_real() {}
    
    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, 
			double &x1, double &x2, 
			double &x3, double &x4)=0;

    /** \brief Compute the discriminant

	The discriminant is zero if and only if at least two roots are
	equal. If the discriminant is non-zero, the discriminant is
	negative if there are two real roots and two complex conjugate
	roots, and it is positive if the roots are either all real or
	all non-real.
    */
    virtual double disc_r(const double a, const double b, const double c, 
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
    
    /// Return a string denoting the type ("quartic_real")
    const char *type() { return "quartic_real"; }
  };

  /** \brief Solve a quartic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  class quartic_real_coeff : public quartic_real {

  public:

    virtual ~quartic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, double &x1, 
			double &x2, double &x3, double &x4) {
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff::solve_r().",
           exc_einval);
      }
      
      std::complex<double> r1,r2,r3,r4;
      int ret=solve_rc(a4,b4,c4,d4,e4,r1,r2,r3,r4);
      x1=r1.real();
      x2=r2.real();
      x3=r3.real();
      x4=r4.real();
      return ret;
    }

    /**  \brief Solves the polynomial 
	 \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	 giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	 \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const double a4, const double b4, const double c4, 
			 const double d4, const double e4, 
			 std::complex<double> &x1, std::complex<double> &x2, 
			 std::complex<double> &x3, 
			 std::complex<double> &x4)=0;

    /// Return a string denoting the type ("quartic_real_coeff")
    const char *type() { return "quartic_real_coeff"; }
  };

  /** \brief Solve a quartic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  class quartic_complex : public quartic_real_coeff {

  public:

    virtual ~quartic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, double &x1, 
			double &x2, double &x3, double &x4) {
			      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_complex::solve_r().",
           exc_einval);
      }
      
      std::complex<double> r1,r2,r3,r4;
      int ret=solve_c(a4,b4,c4,d4,e4,r1,r2,r3,r4);
      x1=r1.real();
      x2=r2.real();
      x3=r3.real();
      x4=r4.real();
      return ret;
    }

    /**  \brief Solves the polynomial 
	 \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	 giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	 \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const double a4, const double b4, const double c4, 
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

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const std::complex<double> a4, 
			const std::complex<double> b4, 
			const std::complex<double> c4, 
			const std::complex<double> d4, 
			const std::complex<double> e4, 
			std::complex<double> &x1, 
			std::complex<double> &x2, std::complex<double> &x3,
			std::complex<double> &x4)=0;

    /// Return a string denoting the type ("quartic_complex")
    const char *type() { return "quartic_complex"; }
  };

  /** \brief Solve a general polynomial with real
      coefficients and complex roots [abstract base]
  */
  class poly_real_coeff : public quadratic_real_coeff,
    public cubic_real_coeff, public quartic_real_coeff {

  public:
    
    virtual ~poly_real_coeff() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_rc_arr(int n, const double co[], 
			     std::complex<double> ro[])=0;

    /// Return a string denoting the type ("poly_real_coeff")
    const char *type() { return "poly_real_coeff"; }
    
  };

  /** \brief Solve a general polynomial with complex
      coefficients [abstract base]
  */
  class poly_complex : public quadratic_complex,
    public cubic_complex, public quartic_complex {

  public:

    virtual ~poly_complex() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_c_arr(int n, const std::complex<double> co[], 
			    std::complex<double> ro[])=0;
    
    /// Polish the roots 
    virtual int polish_c_arr(int n, const std::complex<double> co[],
			     std::complex<double> *ro)=0;

#ifdef O2SCL_NEVER_DEFINED
    {
      mroot_hybrids<> mh;
      ubvector x(2), y(2);
      // Polish all each one of n roots
      for(size_t j=0;j<n;j++) {
        x[0]=ro[0].real();
        x[1]=ro[0].imag();
        mh.msolve(2,x);
      }
      return 0;
    }
    
    virtual int polish_fun(size_t nv, const ubvector &x,
                           double ubvector &y,
                           std::complex<double> *co,
                           std::complex<double> *ro) {
      
      // Using horner's method following, e.g. GSL
      y[0]=co[nv-1];
      y[1]=0.0;
      for(int i=nv-1; i>0; i--) {
        double tmp=co[i-1]+x[0]*y[0]-x[1]*y[1];
        y[1]=x[1]*y[0]+x[0]*y[1];
        y[0]=tmp;
      }
      return 0;
    }
    
#endif

    /// Return a string denoting the type ("poly_complex")
    const char *type() { return "poly_complex"; }
    
  };

  /** \brief Solve a cubic with real coefficients and complex roots 
      (CERNLIB)

      \note The function rrteq3() is based on the CERNLIB routine of
      the same name, but differs slightly. See the documentation of
      that function for details.
  */
  template<class fp_t=double> class cubic_real_coeff_cern :
    public cubic_real_coeff {

  public:

    cubic_real_coeff_cern() {
      eps=1.0e-6;
      delta=1.0e-15;
      improve_scale=true;
    }

    /// Numerical tolerance (default \f$ 10^{-6} \f$)
    fp_t eps;

    /// Numerical tolerance (default \f$ 10^{-15} \f$)
    fp_t delta;

    /// Improve algorithm for poorly-scaled roots (default true)
    bool improve_scale;

    virtual ~cubic_real_coeff_cern() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &r1, 
			 std::complex<fp_t> &r2,
                         std::complex<fp_t> &r3) {
      if (a3==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cubic_real_coeff_cern::solve_rc().",
           exc_einval);
      }
      
      fp_t x[3],d;
      std::complex<fp_t> i(0.0,1.0);
      
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
    
    inline double sign(double a, double b) {
      if (b>=0.0) return fabs(a);
      return -fabs(a);
    }

    /** \brief The CERNLIB-like interface

	This function computes the roots of the cubic equation
	\f[
	x^3 + r x^2 + s x + t =0
	\f]
	returning the value of the discriminant in \c d and the roots
	in the array \c x. If the discriminant is negative, then all
	three real roots are stored in \c x. Otherwise, the real root
	is stored in <tt>x[0]</tt> and the real and imaginary parts of
	the complex conjugate roots are stored in <tt>x[1]</tt> and
	<tt>x[2]</tt>, respectively. This differs from the CERNLIB
	routine where the results were stored in <tt>x[1]</tt>,
	<tt>x[2]</tt>, and <tt>x[3]</tt> instead.
	
	Another small change is that the discriminant for the
	resolvent cubic is evaluated slightly differently in order to
	improve the properties in the case where the roots are not all
	of order unity. The default CERNLIB behavior can be restored
	by setting improve_scale to \c false.
    */	
    virtual int rrteq3(fp_t r, fp_t s, fp_t t, fp_t x[],
                       fp_t &d) {
      
      fp_t delta2=delta;
      fp_t r1=2.0/27.0, r2=0.5, r3=1.0/3.0;
      fp_t w3=sqrt(3.0), r4=w3/2.0;
      fp_t q1=2.0/27.0, q2=0.5, q3=1.0/3.0;
      fp_t y[3];
      std::complex<fp_t> z[3], i(0.0,1.0);
      fp_t h2, h3;
      int j,k;
      
      if (s==0.0 && t==0.0) {
        x[0]=-r;
        x[1]=0.0;
        x[2]=0.0;
        d=0;
        return success;
      }
      fp_t p=s-r3*r*r;
      fp_t q=(r1*r*r-r3*s)*r+t;
      d=r2*r2*q*q+r3*p*r3*p*r3*p;
      if (fabs(d)<=eps) {
        fp_t pp=s-q3*r*r;
        fp_t qq=(q1*r*r-q3*s)*r+t;
        d=q2*q2*qq*qq+q3*pp*q3*pp*q3*pp;
        p=pp;
        q=qq;
      }
      fp_t h=r3*r;
      fp_t h1=r2*q;
      fp_t u,v,d_new;
      
      // The discriminant in 'd' has units of [x]^6 so it is very
      // sensitive to the absolute magnitude of the roots. We attempt to
      // fix this by using the ratio instead of the sum.
      if (improve_scale) {
        fp_t da=r2*r2*q*q;
        fp_t db=r3*p*r3*p*r3*p;
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
        fp_t u0=-h1+h2;
        fp_t v0=-h1-h2;
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
    
    /// Return a string denoting the type ("cubic_real_coeff_cern")
    const char *type() { return "cubic_real_coeff_cern"; }
  };

  /** \brief Solve a quartic with real coefficients and complex 
      roots (CERNLIB)
  */
  template<class fp_t=double> class quartic_real_coeff_cern :
    public quartic_real_coeff {

  public:

    virtual ~quartic_real_coeff_cern() {}

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four complex solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 std::complex<fp_t> &x1, std::complex<fp_t> &x2, 
			 std::complex<fp_t> &x3, std::complex<fp_t> &x4) {
      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff_cern::solve_rc().",
           exc_einval);
      }
      
      int mt;
      fp_t dc;
      std::complex<fp_t> x[4];
      
      rrteq4(b4/a4,c4/a4,d4/a4,e4/a4,x,dc,mt);
      x1=x[0];
      x2=x[1];
      x3=x[2];
      x4=x[3];
      
      return success;
    }      

    /** \brief The CERNLIB-like interface

        There are a couple differences with the original routine.
        The arrays z[] and u[] are now zero-indexed.
    */
    virtual int rrteq4(fp_t a, fp_t b, fp_t c, fp_t d, 
		       std::complex<fp_t> z[], fp_t &dc, 
		       int &mt) {
      
      std::complex<fp_t> i(0.0,1.0), z0[5];
      std::complex<fp_t> w1(0.0,0.0), w2(0.0,0.0), w3;
      fp_t r4=1.0/4.0, r12=1.0/12.0;
      fp_t q2=1.0/2.0, q4=1.0/4.0, q8=1.0/8.0;
      fp_t q1=3.0/8.0, q3=3.0/16.0;
      fp_t u[3], v[4], v1, v2;
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
      fp_t aa=a*a;
      fp_t pp=b-q1*aa;
      fp_t qq=c-q2*a*(b-q4*aa);
      fp_t rr=d-q4*(a*c-q4*aa*(b-q3*aa));
      fp_t rc=q2*pp;
      fp_t sc=q4*(q4*pp*pp-rr);
      fp_t tc=-(q8*qq*q8*qq);
      
      cub_obj.rrteq3(rc,sc,tc,u,dc);
      
      fp_t q=qq;
      fp_t h=r4*a;
      if (dc==0) u[2]=u[1];
      if (dc<=0) {
        mt=2;
        v[1]=fabs(u[0]);
        v[2]=fabs(u[1]);
        v[3]=fabs(u[2]);
        v1=std::max(std::max(v[1],v[2]),v[3]);
        if (v1==v[1]) {
          k1=0;
          v2=std::max(v[2],v[3]);
        } else if (v1==v[2]) {
          k1=1;
          v2=std::max(v[1],v[3]);
        } else {
          k1=2;
          v2=std::max(v[1],v[2]);
        }
        if (v2==v[1]) {
          k2=0;
        } else if (v2==v[2]) {
          k2=1;
        } else {
          k2=2;
        }
        w1=sqrt(((std::complex<fp_t>)(u[k1])));
        w2=sqrt(((std::complex<fp_t>)(u[k2])));
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
    
    /// Return a string denoting the type ("quartic_real_coeff_cern")
    const char *type() { return "quartic_real_coeff_cern"; }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The object to solve for the associated cubic
    cubic_real_coeff_cern<fp_t> cub_obj;

#endif

  };

  /** \brief Solve a quadratic with real coefficients and complex roots (GSL)
   */
  class quadratic_real_coeff_gsl : public quadratic_real_coeff {

  public:

    virtual ~quadratic_real_coeff_gsl() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const double a2, const double b2, const double c2, 
			 std::complex<double> &x1, std::complex<double> &x2);

    /// Return a string denoting the type ("quadratic_real_coeff_gsl")
    const char *type() { return "quadratic_real_coeff_gsl"; }

  };

  /** \brief Solve a cubic with real coefficients and complex roots (GSL)
   */
  class cubic_real_coeff_gsl : public cubic_real_coeff {

  public:

    virtual ~cubic_real_coeff_gsl() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, 
			 std::complex<double> &x2, std::complex<double> &x3);

    /// Return a string denoting the type ("cubic_real_coeff_gsl")
    const char *type() { return "cubic_real_coeff_gsl"; }

    /** \brief An alternative to \c gsl_poly_complex_solve_cubic()
	
	This is an alternative to the function
	<tt>gsl_poly_complex_solve_cubic()</tt> with some small
	corrections to ensure finite values for some cubics. See
	<tt>src/other/poly_ts.cpp</tt> for more.

	\future I think the GSL function is now fixed, so we
	can fall back to the original GSL function here. 
    */
    int gsl_poly_complex_solve_cubic2(double a, double b, double c, 
				      gsl_complex *z0, gsl_complex *z1, 
				      gsl_complex *z2);
    
  };

  /** \brief Solve a quartic with real coefficients and real roots (GSL)

      This class internally uses the GSL functions to solve the
      resolvent cubic and associated quadratics, while
      \ref quartic_real_gsl2 contains explicit code to solve
      them instead.

      \future Optimize value of \c cube_root_tol and compare
      more clearly to \ref o2scl::quartic_real_gsl2
  */
  class quartic_real_gsl : public quartic_real {
    
  public:
    
    quartic_real_gsl() {
      cube_root_tol=1.0e-7;
    }

    virtual ~quartic_real_gsl() {}

    /** \brief A tolerance for determining the proper cube root 
	(default \f$ 10^{-4} \f$ )
    */
    double cube_root_tol;

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four real solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, double &x1, 
			double &x2, double &x3, double &x4);

    /// Return a string denoting the type ("quartic_real_gsl")
    const char *type() { return "quartic_real_gsl"; }

  };

  /** \brief Solve a quartic with real coefficients and real roots (GSL)

      This class directly solves 
      resolvent cubic and associated quadratics without using 
      the GSL functions (as done in \ref quartic_real_gsl).
      
      \future Optimize value of \c cube_root_tol and compare
      more clearly to \ref o2scl::quartic_real_gsl
  */
  class quartic_real_gsl2 : public quartic_real {

  public:

    quartic_real_gsl2() {
      cube_root_tol=1.0e-7;
    }

    virtual ~quartic_real_gsl2() {}

    /** \brief A tolerance for determining the proper cube root 
	(default \f$ 10^{-7} \f$ )
    */
    double cube_root_tol;

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four real solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, double &x1, 
			double &x2, 
			double &x3, double &x4);

    /// Return a string denoting the type ("quartic_real_gsl2")
    const char *type() { return "quartic_real_gsl2"; }
  };

  /** \brief Solve a general polynomial with real coefficients (GSL)
   */
  class poly_real_coeff_gsl : public poly_real_coeff {

  public:

    poly_real_coeff_gsl();

    virtual ~poly_real_coeff_gsl();

    /** \brief Solve a generic polynomial given <tt>n+1</tt> coefficients

	\note In order to be consistent with the other solve_rc()
	functions, the ordering of the coefficients is reversed with
	respect to gsl_poly_complex_solve(). The leading coefficient
	is stored in <tt>co[0]</tt> and the constant term is stored in
	<tt>co[n]</tt>.
     */
    virtual int solve_rc_arr(int n, const double co[], 
			     std::complex<double> ro[]);

    /** \brief Solve a cubic polynomial with real coefficients
     */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, 
			 std::complex<double> &x2, 
			 std::complex<double> &x3);

    /** \brief Solve a quadratic polynomial with real coefficients
     */
    virtual int solve_rc(const double a2, const double b2, const double c2, 
			 std::complex<double> &x1, 
			 std::complex<double> &x2);

    /** \brief Solve a quartic polynomial with real coefficients
     */
    virtual int solve_rc(const double a4, const double b4, const double c4, 
			 const double d4, const double e4, 
			 std::complex<double> &x1, std::complex<double> &x2, 
			 std::complex<double> &x3, std::complex<double> &x4);

    /// Return a string denoting the type ("poly_real_coeff_gsl")
    const char *type() { return "poly_real_coeff_gsl"; }

  protected:

#ifndef DOXYGEN_INTERNAL

    /// Workspace for quadratic polynomials
    gsl_poly_complex_workspace *w2;

    /// Workspace for cubic polynomials
    gsl_poly_complex_workspace *w3;

    /// Workspace for quartic polynomials
    gsl_poly_complex_workspace *w4;

    /// Workspace for general polynomials
    gsl_poly_complex_workspace *wgen;

    /// The size of the workspace \ref wgen
    int gen_size;

#endif

  };

  /** \brief Solve a quadratic with complex coefficients and complex roots
   */
  template<class fp_t=double> class quadratic_complex_std :
    public quadratic_complex {

  public:

    virtual ~quadratic_complex_std() {}

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const std::complex<fp_t> a2, 
			const std::complex<fp_t> b2, 
			const std::complex<fp_t> c2, 
			std::complex<fp_t> &x1,
                        std::complex<fp_t> &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex_std::solve_c().",
           exc_einval);
      }
      x1=(-b2+sqrt(b2*b2-4.0*a2*c2))/2.0/a2;
      x2=(-b2-sqrt(b2*b2-4.0*a2*c2))/2.0/a2;
      return success;
    }

    /// Return a string denoting the type ("quadratic_complex_std")
    const char *type() { return "quadratic_complex_std"; }
  };

  /** \brief Solve a cubic with complex coefficients and complex roots
   */
  template<class fp_t=double> class cubic_complex_std :
    public cubic_complex {

  public:

    virtual ~cubic_complex_std() {}
    
    /** \brief Solves the complex polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the three complex solutions \f$ x=x_1 \f$ , 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_c(const std::complex<fp_t> a3, 
			const std::complex<fp_t> b3, 
			const std::complex<fp_t> c3, 
			const std::complex<fp_t> d3, 
			std::complex<fp_t> &x1, std::complex<fp_t> &x2, 
			std::complex<fp_t> &x3) {
      
      if (a3==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cubic_std::complex::solve_c().",
           exc_einval);
      }
      
      std::complex<fp_t> p3, q3, mo;
      std::complex<fp_t> alpha, beta, cbrta, cbrtb;
      std::complex<fp_t> e2, e4, cacb;
      fp_t test, re_p3;
      
      if (a3==0.0) {
        quadratic_complex_std<> qsc;
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

      fp_t pi=boost::math::constants::pi<fp_t>();
      
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
    
    /// Return a string denoting the type ("cubic_complex_std")
    const char *type() { return "cubic_complex_std"; }
  };

  /** \brief Solve a quartic with real coefficients and real roots
   */
  class quartic_real_simple : public quartic_real {

  public:

    quartic_real_simple() {
      cube_root_tol=1.0e-6;
    }

    virtual ~quartic_real_simple() {}

    virtual int solve_r(const double a4, const double b4, const double c4, 
			const double d4, const double e4, double &x1, 
			double &x2, double &x3, double &x4);

    /// Return a string denoting the type ("quartic_real_simple")
    const char *type() { return "quartic_real_simple"; }

    /** \brief A tolerance for determining the proper cube root 
	(default \f$ 10^{-6} \f$ )
    */
    double cube_root_tol;
  };
  
  /** \brief Solve a quartic with complex coefficients and complex roots
   */
  template <class fp_t=double> class quartic_complex_simple :
    public quartic_complex {

  public:

    virtual ~quartic_complex_simple() {}

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const std::complex<fp_t> a4, 
			const std::complex<fp_t> b4, 
			const std::complex<fp_t> c4, 
			const std::complex<fp_t> d4, 
			const std::complex<fp_t> e4, 
			std::complex<fp_t> &x1, 
			std::complex<fp_t> &x2, 
			std::complex<fp_t> &x3,
			std::complex<fp_t> &x4) {
      
      std::complex<fp_t> p4, q4, r4;
      std::complex<fp_t> a3, b3, c3, d3;
      std::complex<fp_t> b2a, c2a, b2b, c2b;
      std::complex<fp_t> u4, u41, u42;
      
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

    /// Return a string denoting the type ("quartic_complex_simple")
    const char *type() { return "quartic_complex_simple"; }

#ifndef DOXYGEN_NO_O2NS

  protected:

    /// The object to solve for the associated cubic
    cubic_complex_std<> cub_obj;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
