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
#include <o2scl/misc.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Solve a quadratic polynomial with real coefficients and 
      real roots [abstract base]
  */
  template<class fp_t=double> class quadratic_real {
    
  public:

    virtual ~quadratic_real() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2)=0;

    /** \brief Compute the quadratic discriminant, \f$ b^2-4ac \f$
     */
    virtual fp_t disc_r(fp_t a2, fp_t b2, fp_t c2) {
      return b2*b2-4.0*a2*c2;
    }
    
    /// Return a string denoting the type ("quadratic_real")
    const char *type() { return "quadratic_real"; }
  };

  /** \brief Solve a quadratic polynomial with real coefficients to
      get the real roots (GSL)

      This is just a wrapper for gsl_poly_complex_solve_quadratic().
  */
  class quadratic_real_gsl : public quadratic_real<> {
    
  public:

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const double a2, const double b2, const double c2, 
			double &x1, double &x2);
    
  };
  
  /** \brief Solve a quadratic polynomial with real coefficients and 
      complex roots [abstract base]

      This class is designed to be identical to the GSL code in
      gsl_poly_complex_solve_quadratic(), but generalized to generic
      floating point types.
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_real_coeff : public quadratic_real<fp_t> {

  public:

    virtual ~quadratic_real_coeff() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_real_coeff::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2;
      int ret=solve_rc(a2,b2,c2,r1,r2);
      x1=r1.real();
      x2=r2.real();
      return ret;
    }

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &x1, cx_t &x2)=0;
    
    /** \brief Desc
     */
    void test_quadratic_real_coeff_base
    (fp_t alpha, fp_t &s1, fp_t &s2, fp_t &m1,
     fp_t &m2, clock_t &lt1, clock_t &lt2, size_t n=40) {
      
      s1=0.0;
      s2=0.0;
      m1=0.0;
      m2=0.0;
      lt1=clock();
      
      o2scl::gen_test_number<40> ga, gb, gc;
      o2scl::gen_test_number<80> gd, ge;
      
      size_t count=0;
      
      // First pick random coefficients
      for(int j1=0;j1<n;j1++) {
        fp_t ca=ga.gen();
        for(int j2=0;j2<n;j2++) {
          fp_t cb=gb.gen()*alpha;
          for(int j3=0;j3<n;j3++) {
            
            // Ensure that several quadratics near b^2=4*a*c are tested
            fp_t cc=cb*cb/4.0/ca+gc.gen();
            
            // Ensure there is a solution
            if (fabs(ca)>0.0) {

              cx_t cr1, cr2;
              solve_rc(ca,cb,cc,cr1,cr2);
              
              cx_t cbp=-(cr1+cr2)*ca;
              cx_t ccp=(cr1*cr2)*ca;
              
              cx_t czo1=(ca*cr1+cb)*cr1+cc;
              cx_t czo2=(ca*cr2+cb)*cr2+cc;
              fp_t q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
                           fabs(cc-ccp.real())*fabs(cc-ccp.real()));
              fp_t q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
              
              s1+=q1;
              if (q1>m1) m1=q1;
              s2+=q2;
              if (q2>m2) m2=q2;
              count++;
              
              if (!o2isfinite(q1) || !o2isfinite(q2) || 
                  fabs(q1)>1.0e-10 || fabs(q2)>1.0e-10) {
                O2SCL_ERR("Failure in test_quadratic_real_coeff().",
                          exc_esanity);
              }
              
            }
            
          }
        }
      }
      
      // Next, pick random roots which are complex conjugates
      for(int j1=0;j1<n*2;j1++) {
        cx_t cr1, cr2;
        cr1.real(gd.gen());
        for(int j2=0;j2<n*2;j2++) {
          cr1.imag(ge.gen());
          cr2.real(cr1.real());
          cr2.imag(-cr1.imag());
          
          fp_t ca=1.0;
          fp_t cb=(-cr1-cr2).real();
          fp_t cc=(cr1*cr2).real();

          cx_t cr1p, cr2p;
          solve_rc(ca,cb,cc,cr1p,cr2p);
          
          // If the roots are flipped
          if (fabs(cr1.imag()-cr2p.imag())<fabs(cr1.imag()-cr1p.imag())) {
            cx_t temp=cr1p;
            cr1p=cr2p;
            cr2p=temp;
          } 
          
          fp_t q1=sqrt(fabs(cr1.real()-cr1p.real())*
                       fabs(cr1.real()-cr1p.real())+
                       fabs(cr2.real()-cr2p.real())*
                       fabs(cr2.real()-cr2p.real()));
          fp_t q2=sqrt(fabs(cr1.imag()-cr1p.imag())*
                       fabs(cr1.imag()-cr1p.imag())+
                       fabs(cr2.imag()-cr2p.imag())*
                       fabs(cr2.imag()-cr2p.imag()));
          
          s1+=q1;
          if (q1>m1) m1=q1;
          s2+=q2;
          if (q2>m2) m2=q2;
          count++;
          
          if (!o2isfinite(q1) || !o2isfinite(q2) ||
              fabs(q1)>1.0e-10 || fabs(q2)>1.0e-10) {
            O2SCL_ERR("Failure in test_quadratic_real_coeff().",
                      exc_esanity);
          }
        }
      }
      
      // Next, pick random roots which are both real
      for(int j1=0;j1<n*2;j1++) {
        cx_t cr1, cr2;
        cr1.real(gd.gen());
        for(int j2=0;j2<n*2;j2++) {
          cr2.real(ge.gen());
          cr1.imag(0.0);
          cr2.imag(0.0);
          
          fp_t ca=1.0;
          fp_t cb=(-cr1-cr2).real();
          fp_t cc=(cr1*cr2).real();

          cx_t cr1p, cr2p;
          solve_rc(ca,cb,cc,cr1p,cr2p);
          
          // If the roots are flipped
          if (fabs(cr1.real()-cr2p.real())<fabs(cr1.real()-cr1p.real())) {
            cx_t temp=cr1p;
            cr1p=cr2p;
            cr2p=temp;
          } 
          
          fp_t q1=sqrt(fabs(cr1.real()-cr1p.real())*
                       fabs(cr1.real()-cr1p.real())+
                       fabs(cr2.real()-cr2p.real())*
                       fabs(cr2.real()-cr2p.real()));
          fp_t q2=sqrt(fabs(cr1.imag()-cr1p.imag())*
                       fabs(cr1.imag()-cr1p.imag())+
                       fabs(cr2.imag()-cr2p.imag())*
                       fabs(cr2.imag()-cr2p.imag()));
          
          s1+=q1;
          if (q1>m1) m1=q1;
          s2+=q2;
          if (q2>m2) m2=q2;
          count++;
          
          if (!o2isfinite(q1) || !o2isfinite(q2) ||
              fabs(q1)>1.0e-10 || fabs(q2)>1.0e-10) {
            O2SCL_ERR("Failure in test_quadratic_real_coeff().",
                      exc_esanity);
          }
        }
      }
      
      lt2=clock();
      s1/=count;
      s2/=count;
      
      return;
    }
    
    /// Return a string denoting the type ("quadratic_real_coeff")
    const char *type() { return "quadratic_real_coeff"; }
  };

  /** \brief Solve a quadratic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_complex :  public quadratic_real_coeff<fp_t,cx_t> {
    
  public:

    virtual ~quadratic_complex() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2;
      int ret=solve_c(a2,b2,c2,r1,r2);
      x1=r1.real();
      x2=r2.real();
      return ret;
    }

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &x1, cx_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex::solve_rc().",
           exc_einval);
      }
      
      int ret=solve_c(a2,b2,c2,x1,x2);
      return ret;
    }

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const cx_t a2, const cx_t b2, 
			const cx_t c2, cx_t &x1, cx_t &x2)=0;

    /// Return a string denoting the type ("quadratic_complex")
    const char *type() { return "quadratic_complex"; }
  };  

  /** \brief Solve a cubic polynomial with real coefficients and real roots
      [abstract base]
  */
  template<class fp_t=double> class cubic_real {
  public:

    virtual ~cubic_real() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, 
			fp_t &x3)=0;

    /** \brief Compute the cubic discriminant, 
	\f$ b^2 c^2 - 4 a c^3 - 4 b^3 d - 27 a^2 d^2 + 18 a b c d \f$
     */
    virtual fp_t disc_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			  const fp_t d3) {
      return b3*b3*c3*c3-4.0*a3*c3*c3*c3-4.0*b3*b3*b3*d3-
	27.0*a3*a3*d3*d3+18.0*a3*b3*c3*d3;
    }

    /// Return a string denoting the type ("cubic_real")
    const char *type() { return "cubic_real"; }
  };

  /** \brief Solve a cubic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff : public cubic_real<fp_t> {

  public:

    virtual ~cubic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, fp_t &x3) {
      
      cx_t r2,r3;
  
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_real_coeff::solve_r().",exc_einval);
      }
      
      int ret=solve_rc(a3,b3,c3,d3,x1,r2,r3);
      x2=r2.real();
      x3=r3.real();
      return ret;
    }

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, cx_t &x2,
			 cx_t &x3)=0;

    /// Return a string denoting the type ("cubic_real_coeff")
    const char *type() { return "cubic_real_coeff"; }
  };

  /** \brief Solve a cubic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_complex : public cubic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~cubic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, fp_t &x3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_complex::solve_r().",exc_einval);
      }
      
      cx_t r1,r2,r3;
      int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
      x1=r1.real();
      x2=r2.real();
      x3=r3.real();
      return ret;
    }

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, cx_t &x2,
			 cx_t &x3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_complex::solve_rc().",exc_einval);
      }
      
      cx_t r1,r2,r3;
      int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
      fp_t s1,s2,s3;
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
    virtual int solve_c(const cx_t a3, 
			const cx_t b3, 
			const cx_t c3, 
			const cx_t d3, 
			cx_t &x1, cx_t &x2, 
			cx_t &x3)=0;

    /// Return a string denoting the type ("cubic_complex")
    const char *type() { return "cubic_complex"; }
  };

  /** \brief Solve a quartic polynomial with real coefficients and 
      real roots [abstract base]
  */
  template<class fp_t=double> class quartic_real {

  public:

    virtual ~quartic_real() {}
    
    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, 
			fp_t &x1, fp_t &x2, 
			fp_t &x3, fp_t &x4)=0;

    /** \brief Compute the discriminant

	The discriminant is zero if and only if at least two roots are
	equal. If the discriminant is non-zero, the discriminant is
	negative if there are two real roots and two complex conjugate
	roots, and it is positive if the roots are either all real or
	all non-real.
    */
    virtual fp_t disc_r(const fp_t a, const fp_t b, const fp_t c, 
			  const fp_t d, const fp_t e) {
      fp_t a2=a*a;
      fp_t b2=b*b;
      fp_t c2=c*c;
      fp_t d2=d*d;
      fp_t e2=e*e;
      
      fp_t a3=a2*a;
      fp_t b3=b2*b;
      fp_t c3=c2*c;
      fp_t d3=d2*d;
      fp_t e3=e2*e;
      
      fp_t b4=b2*b2;
      fp_t c4=c2*c2;
      fp_t d4=d2*d2;
      
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
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_real_coeff : public quartic_real<fp_t> {

  public:

    virtual ~quartic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2,r3,r4;
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
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, 
			 cx_t &x3, 
			 cx_t &x4)=0;

    /// Return a string denoting the type ("quartic_real_coeff")
    const char *type() { return "quartic_real_coeff"; }
  };

  /** \brief Solve a quartic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_complex : public quartic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~quartic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
			      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_complex::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2,r3,r4;
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
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, cx_t &x3, cx_t &x4) {
			 
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cx_t &x4) {().",
           exc_einval);
      }
      
      return solve_c(a4,b4,c4,d4,e4,x1,x2,x3,x4);
    }      

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const cx_t a4, const cx_t b4, 
			const cx_t c4, const cx_t d4, 
			const cx_t e4, cx_t &x1, 
			cx_t &x2, cx_t &x3, cx_t &x4)=0;

    /// Return a string denoting the type ("quartic_complex")
    const char *type() { return "quartic_complex"; }
  };

  /** \brief Solve a general polynomial with real
      coefficients and complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<fp_t>,
           class root_vec_t=std::vector<cx_t> >
  class poly_real_coeff : public quadratic_real_coeff<fp_t,cx_t>,
                          public cubic_real_coeff<fp_t,cx_t>,
                          public quartic_real_coeff<fp_t,cx_t> {
    
  public:
    
    virtual ~poly_real_coeff() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_rc_arr(int n, const coeff_vec_t &co, 
                             root_vec_t &ro)=0;

    virtual int polish_rc_arr(int n, const coeff_vec_t &co, 
                              root_vec_t &ro) {
      return 0;
    }

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

    /** \brief Desc
     */
    virtual int polish_fun(size_t nv, const ubvector &x,
                           fp_t ubvector &y, const coeff_vec_t &co,
                           root_vec_t &ro) {
      
      // Using horner's method following, e.g. GSL
      y[0]=co[nv-1];
      y[1]=0.0;
      for(int i=nv-1; i>0; i--) {
        fp_t tmp=co[i-1]+x[0]*y[0]-x[1]*y[1];
        y[1]=x[1]*y[0]+x[0]*y[1];
        y[0]=tmp;
      }
      return 0;
    }
    
#endif

    
    /// Return a string denoting the type ("poly_real_coeff")
    const char *type() { return "poly_real_coeff"; }
    
  };

  /** \brief Solve a general polynomial with complex
      coefficients [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<fp_t>,
           class root_vec_t=std::vector<cx_t> >
  class poly_complex : public quadratic_complex<fp_t, cx_t>,
                       public cubic_complex<fp_t, cx_t>,
                       public quartic_complex<fp_t, cx_t> {

  public:

    virtual ~poly_complex() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_c_arr(int n, const coeff_vec_t &co,
                            root_vec_t &ro)=0;
    
    /// Polish the roots 
    virtual int polish_c_arr(int n, const coeff_vec_t &co,
                             root_vec_t &ro)=0;

    /// Return a string denoting the type ("poly_complex")
    const char *type() { return "poly_complex"; }
    
  };

  /** \brief Solve a cubic with real coefficients and complex roots 
      (CERNLIB)

      \note The function rrteq3() is based on the CERNLIB routine of
      the same name, but differs slightly. See the documentation of
      that function for details.
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff_cern : public cubic_real_coeff<fp_t,cx_t> {

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
			 cx_t &r2,
                         cx_t &r3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
                   "cubic_real_coeff_cern::solve_rc().",exc_einval);
      }
      
      fp_t x[3],d;
      cx_t i(0.0,1.0);
      
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

    // Generic sign function
    inline fp_t sign(fp_t a, fp_t b) {
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
      fp_t r1=2.0;
      r1/=27.0;
      fp_t r2=1.0;
      r2/=2.0;
      fp_t r3=1.0;
      r3/=3.0;
      fp_t w3=3.0;
      w3=sqrt(w3);
      fp_t r4=w3/2.0;
      fp_t q1=2.0;
      q1/=27.0;
      fp_t q2=r2;
      fp_t q3=r3;
      fp_t y[3];
      cx_t z[3], i(0.0,1.0);
      fp_t h2, h3;
      int j,k;

      // AWS 4/1/21: These temporaries fix compiling for long double
      // and complex<long double> types
      fp_t one=1.0;
      fp_t two=2.0;
      fp_t three=3.0;
      
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
      if (o2abs(d)<=eps) {
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
          d_new=da/db+one;
        } else {
          d_new=-da/db-one;
        }
      } else {
        d_new=d;
      }
      
      if (d_new>delta2) {
        h2=sqrt(d);
        fp_t u0=-h1+h2;
        fp_t v0=-h1-h2;
        if (o2abs(u0)==0.0) {
          u=sign(0.0,u0);
        }
        else u=sign(o2pow(o2abs(u0),r3),u0);
        if (o2abs(v0)==0.0) v=sign(0.0,v0);
        else v=sign(o2pow(o2abs(v0),r3),v0);
        x[0]=u+v-h;
        x[1]=-r2*(u+v)-h;
        x[2]=r4*o2abs(u-v);
        if (o2abs(u0)<=eps || o2abs(v0)<=eps) {
          y[0]=x[0];
          for(k=0;k<=1;k++) {
            y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/
              ((three*y[k]+two*r)*y[k]+s);
          }
          x[0]=y[2];
          z[0]=x[1]+i*x[2];
          for(k=0;k<=1;k++) {
            z[k+1]=z[k]-(((z[k]+r)*z[k]+s)*z[k]+t)/
              ((three*z[k]+two*r)*z[k]+s);
          }
          x[1]=z[2].real();
          x[2]=z[2].imag();
        }
        
      } else if (o2abs(d_new)<=delta2) {
        
        d=0.0;
        if (o2abs(h1)==0.0) u=sign(0.0,-h1);
        else u=sign(o2pow(o2abs(h1),r3),-h1);
        x[0]=u+u-h;
        x[1]=-u-h;
        x[2]=x[1];
        if (o2abs(h1)<=eps) {
          y[0]=x[0];
          for(k=0;k<=1;k++) {
            h1=(three*y[k]+two*r)*y[k]+s;
            if (o2abs(h1)>delta2) {
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
        
        h3=o2abs(r3*p);
        h3=o2sqrt(h3*h3*h3);
        h2=r3*acos(-h1/h3);
        if (h3==0.0) h1=0.0;
        else h1=o2pow(h3,r3);
        u=h1*cos(h2);
        v=w3*h1*sin(h2);
        x[0]=u+u-h;
        x[1]=-u-v-h;
        x[2]=-u+v-h;
        if (h3<=eps || x[0]<=eps || x[1]<=eps || x[2]<=eps) {
          for(j=0;j<3;j++) {
            y[0]=x[j];
            for(k=0;k<=1;k++) {
              y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/
                ((three*y[k]+two*r)*y[k]+s);
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
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_real_coeff_cern : public quartic_real_coeff<fp_t,cx_t> {
    
  public:

    virtual ~quartic_real_coeff_cern() {}

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four complex solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, cx_t &x3, cx_t &x4) {
      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff_cern::solve_rc().",
           exc_einval);
      }
      
      int mt;
      fp_t dc;
      cx_t x[4];
      
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
		       cx_t z[], fp_t &dc, 
		       int &mt) {
      
      cx_t i(0.0,1.0), z0[5];
      cx_t w1(0.0,0.0), w2(0.0,0.0), w3;

      fp_t q2=1.0;
      q2/=2.0;
      fp_t r4=q2/2.0;
      fp_t q4=r4;
      fp_t q8=r4/2.0;
      fp_t r12=r4/3.0;
      fp_t q1=q8*3.0;
      fp_t q3=q1/2.0;
      fp_t eight=8.0;
        
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
      fp_t zero=0.0;
      
      cub_obj.rrteq3(rc,sc,tc,u,dc);
      
      fp_t q=qq;
      fp_t h=r4*a;
      if (dc==0) u[2]=u[1];
      if (dc<=0) {
        mt=2;
        v[1]=o2abs(u[0]);
        v[2]=o2abs(u[1]);
        v[3]=o2abs(u[2]);
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
        w1=sqrt(((cx_t)(u[k1])));
        w2=sqrt(((cx_t)(u[k2])));
      } else {
        mt=3;
        w1=sqrt(u[1]+i*u[2]);
        w2=sqrt(u[1]-i*u[2]);
      }
      w3=0;
      if (w1*w2!=zero) w3=-q/(eight*w1*w2);
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
    cubic_real_coeff_cern<fp_t,cx_t> cub_obj;

#endif

  };

  /** \brief Solve a quadratic with real coefficients and complex roots (GSL)
   */
  class quadratic_real_coeff_gsl : public quadratic_real_coeff<> {
    
  public:
    
    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const double a, const double b, const double c,
                         std::complex<double> &x1, std::complex<double> &x2);
    
  };

  /** \brief Solve a quadratic with real coefficients and complex roots 
      (C++ rewrite of GSL algorithm)
      
      This code was based on GSL v2.6's function
      gsl_poly_complex_solve_quadratic().
   */
  template<class fp_t=double, class cx_t=std::complex<double> >
  class quadratic_real_coeff_gsl2 : public quadratic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~quadratic_real_coeff_gsl2() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 

        This function returns the number of complex roots (either
        1 or 2). This function calls the error handler if 
        \f$ a=b=0 \f$.
    */
    virtual int solve_rc(const fp_t a, const fp_t b, const fp_t c,
                         cx_t &x1, cx_t &x2) {

      // AWS, 3/31/21: This is equivalent to the GSL code from v2.6,
      // but rewritten in a C++/multiprecision compatible form
      
      fp_t disc=b*b-4*a*c;
      if (a == 0) {
        if (b == 0) {
          O2SCL_ERR2("No solution because a=b=0 in ",
                     "quadratic_real_coeff::solve_rc().",o2scl::exc_einval);
        } else {
          x1.real(-c/b);
          x1.imag(0);
          return 1;
        };
      }
      
      if (disc > 0) {
        if (b == 0) {
          fp_t s=fabs(sqrt(disc)/a/2);
          x1.real(-s);
          x1.imag(0);
          x2.real(s);
          x2.imag(0);
        } else {
          fp_t sgnb=(b > 0 ? 1 : -1);
          fp_t temp=-(b+sgnb*sqrt(disc))/2;
          fp_t r1=temp/a;
          fp_t r2=c/temp;
          
          if (r1 < r2) {
            x1.real(r1);
            x1.imag(0);
            x2.real(r2);
            x2.imag(0);
          } else {
            x1.real(r2);
            x1.imag(0);
            x2.real(r1);
            x2.imag(0);
          }
        }
        return 2;
        
      } else if (disc == 0) {
        
        x1.real(-b/a/2);
        x1.imag(0);
        x2.real(-b/a/2);
        x2.imag(0);
        
        return 2;
        
      }
      
      fp_t s=fabs(sqrt(-disc)/a/2);
      x1.real(-b/a/2);
      x1.imag(-s);
      x2.real(-b/a/2);
      x2.imag(s);
      
      return 2;
    }

    /// Return a string denoting the type ("quadratic_real_coeff_gsl2")
    const char *type() { return "quadratic_real_coeff_gsl2"; }

  };

  /** \brief Solve a cubic with real coefficients and complex roots (GSL)

      This is just a wrapper for gsl_poly_complex_solve_cubic().
   */
  class cubic_real_coeff_gsl : public cubic_real_coeff<double> {

  public:

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, 
			 std::complex<double> &x2, std::complex<double> &x3);
    
  };

  /** \brief Solve a cubic with real coefficients and complex roots
      (C++ rewrite of GSL algorithm)

      This class is designed to be identical to the GSL code in
      gsl_poly_complex_solve_cubic(), but generalized to generic
      floating point types.
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff_gsl2 : public cubic_real_coeff<fp_t,cx_t> {

  protected:
    
    inline void swap(double &a, double &b) {
      double tmp=b;
      b=a;
      a=tmp;
      return;
    }
    
  public:

    virtual ~cubic_real_coeff_gsl2() {}

    /** \brief solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, 
			 std::complex<fp_t> &x2, std::complex<fp_t> &x3) {
      
      if (a3==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cubic_real_coeff_gsl2::solve_rc().",
           exc_einval);
      }
      
      fp_t a=b3/a3;
      fp_t b=c3/a3;
      fp_t c=d3/a3;
      
      fp_t q=(a*a-3*b);
      fp_t r=(2*a*a*a-9*a*b+27*c);
      
      fp_t Q=q/9;
      fp_t R=r/54;
      
      fp_t Q3=Q*Q*Q;
      fp_t R2=R*R;
      
      fp_t CR2=729*r*r;
      fp_t CQ3=2916*q*q*q;
      
      if (R == 0 && Q == 0) {
        
        x1=-a/3;
        x2.real(-a/3);
        x2.imag(0);
        x3.real(-a/3);
        x3.imag(0);

        if (!std::isfinite(x1) || !std::isfinite(x2.real()) ||
            !std::isfinite(x2.imag())) {
          std::cout << "1. " << x1 << " " << x2 << " " << x3 << std::endl;
          exit(-1);
        }
        return 3;
        
      } else if (CR2 == CQ3)  {

        /* [GSL] This test is actually r2 == q3, written in a form suitable
           for exact computation with integers */
        
        /* [GSL] Due to finite precision some double roots may be
           missed, and will be considered to be a pair of complex
           roots z=x +/- epsilon i close to the real axis. */
        
        fp_t sqrtQ=sqrt(Q);
        
        if (R > 0) {
          x1=-2*sqrtQ-a/3;
          x2.real(sqrtQ-a/3);
          x2.imag(0);
          x3.real(sqrtQ-a/3);
          x3.imag(0);
        } else {
          x1=-sqrtQ-a/3;
          x2.real(-sqrtQ-a/3);
          x2.imag(0);
          x3.real(2*sqrtQ-a/3);
          x3.imag(0);
        }
        
        if (!std::isfinite(x1) || !std::isfinite(x2.real()) ||
            !std::isfinite(x2.imag())) {
          std::cout << "2. " << x1 << " " << x2 << " " << x3 << std::endl;
          exit(-1);
        }
        return 3;
        
      } else if (CR2 < CQ3)  {

        /* [GSL] equivalent to R2 < Q3 */
        
        fp_t sqrtQ=sqrt(Q);
        fp_t sqrtQ3=sqrtQ*sqrtQ*sqrtQ;
        fp_t theta=acos(R/sqrtQ3);
        
        // [AWS] Modified from the original GSL routine
        // sometimes R/sqrtQ3 is slightly larger than one
        // when the coefficients are arranged just right
        // so theta becomes not finite.
        
        if (R/sqrtQ3>=1.0) theta=0.0;
        if (R/sqrtQ3<=-1.0) theta=o2scl_const::pi;
        
        fp_t norm=-2*sqrtQ;
        
        fp_t r0=norm*cos(theta/3)-a/3;
        fp_t r1=norm*cos((theta+2.0*M_PI)/3)-a/3;
        fp_t r2=norm*cos((theta-2.0*M_PI)/3)-a/3;
        
        x1=r0;
        x2.real(r1);
        x2.imag(0);
        x3.real(r2);
        x3.imag(0);
        
        if (!std::isfinite(x1) || !std::isfinite(x2.real()) ||
            !std::isfinite(x2.imag())) {
          std::cout << "3. " << x1 << " " << x2 << " " << x3 << std::endl;
          std::cout << Q << std::endl;
          std::cout << norm << " " << sqrtQ << " " << theta << std::endl;
          std::cout << r0 << " " << r1 << " " << r2 << std::endl;
          exit(-1);
        }
        return 3;
        
      }
        
      fp_t sgnR=(R>=0?1:-1);
      fp_t A;
      
      // [AWS] Modification from original GSL behavior: Just in case
      // R2=Q3, finite precision can cause the argument of the sqrt()
      // function to be negative. We correct for this here.
      if (R2<=Q3) {
        A=-sgnR*pow(fabs(R),1.0/3.0);
      } else {
        A=-sgnR*pow(fabs(R)+sqrt(R2-Q3),1.0/3.0);
      }
      
      fp_t B=Q/A;
      
      x1=A+B-a/3;
      x2.real(-0.5*(A+B)-a/3);
      x2.imag(-(sqrt(3.0)/2.0)*fabs(A-B));
      x3.real(-0.5*(A+B)-a/3);
      x3.imag((sqrt(3.0)/2.0)*fabs(A-B));
      
      if (!std::isfinite(x1) || !std::isfinite(x2.real()) ||
          !std::isfinite(x2.imag())) {
        std::cout << R2-Q3 << std::endl;
        std::cout << fabs(R)+sqrt(R2-Q3) << std::endl;
        std::cout << "4. " << x1 << " " << x2 << " " << x3 << std::endl;
        std::cout << sgnR << " " << R2 << " " << Q3 << " " << R << " "
                  << B << " " << Q << " " << A << " " << a << std::endl;
        exit(-1);
      }
        
      return 3;
    }

    /// Return a string denoting the type ("cubic_real_coeff_gsl2")
    const char *type() { return "cubic_real_coeff_gsl2"; }

  };

  /** \brief Solve a quartic with real coefficients and real roots (GSL)

      This class internally uses the GSL functions to solve the
      resolvent cubic and associated quadratics, while
      \ref quartic_real_gsl2 contains explicit code to solve
      them instead.

      \future Optimize value of \c cube_root_tol and compare
      more clearly to \ref o2scl::quartic_real_gsl2
  */
  class quartic_real_gsl : public quartic_real<double> {
    
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
  class quartic_real_gsl2 : public quartic_real<double> {

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
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<double>,
           class root_vec_t=std::vector<std::complex<double> > >
  class poly_real_coeff_gsl : public poly_real_coeff<fp_t,cx_t> {

  public:

    poly_real_coeff_gsl() {
      w2=gsl_poly_complex_workspace_alloc(3);
      w3=gsl_poly_complex_workspace_alloc(4);
      w4=gsl_poly_complex_workspace_alloc(5);
      gen_size=0;
    }

    virtual ~poly_real_coeff_gsl() {
      gsl_poly_complex_workspace_free(w2);
      gsl_poly_complex_workspace_free(w3);
      gsl_poly_complex_workspace_free(w4);
      if (gen_size>0) {
        gsl_poly_complex_workspace_free(wgen);
      }
    }

    /** \brief Solve a generic polynomial given <tt>n+1</tt> coefficients

	\note In order to be consistent with the other solve_rc()
	functions, the ordering of the coefficients is reversed with
	respect to gsl_poly_complex_solve(). The leading coefficient
	is stored in <tt>co[0]</tt> and the constant term is stored in
	<tt>co[n]</tt>.
     */
    virtual int solve_rc_arr(int n, const coeff_vec_t &co,
                             root_vec_t &ro) {
      
      int j;
      typedef boost::numeric::ublas::vector<fp_t> ubvector;
      ubvector a(n+1), z(2*n);
      cx_t i(0.0,1.0);
      
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

    /** \brief Solve a quadratic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &r1, cx_t &r2) {
			 
      if (a2==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
      }
      
      fp_t a[3]={c2,b2,a2};
      fp_t z[4];
      cx_t i(0.0,1.0);

      gsl_poly_complex_solve(a,3,w2,z);
      
      r1=z[0]+i*z[1];
      r2=z[2]+i*z[3];
      
      return success;
    }
    
    /** \brief Solve a cubic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &r1, cx_t &r2, 
			 cx_t &r3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
      }
      
      fp_t a[4]={d3,c3,b3,a3};  
      fp_t z[6],s1,s2,s3;
      cx_t i(0.0,1.0);
      
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
    
    /** \brief Solve a quartic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &r1, cx_t &r2, 
			 cx_t &r3, cx_t &r4) {
      if (a4==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
                   "poly_real_coeff_gsl::solve_rc().",
                   exc_einval);
      }
      
      fp_t a[5]={e4,d4,c4,b4,a4};
      fp_t z[8];
      cx_t i(0.0,1.0);
      
      gsl_poly_complex_solve(a,5,w4,z);
      
      r1=z[0]+i*z[1];
      r2=z[2]+i*z[3];
      r3=z[4]+i*z[5];
      r4=z[6]+i*z[7];
      
      return success;
    }

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
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_complex_std : public quadratic_complex<fp_t,cx_t> {

  public:

    virtual ~quadratic_complex_std() {}

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const cx_t a2, const cx_t b2, 
			const cx_t c2, cx_t &x1, cx_t &x2) {
      
      if (a2.real()==0.0 && a2.imag()==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex_std::solve_c().",
           exc_einval);
      }
      // AWS 4/1/21: These temporaries fix compiling for long double
      // and complex<long double> types
      fp_t two=2.0;
      fp_t four=4.0;
      x1=(-b2+sqrt(b2*b2-four*a2*c2))/two/a2;
      x2=(-b2-sqrt(b2*b2-four*a2*c2))/two/a2;
      return success;
    }

    /// Return a string denoting the type ("quadratic_complex_std")
    const char *type() { return "quadratic_complex_std"; }
  };

  /** \brief Solve a cubic with complex coefficients and complex roots
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_complex_std : public cubic_complex<fp_t,cx_t> {

  public:

    virtual ~cubic_complex_std() {}
    
    /** \brief Solves the complex polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the three complex solutions \f$ x=x_1 \f$ , 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_c(const cx_t a3, const cx_t b3, 
			const cx_t c3, const cx_t d3, 
			cx_t &x1, cx_t &x2, cx_t &x3) {

      if (a3.real()==0.0 && a3.imag()==0.0) {
        quadratic_complex_std<fp_t,cx_t> qsc;
        qsc.solve_c(b3,c3,d3,x1,x2);
        x3=0.0;
        return success;
      }
      
      cx_t p3, q3, mo(0,1);
      cx_t alpha, beta, cbrta, cbrtb;
      cx_t e2, e4, cacb;
      fp_t test, re_p3;

      fp_t two=2.0;
      fp_t three=3.0;
      fp_t four=4.0;
      fp_t nine=9.0;
      fp_t twoseven=27.0;
      fp_t onethird=nine/twoseven;
      
      fp_t pi=boost::math::constants::pi<fp_t>();

      if (false) {

        // AWS 4/4/21: This code appeared to be slightly more accurate
        // on average than the code below, but it occasionally gave
        // incorrect roots. 
        
        p3=(three*a3*c3-b3*b3)/nine/a3/a3;
        q3=(two*b3*b3*b3-nine*a3*b3*c3+twoseven*a3*a3*d3)/twoseven/a3/a3/a3;
        std::cout << "p,q: " << p3 << " " << q3 << std::endl;
        
        alpha=(-q3+sqrt(q3*q3+four*p3*p3*p3))/two;
        beta=(q3+sqrt(q3*q3+four*p3*p3*p3))/two;
        std::cout << "alpha,beta: " << alpha << " " << beta << std::endl;
        
        if (alpha.real()==0.0) cbrta=0.0;
        else cbrta=pow(alpha,onethird);
        if (beta.real()==0.0) cbrtb=0.0;
        else cbrtb=pow(beta,onethird);
        
        cx_t exp_i_pi_three=exp(mo*pi/three);
        
        // It seems that if the real part of alpha is < 0 and the imaginary
        // part is zero, then cbrta is NaN (esp. w/Cygwin). We fix this
        // here:
        if (!o2isfinite(cbrta.real())) {
          cbrta=pow(-alpha,onethird)*exp_i_pi_three;
        }
        if (!o2isfinite(cbrtb.real())) {
          cbrtb=pow(-beta,onethird)*exp_i_pi_three;
        }
        
        e2=exp_i_pi_three*exp_i_pi_three;
        e4=e2*e2;
        
        std::cout << "alpha,beta: " << cbrta << " " << cbrtb << std::endl;
        
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
            std::cout << "1" << std::endl;
          }
          std::cout << "2" << std::endl;
        } else {
          if (fabs(((cacb*e2*e2).real()-re_p3)/re_p3)<test) {
            cbrta*=e2*e2;
            std::cout << "3" << std::endl;
          }
        }
        
        x1=cbrta-cbrtb-b3/three/a3;
        x2=cbrta*e2-cbrtb*e4-b3/three/a3;
        x3=cbrta*e4-cbrtb*e2-b3/three/a3;
        std::cout << "x1,x2,x3: " << x1 << " " << x2 << " " << x3
                  << std::endl;
        
      }

      if (true) {

        // Following
        // https://en.wikipedia.org/wiki/Cubic_equation#Vieta's_substitution

        // Write the cubic as a "depressed" cubic: x^3+px+q=0
        cx_t p=(three*a3*c3-b3*b3)/three/a3/a3;
        cx_t q=(two*b3*b3*b3-nine*a3*b3*c3+twoseven*a3*a3*d3)/
          twoseven/a3/a3/a3;

        // Formulate into a quadratic and then W is one of
        // the roots of the quadratic
        cx_t W=-q/two+sqrt(p*p*p/twoseven+q*q/four);
        // This code fails when W=0
        //std::cout << "W: " << W << std::endl;
        
        // Construct the roots of the depressed cubic from the
        // solution of the quadratic
        cx_t w1=pow(W,onethird);
        cx_t exp_i_pi_three=exp(mo*pi*two/three);
        cx_t w2=w1*exp_i_pi_three;
        cx_t w3=w2*exp_i_pi_three;
        //std::cout << "p: " << p << std::endl;

        // Vieta's substitution
        cx_t r1, r2, r3;
        if (W.real()==0.0 && W.imag()==0.0) {
          // If W=0, then there are three identical roots of the
          // cubic, and they are all -b3/a3/3
          r1=0;
          r2=0;
          r3=0;
        } else {
          //std::cout << "Here " << W.real() << " " << W.imag() << std::endl;
          r1=w1-p/three/w1;
          r2=w2-p/three/w2;
          r3=w3-p/three/w3;
        }
        //std::cout << "b3: " << b3 << std::endl;
        //std::cout << "a3: " << a3 << std::endl;

        // Construct the roots of the original cubic from those of the
        // depressed cubic
        x1=r1-b3/three/a3;
        x2=r2-b3/three/a3;
        x3=r3-b3/three/a3;

        //std::cout << "x1,x2,x3: " << x1 << " " << x2 << " " << x3
        //<< std::endl;
        //if (!o2isfinite(x1.real()) || !o2isfinite(x1.imag())) {
        //exit(-1);
        //}
        
      }
      
      return success;
    }      
    
    /// Return a string denoting the type ("cubic_complex_std")
    const char *type() { return "cubic_complex_std"; }
  };

  /** \brief Solve a quartic with real coefficients and real roots
   */
  template<class fp_t=double> 
  class quartic_real_std : public quartic_real<fp_t> {
    
  protected:

    cubic_real_coeff_gsl2<fp_t,std::complex<fp_t> > cub2;
    quadratic_real_coeff_gsl2<fp_t,std::complex<fp_t> > quad2;
    
  public:

    quartic_real_std() {
    }

    virtual ~quartic_real_std() {}

    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_std::solve_r().",
           exc_einval);
      }
      
      fp_t a34=b4/a4;
      fp_t a24=c4/a4;
      fp_t a14=d4/a4;
      fp_t a04=e4/a4;
      
      //---------------------------------------
      // Solve the resolvent cubic:
      
      fp_t a23=-a24;
      fp_t a13=(a14*a34-4.0*a04);
      fp_t a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);
      
      fp_t u1, u2, u3;
      cub2.solve_r(1,a23,a13,a03,u1,u2,u3);
      
      fp_t u4=u2;
      
      //---------------------------------------
      // Now construct the two quadratics:
      
      fp_t t1=u4+a34*a34/4.0-a24;
      if (t1>0.0) {
        t1=sqrt(t1);
      } else {
        t1=0.0;
      }
      
      fp_t b2a=-t1+a34/2.0;
      fp_t b2b=t1+a34/2.0;
      t1=u4*u4/4.0;
      
      // When numerical errors make t1 slightly smaller than a04.
      if (t1>a04) {
        t1=sqrt(t1-a04);
      } else {
        t1=0;
      }
      
      fp_t c2a=u4/2.0-t1;
      fp_t c2b=u4/2.0+t1;
      
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
      
      quad2.solve_r(1,b2a,c2a,x1,x2);
      quad2.solve_r(1,b2b,c2b,x3,x4);

      if (!std::isfinite(x1) || !std::isfinite(x2) ||
          !std::isfinite(x3) || !std::isfinite(x4)) {
        std::cout << u1 << " " << u2 << " " << u3 << std::endl;
        std::cout << t1 << " " << a34 << std::endl;
        std::cout << b2a << " " << b2b << std::endl;
        std::cout << c2a << " " << c2b << std::endl;
        std::cout << x1 << " " << x2 << std::endl;
        std::cout << x3 << " " << x4 << std::endl;
        exit(-1);
      }
      
      return success;
    }

    /// Return a string denoting the type ("quartic_real_std")
    const char *type() { return "quartic_real_std"; }

  };
  
  /** \brief Solve a quartic with complex coefficients and complex roots
   */
  template <class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_complex_std : public quartic_complex<fp_t,cx_t> {

  public:

    virtual ~quartic_complex_std() {}

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const cx_t a4, const cx_t b4, 
			const cx_t c4, const cx_t d4, 
			const cx_t e4, cx_t &x1, cx_t &x2, 
			cx_t &x3, cx_t &x4) {
      
      cx_t p4, q4, r4;
      cx_t a3, b3, c3, d3;
      cx_t b2a, c2a, b2b, c2b;
      cx_t u4, u41, u42;
      
      if (a4.real()==0.0 && a4.imag()==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_complex_std::solve_c().",
           exc_einval);
      }

      fp_t two=2.0;
      fp_t eight=8.0;
      fp_t four=4.0;
      fp_t sixteen=16.0;
      fp_t sixfour=64.0;
      fp_t twofivesix=256.0;
      fp_t three=3.0;
      
      p4=(eight*a4*c4-three*b4*b4)/eight/a4/a4;
      q4=(b4*b4*b4-four*a4*b4*c4+eight*a4*a4*d4)/eight/(a4*a4*a4);
      r4=(sixteen*a4*b4*b4*c4+twofivesix*a4*a4*a4*e4-three*b4*b4*b4*b4-
          sixfour*a4*a4*b4*d4)/twofivesix/(a4*a4*a4*a4);
      
      //---------------------------------------
      // Solve the resolvent cubic:
      
      a3=1.0;
      b3=-p4;
      c3=-four*r4;
      d3=four*p4*r4-q4*q4;
      
      cub_obj.solve_c(a3,b3,c3,d3,u4,u41,u42);
      
      //---------------------------------------
      
      // What to do when u4==p4?
      // Temporary hack:
      if (u4==p4) {
        b2a=0.0;
        b2b=0.0;
        c2a=u4/two;
        c2b=u4/two;
      } else {
        b2a=sqrt(u4-p4);
        b2b=-sqrt(u4-p4);
        c2a=-sqrt(u4-p4)*q4/two/(u4-p4)+u4/two;
        c2b=sqrt(u4-p4)*q4/two/(u4-p4)+u4/two;
      }
      
      x1=(-b2a+sqrt(b2a*b2a-four*c2a))/two-b4/four/a4;
      x2=(-b2a-sqrt(b2a*b2a-four*c2a))/two-b4/four/a4;
      x3=(-b2b+sqrt(b2b*b2b-four*c2b))/two-b4/four/a4;
      x4=(-b2b-sqrt(b2b*b2b-four*c2b))/two-b4/four/a4;
      
      return success;
    }

    /// Return a string denoting the type ("quartic_complex_std")
    const char *type() { return "quartic_complex_std"; }

#ifndef DOXYGEN_NO_O2NS

  protected:

    /// The object to solve for the associated cubic
    cubic_complex_std<fp_t,cx_t> cub_obj;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
