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
			 std::complex<double> &x3);
    
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
			  const double d, const double e);
    
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
			double &x2, double &x3, double &x4);

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
			double &x2, 
			double &x3, double &x4);

    /**  \brief Solves the polynomial 
	 \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	 giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	 \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const double a4, const double b4, const double c4, 
			 const double d4, const double e4, 
			 std::complex<double> &x1, std::complex<double> &x2, 
			 std::complex<double> &x3, std::complex<double> &x4);

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

    /// Return a string denoting the type ("poly_complex")
    const char *type() { return "poly_complex"; }
  };

  /** \brief Solve a cubic with real coefficients and complex roots 
      (CERNLIB)

      \note The function rrteq3() is based on the CERNLIB routine of
      the same name, but differs slightly. See the documentation of
      that function for details.
  */
  class cubic_real_coeff_cern : public cubic_real_coeff {

  public:

    cubic_real_coeff_cern() {
      eps=1.0e-6;
      delta=1.0e-15;
      improve_scale=true;
    }

    /// Numerical tolerance (default \f$ 10^{-6} \f$)
    double eps;

    /// Numerical tolerance (default \f$ 10^{-15} \f$)
    double delta;

    /// Improve algorithm for poorly-scaled roots (default true)
    bool improve_scale;

    virtual ~cubic_real_coeff_cern() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, 
			 std::complex<double> &x2, std::complex<double> &x3);

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
    virtual int rrteq3(double r, double s, double t, double x[], double &d);

    /// Return a string denoting the type ("cubic_real_coeff_cern")
    const char *type() { return "cubic_real_coeff_cern"; }
  };

  /** \brief Solve a quartic with real coefficients and complex 
      roots (CERNLIB)
  */
  class quartic_real_coeff_cern : public quartic_real_coeff {

  public:

    virtual ~quartic_real_coeff_cern() {}

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four complex solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const double a4, const double b4, const double c4, 
			 const double d4, const double e4, 
			 std::complex<double> &x1, std::complex<double> &x2, 
			 std::complex<double> &x3, std::complex<double> &x4);

    /// The CERNLIB-like interface
    virtual int rrteq4(double a, double b, double c, double d, 
		       std::complex<double> z[], double &dc, 
		       int &mt);

    /// Return a string denoting the type ("quartic_real_coeff_cern")
    const char *type() { return "quartic_real_coeff_cern"; }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The object to solve for the associated cubic
    cubic_real_coeff_cern cub_obj;

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
  class quadratic_complex_std : public quadratic_complex {

  public:

    virtual ~quadratic_complex_std() {}

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const std::complex<double> a2, 
			const std::complex<double> b2, 
			const std::complex<double> c2, 
			std::complex<double> &x1, std::complex<double> &x2);

    /// Return a string denoting the type ("quadratic_complex_std")
    const char *type() { return "quadratic_complex_std"; }
  };

  /** \brief Solve a cubic with complex coefficients and complex roots
   */
  class cubic_complex_std : public cubic_complex {

  public:

    virtual ~cubic_complex_std() {}

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
			std::complex<double> &x3);

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
  class quartic_complex_simple : public quartic_complex {

  public:

    virtual ~quartic_complex_simple() {}

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
			std::complex<double> &x2, 
			std::complex<double> &x3,
			std::complex<double> &x4);

    /// Return a string denoting the type ("quartic_complex_simple")
    const char *type() { return "quartic_complex_simple"; }

#ifndef DOXYGEN_NO_O2NS

  protected:

    /// The object to solve for the associated cubic
    cubic_complex_std cub_obj;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
