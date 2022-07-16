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
/** \file cheb_approx.h
    \brief File for definition of \ref o2scl::cheb_approx_tl
*/
#ifndef O2SCL_CHEB_APPROX_H
#define O2SCL_CHEB_APPROX_H

#include <cmath>
#include <o2scl/funct.h>
#include <o2scl/err_hnd.h>

#include <boost/numeric/ublas/vector.hpp>

namespace o2scl {

  /** \brief Chebyshev approximation (GSL)
      
      Approximate a function on a finite interval using a Chebyshev series:
      \f[
      f(x) = \sum_n c_n T_n(x)
      \f]
      where \f$ T_n(x)=\cos(n \arccos x) \f$

      \verbatim embed:rst
      See also the :ref:`Chebyshev approximation example`.
      \endverbatim
  */
  template<class fp_t=double> class cheb_approx_tl {

  public:

  typedef boost::numeric::ublas::vector<fp_t> ubvector;

  protected:

  /// Coefficients
  ubvector c;
  /// Order of the approximation
  size_t order;
  /// Lower end of the interval
  fp_t a;
  /// Upper end of the interval
  fp_t b;
  /// Single precision order
  size_t order_sp;
  /// Function evaluated at Chebyshev points
  ubvector f;
  /// True if init has been called
  bool init_called;

  /// Desc
  fp_t pi;
  
  public:
    
  cheb_approx_tl() {
    init_called=false;
    pi=boost::math::constants::pi<fp_t>();
  }

  /// Copy constructor
  cheb_approx_tl(const cheb_approx_tl &gc) {
    order=gc.order;
    init_called=gc.init_called;
    c=gc.c;
    f=gc.f;
    a=gc.a;
    b=gc.b;
    order_sp=gc.order_sp;
  }
    
  /// Copy constructor
  cheb_approx_tl &operator=(const cheb_approx_tl &gc) {

    // Check for self-assignment so that we don't
    // reallocate the vectors and lose the coefficients
    if (this==&gc) return *this;

    order=gc.order;
    init_called=gc.init_called;
    c=gc.c;
    f=gc.f;
    a=gc.a;
    b=gc.b;
    order_sp=gc.order_sp;
      
    return *this;
  }

  /// \name Initialization methods
  //@{
  /** \brief Initialize a Chebyshev approximation of the function
      \c func over the interval from \c a1 to \c b1

      The interval must be specified so that \f$ a < b \f$ , so 
      a and b are swapped if this is not the case.
  */
  template<class func_t> 
  void init(func_t &func, size_t ord, fp_t a1, fp_t b1) {
    size_t k, j;

    if(a1>=b1) {
      b=a1;
      a=b1;
    } else {
      a=a1;
      b=b1;
    }
    order=ord;
    order_sp=ord;
    c.resize(order+1);
    f.resize(order+1);

    fp_t two=2;
    fp_t half=1/two;
    fp_t bma=(b-a)/two;
    fp_t bpa=(b+a)/two;
    fp_t fac=two/((fp_t)(order+1));
    
    for(k=0;k<=order;k++) {
      fp_t y=cos(pi*(k+half)/(order+1));
      f[k]=func(y*bma+bpa);
    }
  
    for(j=0;j<=order;j++) {
      fp_t sum=0.0;
      for(k=0;k<=order;k++) {
	sum+=f[k]*cos(pi*j*(k+half)/((fp_t)(order+1)));
      }
      c[j]=fac*sum;
    }

    init_called=true;

    return;
  }

  /// Create an approximation from a vector of coefficients
  template<class vec_t> void init(fp_t a1, fp_t b1, 
				  size_t ord, vec_t &v) {
    order=ord;
    order_sp=order;
    a=a1;
    b=b1;
    c.resize(order+1);
    for(size_t i=0;i<order+1;i++) c[i]=v[i];

    init_called=true;

    return;
  }

  /// Create an approximation from a vector of function values
  template<class vec_t> void init_func_values(fp_t a1, fp_t b1, 
					      size_t ord, vec_t &fval) {
    size_t k, j;
      
    if(a>=b) {
      b=a1;
      a=b1;
    } else {
      a=a1;
      b=b1;
    }
    order=ord;
    order_sp=ord;
    c.resize(order+1);
    f.resize(order+1);
      
    fp_t two=2;
    fp_t half=1/two;
    fp_t bma=(b-a)/two;
    fp_t bpa=(b+a)/two;
    fp_t fac=two/((fp_t)(order+1));
      
    for(j=0;j<=order;j++) {
      fp_t sum=0.0;
      for(k=0;k<=order; k++) {
	sum+=fval[k]*cos(pi*j*(k+half)/((fp_t)(order+1)));
      }
      c[j]=fac*sum;
    }

    init_called=true;

    return;
  }
  //@}

  /// \name Evaulation methods
  //@{

  /** \brief Evaluate the approximation
   */
  fp_t eval(fp_t x) const {

    if (init_called==false) {
      O2SCL_ERR("Series not initialized in cheb_approx::eval()",
		o2scl::exc_einval);
      return 0.0;
    }

    size_t i;
    fp_t d1 = 0.0;
    fp_t d2 = 0.0;
      
    fp_t y = (2*x-a-b)/(b-a);
    fp_t y2 = 2*y;
      
    for (i=order; i >= 1; i--) {
      fp_t temp = d1;
      d1 = y2*d1-d2+c[i];
      d2 = temp;
    }
      
    return y*d1-d2+c[0]/2;
  }

  /** \brief Evaluate the approximation
   */
  fp_t operator()(fp_t x) const {
    return eval(x);
  }

  /** \brief Evaluate the approximation to a specified order
   */
  fp_t eval_n(size_t n, fp_t x) const {
    size_t i;
    fp_t d1 = 0.0;
    fp_t d2 = 0.0;
      
    size_t eval_order;
    if (n<order) eval_order=n;
    else eval_order=order;
      
    fp_t y = (2*x-a-b)/(b-a);
    fp_t y2 = 2*y;
      
    for (i = eval_order; i >= 1; i--) {
      fp_t temp = d1;
      d1 = y2*d1-d2+c[i];
      d2 = temp;
    }
      
    return y*d1-d2+c[0]/2;
  }
    
  /** \brief Evaluate the approximation and give the uncertainty
   */
  void eval_err(fp_t x, fp_t &result, fp_t &abserr) {

    size_t i;
    fp_t d1 = 0.0;
    fp_t d2 = 0.0;
      
    fp_t y = (2*x-a-b)/(b-a);
    fp_t y2 = 2*y;
      
    fp_t absc = 0.0;
      
    for (i = order; i >= 1; i--) {
      fp_t temp = d1;
      d1 = y2*d1-d2+c[i];
      d2 = temp;
    }
      
    result = y*d1-d2+c[0]/2;
      
    /* Estimate cumulative numerical error */
      
    for (i = 0; i <= order; i++) {
      absc += abs(c[i]);
    }
      
    /* Combine truncation error and numerical error */
      
    fp_t dbl_eps=std::numeric_limits<fp_t>::epsilon();

    abserr = abs(c[order])+absc*dbl_eps;
      
    return;
  }

  /** \brief Evaluate the approximation to a specified order
      and give the uncertainty
  */
  void eval_n_err(size_t n, fp_t x, fp_t &result, fp_t &abserr) {
    size_t i;
    fp_t d1 = 0.0;
    fp_t d2 = 0.0;

    fp_t y = (2*x-a-b)/(b-a);
    fp_t y2 = 2*y;

    fp_t absc = 0.0;

    size_t eval_order;
    if (n<order) eval_order=n;
    else eval_order=order;

    for (i = eval_order; i >= 1; i--) {
      fp_t temp = d1;
      d1 = y2*d1-d2+c[i];
      d2 = temp;
    }

    result = y*d1-d2+c[0]/2;

    /* Estimate cumulative numerical error */

    for (i = 0; i <= eval_order; i++) {
      absc += abs(c[i]);
    }

    fp_t dbl_eps=std::numeric_limits<fp_t>::epsilon();

    /* Combine truncation error and numerical error */
    abserr = abs(c[eval_order])+absc*dbl_eps;

    return;
  }
  //@}

  /// \name Maniupulating coefficients and endpoints
  //@{
  /** \brief Get a coefficient

      Legal values of the argument are 0 to \c order (inclusive)
  */
  fp_t get_coefficient(size_t ix) const {
    if (ix<order+1) {
      return c[ix];
    }
    O2SCL_ERR
    ("Requested invalid coefficient in cheb_approx::get_coefficient()",
     o2scl::exc_einval);
    return 0.0;
  }

  /** \brief Set a coefficient

      Legal values of the argument are 0 to \c order (inclusive)
  */
  void set_coefficient(size_t ix, fp_t co) {
    if (ix<order+1) {
      c[ix]=co;
      return;
    }
    O2SCL_ERR
    ("Requested invalid coefficient in cheb_approx::get_coefficient()",
     o2scl::exc_einval);
    return;
  }

  /// Return the endpoints of the approximation
  void get_endpoints(fp_t &la, fp_t &lb) {
    la=a;
    lb=b;
    return;
  }
    
  /** \brief Get the coefficients
   */
  template<class vec_t> void get_coefficients(size_t n, vec_t &v) const {
    for(size_t i=0;i<order+1 && i<n;i++) {
      v[i]=c[i];
    }
    return;
  }

  /** \brief Set the coefficients
   */
  template<class vec_t> void set_coefficients(size_t n, const vec_t &v) {
    for(size_t i=0;i<order+1 && i<n;i++) {
      c[i]=v[i];
    }
    return;
  }
  //@}

  /// \name Derivatives and integrals
  //@{
  /// Make \c gc an approximation to the derivative
  void deriv(cheb_approx_tl &gc) const {

    size_t n=order+1;
      
    const fp_t con=2/(b-a);
    size_t i;

    gc.init_called=true;
    gc.order=order;
    gc.order_sp=order;
    gc.a=a;
    gc.b=b;
    gc.c.resize(n);
    gc.f.resize(n);
  
    gc.c[n-1]=0.0;
      
    if (n > 1) {

      gc.c[n-2]=2*(n-1)*c[n-1];

      for(i=n;i>=3;i--) {
	gc.c[i-3]=gc.c[i-1]+2*(i-2)*c[i-2];
      }
	
      for(i=0;i<n;i++) {
	gc.c[i]*=con;
      }
    }
      
    return;
  }

  /// Make \c gc an approximation to the integral
  void integ(cheb_approx_tl &gc) const {

    size_t n=order+1;
      
    const fp_t con=(b-a)/4;
      
    gc.init_called=true;
    gc.order=order;
    gc.order_sp=order;
    gc.a=a;
    gc.b=b;
    gc.c.resize(n);
    gc.f.resize(n);
      
    if (n == 1) {

      gc.c[0]=0.0;

    } else if (n == 2) {

      gc.c[1]=con*c[0];
      gc.c[0]=2*gc.c[1];

    } else {

      fp_t sum=0;
      fp_t fac=1;

      for(size_t i=1;i<=n-2;i++) {
	gc.c[i]=con*(c[i-1]-c[i+1])/((fp_t)i);
	sum += fac*gc.c[i];
	fac=-fac;
      }
      gc.c[n-1]=con*c[n-2]/(n-1);
      sum += fac*gc.c[n-1];
      gc.c[0]=2*sum;
    }

    return;
  }
  //@}
    
  };

  /** \brief Double-precision version of \ref o2scl::cheb_approx_tl
   */
  typedef cheb_approx_tl<double> cheb_approx;

}

#endif
