 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Jerry Gagelman and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_QAGIU_H
#define O2SCL_GSL_INTE_QAGIU_H

/** \file inte_qagiu_gsl.h
    \brief File defining \ref o2scl::inte_qagiu_gsl
*/
#include <o2scl/inte.h>
#include <o2scl/inte_qags_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Integrate a function over the interval \f$ [a, \infty) \f$ 
      (GSL)

      The integral on the unbounded interval is rewritten over the
      semi-open interval \f$ (0, 1] \f$ via a variable transformation,
      \f[
      \int_a^{\infty} f(x)~dx =
      \int_0^1 f(a + (1-t)/t)t^{-2}~dt,
      \f]
      and the right hand side is evaluated with \ref o2scl::inte_qags_gsl.

      \verbatim embed:rst
      See :ref:`GSL-based integration details` in the User's 
      guide for general information about the GSL integration classes.
      \endverbatim
  */
  template<class func_t=funct> class inte_qagiu_gsl : 
    public inte_transform_gsl<func_t> {
      
#ifndef DOXGYEN_INTERNAL

  protected:

  /// The lower limit
  double lower_limit;

#endif
      
  public:
  
  /** \brief Integrate a function over the interval \f$ [a, \infty) \f$ 
      giving result \c res and error \c err
  */
  virtual int integ_iu_err(func_t &func, double a,
			   double &res, double &err) {
    lower_limit=a;
    return this->qags(func,0.0,1.0,this->tol_abs,this->tol_rel,&res,&err);
  }

  /** \brief Integrate a function over the interval \f$ [a, \infty) \f$ 
      giving result \c res and error \c err

      The value \c b is ignored.
  */
  virtual int integ_err(func_t &func, double a, double b,
			double &res, double &err) {
    return integ_iu_err(func,a,res,err);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:
    
  /// Transform to \f$ t \in (0,1] \f$
  virtual double transform(double t, func_t &func) {
    double x=lower_limit+(1-t)/t, y=0.0;
    y=func(x);
    return y/t/t;
  }

#endif
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
