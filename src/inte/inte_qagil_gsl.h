/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Jerry Gagelman and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_QAGIL_H
#define O2SCL_GSL_INTE_QAGIL_H

/** \file inte_qagil_gsl.h
    \brief File defining \ref o2scl::inte_qagil_gsl
*/
#include <o2scl/inte.h>
#include <o2scl/inte_qags_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Integrate a function over the interval \f$ 
      (-\infty, b] \f$ (GSL)
      
      The integral on the unbounded interval is rewritten over the
      semi-open interval \f$ (0, 1] \f$ via a variable transformation,
      \f[
      \int_{-\infty}^b f(x)~dx =
      \int_0^1 f\left[b - (1-t)/t\right]t^{-2}~dt,
      \f]
      and the right hand side is evaluated with \ref o2scl::inte_qags_gsl.

      \verbatim embed:rst
      See :ref:`GSL-based integration details` in the User's 
      guide for general information about the GSL integration classes.
      \endverbatim
  */
  template<class func_t=funct > class inte_qagil_gsl : 
    public inte_transform_gsl<func_t> {

#ifndef DOXYGEN_INTERNAL      

  protected:

  /// The upper limit
  double upper_limit;

#endif

  public:
    
  /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
      and place the result in \c res and the error in \c err
  */
  virtual int integ_il_err(func_t &func, double b, 
			double &res, double &err) {
    upper_limit=b;
    return this->qags(func,0.0,1.0,this->tol_abs,this->tol_rel,&res,&err);
  }

  /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
      and place the result in \c res and the error in \c err
      
      The value given is \c a is ignored.
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &res, double &err) {
    return integ_il_err(func,b,res,err);
  }

  /// Return string denoting type ("inte_qagil_gsl")
  const char *type() { return "inte_qagil_gsl"; }

  protected:
      
  /// Transform to \f$ t \in (0,1] \f$
  virtual double transform(double t, func_t &func) {
    double x=upper_limit-(1-t)/t, y;
    y=func(x);
    return y/t/t;
  }
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
