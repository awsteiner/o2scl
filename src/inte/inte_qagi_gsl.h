 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Jerry Gagelman and Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_GSL_INTE_QAGI_H
#define O2SCL_GSL_INTE_QAGI_H

/** \file inte_qagi_gsl.h
    \brief File defining \ref o2scl::inte_qagi_gsl
*/
#include <o2scl/inte.h>
#include <o2scl/inte_qags_gsl.h>

namespace o2scl {

  /** \brief Integrate a function over the interval \f$
      (-\infty,\infty) \f$ (GSL)
      
      \verbatim embed:rst
      See :ref:`GSL-based integration details` in the User's 
      guide for general information about the GSL integration classes.
      \endverbatim
  */
  template<class func_t=funct > class inte_qagi_gsl : 
    public inte_transform_gsl<func_t> {
    
  public:
      
  /** \brief Integrate function \c func from \f$ \infty \f$ 
      to \f$ \infty \f$ giving result \c res and error \c err
  */
  virtual int integ_i_err(func_t &func, double &res, double &err) {
    return this->qags(func,0.0,1.0,this->tol_abs,this->tol_rel,&res,&err);
  }
    
  /** \brief Integrate function \c func from \f$ \infty \f$ 
      to \f$ \infty \f$ giving result \c res and error \c err
	
      The values \c a and \c b are ignored
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &res, double &err) {
    return integ_i_err(func,res,err);
  }
    
  /// Return string denoting type ("inte_qagi_gsl")
  const char *type() { return "inte_qagi_gsl"; }

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// Tranformation to \f$ t \in (0,1] \f$ 
  virtual double transform(double t, func_t &func) {
    double x=(1-t)/t, y1, y2;
    y1=func(x);
    y2=func(-x);
    return (y1+y2)/t/t;
  }
    
#endif
      
  };
  
}

#endif
