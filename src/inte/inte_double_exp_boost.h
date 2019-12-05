/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019, Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_TANH_SINH_BOOST_H
#define O2SCL_INTE_TANH_SINH_BOOST_H

/** \file inte_tanh_sinh_boost.h
    \brief File defining \ref o2scl::inte_tanh_sinh_boost
*/

#include <cmath>

#ifndef O2SCL_OLDER_COMPILER

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>

#include <o2scl/inte.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Tanh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
    class inte_tanh_sinh_boost :
    public inte<func_t, fp_t> {
    
  protected:

  /// The boost integration object
  boost::math::quadrature::tanh_sinh<fp_t> it;
  
  /// L1 norm
  fp_t L1norm;

  public:

  inte_tanh_sinh_boost() : it(max_refine) {
  }
  
  virtual ~inte_tanh_sinh_boost() {
  }
    
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			fp_t &res, fp_t &err) {
    res=it.integrate(func,a,b,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      O2SCL_ERR2("Failed to achieve tolerance in ",
		 "inte_tanh_sinh_boost::integ_err().",o2scl::exc_efailed);
    }
    return 0;
  }
  
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, 
			fp_t &res, fp_t &err) {
    res=it.integrate(func,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      O2SCL_ERR2("Failed to achieve tolerance in ",
		 "inte_tanh_sinh_boost::integ_err().",o2scl::exc_efailed);
    }
    return 0;
  }
  
  };
  
  /** \brief Exp-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
    class inte_exp_sinh_boost :
    public inte<func_t, fp_t> {
    
  protected:

  /// The boost integration object
  boost::math::quadrature::exp_sinh<fp_t> it;
  
  /// L1 norm
  fp_t L1norm;

  public:

  inte_exp_sinh_boost() : it(max_refine) {
  }
  
  virtual ~inte_exp_sinh_boost() {
  }
    
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			fp_t &res, fp_t &err) {
    res=it.integrate(func,a,b,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      O2SCL_ERR2("Failed to achieve tolerance in ",
		 "inte_exp_sinh_boost::integ_err().",o2scl::exc_efailed);
    }
    return 0;
  }
  
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, 
			fp_t &res, fp_t &err) {
    res=it.integrate(func,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      O2SCL_ERR2("Failed to achieve tolerance in ",
		 "inte_exp_sinh_boost::integ_err().",o2scl::exc_efailed);
    }
    return 0;
  }
  
  };
  
  /** \brief Sinh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
    class inte_sinh_sinh_boost {
    
  protected:

  /// The boost integration object
  boost::math::quadrature::sinh_sinh<fp_t> it;
  
  /// L1 norm
  fp_t L1norm;

  public:

  inte_sinh_sinh_boost() : it(max_refine) {
  }
  
  virtual ~inte_sinh_sinh_boost() {
  }
    
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, 
			fp_t &res, fp_t &err) {
    res=it.integrate(func,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      O2SCL_ERR2("Failed to achieve tolerance in ",
		 "inte_sinh_sinh_boost::integ_err().",o2scl::exc_efailed);
    }
    return 0;
  }
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

// End of #ifndef O2SCL_OLDER_COMPILER
#endif

#endif