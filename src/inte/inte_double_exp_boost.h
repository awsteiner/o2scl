/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_DOUBLE_EXP_BOOST_H
#define O2SCL_INTE_DOUBLE_EXP_BOOST_H

/** \file inte_tanh_sinh_boost.h
    \brief File defining \ref o2scl::inte_tanh_sinh_boost
*/

#include <cmath>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>

#include <o2scl/inte.h>

namespace o2scl {

  /** \brief Tanh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      The native range of the integrator is -1 to 1, but supports
      infinite limits on either (or both) sides.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_tanh_sinh_boost : public inte<func_t, fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::tanh_sinh<fp_t> it;
  
  public:

    inte_tanh_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_tanh_sinh_boost() {
    }

    /// Return string denoting type ("inte_tanh_sinh_boost")
    virtual const char *type() { return "inte_tanh_sinh_boost"; }
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,a,b,this->tol_rel/10.0,&err,&L1norm,&this->levels);
      if (err>this->tol_rel) {
        if (this->err_nonconv==true) {
          std::cout << "Function inte_tanh_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << this->levels << " " << max_refine
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_tanh_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \f$ \infty \f$
	and place the result in \c res and the error in \c err
    */
    virtual int integ_iu_err(func_t &func, fp_t a, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,a,std::numeric_limits<double>::infinity(),
		       res,err);
    }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
	and place the result in \c res and the error in \c err
    */
    virtual int integ_il_err(func_t &func, fp_t b, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,-std::numeric_limits<double>::infinity(),
		       b,res,err);
    }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
	\infty \f$ and place the result in \c res and the error in \c
	err
    */
    virtual int integ_i_err(func_t &func, 
			    fp_t &res, fp_t &err) {
      return integ_err(func,std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity(),res,err);
    }
  
    /** \brief Integrate function \c func from -1 to 1 and place
	the result in \c res and the error in \c err
    */
    virtual int integ_moo_err(func_t &func, fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (err>this->tol_rel) {
	std::cout << "Function inte_tanh_sinh_boost::integ_moo_err() failed."
                  << std::endl;
        std::cout << "Values err,tol_rel,L1norm,levels,max: "
		  << err << " " << this->tol_rel << " "
		  << L1norm << " " << levels << " " << max_refine
                  << std::endl;
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_tanh_sinh_boost::integ_moo_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }
  
    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };
  
  /** \brief Exp-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      Native range is 0 to \f$ \infty \f$, but 
      any semi-infinite range is supported.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_exp_sinh_boost : public inte<func_t, fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::exp_sinh<fp_t> it;
  
  public:

    inte_exp_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_exp_sinh_boost() {
    }
    
    /// Return string denoting type ("inte_exp_sinh_boost")
    virtual const char *type() { return "inte_exp_sinh_boost"; }
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,a,b,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (err>this->tol_rel) {
        if (this->err_nonconv==true) {
          std::cout << "Function inte_exp_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << levels << " " << max_refine
                    << std::endl;
        }
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_exp_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }
  
    /** \brief Integrate function \c func from \c a to \f$ \infty \f$
	and place the result in \c res and the error in \c err
    */
    virtual int integ_iu_err(func_t &func, fp_t a, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,a,std::numeric_limits<double>::infinity(),
		       res,err);
    }
  
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
	and place the result in \c res and the error in \c err
    */
    virtual int integ_il_err(func_t &func, fp_t b, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,-std::numeric_limits<double>::infinity(),
		       b,res,err);
    }
  
    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };
  
  /** \brief Sinh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      Only infinite limits (upper and lower) are supported.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_sinh_sinh_boost : public inte<func_t,fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::sinh_sinh<fp_t> it;
  
  public:

    inte_sinh_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_sinh_sinh_boost() {
    }
    
    /// Return string denoting type ("inte_sinh_sinh_boost")
    virtual const char *type() { return "inte_sinh_sinh_boost"; }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
	\infty \f$ and place the result in \c res and the error in \c
	err
    */
    virtual int integ_i_err(func_t &func, 
			    fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (err>this->tol_rel) {
        if (this->err_nonconv==true) {
          std::cout << "Function inte_sinh_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << levels << " " << max_refine
                    << std::endl;
        }
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_sinh_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }

    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };
  
}

#endif
