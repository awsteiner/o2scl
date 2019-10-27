/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner and Jerry Gagelman
  
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
#ifndef O2SCL_INTE_H
#define O2SCL_INTE_H

/** \file inte.h
    \brief File defining \ref o2scl::inte
*/

#include <cmath>
#include <o2scl/err_hnd.h>
#include <o2scl/funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Base integration class [abstract base]

      \note Currently \o2 supports only types \c double and, for some
      integration methods, \c long \c double for the floating point
      type \c fp_t . Also, the default values of \ref tol_rel 
      and \ref tol_abs are designed for double precision and
      likely need to be decreased for long double precision
      integration.

      \future It might be useful to have all of the integration
      classes report the number of function evaluations used
      in addition to the number of iterations which were taken.
  */
  template<class func_t=funct, class fp_t=double> class inte {
    
  public:
  
  inte() {
    tol_rel=1.0e-8;
    tol_abs=1.0e-8;
    verbose=0;
    interror=0.0;
    err_nonconv=true;
  }

  virtual ~inte() {}

  /** \brief Verbosity
   */
  int verbose;

  /// The most recent number of iterations taken
  size_t last_iter;
  
  /** \brief The maximum relative uncertainty 
      in the value of the integral (default \f$ 10^{-8} \f$)
  */
  fp_t tol_rel;
    
  /** \brief The maximum absolute uncertainty 
      in the value of the integral (default \f$ 10^{-8} \f$)
  */
  fp_t tol_abs;
  
  /** \brief If true, call the error handler if the routine does not 
      converge or reach the desired tolerance (default true)

      If this is false, the function proceeds normally and 
      may provide convergence information in the integer return value.
  */
  bool err_nonconv;
  
  /** \brief Integrate function \c func from \c a to \c b.
   */
  virtual fp_t integ(func_t &func, fp_t a, fp_t b) {
    fp_t res;
    int ret=integ_err(func,a,b,res,this->interror);
    if (ret!=0) {
      O2SCL_ERR2("Integration failed in inte::integ(), ",
		 "but cannot return int.",o2scl::exc_efailed);
    }
    return res;
  }

  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			fp_t &res, fp_t &err)=0;
  
  /** \brief Return the numerically estimated error in the result from
      the last call to integ()
      
      This will quietly return zero if no integrations have been
      performed or if the integrator does not estimate the error.
  */
  fp_t get_error() { return interror; }
  
  /// Return string denoting type ("inte")
  virtual const char *type() { return "inte"; }
  
#ifndef DOXYGEN_INTERNAL
  
  protected:
  
  /// The uncertainty for the last integration computation
  fp_t interror;
  
#endif
  
  };
  
  /** \brief Integrate over \f$ (-\infty,b] \f$

      \note This class only works if the base integration
      type <tt>def_inte_t</tt> avoids evaluating the function
      at the left-hand end point.
   */
  template<class func_t, class def_inte_t, class fp_t=double>
    class inte_il {
    
  protected:

  /// A pointer to the user-specified function
  func_t *user_func;
  
  /// The upper limit
  fp_t upper_limit;

  /// Transform from \f$ t \in (0,1] \f$ to \f$ x \in (-\infty,b] \f$
  virtual fp_t transform(fp_t t) {
    if (t==0.0) {
      O2SCL_ERR2("Function called with t=0 in ",
		"inte_il::transform().",o2scl::exc_efailed);
    }
    fp_t x=upper_limit-(1-t)/t, y;
    y=(*user_func)(x);
    return y/t/t;
  }

  public:
  
  inte_il() {
    it=&def_inte;
    fo=std::bind(std::mem_fn<fp_t(fp_t)>(&inte_il::transform),
		 this,std::placeholders::_1);
  }    
  
  /** \brief Internal function type based on floating-point type

      \comment 
      This type must be public so the user can change the
      base integration object
      \endcomment
  */
  typedef std::function<fp_t(fp_t)> internal_funct;

  protected:

  /// The base integration object
  inte<internal_funct,fp_t> *it;
  
  /// Function object
  internal_funct fo;
  
  public:  
  
  /** \brief Integrate function \c func from \c a to \c b
      giving result \c res and error \c err
  */
  virtual int integ_err(func_t &func, fp_t b,
			fp_t &res, fp_t &err) {
    user_func=&func;
    upper_limit=b;
    int ret=it->integ_err(fo,0.0,1.0,res,err);
    return ret;
  }
  
  /// \name Integration object
  //@{
  /// Set the base integration object to use
  int set_inte(inte<internal_funct,fp_t> &i) {
    it=&i;
    return 0;
  }
      
  /// Default integration object
  def_inte_t def_inte;
  //@}

  };
  
  /** \brief Integrate over \f$ [a,\infty) \f$

      \note This class only works if the base integration
      type <tt>def_inte_t</tt> avoids evaluating the function
      at the right-hand end point.
   */
  template<class func_t, class def_inte_t, class fp_t=double>
    class inte_iu {
    
  protected:

  /// A pointer to the user-specified function
  func_t *user_func;
  
  /// The lower limit
  fp_t lower_limit;

  /// Transform from \f$ t \in (0,1] \f$ to \f$ x \in [a,\infty) \f$
  virtual fp_t transform(fp_t t) {
    if (t==0.0) {
      O2SCL_ERR2("Function called with t=0 in ",
		 "inte_ul::transform().",o2scl::exc_efailed);
    }
    fp_t x=lower_limit+(1-t)/t, y;
    y=(*user_func)(x);
    return y/t/t;
  }

  public:
  
  inte_iu() {
    it=&def_inte;
    fo=std::bind(std::mem_fn<fp_t(fp_t)>(&inte_iu::transform),
		 this,std::placeholders::_1);
  }    
  
  /** \brief Internal function type based on floating-point type
      
      \comment 
      This type must be public so the user can change the
      base integration object
      \endcomment
  */
  typedef std::function<fp_t(fp_t)> internal_funct;

  protected:

  /// The base integration object
  inte<internal_funct,fp_t> *it;
  
  /// Function object
  internal_funct fo;
  
  public:  
  
  /** \brief Integrate function \c func from \c a to \c b
      giving result \c res and error \c err
  */
  virtual int integ_err(func_t &func, fp_t a, 
			fp_t &res, fp_t &err) {
    user_func=&func;
    lower_limit=a;
    int ret=it->integ_err(fo,0.0,1.0,res,err);
    return ret;
  }
  
  /// \name Integration object
  //@{
  /// Set the base integration object to use
  int set_inte(inte<internal_funct,fp_t> &i) {
    it=&i;
    return 0;
  }
      
  /// Default integration object
  def_inte_t def_inte;
  //@}

  };
  
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif




