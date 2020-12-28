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
#ifndef O2SCL_DERIV_H
#define O2SCL_DERIV_H

/** \file deriv.h
    \brief File defining \ref o2scl::deriv_base
*/

#include <iostream>
#include <cmath>
#include <o2scl/funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Numerical differentiation base [abstract base]
      
      This base class does not perform any actual differentiation. Use
      one of the children \ref o2scl::deriv_cern, \ref
      o2scl::deriv_gsl, or \ref o2scl::deriv_eqi instead.

      This base class contains some code to automatically apply
      the first derivative routines to compute second or third
      derivatives. The error estimates for these will likely
      be underestimated.
      
      \note Because this class template aims to automatically provide
      second and third derivatives, one must overload either both
      deriv() and deriv_int() or both deriv_err() and deriv_err_int().

      \note If \ref err_nonconv is set to false, and the derivative
      computation fails, then the functions \ref deriv(),
      \ref deriv2() and \ref deriv3() may return the wrong result
      without warning. Similarly, if \ref err_nonconv is set to
      false, it is the user's responsibility to check the
      return value from \ref deriv_err(), \ref deriv2_err(), and
      \ref deriv3_err() to see if an error occurred.
      
      \future Improve the methods for second and third derivatives
  */
  template<class func_t=funct, class fp_t=double> class deriv_base {
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /** \brief A structure for passing the function to second and 
	third derivatives [protected]
    */
    typedef struct {
      
    public:

      /// The pointer to the function
      func_t *func;

    } dpars;
    
    /// Avoids infinite loops in case the user calls the base class version
    bool from_deriv;

    typedef std::function<fp_t(fp_t)> internal_func_t;
    
#endif

  public:

    deriv_base() {
      verbose=0;
      from_deriv=false;
      err_nonconv=true;
    }

    virtual ~deriv_base() {}

    /// If true, call the error handler if the routine does not "converge"
    bool err_nonconv;

    /** \brief Calculate the first derivative of \c func w.r.t. x
	
	After calling deriv(), the error may be obtained from 
	\ref get_err().
    */
    virtual fp_t deriv(fp_t x, func_t &func) {
      // There were some uninitialized variable warnings on OSX, so we
      // prevent those by setting the derivative equal to zero.
      fp_t dx=0.0;
      from_deriv=true;
      deriv_err(x,func,dx,derr);
      from_deriv=false;
      return dx;
    }

    /** \brief Calculate the second derivative of \c func w.r.t. x
     */
    virtual fp_t deriv2(fp_t x, func_t &func) {
      fp_t val;
      
      internal_func_t mf=
	std::bind(std::mem_fn<fp_t(fp_t,func_t *)>(&deriv_base::derivfun),
		  this,std::placeholders::_1,&func);
      
      val=deriv_int(x,mf);
      // The error estimate is unavailable, so we set it to zero 
      derr=0.0;
      return val;
    }
  
    /** \brief Calculate the third derivative of \c func w.r.t. x
     */
    virtual fp_t deriv3(fp_t x, func_t &func) {
      fp_t val;

      internal_func_t mf=
      std::bind(std::mem_fn<fp_t(fp_t,func_t *)>(&deriv_base::derivfun2),
		this,std::placeholders::_1,&func);

      val=deriv_int(x,mf);
      // The error estimate is unavailable, so we set it to zero 
      derr=0.0;
      return val;
    }

    /** \brief Get uncertainty of last calculation
     */
    virtual fp_t get_err() {
      return derr;
    }

    /** \brief Output control
     */
    int verbose;
    
    /** \brief Calculate the first derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv_err(fp_t x, func_t &func, fp_t &dfdx, 
			 fp_t &err)=0;

    /** \brief Calculate the second derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv2_err(fp_t x, func_t &func, 
			  fp_t &d2fdx2, fp_t &err) {
      internal_func_t mf=
      std::bind(std::mem_fn<fp_t(fp_t,func_t *)>(&deriv_base::derivfun),
		this,std::placeholders::_1,&func);

      int ret=deriv_err_int(x,mf,d2fdx2,err);
      // The error estimate is unavailable, so we set it to zero 
      err=0.0;
      return ret;
    }
    
    /** \brief Calculate the third derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv3_err(fp_t x, func_t &func, 
			  fp_t &d3fdx3, fp_t &err) {
      internal_func_t mf=
      std::bind(std::mem_fn<fp_t(fp_t,func_t *)>(&deriv_base::derivfun2),
		this,std::placeholders::_1,&func);
      
      int ret=deriv_err_int(x,mf,d3fdx3,err);
      // The error estimate is unavailable, so we set it to zero 
      err=0.0;
      return ret;
    }
  
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /// Return string denoting type ("deriv")
    virtual const char *type() { return "deriv"; }
    
  protected:
    
#ifndef DOXYGEN_INTERNAL
    
    /** \brief Calculate the first derivative of \c func w.r.t. x
	
	This is an internal version of deriv() which is used in
	computing second and third derivatives
    */
    virtual fp_t deriv_int(fp_t x, internal_func_t &func) {
      fp_t dx;
      from_deriv=true;
      deriv_err_int(x,func,dx,derr);
      from_deriv=false;
      return dx;
    }
    
    /** \brief Calculate the first derivative of \c func w.r.t. x and the
	uncertainty
	
	This is an internal version of deriv_err() which is used in
	computing second and third derivatives
    */
    virtual int deriv_err_int(fp_t x, internal_func_t &func, 
			      fp_t &dfdx, fp_t &err)=0;
    
    /// The uncertainity in the most recent derivative computation
    fp_t derr;
    
    /// The function for the second derivative
    fp_t derivfun(fp_t x, func_t *fp) {
      return deriv(x,*fp);
    }
    
    /// The function for the third derivative
    fp_t derivfun2(fp_t x, func_t *fp) {
      internal_func_t mf=
	std::bind(std::mem_fn<fp_t(fp_t,func_t *)>(&deriv_base::derivfun),
		  this,std::placeholders::_1,fp);
      fp_t val=deriv_int(x,mf);
      return val;
    }
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
