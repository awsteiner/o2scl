/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
      
      This base class does not perform any actual differentiation.
      Use one of the children cern_deriv, gsl_deriv, or eqi_deriv
      instead.

      This base class contains some code to automatically apply
      the first derivative routines to compute second or third
      derivatives. The error estimates for these will likely
      be underestimated.
      
      \note Because this class template aims to automatically provide
      second and third derivatives, one must overload either both
      calc() and calc_int() or both calc_err() and calc_err_int().
      
      \future Improve the methods for second and third derivatives
  */
  template<class func_t=funct> class deriv_base {
    
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
    virtual double deriv(double x, func_t &func) {
      double dx;
      from_deriv=true;
      deriv_err(x,func,dx,derr);
      from_deriv=false;
      return dx;
    }

    /** \brief Calculate the second derivative of \c func w.r.t. x
     */
    virtual double deriv2(double x, func_t &func) {
      double val;
      dpars dp={&func};
      funct_mfptr_param<deriv_base,dpars> mf(this,&deriv_base::derivfun,dp);
      val=deriv_int(x,mf);
      // The error estimate is unavailable, so we set it to zero 
      derr=0.0;
      return val;
    }
  
    /** \brief Calculate the third derivative of \c func w.r.t. x
     */
    virtual double deriv3(double x, func_t &func) {
      double val;
      dpars dp={&func};
      funct_mfptr_param<deriv_base,dpars> mf(this,&deriv_base::derivfun2,dp);
      val=deriv_int(x,mf);
      // The error estimate is unavailable, so we set it to zero 
      derr=0.0;
      return val;
    }

    /** \brief Get uncertainty of last calculation
     */
    virtual double get_err() {
      return derr;
    }

    /** \brief Output control
     */
    int verbose;
    
    /** \brief Calculate the first derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv_err(double x, func_t &func, double &dfdx, 
			 double &err)=0;

    /** \brief Calculate the second derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv2_err(double x, func_t &func, 
			  double &d2fdx2, double &err) {
      int ret;
      dpars dp={&func};
      funct_mfptr_param<deriv_base,dpars> mf(this,&deriv_base::derivfun,dp);
      ret=deriv_err_int(x,mf,d2fdx2,err);
      // The error estimate is unavailable, so we set it to zero 
      err=0.0;
      return 0;
    }
    
    /** \brief Calculate the third derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv3_err(double x, func_t &func, 
			  double &d3fdx3, double &err) {
      int ret;
      dpars dp={&func};
      funct_mfptr_param<deriv_base,dpars> mf(this,&deriv_base::derivfun2,dp);
      ret=deriv_err_int(x,mf,d3fdx3,err);
      // The error estimate is unavailable, so we set it to zero 
      err=0.0;
      return 0;
    }
    
    /// Return string denoting type ("deriv")
    virtual const char *type() { return "deriv"; }
    
  protected:

#ifndef DOXYGEN_INTERNAL
    
    /** \brief Calculate the first derivative of \c func w.r.t. x
	
	This is an internal version of deriv() which is used in
	computing second and third derivatives
    */
    virtual double deriv_int(double x, funct &func) {
      double dx;
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
    virtual int deriv_err_int(double x, funct &func, 
			     double &dfdx, double &err)=0;
    
    /// The uncertainity in the most recent derivative computation
    double derr;

    /// The function for the second derivative
    double derivfun(double x, dpars &dp) {
      return deriv(x,*(dp.func));
    }

    /// The function for the third derivative
    double derivfun2(double x, dpars &dp) {
      funct_mfptr_param<deriv_base,dpars> mf(this,&deriv_base::derivfun,dp);
      double val=deriv_int(x,mf);
      return val;
    }
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
