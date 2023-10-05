/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

#ifndef O2SCL_ASTEP_BOOST_H
#define O2SCL_ASTEP_BOOST_H

/** \file astep_boost.h
    \brief File defining \ref o2scl::astep_boost
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>

#include <o2scl/astep.h>

namespace o2scl {

  /** \brief Adaptive stepper based on Boost

      \note This class is experimental
  */
  template<class step_t,
           class vec_y_t=boost::numeric::ublas::vector<double>,
           class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
	   class func_t=ode_funct, class fp_t=double> class astep_boost : 
    public astep_base<vec_y_t,vec_dydx_t,vec_yerr_t,func_t,fp_t> {
      
    public:

    /** \brief The maximum relative uncertainty 
        (default \f$ 10^{-8} \f$)
    */
    double tol_rel;
    
    /** \brief The maximum absolute uncertainty (default \f$ 10^{-10} \f$)
     */
    double tol_abs;
    
    astep_boost() {
      tol_abs=1.0e-10;
      tol_rel=1.0e-8;
    }
      
    virtual ~astep_boost() {
    }

    /** \brief Make an adaptive integration step of the system 
	\c derivs 
	  
	This attempts to take a step of size \c h from the point \c x of
	an \c n-dimensional system \c derivs starting with \c y. On
	exit, \c x and \c y contain the new values at the end of the
	step, \c h contains the size of the step, \c dydx_out contains
	the derivative at the end of the step, and \c yerr contains the
	estimated error at the end of the step.

	If the base stepper returns a non-zero value, that value is
	passed on as the return value of <tt>astep()</tt>, though the
	input \c y may still have been modified by the base stepper.

    */
    virtual int astep(fp_t &x, fp_t xlimit, fp_t &h, 
		      size_t n, vec_y_t &y, vec_dydx_t &dydx_out,
		      vec_yerr_t &yerr, func_t &derivs) {

      ode_funct_boost<boost::numeric::ublas::vector<double>,
                      ode_funct,double> ofb(derivs,n);

      size_t cnt=
        integrate_adaptive(boost::numeric::odeint::make_controlled<step_t>
                           (tol_abs,tol_rel),ofb,y,x,x+h,h);
      x+=h;
      
      derivs(x,n,y,dydx_out);
      
      for(size_t j=0;j<n;j++) {
        yerr[j]=y[j]*tol_rel+tol_abs;
      }
      
      return 0;
    }
      
    /** \brief Make an adaptive integration step of the system 
	\c derivs with derivatives
	  
	This attempts to take a step of size \c h from the point \c x of
	an \c n-dimensional system \c derivs starting with \c y and
	given the initial derivatives \c dydx. On exit, \c x, \c y and
	\c dydx contain the new values at the end of the step, \c h
	contains the size of the step, \c dydx contains the derivative
	at the end of the step, and \c yerr contains the estimated error
	at the end of the step.

	If the base stepper returns a non-zero value, that value is
	passed on as the return value of <tt>astep()</tt>, though the
	inputs \c y and \c dydx may still have been modified by the
	base stepper.
    */
    virtual int astep_derivs(fp_t &x, fp_t xlimit, fp_t &h, 
			     size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			     vec_yerr_t &yerr, func_t &derivs) {

      ode_funct_boost<boost::numeric::ublas::vector<double>,
                      ode_funct,double> ofb(derivs,n);

      integrate_adaptive(boost::numeric::odeint::make_controlled<step_t>
                         (tol_abs,tol_rel),ofb,y,x,x+h,h);
      x+=h;
      
      derivs(x,n,y,dydx);
      
      for(size_t j=0;j<n;j++) {
        yerr[j]=y[j]*tol_rel+tol_abs;
      }
      
      return 0;
    }

    /** \brief Make an adaptive integration step of the system 
	\c derivs

	This function performs an adaptive integration step with the
	\c n-dimensional system \c derivs and parameter \c pa. It
	Begins at \c x with initial stepsize \c h, ensuring that the
	step goes no farther than \c xlimit. At the end of the step,
	the size of the step taken is \c h and the new value of \c x
	is in \c x_out. Initially, the function values and
	derivatives should be specified in \c y and \c dydx. The
	function values, derivatives, and the error at the end of
	the step are given in \c yout, \c yerr, and \c dydx_out.
	Unlike in \c ode_step objects, the objects \c y, \c yout, 
	\c dydx, and \c dydx_out must all be distinct.

	If the base stepper returns a non-zero value, that value is
	passed on as the return value of <tt>astep()</tt>, though the
	output parameters may still have been modified by the
	base stepper.
    */
    virtual int astep_full(fp_t x, fp_t xlimit, fp_t &x_out, fp_t &h, 
			   size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			   vec_y_t &yout, vec_yerr_t &yerr, 
			   vec_dydx_t &dydx_out, func_t &derivs) {
    
      ode_funct_boost<boost::numeric::ublas::vector<double>,
                      ode_funct,double> ofb(derivs,n);

      if (&yout!=&y) {
        for(size_t j=0;j<n;j++) {
          yout[j]=y[j];
        }
      }
      
      size_t cnt=integrate_adaptive
        (boost::numeric::odeint::make_controlled<step_t>
         (tol_abs,tol_rel),ofb,yout,x,x+h,h);
      x_out=x+h;
      
      derivs(x_out,n,yout,dydx_out);
      
      for(size_t j=0;j<n;j++) {
        yerr[j]=y[j]*tol_rel+tol_abs;
      }
      
      return 0;
    }
  
  };
  
}

#endif
