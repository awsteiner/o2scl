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

#ifndef O2SCL_ASTEP_NONADAPT_H
#define O2SCL_ASTEP_NONADAPT_H

/** \file astep_nonadapt.h
    \brief File defining \ref o2scl::astep_nonadapt
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/astep.h>

namespace o2scl {

  /** \brief An non-adaptive stepper implementation of \ref astep_base.

      This class simply calls the specified ODE stepper without any
      attempt to modify the size of the step. It is primarily useful
      to allow for simple comparisons between adaptive and
      non-adaptive solution or to proceed with a solution in places
      where an adaptive stepper is failing.

      To modify the ODE stepper which is used, use the function
      astep_base::set_step().
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
	   class func_t=ode_funct, class fp_t=double> class astep_nonadapt : 
    public astep_base<vec_y_t,vec_dydx_t,vec_yerr_t,func_t,fp_t> {
      
#ifndef DOXYGEN_INTERNAL
      
    protected:

    /// The allocated vector size
    size_t msize;

    /** \brief Internal storage for dydx
      
	\comment
	This isn't named dydx so that we can name function arguments
	dydx in astep_derivs().
	\endcomment
    */
    vec_dydx_t dydx_int;

#endif
    
    public:
    
    astep_nonadapt()  {
      msize=0;
    }
      
    virtual ~astep_nonadapt() {
      if (msize>0) {
	dydx_int.clear();
      }
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

      if (n!=msize) {
	dydx_int.resize(n);
      }
      
      derivs(x,n,y,dydx_int);
      if (h>0) {
	if (x+h>xlimit) h=xlimit-x;
      } else {
	if (x+h<xlimit) h=xlimit-x;
      }
      int ret=this->stepp->step(x,h,n,y,dydx_int,y,yerr,dydx_out,derivs);
      x+=h;

      return ret;
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
      
      if (h>0) {
	if (x+h>xlimit) h=xlimit-x;
      } else {
	if (x+h<xlimit) h=xlimit-x;
      }
      int ret=this->stepp->step(x,h,n,y,dydx,y,yerr,dydx,derivs);
      x+=h;

      return ret;
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
      
      if (h>0) {
	if (x+h>xlimit) h=xlimit-x;
      } else {
	if (x+h<xlimit) h=xlimit-x;
      }
      int ret=this->stepp->step(x,h,n,y,dydx,yout,yerr,dydx_out,derivs);
      x_out=x+h;
    
      return ret;
    }
  
  };
  
}

#endif
