/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file astep.h
    \brief Adaptive stepper base class
*/
#ifndef O2SCL_ASTEP_H
#define O2SCL_ASTEP_H

#include <o2scl/ode_step.h>
#include <o2scl/ode_rkck_gsl.h>
#include <o2scl/ode_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Adaptive stepper [abstract base]
      
      The adaptive stepper routines are based on one or many
      applications of ordinary ODE steppers (implemented in \ref
      ode_step). Each adaptive stepper (\ref astep_gsl or \ref
      astep_nonadapt) can be used with any of the ODE stepper classes
      (e.g. \ref ode_rkck_gsl). By default, ode_rkck_gsl is used. To modify
      the ODE stepper which is used, use the member function
      \ref set_step() documented below.
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
    class func_t=ode_funct11 > class astep_base {
    
  public:
  
  astep_base() {
    stepp=&def_step;
    verbose=0;
  }

  virtual ~astep_base() {}
      
  /** \brief Make an adaptive integration step of the system 
      \c derivs 
	  
      This attempts to take a step of size \c h from the point \c
      x of an \c n-dimensional system \c derivs starting with \c
      y. On exit, \c x and \c y contain the new values at the end
      of the step, \c h contains the size of the step, \c dydx_out
      contains the derivative at the end of the step, and \c yerr
      contains the estimated error at the end of the step.
  */
  virtual int astep(double &x, double xlimit, double &h, 
		    size_t n, vec_y_t &y, vec_dydx_t &dydx_out,
		    vec_yerr_t &yerr, func_t &derivs)=0;

  /** \brief Make an adaptive integration step of the system 
      \c derivs with derivatives
	  
      This attempts to take a step of size \c h from the point \c
      x of an \c n-dimensional system \c derivs starting with \c y
      and given the initial derivatives \c dydx. On exit, \c x, \c
      y and \c dydx contain the new values at the end of the step,
      \c h contains the size of the step, \c dydx contains the
      derivative at the end of the step, and \c yerr contains the
      estimated error at the end of the step.
  */
  virtual int astep_derivs(double &x, double xlimit, double &h, 
			   size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			   vec_yerr_t &yerr, func_t &derivs)=0;

  /** \brief Make an adaptive integration step of the system 
      \c derivs with derivatives 

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
  */
  virtual int astep_full(double x, double xlimit, double &x_out, 
			 double &h, size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			 vec_y_t &yout, vec_yerr_t &yerr, 
			 vec_dydx_t &dydx_out, func_t &derivs)=0;

  /// Set output level
  int verbose;

  /** \brief Set stepper

      This sets the stepper for use in the adaptive step 
      routine. If no stepper is specified, then the default
      (\ref def_step of type \ref ode_rkck_gsl) is used.
  */
  int set_step(ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t> &step) {
    stepp=&step;
    return 0;
  }
  
  /// The default stepper
  ode_rkck_gsl<vec_y_t,vec_dydx_t,vec_yerr_t,func_t> def_step;
  
#ifndef DOXYGEN_INTERNAL

  protected:
      
  /// Pointer to the stepper being used
  ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t> *stepp;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
