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

#ifndef O2SCL_ODE_STEP_H
#define O2SCL_ODE_STEP_H

/** \file ode_step.h
    \brief File defining \ref o2scl::ode_step 
*/

#include <o2scl/ode_funct.h>

namespace o2scl {
  
  /** \brief ODE stepper base [abstract base]
   */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t,
    class func_t=ode_funct, class fp_t=double> class ode_step {

#ifndef DOXGYENP

    protected:
    
    /** \brief The order of the ODE stepper
     */
    int order;

#endif
    
    public:

    ode_step() {
      order=0;
    }

    virtual ~ode_step() {}
    
    /** \brief Return the order of the ODE stepper
	
        This is used, for example, by \ref astep_gsl to adaptively
        adjust the stepsize.
    */
    virtual int get_order() {
      return order;
    }
    
    /** \brief Perform an integration step
	
	Given initial value of the n-dimensional function in \c y and
	the derivative in \c dydx (which must generally be computed
	beforehand) at the point \c x, take a step of size \c h giving
	the result in \c yout, the uncertainty at \f$ x+h \f$ in \c
	yerr, and the new derivative at \f$ x+h \f$ in \c dydx_out
	using function \c derivs to calculate derivatives.
	Implementations which do not calculate \c yerr and/or \c
	dydx_out do not reference these variables so that a blank \c
	vec_t can be given. All of the current \o2 implementations
	allow \c yout=y and \c dydx_out=dydx if necessary
    */
    virtual int step(fp_t x, fp_t h, size_t n, vec_y_t &y, 
		     vec_dydx_t &dydx, vec_y_t &yout, vec_yerr_t &yerr, 
		     vec_dydx_t &dydx_out, func_t &derivs)=0;
    
  };
  
}

#endif
