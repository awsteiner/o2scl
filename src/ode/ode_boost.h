/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_BOOST_H
#define O2SCL_ODE_BOOST_H

/** \file ode_boost.h
    \brief File defining \ref o2scl::ode_boost 
*/

#include <o2scl/err_hnd.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_step.h>
#include <boost/numeric/odeint.hpp>

namespace o2scl {

  /** \brief Simple ODE stepper from boost

      This is a simple implementation of the Boost "Error Steppers"
      into the O2scl framework. The slightly different interface
      requires a function wrapper (\ref ode_funct_boost) and an extra
      copy of the state vector. This works with the steppers \c
      runge_kutta_cash_karp54, \c runge_kutta_dopri5, and \c
      runge_kutta_fehlberg78. The Boost stepper \c
      runge_kutta_cash_karp54 gives identical results to \ref
      ode_rkck_gsl.
   */
  template<class step_t,
           class vec_y_t=boost::numeric::ublas::vector<double>,
	   class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
           class func_t=ode_funct,
           class fp_t=double> class ode_boost :
    public ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t,fp_t> {
    
  protected:

    /** \brief Stepper object
     */
    step_t stepper;
    
  public:

    ode_boost() {
    }
      
    virtual ~ode_boost() {
    }

    /** \brief Perform an integration step
     */
    virtual int step(fp_t x, fp_t h, size_t n, vec_y_t &y, 
                     vec_dydx_t &dydx, vec_y_t &yout, vec_yerr_t &yerr, 
                     vec_dydx_t &dydx_out, func_t &derivs) {

      ode_funct_boost<boost::numeric::ublas::vector<double>,
                      ode_funct,double> ofb(derivs,n);
      if (&yout!=&y) {
        for(size_t j=0;j<n;j++) {
          yout[j]=y[j];
        }
      }
      stepper.do_step(ofb,yout,x,h,yerr);
      derivs(x+h,n,yout,dydx_out);
    
      return 0;
    }
    
  };
  
}

#endif
