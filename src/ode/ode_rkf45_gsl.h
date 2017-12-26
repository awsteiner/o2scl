/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
/* ode-initval/rkf45.c
 * 
 * Copyright (C) 2001, 2004, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_GSL_RKF45_H
#define O2SCL_GSL_RKF45_H

/** \file ode_rkf45_gsl.h
    \brief File defining \ref o2scl::ode_rkf45_gsl 
*/

#include <o2scl/err_hnd.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_step.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Runge-Kutta-Fehlberg embedded Runge-Kutta ODE stepper (GSL)

      Based on \ref Hairer09 .

      \todo Check this because it may not give exact dydt_out.
   */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
    class func_t=ode_funct>
    class ode_rkf45_gsl : public ode_step<vec_y_t,
    vec_dydx_t,vec_yerr_t,func_t> {
    
  protected:
  
  /// \name Storage for the intermediate steps
  //@{
  vec_dydx_t k2, k3, k4, k5, k6;
  vec_y_t ytmp;
  //@}
  
  /// Size of allocated vectors
  size_t ndim;
  
  /** \name Storage for the coefficients
   */
  //@{
  double ah[5], b3[2], b4[3], b5[4], b6[5];
  double c1, c3, c4, c5, c6;
  double ec[7];
  //@}
  
  public:
  
  ode_rkf45_gsl() {
    this->order=5;
    
    ah[0]=1.0/4.0;
    ah[1]=3.0/8.0;
    ah[2]=12.0/13.0;
    ah[3]=1.0;
    ah[4]=1.0/2.0;

    b3[0]=3.0/32.0;
    b3[1]=9.0/32.0;

    b4[0]=1932.0/2197.0;
    b4[1]=-7200.0/2197.0;
    b4[2]=7296.0/2197.0;

    b5[0]=8341.0/4104.0;
    b5[1]=-32832.0/4104.0;
    b5[2]=29440.0/4104.0;
    b5[3]=-845.0/4104.0;

    b6[0]=-6080.0/20520.0;
    b6[1]=41040.0/20520.0;
    b6[2]=-28352.0/20520.0;
    b6[3]=9295.0/20520.0;
    b6[4]=-5643.0/20520.0;

    c1=902880.0/7618050.0;
    c3=3953664.0/7618050.0;
    c4=3855735.0/7618050.0;
    c5=-1371249.0/7618050.0;
    c6=277020.0/7618050.0;

    ec[0]=0.0;
    ec[1]=1.0/360.0;
    ec[2]=0.0;
    ec[3]=-128.0/4275.0;
    ec[4]=-2197.0/75240.0;
    ec[5]=1.0/50.0;
    ec[6]=2.0/55.0;

    ndim=0;
  }
      
  virtual ~ode_rkf45_gsl() {
  }

  /** \brief Perform an integration step

      Given initial value of the n-dimensional function in \c y and
      the derivative in \c dydx (which must be computed beforehand) at
      the point \c x, take a step of size \c h giving the result in \c
      yout, the uncertainty in \c yerr, and the new derivative in \c
      dydx_out using function \c derivs to calculate derivatives. The
      parameters \c yout and \c y and the parameters \c dydx_out and
      \c dydx may refer to the same object.

      If \c derivs always returns zero, then this function will
      also return zero. If not, <tt>step()</tt> will return the first
      non-zero value which was obtained in a call to \c derivs .
      The error handler is never called.
  */
  virtual int step(double x, double h, size_t n, vec_y_t &y, vec_dydx_t &dydx, 
		   vec_y_t &yout, vec_yerr_t &yerr, vec_dydx_t &dydx_out, 
		   func_t &derivs) {
    
    int ret=0;
    size_t i;
      
    if (ndim!=n) {
      k2.resize(n);
      k3.resize(n);
      k4.resize(n);
      k5.resize(n);
      k6.resize(n);
      ytmp.resize(n);
      
      ndim=n;
    }

    // k1 step
    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+ah[0]*h*dydx[i];
    }

    // k2 step
    o2scl::error_update(ret,derivs(x+ah[0]*h,n,ytmp,k2));

    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b3[0]*dydx[i]+b3[1]*k2[i]);
    }
      
    // k3 step
    o2scl::error_update(ret,derivs(x+ah[1]*h,n,ytmp,k3));
      
    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b4[0]*dydx[i]+b4[1]*k2[i]+b4[2]*k3[i]);
    }

    // k4 step
    o2scl::error_update(ret,derivs(x+ah[2]*h,n,ytmp,k4));

    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b5[0]*dydx[i]+b5[1]*k2[i]+b5[2]*k3[i]+
		      b5[3]*k4[i]);
    }
	
    // k5 step
    o2scl::error_update(ret,derivs(x+ah[3]*h,n,ytmp,k5));
      
    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b6[0]*dydx[i]+b6[1]*k2[i]+b6[2]*k3[i]+
		      b6[3]*k4[i]+b6[4]*k5[i]);
    }
      
    // k6 step and final sum
    o2scl::error_update(ret,derivs(x+ah[4]*h,n,ytmp,k6));
      
    for (i=0;i<n;i++) {
      yout[i]=y[i]+h*(c1*dydx[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]);
    }
      
    // We put this before the last function evaluation, in contrast
    // to the GSL version, so that the dydx[i] that appears in the
    // for loop below isn't modified by the subsequent derivative
    // evaluation using dydx_out. (The user could have given the
    // same vector for both)
    for (i=0;i<n;i++) {
      yerr[i]=h*(ec[1]*dydx[i]+ec[3]*k3[i]+ec[4]*k4[i]+ec[5]*k5[i]+
		 ec[6]*k6[i]);
    }

    o2scl::error_update(ret,derivs(x+h,n,yout,dydx_out));
      
    return ret;
  }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
