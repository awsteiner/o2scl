/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

namespace o2scl {

  /** \brief Runge-Kutta-Fehlberg embedded Runge-Kutta ODE stepper (GSL)

      \verbatim embed:rst
      Based on [Hairer09]_.
      \endverbatim

      \todo Check this because it may not give exact dydt_out.
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
           class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
           class func_t=ode_funct, class fp_t=double>
  class ode_rkf45_gsl :
    public ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t,fp_t> {
    
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
  fp_t ah[5], b3[2], b4[3], b5[4], b6[5];
  fp_t c1, c3, c4, c5, c6;
  fp_t ec[7];
  //@}
  
  public:
  
  ode_rkf45_gsl() {
    this->order=5;

    fp_t num=1, den=4;
    ah[0]=num/den;
    num=3; den=8;
    ah[1]=num/den;
    num=12; den=13;
    ah[2]=num/den;
    ah[3]=1;
    num=1; den=2;
    ah[4]=num/den;

    num=3; den=32;
    b3[0]=num/den;
    num=9; den=32;
    b3[1]=num/den;

    num=1932; den=2197;
    b4[0]=num/den;
    num=-7200; den=2197;
    b4[1]=num/den;
    num=7296; den=2197;
    b4[2]=num/den;

    num=8341; den=4104;
    b5[0]=num/den;
    num=-32832; den=4104;
    b5[1]=num/den;
    num=29440; den=4104;
    b5[2]=num/den;
    num=-845; den=4104;
    b5[3]=num/den;

    num=-6080; den=20520;
    b6[0]=num/den;
    num=41040; den=20520;
    b6[1]=num/den;
    num=-28352; den=20520;
    b6[2]=num/den;
    num=9295; den=20520;
    b6[3]=num/den;
    num=-5643; den=20520;
    b6[4]=num/den;

    num=902880; den=7618050;
    c1=num/den;
    num=3953664; den=7618050;
    c3=num/den;
    num=3855735; den=7618050;
    c4=num/den;
    num=-1371249; den=7618050;
    c5=num/den;
    num=277020; den=7618050;
    c6=num/den;

    ec[0]=0;
    num=1; den=360;
    ec[1]=num/den;
    num=0.0;
    ec[2]=0;
    num=-128; den=4275;
    ec[3]=num/den;
    num=-2197; den=75240;
    ec[4]=num/den;
    num=1; den=50;
    ec[5]=num/den;
    num=2; den=55;
    ec[6]=num/den;
    
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
  virtual int step(fp_t x, fp_t h, size_t n, vec_y_t &y, vec_dydx_t &dydx, 
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
  
}

#endif
