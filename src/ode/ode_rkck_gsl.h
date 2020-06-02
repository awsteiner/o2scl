/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
/* ode-initval/rkck.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
#ifndef O2SCL_GSL_RKCK_H
#define O2SCL_GSL_RKCK_H

/** \file ode_rkck_gsl.h
    \brief File defining \ref o2scl::ode_rkck_gsl 
*/

#include <o2scl/err_hnd.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_step.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Cash-Karp embedded Runge-Kutta ODE stepper (GSL)

      \verbatim embed:rst
      Based on [Cash90]_ .
      \endverbatim

      There is an example for the usage of this class in
      <tt>examples/ex_ode.cpp</tt> documented in the \ref ex_ode_sect
      section.
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
    class func_t=ode_funct> class ode_rkck_gsl :
    public ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t> {
    
  protected:
  
  /// \name Storage for the intermediate steps
  //@{
  vec_y_t ytmp;
  vec_dydx_t k2, k3, k4, k5, k6;
  //@}
      
  /// Size of allocated vectors
  size_t ndim;

  /** \name Storage for the coefficients
   */
  //@{
  double ah[5], b3[2], b4[3], b5[4], b6[5], ec[7];
  double b21, c1, c3, c4, c6;
  //@}
      
  public:

  ode_rkck_gsl() {
    this->order=5;

    ah[0]=1.0/5.0;
    ah[1]=3.0/10.0;
    ah[2]=3.0/5.0;
    ah[3]=1.0;
    ah[4]=7.0/8.0;

    b3[0]=3.0/40.0;
    b3[1]=9.0/40.0;

    b4[0]=3.0/10.0;
    b4[1]=-9.0/10.0;
    b4[2]=12.0/10.0;

    b5[0]=-11.0/54.0;
    b5[1]=5.0/2.0;
    b5[2]=-70.0/27.0;
    b5[3]=35.0/27.0;

    b6[0]=1631.0/55296.0;
    b6[1]=175.0/512.0;
    b6[2]=575.0/13824.0;
    b6[3]=44275.0/110592.0;
    b6[4]=253.0/4096.0;

    ec[0]=0.0;
    ec[1]=37.0/378.0-2825.0/27648.0;
    ec[2]=0.0;
    ec[3]=250.0/621.0-18575.0/48384.0;
    ec[4]=125.0/594.0-13525.0/55296.0;
    ec[5]=-277.0/14336.0;
    ec[6]=512.0/1771.0-1.0/4.0;

    b21=1.0/5.0;

    c1=37.0/378.0;
    c3=250.0/621.0;
    c4=125.0/594.0;
    c6=512.0/1771.0;
	
    ndim=0;
  }
      
  virtual ~ode_rkck_gsl() {
  }

  /** \brief Perform an integration step

      Given initial value of the n-dimensional function in \c y and
      the derivative in \c dydx (which must be computed beforehand)
      at the point \c x, take a step of size \c h giving the result
      in \c yout, the uncertainty in \c yerr, and the new derivative
      in \c dydx_out using function \c derivs to calculate
      derivatives. The parameters \c yout and \c y and the
      parameters \c dydx_out and \c dydx may refer to the same
      object.

      If \c derivs always returns zero, then this function will
      also return zero. If not, <tt>step()</tt> will return the first
      non-zero value which was obtained in a call to \c derivs .
      The error handler is never called.
  */
  virtual int step(double x, double h, size_t n, vec_y_t &y, 
		   vec_dydx_t &dydx, vec_y_t &yout, vec_yerr_t &yerr, 
		   vec_dydx_t &dydx_out, func_t &derivs) {
	
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

    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+b21*h*dydx[i];
    }
	
    o2scl::error_update(ret,derivs(x+ah[0]*h,n,ytmp,k2));

    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b3[0]*dydx[i]+b3[1]*k2[i]);
    }
      
    o2scl::error_update(ret,derivs(x+ah[1]*h,n,ytmp,k3));
      
    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b4[0]*dydx[i]+b4[1]*k2[i]+b4[2]*k3[i]);
    }

    o2scl::error_update(ret,derivs(x+ah[2]*h,n,ytmp,k4));

    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b5[0]*dydx[i]+b5[1]*k2[i]+b5[2]*k3[i]+
		      b5[3]*k4[i]);
    }
	
    o2scl::error_update(ret,derivs(x+ah[3]*h,n,ytmp,k5));
      
    for (i=0;i<n;i++) {
      ytmp[i]=y[i]+h*(b6[0]*dydx[i]+b6[1]*k2[i]+b6[2]*k3[i]+
		      b6[3]*k4[i]+b6[4]*k5[i]);
    }
      
    o2scl::error_update(ret,derivs(x+ah[4]*h,n,ytmp,k6));
      
    for (i=0;i<n;i++) {
      yout[i]=y[i]+h*(c1*dydx[i]+c3*k3[i]+c4*k4[i]+c6*k6[i]);
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
