/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
/* ode-initval/rk8pd.c
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
#ifndef O2SCL_GSL_RK8PD_H
#define O2SCL_GSL_RK8PD_H

/** \file ode_rk8pd_gsl.h
    \brief File defining \ref o2scl::ode_rk8pd_gsl 
*/

#include <o2scl/ode_step.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Embedded Runge-Kutta Prince-Dormand ODE stepper (GSL)

      Based on \ref Prince81 .

      There is an example for the usage of this class in
      <tt>examples/ex_ode.cpp</tt> documented in the \ref ex_ode_sect
      section.
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
    class func_t=ode_funct11>
    class ode_rk8pd_gsl : public ode_step<vec_y_t,
    vec_dydx_t,vec_yerr_t,func_t> {

  protected:

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  
    
    /// \name Storage for the intermediate steps
    //@{
    vec_dydx_t k2, k3, k4, k5, k6, k7;
    vec_dydx_t k8, k9, k10, k11, k12, k13;
    vec_y_t ytmp;
    //@}
    
    /// Size of allocated vectors
    size_t ndim;

    /** \name Storage for the coefficients
     */
    //@{
    double Abar[13], A[12], ah[10], b21, b3[2], b4[3], b5[4], b6[5];
    double b7[6], b8[7], b9[8], b10[9], b11[10], b12[11], b13[12];
    //@}
      
  public:

    ode_rk8pd_gsl() {
      this->order=8;

      Abar[0]=14005451.0/335480064.0;
      Abar[1]=0.0;
      Abar[2]=0.0;
      Abar[3]=0.0;
      Abar[4]=0.0;
      Abar[5]=-59238493.0/1068277825.0;
      Abar[6]=181606767.0/758867731.0;
      Abar[7]=561292985.0/797845732.0;
      Abar[8]=-1041891430.0/1371343529.0;
      Abar[9]=760417239.0/1151165299.0;
      Abar[10]=118820643.0/751138087.0;
      Abar[11]=-528747749.0/2220607170.0;
      Abar[12]=1.0/4.0;

      A[0]=13451932.0/455176623.0;
      A[1]=0.0;
      A[2]=0.0;
      A[3]=0.0;
      A[4]=0.0;
      A[5]=-808719846.0/976000145.0;
      A[6]=1757004468.0/5645159321.0;
      A[7]=656045339.0/265891186.0;
      A[8]=-3867574721.0/1518517206.0;
      A[9]=465885868.0/322736535.0;
      A[10]=53011238.0/667516719.0;
      A[11]=2.0/45.0;

      ah[0]=1.0/18.0;
      ah[1]=1.0/12.0;
      ah[2]=1.0/8.0;
      ah[3]=5.0/16.0;
      ah[4]=3.0/8.0;
      ah[5]=59.0/400.0;
      ah[6]=93.0/200.0;
      ah[7]=5490023248.0/9719169821.0;
      ah[8]=13.0/20.0;
      ah[9]=1201146811.0/1299019798.0;
      
      b21=1.0/18.0;

      b3[0]=1.0/48.0;
      b3[1]=1.0/16.0;

      b4[0]=1.0/32.0;
      b4[1]=0.0;
      b4[2]=3.0/32.0;

      b5[0]=5.0/16.0;
      b5[1]=0.0;
      b5[2]=-75.0/64.0;
      b5[3]=75.0/64.0;

      b6[0]=3.0/80.0;
      b6[1]=0.0;
      b6[2]=0.0;
      b6[3]=3.0/16.0;
      b6[4]=3.0/20.0;

      b7[0]=29443841.0/614563906.0;
      b7[1]=0.0;
      b7[2]=0.0;
      b7[3]=77736538.0/692538347.0;
      b7[4]=-28693883.0/1125000000.0;
      b7[5]=23124283.0/1800000000.0;

      b8[0]=16016141.0/946692911.0;
      b8[1]=0.0;
      b8[2]=0.0;
      b8[3]=61564180.0/158732637.0;
      b8[4]=22789713.0/633445777.0;
      b8[5]=545815736.0/2771057229.0;
      b8[6]=-180193667.0/1043307555.0;

      b9[0]=39632708.0/573591083.0;
      b9[1]=0.0;
      b9[2]=0.0;
      b9[3]=-433636366.0/683701615.0;
      b9[4]=-421739975.0/2616292301.0;
      b9[5]=100302831.0/723423059.0;
      b9[6]=790204164.0/839813087.0;
      b9[7]=800635310.0/3783071287.0;

      b10[0]=246121993.0/1340847787.0;
      b10[1]=0.0;
      b10[2]=0.0;
      b10[3]=-37695042795.0/15268766246.0;
      b10[4]=-309121744.0/1061227803.0;
      b10[5]=-12992083.0/490766935.0;
      b10[6]=6005943493.0/2108947869.0;
      b10[7]=393006217.0/1396673457.0;
      b10[8]=123872331.0/1001029789.0;

      b11[0]=-1028468189.0/846180014.0;
      b11[1]=0.0;
      b11[2]=0.0;
      b11[3]=8478235783.0/508512852.0;
      b11[4]=1311729495.0/1432422823.0;
      b11[5]=-10304129995.0/1701304382.0;
      b11[6]=-48777925059.0/3047939560.0;
      b11[7]=15336726248.0/1032824649.0;
      b11[8]=-45442868181.0/3398467696.0;
      b11[9]=3065993473.0/597172653.0;

      b12[0]=185892177.0/718116043.0;
      b12[1]=0.0;
      b12[2]=0.0;
      b12[3]=-3185094517.0/667107341.0;
      b12[4]=-477755414.0/1098053517.0;
      b12[5]=-703635378.0/230739211.0;
      b12[6]=5731566787.0/1027545527.0;
      b12[7]=5232866602.0/850066563.0;
      b12[8]=-4093664535.0/808688257.0;
      b12[9]=3962137247.0/1805957418.0;
      b12[10]=65686358.0/487910083.0;

      b13[0]=403863854.0/491063109.0;
      b13[1]=0.0;
      b13[2]=0.0;
      b13[3]=-5068492393.0/434740067.0;
      b13[4]=-411421997.0/543043805.0;
      b13[5]=652783627.0/914296604.0;
      b13[6]=11173962825.0/925320556.0;
      b13[7]=-13158990841.0/6184727034.0;
      b13[8]=3936647629.0/1978049680.0;
      b13[9]=-160528059.0/685178525.0;
      b13[10]=248638103.0/1413531060.0;
      b13[11]=0.0;

      ndim=0;
    }
      
    virtual ~ode_rk8pd_gsl() {
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
	k7.resize(n);
	k8.resize(n);
	k9.resize(n);
	k10.resize(n);
	k11.resize(n);
	k12.resize(n);
	k13.resize(n);
	ytmp.resize(n);
	
	ndim=n;
      }

      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+b21*h*dydx[i];
      }
	
      error_update(ret,derivs(x+ah[0]*h,n,ytmp,k2));

      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b3[0]*dydx[i]+b3[1]*k2[i]);
      }
      
      error_update(ret,derivs(x+ah[1]*h,n,ytmp,k3));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b4[0]*dydx[i]+b4[2]*k3[i]);
      }

      error_update(ret,derivs(x+ah[2]*h,n,ytmp,k4));

      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b5[0]*dydx[i]+b5[2]*k3[i]+b5[3]*k4[i]);
      }
      
      error_update(ret,derivs(x+ah[3]*h,n,ytmp,k5));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b6[0]*dydx[i]+b6[3]*k4[i]+b6[4]*k5[i]);
      }
      
      error_update(ret,derivs(x+ah[4]*h,n,ytmp,k6));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b7[0]*dydx[i]+b7[3]*k4[i]+b7[4]*k5[i]+b7[5]*k6[i]);
      }
      
      error_update(ret,derivs(x+ah[5]*h,n,ytmp,k7));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b8[0]*dydx[i]+b8[3]*k4[i]+b8[4]*k5[i]+b8[5]*k6[i]+
		   b8[6]*k7[i]);
      }

      error_update(ret,derivs(x+ah[6]*h,n,ytmp,k8));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b9[0]*dydx[i]+b9[3]*k4[i]+b9[4]*k5[i]+b9[5]*k6[i]+
		   b9[6]*k7[i]+b9[7]*k8[i]);
      }

      error_update(ret,derivs(x+ah[7]*h,n,ytmp,k9));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b10[0]*dydx[i]+b10[3]*k4[i]+b10[4]*k5[i]+
			b10[5]*k6[i]+b10[6]*k7[i]+b10[7]*k8[i]+
			b10[8]*k9[i]);
      }

      error_update(ret,derivs(x+ah[8]*h,n,ytmp,k10));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b11[0]*dydx[i]+b11[3]*k4[i]+b11[4]*k5[i]+
			b11[5]*k6[i]+b11[6]*k7[i]+b11[7]*k8[i]+
			b11[8]*k9[i]+b11[9]*k10[i]);
      }

      error_update(ret,derivs(x+ah[9]*h,n,ytmp,k11));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b12[0]*dydx[i]+b12[3]*k4[i]+b12[4]*k5[i]+
			b12[5]*k6[i]+b12[6]*k7[i]+b12[7]*k8[i]+
			b12[8]*k9[i]+b12[9]*k10[i]+b12[10]*k11[i]);
      }

      error_update(ret,derivs(x+h,n,ytmp,k12));
      
      for (i=0;i<n;i++) {
	ytmp[i]=y[i]+h*(b13[0]*dydx[i]+b13[3]*k4[i]+b13[4]*k5[i]+
			b13[5]*k6[i]+b13[6]*k7[i]+b13[7]*k8[i]+
			b13[8]*k9[i]+b13[9]*k10[i]+b13[10]*k11[i]+
			b13[11]*k12[i]);
      }

      error_update(ret,derivs(x+h,n,ytmp,k13));

      // final sum

      for (i=0;i<n;i++) {
	double ksum8=Abar[0]*dydx[i]+Abar[5]*k6[i]+Abar[6]*k7[i]+
	  Abar[7]*k8[i]+Abar[8]*k9[i]+Abar[9]*k10[i]+
	  Abar[10]*k11[i]+Abar[11]*k12[i]+Abar[12]*k13[i];

	yout[i]=y[i]+h*ksum8;
      }
      
      // We put this before the last function evaluation, in contrast
      // to the GSL version, so that the dydx[i] that appears in the
      // for loop below isn't modified by the subsequent derivative
      // evaluation using dydx_out. (The user could have given the
      // same vector for both)
      for (i=0;i<n;i++) {

	double ksum8=Abar[0]*dydx[i]+Abar[5]*k6[i]+Abar[6]*k7[i]+
	  Abar[7]*k8[i]+Abar[8]*k9[i]+Abar[9]*k10[i]+
	  Abar[10]*k11[i]+Abar[11]*k12[i]+Abar[12]*k13[i];
	double ksum7=A[0]*dydx[i]+A[5]*k6[i]+A[6]*k7[i]+A[7]*k8[i]+
	  A[8]*k9[i]+A[9]*k10[i]+A[10]*k11[i]+A[11]*k12[i];
	
	yerr[i]=h*(ksum7-ksum8);
      }
      
      error_update(ret,derivs(x+h,n,yout,dydx_out));

      return ret;
    }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
