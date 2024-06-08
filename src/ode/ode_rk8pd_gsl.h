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

namespace o2scl {

  /** \brief Embedded Runge-Kutta Prince-Dormand ODE stepper (GSL)

      \verbatim embed:rst
      Based on [Prince81]_.
      \endverbatim

      \verbatim embed:rst
      There is an example for the usage of this class in
      ``examples/ex_ode.cpp<`` documented in the
      :ref:`Ordinary differential equations example`.
      \endverbatim
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
           class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
           class func_t=ode_funct, class fp_t=double>
  class ode_rk8pd_gsl :
    public ode_step<vec_y_t,vec_dydx_t,vec_yerr_t,func_t,fp_t> {
    
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

      fp_t num=14005451, den=335480064;
      Abar[0]=num/den;
      Abar[1]=0;
      Abar[2]=0;
      Abar[3]=0;
      Abar[4]=0;
      num=-59238493; den=1068277825;
      Abar[5]=num/den;
      num=181606767; den=758867731;
      Abar[6]=num/den;
      num=561292985; den=797845732;
      Abar[7]=num/den;
      num=-1041891430; den=1371343529;
      Abar[8]=num/den;
      num=760417239; den=1151165299;
      Abar[9]=num/den;
      num=118820643; den=751138087;
      Abar[10]=num/den;
      num=-528747749; den=2220607170;
      Abar[11]=num/den;
      num=1; den=4;
      Abar[12]=num/den;

      num=13451932; den=455176623;
      A[0]=num/den;
      A[1]=0;
      A[2]=0;
      A[3]=0;
      A[4]=0;
      num=-808719846; den=976000145;
      A[5]=num/den;
      num=1757004468; den=5645159321;
      A[6]=num/den;
      num=656045339; den=265891186;
      A[7]=num/den;
      num=-3867574721; den=1518517206;
      A[8]=num/den;
      num=465885868; den=322736535;
      A[9]=num/den;
      num=53011238; den=667516719;
      A[10]=num/den;
      num=2; den=45;
      A[11]=num/den;

      num=1; den=18;
      ah[0]=num/den;
      num=1; den=12;
      ah[1]=num/den;
      num=1; den=8;
      ah[2]=num/den;
      num=5; den=16;
      ah[3]=num/den;
      num=3; den=8;
      ah[4]=num/den;
      num=59; den=400;
      ah[5]=num/den;
      num=93; den=200;
      ah[6]=num/den;
      num=5490023248; den=9719169821;
      ah[7]=num/den;
      num=13; den=20;
      ah[8]=num/den;
      num=1201146811; den=1299019798;
      ah[9]=num/den;

      num=1; den=18;
      b21=num/den;

      num=1; den=48;
      b3[0]=num/den;
      num=1; den=16;
      b3[1]=num/den;

      num=1; den=32;
      b4[0]=num/den;
      b4[1]=0;
      num=3; den=32;
      b4[2]=num/den;

      num=5; den=16;
      b5[0]=num/den;
      b5[1]=0;
      num=-75; den=64;
      b5[2]=num/den;
      num=75; den=64;
      b5[3]=num/den;

      num=3; den=80;
      b6[0]=num/den;
      b6[1]=0;
      b6[2]=0;
      num=3; den=16;
      b6[3]=num/den;
      num=3; den=20;
      b6[4]=num/den;

      num=29443841; den=614563906;
      b7[0]=num/den;
      b7[1]=0;
      b7[2]=0;
      num=77736538; den=692538347;
      b7[3]=num/den;
      num=-28693883; den=1125000000;
      b7[4]=num/den;
      num=23124283; den=1800000000;
      b7[5]=num/den;

      num=16016141; den=946692911;
      b8[0]=num/den;
      b8[1]=0;
      b8[2]=0;
      num=61564180; den=158732637;
      b8[3]=num/den;
      num=22789713; den=633445777;
      b8[4]=num/den;
      num=545815736; den=2771057229;
      b8[5]=num/den;
      num=-180193667; den=1043307555;
      b8[6]=num/den;

      num=39632708; den=573591083;
      b9[0]=num/den;
      b9[1]=0;
      b9[2]=0;
      num=-433636366; den=683701615;
      b9[3]=num/den;
      num=-421739975; den=2616292301;
      b9[4]=num/den;
      num=100302831; den=723423059;
      b9[5]=num/den;
      num=790204164; den=839813087;
      b9[6]=num/den;
      b9[7]=800635310/3783071287.0;

      num=246121993; den=1340847787;
      b10[0]=num/den;
      b10[1]=0;
      b10[2]=0;
      num=-37695042795; den=15268766246;
      b10[3]=num/den;
      num=-309121744; den=1061227803;
      b10[4]=num/den;
      num=-12992083; den=490766935;
      b10[5]=num/den;
      num=6005943493; den=2108947869;
      b10[6]=num/den;
      num=393006217; den=1396673457;
      b10[7]=num/den;
      num=123872331; den=1001029789;
      b10[8]=num/den;

      num=-1028468189; den=846180014;
      b11[0]=num/den;
      b11[1]=0;
      b11[2]=0;
      num=8478235783; den=508512852;
      b11[3]=num/den;
      num=1311729495; den=1432422823;
      b11[4]=num/den;
      num=-10304129995; den=1701304382;
      b11[5]=num/den;
      num=-48777925059; den=3047939560;
      b11[6]=num/den;
      num=15336726248; den=1032824649;
      b11[7]=num/den;
      num=-45442868181; den=3398467696;
      b11[8]=num/den;
      num=3065993473; den=597172653;
      b11[9]=num/den;

      num=185892177; den=718116043;
      b12[0]=num/den;
      b12[1]=0;
      b12[2]=0;
      num=-3185094517; den=667107341;
      b12[3]=num/den;
      num=-477755414; den=1098053517;
      b12[4]=num/den;
      num=-703635378; den=230739211;
      b12[5]=num/den;
      num=5731566787; den=1027545527;
      b12[6]=num/den;
      num=5232866602; den=850066563;
      b12[7]=num/den;
      num=-4093664535; den=808688257;
      b12[8]=num/den;
      num=3962137247; den=1805957418;
      b12[9]=num/den;
      num=65686358; den=487910083;
      b12[10]=num/den;

      num=403863854; den=491063109;
      b13[0]=num/den;
      b13[1]=0;
      b13[2]=0;
      num=-5068492393; den=434740067;
      b13[3]=num/den;
      num=-411421997; den=543043805;
      b13[4]=num/den;
      num=652783627; den=914296604;
      b13[5]=num/den;
      num=11173962825; den=925320556;
      b13[6]=num/den;
      num=-13158990841; den=6184727034;
      b13[7]=num/den;
      num=3936647629; den=1978049680;
      b13[8]=num/den;
      num=-160528059; den=685178525;
      b13[9]=num/den;
      num=248638103; den=1413531060;
      b13[10]=num/den;
      b13[11]=0;

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
  
}

#endif
