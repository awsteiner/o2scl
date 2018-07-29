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
/** \file inte_ld_gauss_cern.h
    \brief File defining \ref o2scl::inte_ld_gauss_cern
*/
#ifndef O2SCL_INTE_LD_GAUSS_CERN_H
#define O2SCL_INTE_LD_GAUSS_CERN_H

#include <o2scl/inte_ld.h>
 
#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Gaussian quadrature (CERNLIB)
 
      For any interval \f$ (a,b) \f$, we define \f$ g_8(a,b) \f$ and
      \f$ g_{16}(a,b) \f$ to be the 8- and 16-point Gaussian
      quadrature approximations to
      \f[
      I=\int_a^b f(x)~dx
      \f]
      and define
      \f[
      r(a,b)=\frac{|g_{16}(a,b)-g_{8}(a,b)|}{1+g_{16}(a,b)}
      \f]
      The function integ() returns \f$ G \f$ given by
      \f[
      G=\sum_{i=1}^{k} g_{16}(x_{i-1},x_i)
      \f]
      where \f$ x_0=a \f$ and \f$ x_k=b \f$ and the subdivision
      points \f$ x_i \f$ are given by
      \f[
      x_i=x_{i-1}+\lambda(B-x_{i-1})
      \f]
      where \f$ \lambda \f$ is the first number in the sequence \f$
      1,\frac{1}{2},\frac{1}{4},... \f$ for which
      \f[
      r(x_{i-1},x_i)<\mathrm{eps}.
      \f]
      If, at any stage, the ratio
      \f[
      q=\left| \frac{x_i-x_{i-1}}{b-a} \right|
      \f]
      is so small so that \f$ 1+0.005 q \f$ is indistinguishable from
      unity, then the accuracy is required is not reachable and 
      the error handler is called.

      Unless there is severe cancellation, inte::tol_rel may be
      considered as specifying a bound on the relative error of the
      integral in the case that \f$ |I|>1 \f$ and an absolute error if
      \f$ |I|<1 \f$. More precisely, if \f$ k \f$ is the number of
      subintervals from above, and if
      \f[
      I_{abs} = \int_a^b |f(x)|~dx
      \f]
      then
      \f[
      \frac{|G-I|}{I_{abs}+k}<\mathrm{tol_rel}
      \f]
      will nearly always be true when no error is returned.  For
      functions with no singualarities in the interval, the accuracy
      will usually be higher than this.

      If the desired accuracy is not achieved, the integration
      functions will call the error handler and return the best guess,
      unless \ref inte::err_nonconv is false, in which case the error
      handler is not called. 
    
      This function is based on the CERNLIB routines GAUSS and
      DGAUSS which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d103/top.html

      \future Allow user to change \c cst?
  */
  template<class func_t=funct_ld> class inte_ld_gauss_cern :
    public inte_ld<func_t> {

    public:
  
    inte_ld_gauss_cern() {

      x[0]=0.96028985649753623168356086856947299L;
      x[1]=0.79666647741362673959155393647583044L;
      x[2]=0.52553240991632898581773904918924635L;
      x[3]=0.18343464249564980493947614236018398L;

      w[0]=0.10122853629037625915253135430996219L;
      w[1]=0.22238103445337447054435599442624088L;
      w[2]=0.31370664587788728733796220198660131L;
      w[3]=0.36268378337836198296515044927719561L;

      x[4]=0.98940093499164993259615417345033263L;
      x[5]=0.94457502307323257607798841553460835L;
      x[6]=0.86563120238783174388046789771239313L;
      x[7]=0.75540440835500303389510119484744227L;
      x[8]=0.61787624440264374844667176404879102L;
      x[9]=0.45801677765722738634241944298357757L;
      x[10]=0.28160355077925891323046050146049611L;
      x[11]=0.095012509837637440185319335424958063L;
      
      w[4]=0.027152459411754094851780572456018104L;
      w[5]=0.062253523938647892862843836994377694L;
      w[6]=0.095158511682492784809925107602246226L;
      w[7]=0.12462897125553387205247628219201642L;
      w[8]=0.14959598881657673208150173054747855L;
      w[9]=0.16915651939500253818931207903035996L;
      w[10]=0.18260341504492358886676366796921994L;
      w[11]=0.18945061045506849628539672320828311L;
    }

    /** \brief Integrate function \c func from \c a to \c b.
    */
    virtual int integ_err(func_t &func, long double a, long double b, 
			  long double &res, long double &err) {

      long double y1, y2;
      err=0.0;

      size_t itx=0;

      int i;
      bool loop=true, loop2=false;
      static const long double cst=0.005L;
      long double h=0.0L;
      if (b==a) {
	res=0.0L;
	return o2scl::success;
      }
      long double cnst=cst/(b-a);
      long double aa=0.0L, bb=a;
      while (loop==true || loop2==true) {
	itx++;
	if (loop==true) {
	  aa=bb;
	  bb=b;
	}
	long double c1=(bb+aa)/2.0L;
	long double c2=(bb-aa)/2.0L;
	long double s8=0.0L;
	for(i=0;i<4;i++) {
	  long double u=c2*x[i];
	  y1=func(c1+u);
	  y2=func(c1-u);
	  s8+=w[i]*(y1+y2);
	}
	long double s16=0.0L;
	for(i=4;i<12;i++) {
	  long double u=c2*x[i];
	  y1=func(c1+u);
	  y2=func(c1-u);
	  s16+=w[i]*(y1+y2);
	}
	s16*=c2;
 
	loop=false;
	loop2=false;
	if (fabs(s16-c2*s8)<this->tol_rel*(1.0L+fabs(s16))) {
	  h+=s16;
	  if (bb!=b) loop=true;
	} else {
	  bb=c1;
	  if (1.0L+cnst*fabs(c2)!=1.0L) {
	    loop2=true;
	  } else {
	    this->last_iter=itx;
	    O2SCL_CONV2_RET("Failed to reach required accuracy in cern_",
			    "gauss::integ_ld().",
			    exc_efailed,this->err_nonconv);
	  }
	}
      }
      this->last_iter=itx;
      res=h;
      return o2scl::success;
    }

    protected:

#ifndef DOXYGEN_INTERNAL

    /** \name Integration constants (Set in the constructor)
    */
    //@{
    long double x[12], w[12];
    //@}

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
