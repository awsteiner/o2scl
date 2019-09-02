/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/** \file inte_gauss_cern.h
    \brief File defining \ref o2scl::inte_gauss_cern
*/
#ifndef O2SCL_CERN_GAUSS_H
#define O2SCL_CERN_GAUSS_H

#include <o2scl/inte.h>
 
#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  class inte_gauss_coeffs_double {
    
  public:
    
    /** \brief Integration abscissas for \ref o2scl::inte_gauss_cern and 
	\ref o2scl::inte_cauchy_cern in double precision
    */
    double x[12];
    
    /** \brief Integration weights for \ref o2scl::inte_gauss_cern and 
	\ref o2scl::inte_cauchy_cern in double precision
    */
    double w[12];
    
    inte_gauss_coeffs_double() {
      x[0]=0.96028985649753623;
      x[1]=0.79666647741362674;
      x[2]=0.52553240991632899;
      x[3]=0.18343464249564980;
      x[4]=0.98940093499164993;
      x[5]=0.94457502307323258;
      x[6]=0.86563120238783175;
      x[7]=0.75540440835500303;
      x[8]=0.61787624440264375;
      x[9]=0.45801677765722739;
      x[10]=0.28160355077925891;
      x[11]=0.95012509837637440e-1;
      
      w[0]=0.10122853629037626;
      w[1]=0.22238103445337447;
      w[2]=0.31370664587788729;
      w[3]=0.36268378337836198;
      w[4]=0.27152459411754095e-1;
      w[5]=0.62253523938647893e-1;
      w[6]=0.95158511682492785e-1;
      w[7]=0.12462897125553387;
      w[8]=0.14959598881657673;
      w[9]=0.16915651939500254;
      w[10]=0.18260341504492359;
      w[11]=0.18945061045506850;
    }
  };
  
  class inte_gauss_coeffs_long_double {
    
  public:
    
    /** \brief Integration abscissas for \ref o2scl::inte_gauss_cern and 
	\ref o2scl::inte_cauchy_cern in double precision
    */
    long double x[12];
    
    /** \brief Integration weights for \ref o2scl::inte_gauss_cern and 
	\ref o2scl::inte_cauchy_cern in double precision
    */
    long double w[12];

    inte_gauss_coeffs_long_double() {
      
      x[0]=0.96028985649753623168356086856947299L;
      x[1]=0.79666647741362673959155393647583044L;
      x[2]=0.52553240991632898581773904918924635L;
      x[3]=0.18343464249564980493947614236018398L;
      x[4]=0.98940093499164993259615417345033263L;
      x[5]=0.94457502307323257607798841553460835L;
      x[6]=0.86563120238783174388046789771239313L;
      x[7]=0.75540440835500303389510119484744227L;
      x[8]=0.61787624440264374844667176404879102L;
      x[9]=0.45801677765722738634241944298357757L;
      x[10]=0.28160355077925891323046050146049611L;
      x[11]=0.095012509837637440185319335424958063L;
      
      w[0]=0.10122853629037625915253135430996219L;
      w[1]=0.22238103445337447054435599442624088L;
      w[2]=0.31370664587788728733796220198660131L;
      w[3]=0.36268378337836198296515044927719561L;
      w[4]=0.027152459411754094851780572456018104L;
      w[5]=0.062253523938647892862843836994377694L;
      w[6]=0.095158511682492784809925107602246226L;
      w[7]=0.12462897125553387205247628219201642L;
      w[8]=0.14959598881657673208150173054747855L;
      w[9]=0.16915651939500253818931207903035996L;
      w[10]=0.18260341504492358886676366796921994L;
      w[11]=0.18945061045506849628539672320828311L;
    }
    
  };

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
      
      \note Currently \o2 supports only types \c double and
      \c long \c double for the floating point type \c fp_t .

      \future Allow user to change \c cst?
  */
  template<class func_t=funct, class fp_t=double,
    class weights_t=inte_gauss_coeffs_double>
    class inte_gauss_cern : public inte<func_t,fp_t> {

    public:
  
    const fp_t *w, *x;
    weights_t wgts;
    
    inte_gauss_cern() {
      w=&(wgts.w[0]);
      x=&(wgts.x[0]);
    }

    virtual ~inte_gauss_cern() {
    }

    /** \brief Integrate function \c func from \c a to \c b.
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {

      fp_t y1, y2;
      err=0.0;

      size_t itx=0;

      int i;
      bool loop=true, loop2=false;
      static const fp_t cst=0.005;
      fp_t h=0.0;
      if (b==a) {
	res=0.0;
	return o2scl::success;
      }
      fp_t cnst=cst/(b-a);
      fp_t aa=0.0, bb=a;
      while (loop==true || loop2==true) {
	itx++;
	if (loop==true) {
	  aa=bb;
	  bb=b;
	}
	fp_t c1=(bb+aa)/2.0;
	fp_t c2=(bb-aa)/2.0;
	fp_t s8=0.0;
	for(i=0;i<4;i++) {
	  fp_t u=c2*x[i];
	  y1=func(c1+u);
	  y2=func(c1-u);
	  s8+=w[i]*(y1+y2);
	}
	fp_t s16=0.0;
	for(i=4;i<12;i++) {
	  fp_t u=c2*x[i];
	  y1=func(c1+u);
	  y2=func(c1-u);
	  s16+=w[i]*(y1+y2);
	}
	s16*=c2;
 
	loop=false;
	loop2=false;
	if (std::abs(s16-c2*s8)<this->tol_rel*(1.0+std::abs(s16))) {
	  h+=s16;
	  if (bb!=b) loop=true;
	} else {
	  bb=c1;
	  if (1.0+cnst*std::abs(c2)!=1.0) {
	    loop2=true;
	  } else {
	    this->last_iter=itx;
	    O2SCL_CONV2_RET("Failed to reach required accuracy in cern_",
			    "gauss::integ().",exc_efailed,this->err_nonconv);
	  }
	}
      }
      this->last_iter=itx;
      res=h;
      return o2scl::success;
    }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
