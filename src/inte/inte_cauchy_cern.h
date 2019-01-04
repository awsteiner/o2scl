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
/** \file inte_cauchy_cern.h
    \brief File defining \ref o2scl::inte_cauchy_cern
*/
#ifndef O2SCL_CERN_CAUCHY_H
#define O2SCL_CERN_CAUCHY_H

#include <o2scl/inte.h>
#include <o2scl/inte_gauss_cern.h>
 
#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
 
  /** \brief Cauchy principal value integration (CERNLIB)
 
      The location of the singularity must be specified before-hand in
      inte_cauchy_cern::s, and the singularity must not be at one of the
      endpoints. Note that when integrating a function of the form \f$
      \frac{f(x)}{(x-s)} \f$, the denominator \f$ (x-s) \f$ must be
      specified in the argument \c func to integ(). This is different
      from how the \ref inte_qawc_gsl operates.
      
      The method from \ref Longman58 is used for the decomposition of
      the integral, and the resulting integrals are computed using
      a user-specified base integration object.
      
      The uncertainty in the integral is not calculated, and is always
      given as zero. The default base integration object is of type
      \ref inte_gauss_cern.  This is the CERNLIB default, but can be
      modified by calling set_inte(). If the singularity is outside
      the region of integration, then the result from the base
      integration object is returned without calling the error
      handler.
    
      Possible errors for integ() and integ_err():
      - exc_einval - Singularity is on an endpoint
      - exc_efailed - Could not reach requested accuracy 
      (occurs only if inte::err_nonconv is true)
      
      \note Currently \o2 supports only types \c double and
      \c long \c double for the floating point type \c fp_t .

      This function is based on the CERNLIB routines RCAUCH and
      DCAUCH which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d104/top.html
  */
  template<class func_t, class fp_t=double,
    const fp_t x[]=inte_gauss_cern_x_double,
    const fp_t w[]=inte_gauss_cern_w_double>
    class inte_cauchy_cern : public inte<func_t,fp_t> {

    public:
  
    inte_cauchy_cern() {
      it=&def_inte;
    }

    /// The singularity (must be set before calling integ() or integ_err())
    fp_t s;

    /** \brief Set the base integration object to use (default is \ref
	inte_cauchy_cern::def_inte of type \ref inte_gauss_cern)
    */
    int set_inte(inte<func_t,fp_t> &i) {
      it=&i;
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b
     */
    virtual int integ_err(func_t &func, fp_t a, fp_t b,
			  fp_t &res, fp_t &err) {
      fp_t y1, y2, y3, y4;
      size_t itx=0;
      
      err=0.0;
      
      if (s==a || s==b) {
	O2SCL_CONV2_RET("Singularity on boundary in ",
			"inte_cauchy_cern::integ_err().",
			exc_einval,this->err_nonconv);
      } else if ((s<a && s<b) || (s>a && s>b)) {
	return it->integ_err(func,a,b,res,err);
      }
      fp_t h, b0;
      if (2.0*s<a+b) {
	h=it->integ(func,2.0*s-a,b);
	b0=s-a;
      } else {
	h=it->integ(func,a,2.0*s-b);
	b0=b-s;
      }
      fp_t c=0.005/b0;
      fp_t bb=0.0;
      bool loop1=true;
      while(loop1==true) {
	itx++;
	fp_t s8, s16;
	fp_t aa=bb;
	bb=b0;
	bool loop2=true;
	while(loop2==true) {
	  fp_t c1=(bb+aa)/2.0;
	  fp_t c2=(bb-aa)/2.0;
	  fp_t c3=s+c1;
	  fp_t c4=s-c1;
	  s8=0.0;
	  s16=0.0;
	  for(int i=0;i<4;i++) {
	    fp_t u=c2*x[i];
	    y1=func(c3+u);
	    y2=func(c4-u);
	    y3=func(c3-u);
	    y4=func(c4+u);
	    s8+=w[i]*((y1+y2)+(y3+y4));
	  }
	  s8*=c2;
	  for(int i=4;i<12;i++) {
	    fp_t u=c2*x[i];
	    y1=func(c3+u);
	    y2=func(c4-u);
	    y3=func(c3-u);
	    y4=func(c4+u);
	    s16+=w[i]*((y1+y2)+(y3+y4));
	  }
	  s16*=c2;

	  if (std::abs(s16-s8)<=this->tol_rel*(1.0+std::abs(s16))) {
	    loop2=false;
	  } else {
	    bb=c1;
	    if ((1.0+std::abs(c*c2))==1.0) {
	      loop2=false;
	    } else {
	      this->last_iter=itx;
	      O2SCL_CONV2_RET("Could not reach required accuracy in cern_",
			      "cauchy::integ()",exc_efailed,this->err_nonconv);
	    }
	  }
	}
	h+=s16;
	if (bb==b0) loop1=false;
      }
      this->last_iter=itx;
      res=h;
      return o2scl::success;
    }

    /// Default integration object
    inte_gauss_cern<func_t,fp_t> def_inte;

    protected:

#ifndef DOXYGEN_INTERNAL

    /// The base integration object
    inte<func_t,fp_t> *it;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
