/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
      - exc_efailed - Couldn't reach requested accuracy 
      (occurs only if inte::err_nonconv is true)
      
      This function is based on the CERNLIB routines RCAUCH and
      DCAUCH which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d104/top.html
  */
  template<class func_t> class inte_cauchy_cern : 
    public inte<func_t> {

    public:
  
    inte_cauchy_cern() {
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

      it=&def_inte;
    }

    /// The singularity (must be set before calling integ() or integ_err())
    double s;

    /** \brief Set the base integration object to use (default is \ref
	inte_cauchy_cern::def_inte of type \ref inte_gauss_cern)
    */
    int set_inte(inte<func_t> &i) {
      it=&i;
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b
	giving result \c res and error \c err

	\todo Fix converge error issue here.
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {
      if (s==a || s==b) {
	O2SCL_ERR("Singularity on boundary in inte_cauchy_cern::integ().",
		      exc_einval);
      }
      res=integ(func,a,b);
      err=0.0;
      //if (this->err_nonconv) return this->last_conv;
      return success;
    }
    
    /** \brief Integrate function \c func from \c a to \c b
    */
    virtual double integ(func_t &func, double a, double b) {
      double y1, y2, y3, y4;
      size_t itx=0;
      
      if (s==a || s==b) {
	O2SCL_ERR("Singularity on boundary in inte_cauchy_cern::integ().",
		  exc_einval);
	return 0.0;
      } else if ((s<a && s<b) || (s>a && s>b)) {
	return it->integ(func,a,b );
      }
      double h, b0;
      if (2.0*s<a+b) {
	h=it->integ(func,2.0*s-a,b );
	b0=s-a;
      } else {
	h=it->integ(func,a,2.0*s-b );
	b0=b-s;
      }
      double c=0.005/b0;
      double bb=0.0;
      bool loop1=true;
      while(loop1==true) {
	itx++;
	double s8, s16;
	double aa=bb;
	bb=b0;
	bool loop2=true;
	while(loop2==true) {
	  double c1=(bb+aa)/2.0;
	  double c2=(bb-aa)/2.0;
	  double c3=s+c1;
	  double c4=s-c1;
	  s8=0.0;
	  s16=0.0;
	  for(int i=0;i<4;i++) {
	    double u=c2*x[i];
	    y1=func(c3+u);
	    y2=func(c4-u);
	    y3=func(c3-u);
	    y4=func(c4+u);
	    s8+=w[i]*((y1+y2)+(y3+y4));
	  }
	  s8*=c2;
	  for(int i=4;i<12;i++) {
	    double u=c2*x[i];
	    y1=func(c3+u);
	    y2=func(c4-u);
	    y3=func(c3-u);
	    y4=func(c4+u);
	    s16+=w[i]*((y1+y2)+(y3+y4));
	  }
	  s16*=c2;

	  if (fabs(s16-s8)<=this->tol_rel*(1.0+fabs(s16))) {
	    loop2=false;
	  } else {
	    bb=c1;
	    if ((1.0+fabs(c*c2))==1.0) {
	      loop2=false;
	    } else {
	      this->last_iter=itx;
	      O2SCL_CONV2("Couldn't reach required accuracy in cern_",
			  "cauchy::integ()",exc_efailed,this->err_nonconv);
	      return 0.0;
	    }
	  }
	}
	h+=s16;
	if (bb==b0) loop1=false;
      }
      this->last_iter=itx;
      return h;
    }

    /// Default integration object
    inte_gauss_cern<func_t> def_inte;

    protected:

#ifndef DOXYGEN_INTERNAL

    /** \name Integration constants
    */
    //@{
    double x[12], w[12];
    //@}

    /// The base integration object
    inte<func_t> *it;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
