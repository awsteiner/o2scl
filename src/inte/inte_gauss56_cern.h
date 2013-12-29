/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
/** \file inte_gauss56_cern.h
    \brief File defining \ref o2scl::inte_gauss56_cern
*/
#ifndef O2SCL_CERN_GAUSS56_H
#define O2SCL_CERN_GAUSS56_H

#include <o2scl/inte.h>
#include <o2scl/funct.h>
 
#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief 5,6-point Gaussian quadrature (CERNLIB)
      
      If \f$ I_5 \f$ is the 5-point approximation, and \f$ I_6 \f$ is the
      6-point approximation to the integral, then integ_err() returns
      the result \f$ \frac{1}{2}(I_5+I_6) \f$ with uncertainty
      \f$ |I_5-I_6| \f$.

      This class is based on the CERNLIB routines RGS56P and
      DGS56P which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d106/top.html

  */
  template<class func_t=funct> class inte_gauss56_cern : public inte<func_t>  {

  public:
  
    inte_gauss56_cern() {
      x5[0]=4.6910077030668004e-02;
      w5[0]=1.1846344252809454e-01;
      x5[1]=2.3076534494715846e-01;
      w5[1]=2.3931433524968324e-01;
      x5[2]=5.0000000000000000e-01;
      w5[2]=2.8444444444444444e-01;
      x5[3]=7.6923465505284154e-01;
      w5[3]=2.3931433524968324e-01;
      x5[4]=9.5308992296933200e-01;
      w5[4]=1.1846344252809454e-01;
 
      x6[0]=3.3765242898423989e-02;
      w6[0]=8.5662246189585178e-02;
      x6[1]=1.6939530676686775e-01;
      w6[1]=1.8038078652406930e-01;
      x6[2]=3.8069040695840155e-01;
      w6[2]=2.3395696728634552e-01;
      x6[3]=6.1930959304159845e-01;
      w6[3]=2.3395696728634552e-01;
      x6[4]=8.3060469323313225e-01;
      w6[4]=1.8038078652406930e-01;
      x6[5]=9.6623475710157601e-01;
      w6[5]=8.5662246189585178e-02;
    }

    /** \brief Integrate function \c func from \c a to \c b
	giving result \c res and error \c err

	This function always returns \ref success.
    */
    virtual int integ_err(func_t &func, double a, double b,
			  double &res, double &err) {

      double rang=b-a, e5=0.0, e6=0.0, ytmp;

      for(int i=0;i<5;i++) {
	ytmp=func(a+rang*x5[i]);
	e5+=w5[i]*ytmp;
	ytmp=func(a+rang*x6[i]);
	e6+=w6[i]*ytmp;
      }
      ytmp=func(a+rang*x6[5]);
      e6+=w6[5]*ytmp;
      res=(e6+e5)*rang/2.0;
      err=fabs((e6-e5)*rang);

      return success;
    }

  protected:

#ifndef DOXYGEN_INTERNAL

    /** \name Integration constants
    */
    //@{
    double x5[5], w5[5];
    double x6[6], w6[6];
    //@}

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
