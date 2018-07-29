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
/** \file inte_ld_gauss56_cern.h
    \brief File defining \ref o2scl::inte_ld_gauss56_cern
*/
#ifndef O2SCL_INTE_LD_GAUSS56_CERN_H
#define O2SCL_INTE_LD_GAUSS56_CERN_H

#include <o2scl/inte_ld.h>
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
  template<class func_t=funct_ld> class inte_ld_gauss56_cern :
    public inte_ld<func_t>  {

  public:
  
    inte_ld_gauss56_cern() {
      
      x5[0]=0.04691007703066800360118656085030352L;
      w5[0]=0.11846344252809454375713202035995868L;
      x5[1]=0.23076534494715845448184278964989560L;
      w5[1]=0.23931433524968323402064575741781910L;
      x5[2]=0.5L;
      w5[2]=0.28444444444444444444444444444444444L;
      x5[3]=0.76923465505284154551815721035010440L;
      w5[3]=0.23931433524968323402064575741781910L;
      x5[4]=0.95308992296933199639881343914969648L;
      w5[4]=0.11846344252809454375713202035995868L;
 
      x6[0]=0.03376524289842398609384922275300270L;
      w6[0]=0.08566224618958517252014807108636645L;
      x6[1]=0.16939530676686774316930020249004733L;
      w6[1]=0.18038078652406930378491675691885806L;
      x6[2]=0.38069040695840154568474913915964403L;
      w6[2]=0.23395696728634552369493517199477550L;
      x6[3]=0.61930959304159845431525086084035597L;
      w6[3]=0.23395696728634552369493517199477550L;
      x6[4]=0.83060469323313225683069979750995267L;
      w6[4]=0.18038078652406930378491675691885806L;
      x6[5]=0.96623475710157601390615077724699730L;
      w6[5]=0.08566224618958517252014807108636645L;
    }

    /** \brief Integrate function \c func from \c a to \c b
	giving result \c res and error \c err

	This function always returns \ref success.
    */
    virtual int integ_err(func_t &func, long double a, long double b,
			  long double &res, long double &err) {

      long double rang=b-a, e5=0.0L, e6=0.0L, ytmp;

      for(int i=0;i<5;i++) {
	ytmp=func(a+rang*x5[i]);
	e5+=w5[i]*ytmp;
	ytmp=func(a+rang*x6[i]);
	e6+=w6[i]*ytmp;
      }
      ytmp=func(a+rang*x6[5]);
      e6+=w6[5]*ytmp;
      res=(e6+e5)*rang/2.0L;
      err=fabs((e6-e5)*rang);

      return success;
    }

  protected:

#ifndef DOXYGEN_INTERNAL

    /** \name Integration constants
    */
    //@{
    long double x5[5], w5[5];
    long double x6[6], w6[6];
    //@}

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
