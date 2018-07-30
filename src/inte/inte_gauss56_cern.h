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

  /** \brief Fifth order integration abscissas for 
      \ref o2scl::inte_gauss56_cern in double precision
   */
  static const double inte_gauss56_cern_x5_double[5]=
    {4.6910077030668004e-02,2.3076534494715846e-01,
     5.0000000000000000e-01,7.6923465505284154e-01,
     9.5308992296933200e-01};

  /** \brief Fifth order integration weights for 
      \ref o2scl::inte_gauss56_cern in double precision
   */
  static const double inte_gauss56_cern_w5_double[5]=
    {1.1846344252809454e-01,2.3931433524968324e-01,
     2.8444444444444444e-01,2.3931433524968324e-01,
     1.1846344252809454e-01};
     
  /** \brief Sixth order integration abscissas for 
      \ref o2scl::inte_gauss56_cern in double precision
   */
  static const double inte_gauss56_cern_x6_double[6]=
    {3.3765242898423989e-02,1.6939530676686775e-01,
     3.8069040695840155e-01,6.1930959304159845e-01,
     8.3060469323313225e-01,9.6623475710157601e-01};
      
  /** \brief Sixth order integration weights for 
      \ref o2scl::inte_gauss56_cern in double precision
   */
  static const double inte_gauss56_cern_w6_double[6]=
    {8.5662246189585178e-02,1.8038078652406930e-01,
     2.3395696728634552e-01,2.3395696728634552e-01,
     1.8038078652406930e-01,8.5662246189585178e-02};
    
  /** \brief Fifth order integration abscissas for 
      \ref o2scl::inte_gauss56_cern in long double precision
   */
  static const long double inte_gauss56_cern_x5_long_double[5]=
    {0.04691007703066800360118656085030352L,
     0.23076534494715845448184278964989560L,
     0.5L,
     0.76923465505284154551815721035010440L,
     0.95308992296933199639881343914969648L};
  
  /** \brief Fifth order integration weights for 
      \ref o2scl::inte_gauss56_cern in long double precision
   */
  static const long double inte_gauss56_cern_w5_long_double[5]=
    {0.11846344252809454375713202035995868L,
     0.23931433524968323402064575741781910L,
     0.28444444444444444444444444444444444L,
     0.23931433524968323402064575741781910L,
     0.11846344252809454375713202035995868L};
  
  /** \brief Sixth order integration abscissas for 
      \ref o2scl::inte_gauss56_cern in long double precision
   */
  static const long double inte_gauss56_cern_x6_long_double[6]=
    {0.03376524289842398609384922275300270L,
     0.16939530676686774316930020249004733L,
     0.38069040695840154568474913915964403L,
     0.61930959304159845431525086084035597L,
     0.83060469323313225683069979750995267L,
     0.96623475710157601390615077724699730L};
  
  /** \brief Sixth order integration weights for 
      \ref o2scl::inte_gauss56_cern in long double precision
   */
  static const long double inte_gauss56_cern_w6_long_double[6]=
    {0.08566224618958517252014807108636645L,
     0.18038078652406930378491675691885806L,
     0.23395696728634552369493517199477550L,
     0.23395696728634552369493517199477550L,
     0.18038078652406930378491675691885806L,
     0.08566224618958517252014807108636645L};
  
  /** \brief 5,6-point Gaussian quadrature (CERNLIB)
      
      If \f$ I_5 \f$ is the 5-point approximation, and \f$ I_6 \f$ is the
      6-point approximation to the integral, then integ_err() returns
      the result \f$ \frac{1}{2}(I_5+I_6) \f$ with uncertainty
      \f$ |I_5-I_6| \f$.

      This class is based on the CERNLIB routines RGS56P and
      DGS56P which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d106/top.html

      \note Currently \o2 supports only types \c double and
      \c long \c double for the floating point type \c fp_t .
  */
  template<class func_t=funct, class fp_t=double,
    const fp_t x5[]=inte_gauss56_cern_x5_double,
    const fp_t w5[]=inte_gauss56_cern_w5_double,
    const fp_t x6[]=inte_gauss56_cern_x6_double,
    const fp_t w6[]=inte_gauss56_cern_w6_double>
    class inte_gauss56_cern : public inte<func_t,fp_t> {

  public:
  
    inte_gauss56_cern() {
    }

    /** \brief Integrate function \c func from \c a to \c b
	giving result \c res and error \c err

	This function always returns \ref success.
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b,
			  fp_t &res, fp_t &err) {

      fp_t rang=b-a, e5=0.0, e6=0.0, ytmp;

      for(int i=0;i<5;i++) {
	ytmp=func(a+rang*x5[i]);
	e5+=w5[i]*ytmp;
	ytmp=func(a+rang*x6[i]);
	e6+=w6[i]*ytmp;
      }
      ytmp=func(a+rang*x6[5]);
      e6+=w6[5]*ytmp;
      res=(e6+e5)*rang/2.0;
      err=std::abs((e6-e5)*rang);

      return success;
    }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
