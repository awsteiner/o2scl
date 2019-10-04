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
/** \file inte_gauss56_cern.h
    \brief File defining \ref o2scl::inte_gauss56_cern
*/
#ifndef O2SCL_CERN_GAUSS56_H
#define O2SCL_CERN_GAUSS56_H

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

#include <o2scl/inte.h>
#include <o2scl/funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  class inte_gauss56_coeffs_double {

  public:
    
    /** \brief Fifth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in double precision
    */
    double x5[5];
    
    /** \brief Fifth order integration weights for 
	\ref o2scl::inte_gauss56_cern in double precision
    */
    double w5[5];
    
    /** \brief Sixth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in double precision
    */
    double x6[6];
    
    /** \brief Sixth order integration weights for 
	\ref o2scl::inte_gauss56_cern in double precision
    */
    double w6[6];

    inte_gauss56_coeffs_double() {
      x5[0]=4.6910077030668004e-02;
      x5[1]=2.3076534494715846e-01;
      x5[2]=5.0000000000000000e-01;
      x5[3]=7.6923465505284154e-01;
      x5[4]=9.5308992296933200e-01;
    
      w5[0]=1.1846344252809454e-01;
      w5[1]=2.3931433524968324e-01;
      w5[2]=2.8444444444444444e-01;
      w5[3]=2.3931433524968324e-01;
      w5[4]=1.1846344252809454e-01;
    
      x6[0]=3.3765242898423989e-02;
      x6[1]=1.6939530676686775e-01;
      x6[2]=3.8069040695840155e-01;
      x6[3]=6.1930959304159845e-01;
      x6[4]=8.3060469323313225e-01;
      x6[5]=9.6623475710157601e-01;
    
      w6[0]=8.5662246189585178e-02;
      w6[1]=1.8038078652406930e-01;
      w6[2]=2.3395696728634552e-01;
      w6[3]=2.3395696728634552e-01;
      w6[4]=1.8038078652406930e-01;
      w6[5]=8.5662246189585178e-02;
    }      
  };
  
  class inte_gauss56_coeffs_long_double {

  public:

    /** \brief Fifth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in long double precision
    */
    long double x5[5];
  
    /** \brief Fifth order integration weights for 
	\ref o2scl::inte_gauss56_cern in long double precision
    */
    long double w5[5];
  
    /** \brief Sixth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in long double precision
    */
    long double x6[6];
  
    /** \brief Sixth order integration weights for 
	\ref o2scl::inte_gauss56_cern in long double precision
    */
    long double w6[6];

    inte_gauss56_coeffs_long_double() {
    
      x5[0]=0.04691007703066800360118656085030352L;
      x5[1]=0.23076534494715845448184278964989560L;
      x5[2]=0.5L;
      x5[3]=0.76923465505284154551815721035010440L;
      x5[4]=0.95308992296933199639881343914969648L;
  
      w5[0]=0.11846344252809454375713202035995868L;
      w5[1]=0.23931433524968323402064575741781910L;
      w5[2]=0.28444444444444444444444444444444444L;
      w5[3]=0.23931433524968323402064575741781910L;
      w5[4]=0.11846344252809454375713202035995868L;
  
      x6[0]=0.03376524289842398609384922275300270L;
      x6[1]=0.16939530676686774316930020249004733L;
      x6[2]=0.38069040695840154568474913915964403L;
      x6[3]=0.61930959304159845431525086084035597L;
      x6[4]=0.83060469323313225683069979750995267L;
      x6[5]=0.96623475710157601390615077724699730L;
  
      w6[0]=0.08566224618958517252014807108636645L;
      w6[1]=0.18038078652406930378491675691885806L;
      w6[2]=0.23395696728634552369493517199477550L;
      w6[3]=0.23395696728634552369493517199477550L;
      w6[4]=0.18038078652406930378491675691885806L;
      w6[5]=0.08566224618958517252014807108636645L;

    }
  
  };

#ifdef O2SCL_LD_TYPES
  
  class inte_gauss56_coeffs_float128 {

  public:

    /** \brief Fifth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in __float128 precision
    */
    __float128 x5[5];
  
    /** \brief Fifth order integration weights for 
	\ref o2scl::inte_gauss56_cern in __float128 precision
    */
    __float128 w5[5];
  
    /** \brief Sixth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in __float128 precision
    */
    __float128 x6[6];
  
    /** \brief Sixth order integration weights for 
	\ref o2scl::inte_gauss56_cern in __float128 precision
    */
    __float128 w6[6];

    inte_gauss56_coeffs_float128() {
    
      x5[0]=0.04691007703066800360118656085030352L;
      x5[1]=0.23076534494715845448184278964989560L;
      x5[2]=0.5L;
      x5[3]=0.76923465505284154551815721035010440L;
      x5[4]=0.95308992296933199639881343914969648L;
  
      w5[0]=0.11846344252809454375713202035995868L;
      w5[1]=0.23931433524968323402064575741781910L;
      w5[2]=0.28444444444444444444444444444444444L;
      w5[3]=0.23931433524968323402064575741781910L;
      w5[4]=0.11846344252809454375713202035995868L;
  
      x6[0]=0.03376524289842398609384922275300270L;
      x6[1]=0.16939530676686774316930020249004733L;
      x6[2]=0.38069040695840154568474913915964403L;
      x6[3]=0.61930959304159845431525086084035597L;
      x6[4]=0.83060469323313225683069979750995267L;
      x6[5]=0.96623475710157601390615077724699730L;
  
      w6[0]=0.08566224618958517252014807108636645L;
      w6[1]=0.18038078652406930378491675691885806L;
      w6[2]=0.23395696728634552369493517199477550L;
      w6[3]=0.23395696728634552369493517199477550L;
      w6[4]=0.18038078652406930378491675691885806L;
      w6[5]=0.08566224618958517252014807108636645L;

    }
  
  };

  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<50> > cpp_dec_float_50;
  
  class inte_gauss56_coeffs_cpp_dec_float_50 {

  public:

    /** \brief Fifth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in cpp_dec_float_50 precision
    */
    cpp_dec_float_50 x5[5];
  
    /** \brief Fifth order integration weights for 
	\ref o2scl::inte_gauss56_cern in cpp_dec_float_50 precision
    */
    cpp_dec_float_50 w5[5];
  
    /** \brief Sixth order integration abscissas for 
	\ref o2scl::inte_gauss56_cern in cpp_dec_float_50 precision
    */
    cpp_dec_float_50 x6[6];
  
    /** \brief Sixth order integration weights for 
	\ref o2scl::inte_gauss56_cern in cpp_dec_float_50 precision
    */
    cpp_dec_float_50 w6[6];

    inte_gauss56_coeffs_cpp_dec_float_50() {
    
      x5[0]=0.04691007703066800360118656085030352L;
      x5[1]=0.23076534494715845448184278964989560L;
      x5[2]=0.5L;
      x5[3]=0.76923465505284154551815721035010440L;
      x5[4]=0.95308992296933199639881343914969648L;
  
      w5[0]=0.11846344252809454375713202035995868L;
      w5[1]=0.23931433524968323402064575741781910L;
      w5[2]=0.28444444444444444444444444444444444L;
      w5[3]=0.23931433524968323402064575741781910L;
      w5[4]=0.11846344252809454375713202035995868L;
  
      x6[0]=0.03376524289842398609384922275300270L;
      x6[1]=0.16939530676686774316930020249004733L;
      x6[2]=0.38069040695840154568474913915964403L;
      x6[3]=0.61930959304159845431525086084035597L;
      x6[4]=0.83060469323313225683069979750995267L;
      x6[5]=0.96623475710157601390615077724699730L;
  
      w6[0]=0.08566224618958517252014807108636645L;
      w6[1]=0.18038078652406930378491675691885806L;
      w6[2]=0.23395696728634552369493517199477550L;
      w6[3]=0.23395696728634552369493517199477550L;
      w6[4]=0.18038078652406930378491675691885806L;
      w6[5]=0.08566224618958517252014807108636645L;

    }
  
  };

#endif
  
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
    class weights_t=inte_gauss56_coeffs_double>
    class inte_gauss56_cern : public inte<func_t,fp_t> {

  protected:

  const fp_t *w5, *x5, *w6, *x6;
  weights_t wgts;
    
  public:
  
  inte_gauss56_cern() {
    w5=&(wgts.w5[0]);
    x5=&(wgts.x5[0]);
    w6=&(wgts.w6[0]);
    x6=&(wgts.x6[0]);
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
    
    // AWS 10/4/19: I replaced std::abs here with a
    // more type-agnostic code
    err=(e6-e5)*rang;
    if (err<0.0) err=-err;

    return success;
  }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
