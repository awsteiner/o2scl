/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file polylog.h
    \brief File defining \ref o2scl::polylog
*/
#ifndef O2SCL_POLYLOG_H
#define O2SCL_POLYLOG_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/lib_settings.h>
#include <gsl/gsl_sf_dilog.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Polylogarithms (approximate) \f$ Li_n(x)\f$ 
    
      This class is experimental.

      This gives an approximation to the polylogarithm functions.

      Only works at present for \f$n=0,1,...,6\f$. Uses GSL library
      for n=2.

      Uses linear interpolation for \f$-1<x<0\f$
      and a series expansion for \f$x<-1\f$

      \future
      - Give error estimate? 
      - Improve accuracy?
      - Use more sophisticated interpolation?
      - Add the series \f$Li(n,x)=x+2^{-n} x^2+3^{-n} x^3+...\f$ 
      for \f$ x \rightarrow 0\f$?
      - Implement for positive arguments < 1.0
      - Make another polylog class which implements series acceleration?

      For reference, there are exact relations
      \f[
      \mathrm{Li}_2 \left(\frac{1}{2}\right) =
      \frac{1}{12}\left[\pi^2-6\left(\ln 2\right)^2\right]
      \f]
      \f[
      \mathrm{Li}_3 \left(\frac{1}{2}\right) =
      \frac{1}{24}\left[ 4\left(\ln 2\right)^3 - 2 \pi^2 \ln 2 +
      21 \zeta (3) \right]
      \f]
      \f[
      \mathrm{Li}_{-1} (x) = \frac{x}{\left(1-x\right)^2}
      \f]
      \f[
      \mathrm{Li}_{-2} (x) = \frac{x\left(x+1\right)}{\left(1-x\right)^3}
      \f]

  */
  class polylog {

  public:

    /// 0-th order polylogarithm = \f$ x/(1-x)\f$
    double li0(double x);

    /// 1-st order polylogarithm = \f$ -\ln(1-x) \f$
    double li1(double x);

    /// 2-nd order polylogarithm
    double li2(double x);

    /// 3-rd order polylogarithm
    double li3(double x);

    /// 4-th order polylogarithm
    double li4(double x);

    /// 5-th order polylogarithm
    double li5(double x);

    /// 6-th order polylogarithm
    double li6(double x);

    polylog(); 
    ~polylog();

  protected:

#ifndef DOXYGEN_NO_O2NS

    double *arg;
    double *two;
    double *three;
    double *four;
    double *five;
    double *six;
    double li2neg1;
    double li4neg1;
    double li6neg1;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
