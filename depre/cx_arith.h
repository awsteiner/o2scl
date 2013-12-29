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
#ifndef O2SCL_CX_ARITH_H
#define O2SCL_CX_ARITH_H

/** \file cx_arith.h
    \brief Complex arithmetic

    \note One should be careful about including this header file,
    especially in other header files as it overloads some of the
    standard mathematical functions, i.e. <tt>sqrt()</tt>. If you need
    access to both these functions and the standard functions for
    <tt>double</tt> objects, the easiest way is probably to include
    this file in a source (not header file) and use <tt>using
    namespace std</tt>.

    \note This used to be in a separate namespace, called
    <tt>o2scl_arith</tt>, but this causes problems with Koenig lookup
    in template classes for operator*() when defined for vector
    addition (for example).

    \future Define operators with assignment for complex + double?
*/
#include <iostream>
#include <complex>
#include <cmath>

// For M_PI and M_PI_2
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Convert a complex number to GSL form
  gsl_complex complex_to_gsl(std::complex<double> &d);
  
  /// Convert a complex number to STL form
  std::complex<double> gsl_to_complex(gsl_complex &g);
  
  /// \name Binary operators for two complex numbers
  //@{
  /// Add two complex numbers
  gsl_complex operator+(gsl_complex x, gsl_complex y);

  /// Subtract two complex numbers
  gsl_complex operator-(gsl_complex x, gsl_complex y);

  /// Multiply two complex numbers
  gsl_complex operator*(gsl_complex x, gsl_complex y);

  /// Divide two complex numbers
  gsl_complex operator/(gsl_complex x, gsl_complex y);
  //@}

  /// \name Binary operators with assignment for two complex numbers
  //@{
  /// Add a complex number
  gsl_complex operator+=(gsl_complex &x, gsl_complex y);

  /// Subtract a complex number
  gsl_complex operator-=(gsl_complex &x, gsl_complex y);

  /// Multiply a complex number
  gsl_complex operator*=(gsl_complex &x, gsl_complex y);

  /// Divide a complex number
  gsl_complex operator/=(gsl_complex &x, gsl_complex y);
  //@}

  /// \name Binary operators with assignment for a complex and real
  //@{
  /// Add a complex and real number
  gsl_complex operator+(gsl_complex x, double y);

  /// Add a complex and real number
  gsl_complex operator+(double y, gsl_complex x);

  /// Subtract a complex and real number
  gsl_complex operator-(gsl_complex x, double y);

  /// Subtract a complex and real number
  gsl_complex operator-(double y, gsl_complex x);

  /// Multiply a complex and real number
  gsl_complex operator*(gsl_complex x, double y);

  /// Multiply a complex and real number
  gsl_complex operator*(double y, gsl_complex x);

  /// Divide a complex and real number
  gsl_complex operator/(gsl_complex x, double y);
  //@}

  /// \name Miscellaneous functions
  //@{
  double arg(gsl_complex x);
  double abs(gsl_complex x);
  double abs2(gsl_complex z);
  gsl_complex conjugate(gsl_complex a);
  //@}

  /// \name Square root and exponent functions
  //@{
  gsl_complex sqrt(gsl_complex a);
  gsl_complex sqrt_real(double x);
  gsl_complex pow(gsl_complex a, gsl_complex b);
  gsl_complex pow_real(gsl_complex a, double b);
  //@}
  
  /// \name Logarithmic and exponential functions
  //@{
  double logabs(gsl_complex z);
  gsl_complex exp(gsl_complex a);
  gsl_complex log(gsl_complex a);
  gsl_complex log10(gsl_complex a);
  gsl_complex log_b(gsl_complex a, gsl_complex b);
  //@}

  /// \name Trigonometric functions
  //@{
  gsl_complex sin(gsl_complex a);
  gsl_complex cos(gsl_complex a);
  gsl_complex tan(gsl_complex a);
  gsl_complex sec(gsl_complex a);
  gsl_complex csc(gsl_complex a);
  gsl_complex cot(gsl_complex a);
  gsl_complex asin(gsl_complex a);
  gsl_complex asin_real(double a);
  gsl_complex acos(gsl_complex a);
  gsl_complex acos_real(double a);
  gsl_complex atan(gsl_complex a);
  gsl_complex asec(gsl_complex a);
  gsl_complex asec_real(double a);
  gsl_complex acsc(gsl_complex a);
  gsl_complex acsc_real(double a);
  gsl_complex acot(gsl_complex a);
  //@}

  /// \name Hyperbolic trigonometric functions
  //@{
  gsl_complex sinh(gsl_complex a);
  gsl_complex cosh(gsl_complex a);
  gsl_complex tanh(gsl_complex a);
  gsl_complex sech(gsl_complex a);
  gsl_complex csch(gsl_complex a);
  gsl_complex coth(gsl_complex a);
  gsl_complex asinh(gsl_complex a);
  gsl_complex acosh(gsl_complex a);
  gsl_complex acosh_real(double a);
  gsl_complex atanh(gsl_complex a);
  gsl_complex atanh_real(double a);
  gsl_complex asech(gsl_complex a);
  gsl_complex acsch(gsl_complex a);
  gsl_complex acoth(gsl_complex a);
  //@}

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
