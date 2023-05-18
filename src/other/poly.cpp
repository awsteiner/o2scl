/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
/* poly/solve_cubic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/poly.h>
// For is_finite()
#include <o2scl/misc.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int quadratic_real_gsl::solve_r(const double a2, const double b2,
                                const double c2, 
                                double &x1, double &x2) {
  // This function returns the number of real roots
  return gsl_poly_solve_quadratic(a2,b2,c2,&x1,&x2);
}

int quadratic_real_coeff_gsl::solve_rc(const double a2, const double b2,
                                       const double c2, 
                                       std::complex<double> &x1,
                                       std::complex<double> &x2) {
  gsl_complex z0, z1;
  // This function returns the number of complex roots
  int ret=gsl_poly_complex_solve_quadratic(a2,b2,c2,&z0,&z1);
  x1.real(GSL_REAL(z0));
  x1.imag(GSL_IMAG(z0));
  x2.real(GSL_REAL(z1));
  x2.imag(GSL_IMAG(z1));
  return 2-ret;
}

int cubic_real_coeff_gsl::solve_rc(const double a3, const double b3,
                                   const double c3, const double d3,
                                   double &x1,
                                   std::complex<double> &x2,
                                   std::complex<double> &x3) {
  
  gsl_complex z0, z1, z2;
  // This function always returns 3
  int ret=gsl_poly_complex_solve_cubic(b3/a3,c3/a3,d3/a3,&z0,&z1,&z2);
  if (GSL_IMAG(z0)==0.0) {
    x1=GSL_REAL(z0);
    x2.real(GSL_REAL(z1));
    x2.imag(GSL_IMAG(z1));
    x3.real(GSL_REAL(z2));
    x3.imag(GSL_IMAG(z2));
  } else if (GSL_IMAG(z1)==0.0) {
    x1=GSL_REAL(z1);
    x2.real(GSL_REAL(z0));
    x2.imag(GSL_IMAG(z0));
    x3.real(GSL_REAL(z2));
    x3.imag(GSL_IMAG(z2));
  } else if (GSL_IMAG(z2)==0.0) {
    x1=GSL_REAL(z2);
    x2.real(GSL_REAL(z0));
    x2.imag(GSL_IMAG(z0));
    x3.real(GSL_REAL(z1));
    x3.imag(GSL_IMAG(z1));
  } else {
    O2SCL_ERR("GSL returned three complex roots.",o2scl::exc_einval);
  }

  // Return the number of real roots
  if (x2.imag()==0.0 && x3.imag()==0.0) return 3;
  return 1;
}

