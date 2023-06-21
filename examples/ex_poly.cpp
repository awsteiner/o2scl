/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/

// sphinx-example-start
/* Example: ex_poly.cpp
   -------------------------------------------------------------------
   Demonstrate the solution of the Chebyshev polynomials. See "License 
   Information" section of the documentation for license information.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/poly.h>
// For pi
#include <o2scl/constants.h>
#include <o2scl/vector.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);

  test_mgr t;
  t.set_output_level(1);

  typedef boost::numeric::ublas::vector<double> ubvector;

  // Quadratic solver
  quadratic_real_coeff_gsl2<> quad;
  // Cubic solver
  cubic_real_coeff_cern<> cubic;
  // Quartic solver
  quartic_real_coeff_cern<> quart;
  // Generic polynomial solver
  poly_real_coeff_gsl<> gen;

  // Storage for the roots
  ubvector v(5);
  double d;
  std::vector<std::complex<double> > ca(5);

  // The second order polynomial
  cout << "Second order roots: " << endl;
  quad.solve_r(2.0,0.0,-1.0,v[0],v[1]);

  // Sort the roots and compare with the exact results
  vector_sort<ubvector,double>(2,v);
  for(size_t i=0;i<2;i++) {
    double exact=cos(pi*(((double)(2-i))-0.5)/2.0);
    cout << v[i] << " " << exact << endl;
    t.test_abs(v[i],exact,1.0e-14,"2nd order");
  }
  cout << endl;

  // The third order polynomial
  cout << "Third order roots: " << endl;
  cubic.solve_rc(4.0,0.0,-3.0,0.0,v[0],ca[0],ca[1]);

  // Sort the roots and compare with the exact results
  v[1]=ca[0].real();
  v[2]=ca[1].real();
  vector_sort<ubvector,double>(3,v);
  for(size_t i=0;i<3;i++) {
    double exact=cos(pi*(((double)(3-i))-0.5)/3.0);
    cout << v[i] << " " << exact << endl;
    if (i==1) {
      t.test_abs(v[i],exact,1.0e-14,"3rd order");
    } else {
      t.test_abs(v[i],exact,1.0e-14,"3rd order");
    }
  }
  cout << endl;

  // The fourth order polynomial
  cout << "Fourth order roots: " << endl;
  quart.solve_rc(8.0,0.0,-8.0,0.0,1.0,ca[0],ca[1],ca[2],ca[3]);

  // Sort the roots and compare with the exact results
  for(size_t i=0;i<4;i++) v[i]=ca[i].real();
  vector_sort<ubvector,double>(4,v);
  for(size_t i=0;i<4;i++) {
    double exact=cos(pi*(((double)(4-i))-0.5)/4.0);
    cout << v[i] << " " << exact << endl;
    t.test_abs(v[i],exact,1.0e-14,"4th order");
  }
  cout << endl;

  // The fifth order polynomial
  cout << "Fifth order roots: " << endl;
  vector<double> co={16.0,0.0,-20.0,0.0,5.0,0.0};
  gen.solve_rc_arr(5,co,ca);

  // Sort the roots and compare with the exact results
  for(size_t i=0;i<5;i++) v[i]=ca[i].real();
  vector_sort<ubvector,double>(5,v);
  for(size_t i=0;i<5;i++) {
    double exact=cos(pi*(((double)(5-i))-0.5)/5.0);
    cout << v[i] << " " << exact << endl;
    t.test_abs(v[i],exact,1.0e-14,"5th order");
  }
  cout << endl;

  t.report();
  return 0;
}
