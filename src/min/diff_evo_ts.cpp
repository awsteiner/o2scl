/* 
   -------------------------------------------------------------------
   
   Copyright (C) 2006-2014, Andrew W. Steiner and Edwin van Leeuwen
   
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
#include <iostream>
#include <cmath>

#include <gsl/gsl_sf_bessel.h>

#include <o2scl/multi_funct.h>
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/diff_evo.h>

typedef boost::numeric::ublas::vector<double> ubvector;

using namespace std;
using namespace o2scl;

// A simple function with many local minima. A "greedy" minr
// would likely fail to find the correct minimum.
double func(size_t nvar, const ubvector &x) {
  double a, b;
  a=(x[0]-2.0);
  b=(x[1]+3.0);
  return -gsl_sf_bessel_J0(a)*gsl_sf_bessel_J0(b);
}

rng_gsl gr;

int init_function( size_t dim, const ubvector &x, ubvector &y ) {
  for (size_t i = 0; i < dim; ++i) {
    y[i] = 20*gr.random()-10;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  diff_evo<multi_funct<> > de;
  double result;
  ubvector init(2);
  
  multi_funct_fptr<> fx(func);
  mm_funct_fptr<ubvector> init_f( init_function );

  de.set_init_function( init_f );
  de.verbose = 1;
  de.ntrial=1000;
  
  // Perform the minimization
  de.mmin(2,init,result,fx);
  cout << "x: " << init[0] << " " << init[1] 
       << ", minimum function value: " << result << endl;
  cout << endl;

  // Test that it found the global minimum
  t.test_rel(init[0],2.0,1.0e-2,"another test - value");
  t.test_rel(init[1],-3.0,1.0e-2,"another test - value 2");
  t.test_rel(result,-1.0,1.0e-2,"another test - min");

  t.report();
  
  return 0;
}

