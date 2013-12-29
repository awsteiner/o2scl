/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2012, Andrew W. Steiner
  
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

/* Example: ex_deriv.cpp
   -------------------------------------------------------------------
   An example to demonstrate numerical differentiation
*/

#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/deriv_cern.h>

using namespace std;
using namespace o2scl;

class cl {

public:

  // This is the function we'll take the derivative of
  double function(double x) {
    return sin(2.0*x)+0.5;
  }
};

int main(void) {

  test_mgr t;
  t.set_output_level(2);

  // The class and associated function
  cl acl;
  funct_mfptr<cl> f1(&acl,&cl::function);

  deriv_gsl<> gd;
  // Note that the GSL derivative routine requires an initial stepsize
  gd.h=1.0e-3;
  deriv_cern<> cd;
  
  // Compute the first derivative using the deriv_gsl class and
  // verify that the answer is correct
  double d1=gd.deriv(1.0,f1);
  t.test_rel(d1,2.0*cos(2.0),1.0e-10,"deriv_gsl");
  
  // Compute the first derivative using the deriv_cern class and
  // verify that the answer is correct
  double d2=cd.deriv(1.0,f1);
  t.test_rel(d2,2.0*cos(2.0),1.0e-10,"deriv_cern");

  // Compute the second derivative also
  double d3=gd.deriv2(1.0,f1);
  t.test_rel(d3,-4.0*sin(2.0),5.0e-7,"deriv_gsl");
  
  double d4=cd.deriv2(1.0,f1);
  t.test_rel(d4,-4.0*sin(2.0),1.0e-8,"deriv_cern");

  t.report();
  return 0;
}
// End of example

