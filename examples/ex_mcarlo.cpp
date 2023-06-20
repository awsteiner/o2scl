/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/* Example: ex_mcarlo.cpp
   -------------------------------------------------------------------
   An example to demonstrate multidimensional integration. See "License 
   Information" section of the documentation for license information.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>

/// For M_PI
#include <gsl/gsl_math.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double test_fun(size_t nv, const ubvector &x) {
  double y=1.0/(1.0-cos(x[0])*cos(x[1])*cos(x[2]))/M_PI/M_PI/M_PI;
  return y;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
 
  cout.setf(ios::scientific);

  double exact=1.3932039296;
  double res;

  double err;
  
  mcarlo_vegas<> gm;
  ubvector a(3), b(3);
  a[0]=0.0; b[0]=M_PI;
  a[1]=0.0; b[1]=M_PI;
  a[2]=0.0; b[2]=M_PI;

  multi_funct tf=test_fun;

  gm.n_points=100000;
  gm.minteg_err(tf,3,a,b,res,err);

  cout << res << " " << exact << " " << (res-exact)/err << endl;
  t.test_rel(res,exact,err*10.0,"O2scl");

  t.report();
 
  return 0;
}
// End of example
