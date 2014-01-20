/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

/* Example: ex_minte.cpp
   -------------------------------------------------------------------
   An example to demonstrate multidimensional integration with
   inte_multi_comp. 
*/

#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_multi_comp.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/mcarlo_vegas.h>

/// For pi
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double test_fun(size_t nv, const ubvector &x) {
  return sqrt(x[0]*x[0]*x[0]+x[1]*x[1]*x[1]+x[2]*x[2]*x[2]+
	      x[0]*x[1]*x[1]*x[2]);
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
 
  cout.setf(ios::scientific);

  double exact=0.8477219883;
  double res, err;

  inte_multi_comp<> ci;

  // The integration limits
  ubvector a(3), b(3);
  a[0]=0.0; b[0]=1.0;
  a[1]=0.0; b[1]=1.0;
  a[2]=0.0; b[2]=1.0;

  /// The individual integration objects for inte_multi_comp
  inte_qag_gsl<> gl[3];
  ci.set_oned_inte(gl[0],0);
  ci.set_oned_inte(gl[1],1);
  ci.set_oned_inte(gl[2],2);

  // Specify the function to integrate
  multi_funct_fptr<> tf(test_fun);

  // Integrate with inte_multi_comp and test result
  res=ci.minteg(tf,3,a,b);
  cout << res << " " << fabs(res-exact)/exact << endl;
  t.test_rel(res,exact,1.0e-9,"inte_multi_comp");

  // Compare to Monte Carlo integration
  mcarlo_vegas<multi_funct<>,ubvector,int,rng_gsl> gv;
  gv.n_points=100000;
  gv.minteg_err(tf,3,a,b,res,err);
  cout << res << " " << err << " " << fabs(res-exact)/exact << endl;
  t.test_rel(res,exact,5.0e-5,"mcarlo_vegas");

  t.report();
 
  return 0;
}
// End of example
