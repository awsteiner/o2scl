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

/* Example: ex_inte.cpp
   -------------------------------------------------------------------
   An example to demonstrate numerical integration.
*/

#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagi_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qagil_gsl.h>
#include <o2scl/inte_adapt_cern.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

class cl {

public:

  // We'll use this to count the number of function
  // evaulations required by the integration routines
  int nf;

  // A function to be integrated
  double integrand(double x) {
    nf++;
    return exp(-x*x);
  }

  // Another function to be integrated
  double integrand2(double x) {
    nf++;
    return sin(2.0*x)+0.5;
  }
};

int main(void) {
  cl acl;
  test_mgr t;

  t.set_output_level(1);

  funct f1=std::bind(std::mem_fn<double(double)>
		       (&cl::integrand),&acl,std::placeholders::_1);
  funct f2=std::bind(std::mem_fn<double(double)>
		       (&cl::integrand2),&acl,std::placeholders::_1);

  // We don't need to specify the function type in the integration
  // objects, because we're using the default function type (type
  // funct).
  inte_qag_gsl<> g;
  inte_qagi_gsl<> gi;
  inte_qagiu_gsl<> gu;
  inte_qagil_gsl<> gl;
  inte_adapt_cern<> ca;

  // The result and the uncertainty
  double res, err;
  
  // An integral from -infinity to +infinity (the limits are ignored)
  acl.nf=0;
  int ret1=gi.integ_err(f1,0.0,0.0,res,err);
  cout << "inte_qagi_gsl: " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Result: " << res << " Uncertainty: " << err << endl;
  cout << "Number of iterations: " << gi.last_iter << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(res,sqrt(pi),1.0e-8,"inte 1");
  
  // An integral from 0 to +infinity (the second limit argument is
  // ignored in the line below)
  acl.nf=0;
  gu.integ_err(f1,0.0,0.0,res,err);
  cout << "inte_qagiu_gsl: " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Result: " << res << " Uncertainty: " << err << endl;
  cout << "Number of iterations: " << gu.last_iter << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(res,sqrt(pi)/2.0,1.0e-8,"inte 2");

  // An integral from -infinity to zero (the first limit argument is
  // ignored in the line below)
  acl.nf=0;
  gl.integ_err(f1,0.0,0.0,res,err);
  cout << "inte_qagil_gsl: " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Result: " << res << " Uncertainty: " << err << endl;
  cout << "Number of iterations: " << gl.last_iter << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(res,sqrt(pi)/2.0,1.0e-8,"inte 3");

  // An integral from 0 to 1 with the GSL integrator
  acl.nf=0;
  g.integ_err(f2,0.0,1.0,res,err);
  cout << "inte_qag_gsl: " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Result: " << res << " Uncertainty: " << err << endl;
  cout << "Number of iterations: " << g.last_iter << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(res,0.5+sin(1.0)*sin(1.0),1.0e-8,"inte 4");

  // The same integral with the CERNLIB integrator
  acl.nf=0;
  ca.integ_err(f2,0.0,1.0,res,err);
  cout << "inte_adapt_cern: " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Result: " << res << " Uncertainty: " << err << endl;
  cout << "Number of iterations: " << ca.last_iter << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(res,0.5+sin(1.0)*sin(1.0),1.0e-8,"inte 5");

  t.report();
  return 0;
}
// End of example
