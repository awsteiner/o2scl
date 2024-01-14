/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
  
#include <cmath>
#include <cstdlib>

#include <gsl/gsl_siman.h>

#include <o2scl/multi_funct.h>
#include <o2scl/funct.h>
#include <o2scl/anneal_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double funx(size_t nv, const ubvector &x) {
  double ret, a;
  a=x[0]-2.0;
  // AWS, 10/6/21: changed to add an a*a term. -exp(-a*a) is a bad
  // example because, even when the temperature is small, the routine
  // can escape the minimum to where the function is flat and then the
  // algorithm takes an anomalously long time to converge.
  ret=-exp(-a*a)+a*a;
  return ret;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  anneal_gsl<multi_funct,ubvector> ga;
  double result;
  ubvector init(1);
    
  multi_funct fx=funx;
    
  /// 1d to vectors
    
  init[0]=0.1;
  ga.tol_abs=1.0e-6;
  ga.mmin(1,init,result,fx);
  cout << init[0] << " " << result << endl;
  t.test_rel(init[0],2.0,1.0e-3,"another test - value");
  t.test_rel(result,-1.0,1.0e-3,"another test - min");
    
  // Test verbose=1
    
  init[0]=15.5;
  ga.verbose=1;
  ga.mmin(1,init,result,fx);

  t.report();
  
  return 0;
}

