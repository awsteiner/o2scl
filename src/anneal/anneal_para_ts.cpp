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
  
#include <cmath>
#include <cstdlib>

#include <gsl/gsl_siman.h>

#include <o2scl/multi_funct.h>
#include <o2scl/funct.h>
#include <o2scl/anneal_para.h>
#include <o2scl/test_mgr.h>

#ifdef O2SCL_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double funx(size_t nv, const ubvector &x) {
  double ret, a;
  a=x[0]-2.0;
  ret=-exp(-a*a);
  return ret;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  anneal_para<multi_funct,ubvector> ga;

#ifdef O2SCL_OPENMP
  ga.n_threads=omp_get_max_threads();
#endif
  
  double result;
  ubvector init(1);
    
  multi_funct fx=funx;
    
  init[0]=0.1;
  ga.verbose=1;
  ga.tol_abs=1.0e-6;
  ga.mmin(1,init,result,fx);
  cout << init[0] << " " << result << endl;
  t.test_rel(init[0],2.0,1.0e-3,"another test - value");
  t.test_rel(result,-1.0,1.0e-3,"another test - min");
    
  t.report();
  
  return 0;
}

