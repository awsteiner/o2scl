/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

/* Example: ex_conmin.cpp
   -------------------------------------------------------------------
   This gives an example of the use of a constrained minimizer. This
   code finds the global minimum of a two-dimensional function which
   is not well-defined outside the region of interest.
*/
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mmin_constr_spg.h>
#include <o2scl/mmin_simp2.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double func(size_t nv, const ubvector &x) {
  if (x[0]<0.0 || x[1]<0.0 || x[0]>100.0 || x[1]>100.0) {
    cout << "Outside constraint region." << endl;
  }
  double ret=(log(x[0])*x[0]*x[0]+1.0)*(sqrt(x[1])*(x[1]-1.0)+1.0);
  return ret;
}

int dfunc(size_t nv, ubvector &x, ubvector &g) {
  g[0]=(x[0]+2.0*x[0]*log(x[0]))*(sqrt(x[1])*(x[1]-1.0)+1.0);
  g[1]=(log(x[0])*x[0]*x[0]+1.0)*(sqrt(x[1])+(x[1]-1.0)/2.0/sqrt(x[1]));
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  static const size_t nv=2;
  
  // Specify the function to minimize and its gradient
  multi_funct mff11=func;
  grad_funct gff=dfunc;
  
  // The unconstrained minimizer
  mmin_simp2<> gm1;
  // The constrained minimizer
  mmin_constr_spg<> omp;

  // The constraints and the location of the minimum
  ubvector c1(nv), c2(nv), x(nv);
  double fmin;
    
  cout << "Simple minimizer: " << endl;

  // Initial guess
  for(size_t i=0;i<nv;i++) {
    x[i]=10.0;
  }

  // Minimize
  gm1.mmin(nv,x,fmin,mff11);
  cout << endl;

  cout << "Constrained minimizer: " << endl;

  // Initial guess
  for(size_t i=0;i<nv;i++) {
    x[i]=10.0;
  }

  // Set constraints
  for(size_t i=0;i<nv;i++) {
    c1[i]=1.0e-9;
    c2[i]=100.0;
  }
  omp.set_constraints(nv,c1,c2);
    
  // Minimize
  omp.mmin_de(nv,x,fmin,mff11,gff);

  // Output results
  cout << x[0] << " " << x[1] << " " << fmin << endl;

  // Test the constrained minimizer results
  t.test_rel(x[0],0.60655,1.0e-4,"x0");
  t.test_rel(x[1],1.0/3.0,1.0e-4,"x1");

  t.report();
  return 0;
}
// End of example
