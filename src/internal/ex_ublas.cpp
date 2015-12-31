/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

/* Example: ex_ublas.cpp
   -------------------------------------------------------------------
   This gives an example of the how one can use the uBlas vector
   and matrix types in an O2scl equation solver class
*/

#include <iostream>
#include <o2scl/vec_arith.h>
#include <o2scl/mm_funct.h>
#include <o2scl/gsl_mroot_hybrids.h>
#include <o2scl/test_mgr.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;

// Some convenient typedefs
typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/// An allocation object for ublas vectors 
class ubvector_alloc {
public:
  /// Ensure space in \c u for \c i elements
  void allocate(ubvector &u, size_t i) { u.resize(i); }
  /// Free memory
  void free(ubvector &u) { }
};

/// An allocation object for ublas matrices
class ubmatrix_alloc {
public:
  /// Ensure space in \c u for \c i elements
  void allocate(ubmatrix &u, size_t i, size_t j) { u.resize(i,j); }
  /// Free memory
  void free(ubmatrix &u, size_t i) { }
};

/// A class with a function to solve
class cl {

public:

  int function(size_t nv, const ubvector &x, ubvector &y, int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

};

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  cl acl;

  // Create a new function object and solver with the 
  // ublas vector and matrix types.
  mm_funct_mfptr<cl,int,ubvector> f1(&acl,&cl::function);
  gsl_mroot_hybrids<int,mm_funct<int,ubvector>,ubvector,ubvector,
    ubvector_alloc,ubmatrix,ubmatrix,ubmatrix_alloc> gmr;
  
  ubvector sol(2);
  sol[0]=0.5;
  sol[1]=0.5;
  int pa;
  gmr.msolve(2,sol,pa,f1);
  cout << sol[0] << " " << sol[1] << endl;
  t.test_rel(sol[0],0.25,1.0e-6,"Solution 1.");
  t.test_rel(sol[1],0.2,1.0e-6,"Solution 2.");

  t.report();
  return 0;
}
// End of example
