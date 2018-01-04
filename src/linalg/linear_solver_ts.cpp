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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <o2scl/linear_solver.h>
#include <o2scl/test_mgr.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  
  ubmatrix gm1(5,5), gm2(5,5);
  ubvector gv1(5), gv2(5), gv3(5);
  permutation gp1(5);
  int sig, ret;

  // Test LU solve using O2SCL
  {

    for(size_t i=0;i<5;i++) {
      gv1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gm1(i,j)=1.0/(i+j+1);
	gm2(i,j)=1.0/(i+j+1);
      }
    }

    cout << "A: " << endl;
    //cout << gm1 << endl;
    cout << "b: " << endl;
    //cout << gv1 << endl;

    // -------------------------------------------------

    linear_solver_LU<ubvector,ubmatrix> lus;
    lus.solve(5,gm1,gv1,gv2);

    cout << "x: " << endl;
    //cout << gv2 << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    cout << "A*x: " << endl;
    //cout << gv3 << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"solve 1");
  }

  // Test QR solve using O2SCL
  {

    for(size_t i=0;i<5;i++) {
      gv1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gm1(i,j)=1.0/(i+j+1);
	gm2(i,j)=1.0/(i+j+1);
      }
    }

    cout << "A: " << endl;
    //cout << gm1 << endl;
    cout << "b: " << endl;
    //cout << gv1 << endl;

    // -------------------------------------------------

    linear_solver_QR<ubvector,ubmatrix> lus;
    lus.solve(5,gm1,gv1,gv2);

    cout << "x: " << endl;
    //cout << gv2 << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    cout << "A*x: " << endl;
    //cout << gv3 << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"solve 1");
  }

  // Test HH solve using O2SCL
  {

    for(size_t i=0;i<5;i++) {
      gv1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gm1(i,j)=1.0/(i+j+1);
	gm2(i,j)=1.0/(i+j+1);
      }
    }

    cout << "A: " << endl;
    //cout << gm1 << endl;
    cout << "b: " << endl;
    //cout << gv1 << endl;

    // -------------------------------------------------

    linear_solver_HH<ubvector,ubmatrix> lus;
    lus.solve(5,gm1,gv1,gv2);

    cout << "x: " << endl;
    //cout << gv2 << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    cout << "A*x: " << endl;
    //cout << gv3 << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"solve 1");
  }

  t.report();
  return 0;
}

