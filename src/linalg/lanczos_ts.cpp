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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#include <o2scl/test_mgr.h>
#include <o2scl/lanczos.h>
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  gsl_matrix *oldgsl;
  int i, j, k, size=10;

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  {
    using namespace o2scl_linalg;

    rng_gsl gr;

    ubvector eigen(size);
    ubvector diag(size);
    ubvector off_diag(size);
    ubmatrix mat(size,size);

    cout.setf(ios::scientific);
    cout.precision(5);

    for(int iter=0;iter<100;iter++) {

      // Make random symmetric matrix
      for(i=0;i<size;i++) {
	for(j=i;j<size;j++) {
	  mat(i,j)=gr.random();
	  mat(j,i)=mat(i,j);
	}
      }

      // Compute the Lanczos eigenvalues
    
      lanczos<ubvector,ubmatrix> la;

      gsl_vector *gv=gsl_vector_alloc(size);
      for(int m=1;m<=size;m++) {
	la.eigenvalues(size,mat,m,eigen,diag,off_diag);
      }
      for(i=0;i<size;i++) {
	gsl_vector_set(gv,i,eigen[i]);
      }
      gsl_sort_vector(gv);

      // Compute the eigenvalues using GSL 

      oldgsl=gsl_matrix_alloc(size,size);
      gsl_eigen_symm_workspace *ww=gsl_eigen_symm_alloc(size);
      gsl_vector *eval=gsl_vector_alloc(size);

      for(i=0;i<size;i++) {
	for(j=0;j<size;j++) {
	  gsl_matrix_set(oldgsl,i,j,mat(i,j));
	}
      }

      gsl_eigen_symm(oldgsl,eval,ww);
      gsl_sort_vector(eval);

      // Now compute the eigenvalues of the tridiagonal matrix
  
      // Shift sub-diagonals:
      for(i=size-1;i>=1;i--) off_diag[i]=off_diag[i-1];
  
      ubvector d2(10), od2(10);
      for(i=0;i<size;i++) {
	d2[i]=diag[i];
	od2[i]=off_diag[i];
      }

      la.eigen_tdiag(size,diag,off_diag);
      gsl_vector *gv2=gsl_vector_alloc(size);
      for(i=0;i<size;i++) {
	gsl_vector_set(gv2,i,diag[i]);
      }
      gsl_sort_vector(gv2);

      // Compute with imtql
      la.eigen_tdiag(size,d2,od2);
    
      // Compare the different approaches
      for(i=0;i<size;i++) {
	t.test_rel(gsl_vector_get(gv,i),gsl_vector_get(gv2,i),1.0e-8,"eigen1");
	t.test_rel(gsl_vector_get(gv,i),gsl_vector_get(eval,i),1.0e-8,"eigen2");
	t.test_rel(gsl_vector_get(gv,i),d2[i],1.0e-8,"eigen3");
      }

      gsl_vector_free(eval);
      gsl_eigen_symm_free(ww);
      gsl_matrix_free(oldgsl);

    }

  }

  t.report();

  return 0;
}

