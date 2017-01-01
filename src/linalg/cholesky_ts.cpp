/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2017, Andrew W. Steiner

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

#include <o2scl/cholesky.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  
  if (o2scl_settings.range_check()) {
    cout << "O2scl range checking on." << endl;
  } else {
    cout << "O2scl range checking off." << endl;
  }

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    gsl_vector *gv1=gsl_vector_alloc(5);
    gsl_vector *gv2=gsl_vector_alloc(5);
    gsl_vector *gv3=gsl_vector_alloc(5);
    gsl_matrix *gm1=gsl_matrix_alloc(5,5);
    gsl_matrix *gm2=gsl_matrix_alloc(5,5);
    ubvector ov1(5), ov2(5), ov3(5);
    ubmatrix om1(5,5), om2(5,5);
      
    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      gsl_vector_set(gv2,i,cos(((double)(i))));
      ov2[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,((double)(1))/(1.0+i+j));
	om1(i,j)=((double)(1))/(1.0+i+j);
	gsl_matrix_set(gm2,i,j,((double)(1))/(1.0+i+j));
	om2(i,j)=((double)(1))/(1.0+i+j);
      }
    }

    // Test cholesky decomposition
      
    gsl_linalg_cholesky_decomp(gm1);
    cholesky_decomp(5,om1);
      
#ifdef O2SCL_EIGEN
    
    cout << "Performing eigen test." << endl;
    Eigen::MatrixXd em1(5,5);
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	em1(i,j)=((double)(1))/(1.0+i+j);
      }
    }
    cholesky_decomp(5,em1);
    t.test_rel_mat(5,5,em1,om1,1.0e-11,"ch decomp");
      
#endif
    
    t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),5.0e-12,"cholesky decomp");
      
    // Test solve 

    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      gsl_vector_set(gv2,i,cos(((double)(i))));
      ov2[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,((double)(1))/(1.0+i+j));
	om1(i,j)=((double)(1))/(1.0+i+j);
	gsl_matrix_set(gm2,i,j,((double)(1))/(1.0+i+j));
	om2(i,j)=((double)(1))/(1.0+i+j);
      }
    }

    // -------------------------------------------------

    gsl_linalg_cholesky_decomp(gm1);
    gsl_linalg_cholesky_solve(gm1,gv2,gv1);
    gsl_blas_dgemv(CblasTrans,1.0,gm2,gv1,0.0,gv3);

    t.test_rel_vec(5,gsl_vector_wrap(gv2),
		   gsl_vector_wrap(gv3),1.0e-10,"solve 1");
    
    // -------------------------------------------------

    cholesky_decomp(5,om1);
    cholesky_solve(5,om1,ov2,ov1);
    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,om2,ov1,0.0,ov3);

    t.test_rel_vec(5,ov2,ov3,1.0e-10,"solve 2");

    // -------------------------------------------------

    // Test invert

    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,((double)(1))/(1.0+i+j));
	om1(i,j)=((double)(1))/(1.0+i+j);
	gsl_matrix_set(gm2,i,j,((double)(1))/(1.0+i+j));
	om2(i,j)=((double)(1))/(1.0+i+j);
      }
    }

    gsl_linalg_cholesky_decomp(gm1);
    gsl_linalg_cholesky_invert(gm1);

    cholesky_decomp(5,om1);
    cholesky_invert<ubmatrix>(5,om1);
    
    t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),5.0e-12,"cholesky invert 1");

  }

  t.report();
  return 0;
}

