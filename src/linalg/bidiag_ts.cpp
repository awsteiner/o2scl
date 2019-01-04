/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2019 Andrew W. Steiner

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

#include <gsl/gsl_linalg.h>

#include <o2scl/test_mgr.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "bidiag.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  static const size_t arr_size=10;

  gsl_matrix *gm1=gsl_matrix_alloc(arr_size,arr_size);
  gsl_matrix *gm2=gsl_matrix_alloc(arr_size,arr_size);
  gsl_matrix *gm3=gsl_matrix_alloc(arr_size,arr_size);
  gsl_vector *gv1=gsl_vector_alloc(arr_size);
  gsl_vector *gv2=gsl_vector_alloc(arr_size-1);
  gsl_vector *gv3=gsl_vector_alloc(arr_size);
  gsl_vector *gv4=gsl_vector_alloc(arr_size-1);
  ubmatrix om1(arr_size,arr_size);
  ubmatrix om2(arr_size,arr_size);
  ubmatrix om3(arr_size,arr_size);
  ubvector ov1(arr_size);
  ubvector ov2(arr_size-1);
  ubvector ov3(arr_size);
  ubvector ov4(arr_size-1);

  // Setup original matrix
  for(size_t i=0;i<arr_size;i++) {
    for(size_t j=0;j<arr_size;j++) {
      gsl_matrix_set(gm1,i,j,fabs((double)(i-j)));
      om1(i,j)=fabs((double)(i-j));
    }
  }

  // Test decomposition
  gsl_linalg_bidiag_decomp(gm1,gv1,gv2);
  bidiag_decomp(arr_size,arr_size,om1,ov1,ov2);

  t.test_rel_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),5.0e-12,"m1");
  t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),1.0e-12,"v1");
  t.test_rel_vec(arr_size-1,ov2,gsl_vector_wrap(gv2),1.0e-12,"v2");

  // Test unpack
  gsl_linalg_bidiag_unpack(gm1,gv1,gm2,gv2,gm3,gv3,gv4);
  bidiag_unpack(arr_size,arr_size,om1,ov1,om2,ov2,om3,ov3,ov4);

  t.test_abs_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),1.0e-12,"m2");
  t.test_abs_mat(arr_size,arr_size,om3,gsl_matrix_wrap(gm3),1.0e-12,"m3");
  t.test_rel_vec(arr_size,ov3,gsl_vector_wrap(gv3),1.0e-12,"v3");
  t.test_rel_vec(arr_size-1,ov4,gsl_vector_wrap(gv4),1.0e-12,"v4");

  // Setup original matrix
  for(size_t i=0;i<arr_size;i++) {
    for(size_t j=0;j<arr_size;j++) {
      gsl_matrix_set(gm1,i,j,fabs((double)(i-j)));
      om1(i,j)=fabs((double)(i-j));
    }
  }

  // Perform initial decomposition
  gsl_linalg_bidiag_decomp(gm1,gv1,gv2);
  bidiag_decomp(arr_size,arr_size,om1,ov1,ov2);

  // Test unpack2
  gsl_linalg_bidiag_unpack2(gm1,gv1,gv2,gm2);
  bidiag_unpack2(arr_size,arr_size,om1,ov1,ov2,om2);

  t.test_abs_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),1.0e-12,"m4");

  // Setup original matrix
  for(size_t i=0;i<arr_size;i++) {
    for(size_t j=0;j<arr_size;j++) {
      gsl_matrix_set(gm1,i,j,fabs((double)(i-j)));
      om1(i,j)=fabs((double)(i-j));
    }
  }

  // Perform initial decomposition
  gsl_linalg_bidiag_decomp(gm1,gv1,gv2);
  bidiag_decomp(arr_size,arr_size,om1,ov1,ov2);

  // Test unpack_B
  gsl_linalg_bidiag_unpack_B(gm1,gv3,gv4);
  bidiag_unpack_B(arr_size,arr_size,om1,ov3,ov4);

  t.test_rel_vec(arr_size,ov3,gsl_vector_wrap(gv3),1.0e-12,"v5");
  t.test_rel_vec(arr_size-1,ov4,gsl_vector_wrap(gv4),1.0e-12,"v6");

  t.report();
  return 0;
}

