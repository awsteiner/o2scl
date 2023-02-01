/*
  -------------------------------------------------------------------

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

  -------------------------------------------------------------------
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <o2scl/lu.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  
  if (o2scl_settings.range_check()) {
    cout << "O2scl range checking on." << endl;
  } else {
    cout << "O2scl range checking off." << endl;
  }

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    for(size_t ki=1;ki<10;ki++) {

      // -------------------------------------------------
      // Test LU decomp

      gsl_matrix *gm1=gsl_matrix_alloc(5,5);
      gsl_permutation *gp1=gsl_permutation_alloc(5);
      ubmatrix om1(5,5);
      permutation op1(5);
      int sig;

      for(size_t i=0;i<5;i++) {
	for(size_t j=0;j<5;j++) {
	  double diag=0.0;
	  if (i==j) diag+=3.0;
	  gsl_matrix_set(gm1,i,j,1.0+diag+sin(((double)ki))+
			 sin((double)(i))+tan(((double)(j))));
	  om1(i,j)=1.0+diag+sin(((double)ki))+sin((double)(i))+
	    tan(((double)(j)));
	}
      }

      gsl_linalg_LU_decomp(gm1,gp1,&sig);
      LU_decomp(5,om1,op1,sig);

      matrix_out(cout,5,5,om1);
      cout << endl;
      matrix_out(cout,5,5,gsl_matrix_wrap(gm1));
      cout << endl;
      t.test_rel_nonzero_mat(5,5,om1,gsl_matrix_wrap(gm1),1.0e-8,
			     1.0e-14,"LU decomp");

      // -------------------------------------------------
      // Prepare data for solve and refine

      ubmatrix om2(5,5);
      ubvector ov1(5), ov2(5), ov3(5);
      gsl_vector *gv1=gsl_vector_alloc(5);
      gsl_matrix *gm2=gsl_matrix_alloc(5,5);
      gsl_vector *gv2=gsl_vector_alloc(5);
      gsl_vector *gv3=gsl_vector_alloc(5);

      for(size_t i=0;i<5;i++) {
        gsl_vector_set(gv1,i,cos(((double)(i+ki))));
        ov1[i]=cos(((double)(i+ki)));
        for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,1.0/(i+j+1));
          om1(i,j)=1.0/(i+j+1);
	  gsl_matrix_set(gm2,i,j,1.0/(i+j+1));
          om2(i,j)=1.0/(i+j+1);
        }
      }

      dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,om1,ov1,0.0,ov2);

      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv2,i,ov2[i]);
	gsl_vector_set(gv1,i,0.0);
        ov1[i]=0.0;
      }

      // -------------------------------------------------
      // GSL version of solve and refine

      gsl_linalg_LU_decomp(gm1,gp1,&sig);
      gsl_linalg_LU_solve(gm1,gp1,gv2,gv1);

      gsl_blas_dgemv(CblasNoTrans,1.0,gm2,gv1,0.0,gv3);
      t.test_rel_vec(5,gsl_vector_wrap(gv2),
		     gsl_vector_wrap(gv3),1.0e-10,"solve 1");
    
      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv1,i,gsl_vector_get(gv1,i)*1.2);
      }

      gsl_linalg_LU_refine(gm2,gm1,gp1,gv2,gv1,gv3);

      gsl_blas_dgemv(CblasNoTrans,1.0,gm2,gv1,0.0,gv3);
      t.test_rel_vec(5,gsl_vector_wrap(gv2),
		     gsl_vector_wrap(gv3),1.0e-10,"refine 1");

      // -------------------------------------------------
      // O2scl version of solve and refine

      LU_decomp(5,om1,op1,sig);
      LU_solve(5,om1,op1,ov1,ov2);

      dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,om2,ov1,0.0,ov3);
      t.test_rel_vec(5,ov2,ov3,1.0e-10,"solve 2 (paren)");

      for(size_t i=0;i<5;i++) ov1[i]*=1.2;

      LU_refine(5,om2,om1,op1,ov2,ov1,ov3);

      dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,om2,ov1,0.0,ov3);
      t.test_rel_vec(5,ov2,ov3,1.0e-10,"refine 2 (paren)");

      // -------------------------------------------------
      // Test invert

      for(size_t i=0;i<5;i++) {
	for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,ki/(1.0+i+j));
	  om1(i,j)=ki/(1.0+i+j);
	  gsl_matrix_set(gm2,i,j,ki/(1.0+i+j));
	  om2(i,j)=ki/(1.0+i+j);
	}
      }

      typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
      
      gsl_linalg_LU_decomp(gm1,gp1,&sig);
      gsl_linalg_LU_invert(gm1,gp1,gm2);

      LU_decomp(5,om1,op1,sig);
      LU_invert<ubmatrix,ubmatrix,ubmatrix_column>(5,om1,op1,om2);
      
      t.test_rel_nonzero_mat(5,5,om2,gsl_matrix_wrap(gm2),1.0e-7,1.0e-14,
			     "LU invert 1 (paren)");

    }
  }

  t.report();
  return 0;
}

