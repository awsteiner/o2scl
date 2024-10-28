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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <o2scl/qr.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  bool eigen_tests=false;

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    gsl_matrix *gm1=gsl_matrix_alloc(5,5);
    gsl_matrix *gm2=gsl_matrix_alloc(5,5);
    gsl_matrix *gm3=gsl_matrix_alloc(5,5);
    gsl_vector *gv1=gsl_vector_alloc(5);
    gsl_vector *gv2=gsl_vector_alloc(5);
    gsl_vector *gv3=gsl_vector_alloc(5);
    ubmatrix om1(5,5), om2(5,5), om3(5,5);
    ubvector ov1(5), ov2(5), ov3(5);
    permutation gp1(5), op1(5), mp1(5);

    for(size_t ik=0;ik<10;ik++) {

      // Test decomp

      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv1,i,cos(((double)(i))));
	gsl_vector_set(gv2,i,cos(((double)(i))));
	ov1[i]=cos(((double)(i)));
	ov2[i]=cos(((double)(i)));
	for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,((double)(ik+1))/(1.0+i+j));
	  om1(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }

      gsl_linalg_QR_decomp(gm1,gv1);
      QR_decomp<ubmatrix,ubvector,double>(5,5,om1,ov1);
      t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),1.0e-11,"qr decomp 1");
      t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-11,"qr decomp 2");

#ifdef O2SCL_SET_EIGEN

      Eigen::MatrixXd em1(5,5);
      for(size_t i=0;i<5;i++) {
	for(size_t j=0;j<5;j++) {
	  em1(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }
      Eigen::HouseholderQR<Eigen::MatrixXd> hqr(em1);
      t.test_rel_mat(5,5,hqr.matrixQR(),om1,1.0e-11,"qr decomp 3");
      eigen_tests=true;
      
#endif
      
#ifdef O2SCL_NEVER_DEFINED
      
      // This doesn't work yet
      arma::mat am1(5,5), amq, amr;
      for(size_t i=0;i<5;i++) {
	for(size_t j=0;j<5;j++) {
	  am1(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }
      qr(amq,amr,am1);
      matrix_out(cout,amq,5,5);
      matrix_out(cout,amr,5,5);
      
#endif

      // Test QTvec

      gsl_linalg_QR_QTvec(gm1,gv1,gv2);
      QR_QTvec<ubmatrix,ubvector,ubvector,double>(5,5,om1,ov1,ov2);
      t.test_rel_vec(5,ov2,gsl_vector_wrap(gv2),1.0e-11,"qr qtvec 1");

      // Test solve

      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv1,i,cos(((double)(i))));
	ov1[i]=cos(((double)(i)));
	gsl_vector_set(gv2,i,cos(((double)(i))));
	ov2[i]=cos(((double)(i)));
	gsl_vector_set(gv3,i,cos(((double)(i))));
	ov3[i]=cos(((double)(i)));
	for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,((double)(ik+1))/(1.0+i+j));
	  om1(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }

      gsl_linalg_QR_decomp(gm1,gv2);
      gsl_linalg_QR_solve(gm1,gv2,gv1,gv3);

      QR_decomp<ubmatrix,ubvector,double>(5,5,om1,ov2);
      QR_solve<ubmatrix,ubvector,ubvector,
               ubvector,double>(5,om1,ov2,ov1,ov3);

      t.test_rel_vec(5,ov3,gsl_vector_wrap(gv3),1.0e-11,"qr solve 1");

      // Test unpack

      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv1,i,cos(((double)(i))));
	ov1[i]=cos(((double)(i)));
	gsl_vector_set(gv2,i,cos(((double)(i))));
	ov2[i]=cos(((double)(i)));
	for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,((double)(ik+1))/(1.0+i+j));
	  om1(i,j)=((double)(ik+1))/(1.0+i+j);
	  gsl_matrix_set(gm2,i,j,((double)(ik+1))/(1.0+i+j));
	  om2(i,j)=((double)(ik+1))/(1.0+i+j);
	  gsl_matrix_set(gm3,i,j,((double)(ik+1))/(1.0+i+j));
	  om3(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }
    
      gsl_linalg_QR_decomp(gm1,gv2);
      gsl_linalg_QR_unpack(gm1,gv2,gm2,gm3);

      QR_decomp<ubmatrix,ubvector,double>(5,5,om1,ov2);
      QR_unpack<ubmatrix,ubmatrix,ubmatrix,ubvector,double>
        (5,5,om1,ov2,om2,om3);

      t.test_rel_mat(5,5,om2,gsl_matrix_wrap(gm2),1.0e-11,"qr decomp 5");
      t.test_rel_mat(5,5,om3,gsl_matrix_wrap(gm3),1.0e-11,"qr decomp 6");

    }

  }

  if (eigen_tests) {
    cout << "Included tests with Eigen." << endl;
  }

  t.report();

  return 0;
}

