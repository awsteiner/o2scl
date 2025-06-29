/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2025, Andrew W. Steiner

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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <o2scl/invert.h>
#include <o2scl/test_mgr.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/columnify.h>
#include <o2scl/set_cuda.h>
#include <o2scl/invert_cuda.h>
#include <o2scl/invert_auto.h>
#include <o2scl/tensor.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::vector<long double> ubvector_ld;
typedef boost::numeric::ublas::matrix<long double> ubmatrix_ld;

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

#ifdef O2SCL_SET_CUDA
  
  if (true) {
    matrix_invert_det_cholesky_cuda midcc;
    
    vector<double> A={
      4.0,1.0,1.0,
      1.0,3.0,0.0,
      1.0,0.0,2.0
    };

    vector<double> A2(9);
    midcc.invert(3,A,A2);

    tensor2<> t2(3,3);
    t2.swap_data(A2);
    matrix_out(cout,t2);
  }
  
#endif

  // Create a 5x5 identity matrix for testing
  ubmatrix id(5,5);
  for(size_t i=0;i<5;i++) {
    for(size_t j=0;j<5;j++) {
      if (i==j) id(i,j)=1.0;
      else id(i,j)=0.0;
    }
  }
  
  ubmatrix gm1(5,5), gm2(5,5), gm3(5,5);
    
  {

    cout << "Class matrix_invert_det_LU:" << endl;

    // We choose a nearly diagonal positive symmetric matrix which
    // is easy to invert
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) gm1(i,j)=((double)(i+2));
        else gm1(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
      }
    }
    
    matrix_invert_det_LU<> mi;
    mi.invert(5,gm1,gm2);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,gm1,gm2,0.0,gm3);

    matrix_out(cout,5,5,gm2);
    cout << endl;
    matrix_out(cout,5,5,gm3);
    cout << endl;

  }

  t.test_abs_mat<ubmatrix,ubmatrix,double>(5,5,id,gm3,1.0e-3,
                                           "LU inverse");
  
  ubmatrix gm4(5,5), gm5(5,5), gm6(5,5);
    
  {
    cout << "Class matrix_invert_det_cholesky:" << endl;

    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) gm4(i,j)=((double)(i+2));
        else gm4(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
        gm5(i,j)=0.0;
        gm6(i,j)=0.0;
      }
    }
    
    matrix_invert_det_cholesky<> mi;
    mi.invert(5,gm4,gm5);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,gm4,gm5,0.0,gm6);

    matrix_out(cout,5,5,gm5);
    cout << endl;
    matrix_out(cout,5,5,gm6);
    cout << endl;

  }

  t.test_abs_mat<ubmatrix,ubmatrix,double>(5,5,id,gm6,1.0e-12,
                                           "Cholesky inverse");
  t.test_abs_mat<ubmatrix,ubmatrix,double>(5,5,gm2,gm5,1.0e-12,
                                           "LU vs. Cholesky");

  cout << "Test LU inverse with multiprecision:" << endl;
  {

    cout << "Class matrix_invert_det_LU:" << endl;

    // We choose a nearly diagonal positive symmetric matrix which
    // is easy to invert
    ubmatrix um(5,5), umi(5,5);
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	um(i,j)=1.0/(i+j+1);
      }
    }
    
    matrix_invert_det_LU<> mi;
    mi.invert(5,um,umi);

    matrix_out(cout,5,5,umi);
    cout << dtos(umi(1,0),0) << endl;
    cout << dtos(umi(2,1),0) << endl;
    cout << endl;
  }
  {

    cout << "Class matrix_invert_det_LU (long double):" << endl;

    // We choose a nearly diagonal positive symmetric matrix which
    // is easy to invert
    ubmatrix_ld um(5,5), umi(5,5);
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	um(i,j)=1.0L/(i+j+1);
      }
    }
    
    matrix_invert_det_LU<boost::numeric::ublas::matrix<long double>,
                         boost::numeric::ublas::matrix_column<
                           boost::numeric::ublas::matrix<long double> >,
                         long double> mi;
    mi.invert(5,um,umi);

    matrix_out(cout,5,5,umi);
    cout << dtos(umi(1,0),0) << endl;
    cout << dtos(umi(2,1),0) << endl;
    cout << endl;
  }
  
  
#ifdef O2SCL_SET_ARMA
  
    arma::mat am1(5,5), am2(5,5), am3(5,5);

  {

    cout << "Class matrix_invert_det_arma: " << endl;
    
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) am1(i,j)=((double)(i+2));
        else am1(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
      }
    }
    
    matrix_invert_det_arma<arma::mat> mi;
    mi.invert(5,am1,am2);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,am1,am2,0.0,am3);

    matrix_out(cout,5,5,am2);
    cout << endl;
    matrix_out(cout,5,5,am3);
    cout << endl;

  }
  
#endif

#ifdef O2SCL_SET_EIGEN
  
  Eigen::MatrixXd em1(5,5), em2(5,5), em3(5,5);
    
  if (true) {
    
    cout << "Class matrix_invert_det_eigen: " << endl;

    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) em1(i,j)=((double)(i+2));
        else em1(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
      }
    }
    
    matrix_invert_det_eigen<Eigen::MatrixXd> mi;
    mi.invert(5,em1,em2);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,em1,em2,0.0,em3);

    matrix_out(cout,5,5,em2);
    cout << endl;
    matrix_out(cout,5,5,em3);
    cout << endl;

    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>(5,5,id,em3,1.0e-12,
                                                    "Eigen inverse");
    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>(5,5,gm5,em2,1.0e-12,
                                                    "Cholesky vs. Eigen");
    
  }

  if (true) {
    
    cout << "Class matrix_invert_det_eigen_decomp with PartialPivLU: "
         << endl;

    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) em1(i,j)=((double)(i+2));
        else em1(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
      }
    }
    
    matrix_invert_det_eigen_decomp<Eigen::MatrixXd,
                                   Eigen::PartialPivLU
                                   <Eigen::MatrixXd> > mide;
    mide.invert(5,em1,em2);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,em1,em2,0.0,em3);

    matrix_out(cout,5,5,em2);
    cout << endl;
    matrix_out(cout,5,5,em3);
    cout << endl;

    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>
      (5,5,id,em3,1.0e-12,"Eigen decomp inverse");
    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>
      (5,5,gm5,em2,1.0e-12,"Cholesky vs. Eigen decomp");
    
  }

  if (true) {
    
    cout << "Class matrix_invert_det_eigen_decomp with FullPivLU: "
         << endl;

    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
        if (i==j) em1(i,j)=((double)(i+2));
        else em1(i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
      }
    }
    
    matrix_invert_det_eigen_decomp<Eigen::MatrixXd,
                                   Eigen::FullPivLU
                                   <Eigen::MatrixXd> > mide;
    mide.invert(5,em1,em2);

    dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
          5,5,5,1.0,em1,em2,0.0,em3);

    matrix_out(cout,5,5,em2);
    cout << endl;
    matrix_out(cout,5,5,em3);
    cout << endl;
    
    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>
      (5,5,id,em3,1.0e-12,"Eigen decomp inverse");
    t.test_abs_mat<ubmatrix,Eigen::MatrixXd,double>
      (5,5,gm5,em2,1.0e-12,"Cholesky vs. Eigen decomp");
    
  }

#endif

  t.report();
  return 0;
}

