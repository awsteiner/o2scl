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
  t.set_output_level(2);

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
    
    /*
      cout << "A: " << endl;
      cout << gm1 << endl;
      cout << "b: " << endl;
      cout << gv1 << endl;
    */

    // -------------------------------------------------

    linear_solver_LU<ubvector,ubmatrix,double> lus;
    lus.solve(5,gm1,gv1,gv2);

    //cout << "x: " << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    //cout << "A*x: " << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"linear_solver_LU");
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

    /*
      cout << "A: " << endl;
      cout << gm1 << endl;
      cout << "b: " << endl;
      cout << gv1 << endl;
    */

    // -------------------------------------------------

    linear_solver_QR<ubvector,ubmatrix> lus;
    lus.solve(5,gm1,gv1,gv2);

    //cout << "x: " << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    //cout << "A*x: " << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"linear_solver_QR");
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

    /*
    cout << "A: " << endl;
    cout << gm1 << endl;
    cout << "b: " << endl;
    cout << gv1 << endl;
    */

    // -------------------------------------------------

    linear_solver_HH<ubvector,ubmatrix> lus;
    lus.solve(5,gm1,gv1,gv2);

    //cout << "x: " << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,gm2,gv2,0.0,gv3);

    //cout << "A*x: " << endl;

    t.test_rel_vec(5,gv1,gv3,1.0e-10,"linear_solver_HH");
  }

#ifdef O2SCL_SET_ARMA

  arma::mat am1(5,5), am2(5,5);
  arma::vec av1(5), av2(5), av3(5);
  
  {

    for(size_t i=0;i<5;i++) {
      av1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	am1(i,j)=1.0/(i+j+1);
	am2(i,j)=1.0/(i+j+1);
      }
    }

    // -------------------------------------------------

    linear_solver_arma<arma::vec,arma::mat> lsa;
    lsa.solve(5,am1,av1,av2);

    //cout << "x: " << endl;

    dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,am2,av2,0.0,av3);

    //cout << "A*x: " << endl;

    t.test_rel_vec(5,av1,av3,1.0e-10,"linear_solver_arma");
  }

#endif
  
#ifdef O2SCL_SET_EIGEN

  for(size_t kk=0;kk<7;kk++) {
  
    Eigen::MatrixXd em1(5,5), em2(5,5);
    Eigen::VectorXd ev1(5), ev2(5), ev3(5);
    
    // Test 
    {
      
      for(size_t i=0;i<5;i++) {
        ev1[i]=cos(((double)(i)));
        for(size_t j=0;j<5;j++) {
          em1(i,j)=1.0/(i+j+1);
          em2(i,j)=1.0/(i+j+1);
        }
      }
      
      // -------------------------------------------------

      if (kk==0) {
        linear_solver_eigen_houseQR<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==1) {
        linear_solver_eigen_colQR<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==2) {
        linear_solver_eigen_fullQR<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==3) {
        linear_solver_eigen_partLU<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==4) {
        linear_solver_eigen_fullLU<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==5) {
        linear_solver_eigen_LLT<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      } else if (kk==6) {
        linear_solver_eigen_LDLT<Eigen::VectorXd,Eigen::MatrixXd> ls;
        ls.solve(5,em1,ev1,ev2);
      }
      
      dgemv(o2cblas_RowMajor,o2cblas_NoTrans,5,5,1.0,em2,ev2,0.0,ev3);

      if (kk==0) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_houseQR");
      } else if (kk==1) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_colQR");
      } else if (kk==2) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_fullQR");
      } else if (kk==3) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_partLU");
      } else if (kk==4) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_fullLU");
      } else if (kk==5) {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_LLT");
      } else {
        t.test_rel_vec(5,ev1,ev3,1.0e-10,"linear_solver_eigen_LDLT");
      }
    }
    
  }

#endif
  
  t.report();
  return 0;
}

