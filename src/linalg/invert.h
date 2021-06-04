/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
/** \file matrix_invert.h
    \brief File for linear solvers
*/
#ifndef O2SCL_MATRIX_INVERT_H
#define O2SCL_MATRIX_INVERT_H

#include <gsl/gsl_linalg.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/permutation.h>
#include <o2scl/lu.h>
#include <o2scl/cholesky.h>

namespace o2scl_linalg {

  /** \brief Invert a matrix and compute its determinant
   */
  template<class mat_t=boost::numeric::ublas::matrix<double> >
  class matrix_invert_det {
    
  public:

    bool err_on_fail;

    matrix_invert_det() {
      err_on_fail=true;
    }
    
    virtual ~matrix_invert_det() {}
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const mat_t &A, mat_t &A_inv)=0;
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const mat_t &A, mat_t &A_inv,
                            double &A_det)=0;

    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const mat_t &A)=0;
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        modifying the original matrix A 
    */
    virtual int invert_dest(size_t n, mat_t &A, mat_t &A_inv) {
      return invert(n,A,A_inv);
    }
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, mat_t &A)=0;

    /// Invert a matrix, returning the inverse
    virtual mat_t invert(size_t n, const mat_t &A) {
      mat_t A_inv(n,n);
      int ret=invert(n,A,A_inv);
      if (ret!=0) {
        O2SCL_ERR("Matrix inversion failed.",o2scl::exc_efailed);
      }
      return A_inv;
    }
    
  };

  /** \brief Generic inverse using LU decomposition
   */
  template <class mat_t=boost::numeric::ublas::matrix<double>,
            class mat_col_t=boost::numeric::ublas::matrix_column<
              boost::numeric::ublas::matrix<double> > >
  class matrix_invert_det_LU : public matrix_invert_det<mat_t> {
    
  public:

    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert_dest(size_t n, mat_t &A, mat_t &A_inv) {
      int sig;
      o2scl::permutation p(n);
      LU_decomp(n,A,p,sig);
      if (o2scl_linalg::diagonal_has_zero(n,A)) {
        O2SCL_ERR2("Matrix singular (LU method) ",
                   "in matrix_invert_det_LU::invert().",o2scl::exc_esing);
      }
      LU_invert<mat_t,mat_t,mat_col_t>(n,A,p,A_inv);
      return 0;
    };
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const mat_t &A, mat_t &A_inv) {
      mat_t A2=A;
      invert_dest(n,A2,A_inv);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const mat_t &A) {
      mat_t A2=A;
      int sig;
      o2scl::permutation p(n);
      LU_decomp(n,A2,p,sig);
      if (o2scl_linalg::diagonal_has_zero(n,A)) {
        O2SCL_ERR2("Matrix singular (LU method) ",
                   "in matrix_invert_det_LU::invert().",o2scl::exc_esing);
      }
      return LU_det(n,A2,sig);
    }
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const mat_t &A, mat_t &A_inv,
                            double &A_det) {
      mat_t A2=A;
      int sig;
      o2scl::permutation p(n);
      LU_decomp(n,A2,p,sig);
      if (o2scl_linalg::diagonal_has_zero(n,A)) {
        O2SCL_ERR2("Matrix singular (LU method) ",
                   "in matrix_invert_det_LU::invert().",o2scl::exc_esing);
      }
      A_det=LU_det(n,A2,sig);
      LU_invert<mat_t,mat_t,mat_col_t>(n,A2,p,A_inv);
      return 0;
    }
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, mat_t &A) {
      mat_t A_inv(n,n);
      invert(n,A,A_inv);
      A=A_inv;
      return 0;
    }
    
    virtual ~matrix_invert_det_LU() {}
    
  };

  /** \brief Inverse symmetric positive matrix using Cholesky decomposition
   */
  template <class mat_t=boost::numeric::ublas::matrix<double> > 
  class matrix_invert_det_cholesky : public matrix_invert_det<mat_t> {
    
  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const mat_t &A, mat_t &A_inv) {
      A_inv=A;
      invert_inplace(n,A_inv);
      return 0;
    };

    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const mat_t &A, mat_t &A_inv,
                            double &A_det) {
      A_inv=A;
      cholesky_decomp(n,A_inv,false);
      double sqrt_det=cholesky_det(n,A_inv);
      A_det=sqrt_det*sqrt_det;
      cholesky_invert(n,A_inv);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const mat_t &A) {
      mat_t A_copy=A;
      cholesky_decomp(n,A_copy,false);
      double sqrt_det=cholesky_det(n,A_copy);
      return sqrt_det*sqrt_det;
    }
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, mat_t &A) {
      cholesky_decomp(n,A,false);
      cholesky_invert(n,A);
      return 0;
    }
    
    virtual ~matrix_invert_det_cholesky() {}
    
  };
  
  // End of namespace o2scl_linalg
}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)
#if defined (O2SCL_ARMA) || defined (DOXYGEN)
#include <armadillo>
namespace o2scl_linalg {
  
  /** \brief Armadillo inverse 

      This class is only defined if Armadillo support was enabled
      during installation
  */
  template<class arma_mat_t> class matrix_invert_det_arma : 
    public matrix_invert_det<arma_mat_t> {

  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const arma_mat_t &A, arma_mat_t &A_inv) {
      A_inv=inv(A);
      return 0;
    }

    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const arma_mat_t &A, arma_mat_t &A_inv,
                            double &A_det) {
      A_det=det(A);
      A_inv=inv(A);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const arma_mat_t &A) {
      return det(A);
    }
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, arma_mat_t &A) {
      A=inv(A);
      return 0;
    }

    virtual ~matrix_invert_det_arma() {}
    
  };

  /** \brief Armadillo inverse of symmetric positive definite matrix

      This class is only defined if Armadillo support was enabled
      during installation
  */
  template<class arma_mat_t> class matrix_invert_det_sympd_arma : 
    public matrix_invert_det<arma_mat_t> {
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const arma_mat_t &A, arma_mat_t &A_inv) {
      A_inv=inv_sympd(A);
      return 0;
    }

    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const arma_mat_t &A, arma_mat_t &A_inv,
                            double &A_det) {
      A_det=det(A);
      A_inv=inv_sympd(A);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const arma_mat_t &A) {
      return det(A);
    }
    
    /// Inver matrix \c A in place
    virtual int invert_inplace(size_t n, arma_mat_t &A) {
      A=inv_sympd(A);
      return 0;
    }

    virtual ~matrix_invert_det_sympd_arma() {}
    
  };
}
#endif

#if defined (O2SCL_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  
  /** \brief Eigen inverse using QR decomposition with 
      column pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_mat_t> class matrix_invert_det_eigen : 
    public matrix_invert_det<eigen_mat_t> {
    
  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const eigen_mat_t &A, eigen_mat_t &A_inv) {
      A_inv=A.inverse();
      return 0;
    }
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const eigen_mat_t &A,
                           eigen_mat_t &A_inv, double &A_det) {
      A_inv=A.inverse();
      A_det=A.determinant();
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const eigen_mat_t &A) {
      return A.determinant();
    }
    
    /// Inver matrix \c A in place
    virtual int invert_inplace(size_t n, eigen_mat_t &A) {
      A=A.inverse();
      return 0;
    }
    
  };
  
}
#endif

#else
#include <o2scl/linear_special.h>
#endif

#endif
