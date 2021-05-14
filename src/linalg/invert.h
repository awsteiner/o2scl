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

  /** \brief Invert a matrix
   */
  template<class mat_t=boost::numeric::ublas::matrix<double> >
  class matrix_invert {
    
  public:
    
    virtual ~matrix_invert() {}
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, mat_t &A, mat_t &Ainv)=0;
    
    /// Invert matrix \c A in place
    virtual void invert_inplace(size_t n, mat_t &A)=0;

    /// Invert a matrix, returning the inverse
    virtual mat_t invert(size_t n, mat_t &A) {
      mat_t Ainv(n,n);
      invert(n,A,Ainv);
      return Ainv;
    }
    
  };

  /** \brief Generic inverse using LU decomposition
   */
  template <class mat_t=boost::numeric::ublas::matrix<double>,
            class mat_col_t=boost::numeric::ublas::matrix_column<
              boost::numeric::ublas::matrix<double> > >
  class matrix_invert_LU : public matrix_invert<mat_t> {
    
  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, mat_t &A, mat_t &A_inv) {
      int sig;
      o2scl::permutation p(n);
      LU_decomp(n,A,p,sig);
      if (o2scl_linalg::diagonal_has_zero(n,A)) {
        O2SCL_ERR2("Matrix singular (LU method) ",
                   "in matrix_invert_LU::invert().",o2scl::exc_esing);
      }
      LU_invert<mat_t,mat_t,mat_col_t>(n,A,p,A_inv);
      return;
    };
    
    /// Invert matrix \c A in place
    virtual void invert_inplace(size_t n, mat_t &A) {
      mat_t Ainv(n,n);
      invert(n,A,Ainv);
      A=Ainv;
      return;
    }
    
    virtual ~matrix_invert_LU() {}
    
  };

  /** \brief Inverse symmetric positive matrix using Cholesky decomposition
   */
  template <class mat_t=boost::numeric::ublas::matrix<double> > 
  class matrix_invert_cholesky : public matrix_invert<mat_t> {
    
  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, mat_t &A, mat_t &A_inv) {
      A_inv=A;
      invert_inplace(n,A_inv);
      return;
    };

    /// Invert matrix \c A in place
    virtual void invert_inplace(size_t n, mat_t &A) {
      cholesky_decomp(n,A,false);
      cholesky_invert(n,A);
      return;
    }
    
    virtual ~matrix_invert_cholesky() {}
    
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
  template<class arma_mat_t> class matrix_invert_arma : 
    public matrix_invert<arma_mat_t> {
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, const arma_mat_t &A, arma_mat_t &A_inv) {
      A_inv=inv(A);
      return;
    }

    /// Invert matrix \c A in place
    virtual void invert_inplace(size_t n, arma_mat_t &A) {
      A=inv(A);
      return;
    }

    virtual ~matrix_invert_arma() {}
    
  };

  /** \brief Armadillo inverse of symmetric positive definite matrix

      This class is only defined if Armadillo support was enabled
      during installation
  */
  template<class arma_mat_t> class matrix_invert_sympd_arma : 
    public matrix_invert<arma_mat_t> {
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, const arma_mat_t &A, arma_mat_t &A_inv) {
      A_inv=inv_sympd(A);
      return;
    }

    /// Inver matrix \c A in place
    virtual void invert_inplace(size_t n, arma_mat_t &A) {
      A=inv_sympd(A);
      return;
    }

    virtual ~matrix_invert_sympd_arma() {}
    
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
  template<class eigen_mat_t> class matrix_invert_eigen : 
    public matrix_invert<eigen_mat_t> {
    
  public:
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual void invert(size_t n, eigen_mat_t &A, eigen_mat_t &A_inv) {
      A_inv=A.inverse();
      return;
    }
    
    /// Inver matrix \c A in place
    virtual void invert_inplace(size_t n, eigen_mat_t &A) {
      A=A.inverse();
      return;
    }
    
  };
  
}
#endif

#else
#include <o2scl/linear_special.h>
#endif

#endif
