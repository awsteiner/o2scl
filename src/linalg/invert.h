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
  template<class mat_t=boost::numeric::ublas::matrix<double>,
           class fp_t=double>
  class matrix_invert_det {
    
  public:

    /// If true, call the error handler on failure (default true)
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
                           fp_t &A_det)=0;

    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual fp_t det(size_t n, const mat_t &A)=0;
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        modifying the original matrix A 
    */
    virtual int invert_det(size_t n, mat_t &A, mat_t &A_inv) {
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
              boost::numeric::ublas::matrix<double> >,
            class fp_t=double> 
  class matrix_invert_det_LU : public matrix_invert_det<mat_t,fp_t> {
    
  public:

    /** \brief Invert matrix \c A, destroying \c A and
        returning the inverse in \c A_inv
    */
    virtual int invert_dest(size_t n, mat_t &A, mat_t &A_inv) {
      int sig;
      o2scl::permutation p(n);
      LU_decomp<mat_t,fp_t>(n,A,p,sig);
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
    virtual fp_t det(size_t n, const mat_t &A) {
      mat_t A2=A;
      int sig;
      o2scl::permutation p(n);
      LU_decomp<mat_t,fp_t>(n,A2,p,sig);
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
                            fp_t &A_det) {
      mat_t A2=A;
      int sig;
      o2scl::permutation p(n);
      LU_decomp<mat_t,fp_t>(n,A2,p,sig);
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
  template <class mat_t=boost::numeric::ublas::matrix<double>,
            class fp_t=double> 
  class matrix_invert_det_cholesky : public matrix_invert_det<mat_t,fp_t> {
    
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
                            fp_t &A_det) {
      A_inv=A;
      cholesky_decomp(n,A_inv,false);
      fp_t sqrt_det=cholesky_det(n,A_inv);
      A_det=sqrt_det*sqrt_det;
      cholesky_invert(n,A_inv);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual fp_t det(size_t n, const mat_t &A) {
      mat_t A_copy=A;
      cholesky_decomp(n,A_copy,false);
      fp_t sqrt_det=cholesky_det(n,A_copy);
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
#if defined (O2SCL_SET_ARMA) || defined (DOXYGEN)

#include <armadillo>
namespace o2scl_linalg {
  
  /** \brief Armadillo inverse 

      This class is only defined if Armadillo support was enabled
      during installation
  */
  template<class arma_mat_t=arma::mat> class matrix_invert_det_arma : 
    public matrix_invert_det<arma_mat_t,double> {

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
      A_det=arma::det(A);
      A_inv=inv(A);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const arma_mat_t &A) {
      return arma::det(A);
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
  template<class arma_mat_t=arma::mat> class matrix_invert_det_sympd_arma : 
    public matrix_invert_det<arma_mat_t,double> {

  public:
    
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
      A_det=arma::det(A);
      A_inv=inv_sympd(A);
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const arma_mat_t &A) {
      return arma::det(A);
    }
    
    /// Inver matrix \c A in place
    virtual int invert_inplace(size_t n, arma_mat_t &A) {
      A=inv_sympd(A);
      return 0;
    }

    virtual ~matrix_invert_det_sympd_arma() {}
    
  };
}

// End of #if defined (O2SCL_SET_ARMA) || defined (DOXYGEN)
#endif

#if defined (O2SCL_SET_EIGEN) || defined (DOXYGEN)

#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  
  /** \brief Eigen generic inverse and determinant

      AWS, 7/21/24: I believe this uses Eigen's PartialPivLU method
      except for smaller matrices where it uses faster explicit
      calculations.
  */
  template<class eigen_mat_t=Eigen::Matrix
           <double,Eigen::Dynamic,Eigen::Dynamic>,
           class fp_t=double> class matrix_invert_det_eigen : 
    public matrix_invert_det<eigen_mat_t,fp_t> {
    
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

  /** \brief Eigen generic inverse and determinant using a
      user-specified decomposition
      
      AWS, 10/26/24: Currently, only Eigen's PartialPivLU
      and FullPivLU methods are supported.
   */
  template<class eigen_mat_t, class eigen_decomp_t>
  class matrix_invert_det_eigen_decomp : 
    public matrix_invert_det<eigen_mat_t,double> {

  protected:

    /// Pointer to see if matrix changes
    eigen_mat_t *Aptr;

    /// The decomposition object
    std::shared_ptr<eigen_decomp_t> decomp_sp;
    
  public:

    matrix_invert_det_eigen_decomp() {
      Aptr=0;
    }
    
    /// Perform the decomposition
    void decomp(eigen_mat_t &A) {
      Aptr=&A;
      std::shared_ptr<eigen_decomp_t> decomp_loc=
        new std::shared_ptr<eigen_decomp_t>(A);
      decomp_sp=decomp_loc;
      return;
    }
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert_dest(size_t n, eigen_mat_t &A, eigen_mat_t &A_inv) {
      if (&A!=Aptr) {
        Aptr=&A;
        std::shared_ptr<eigen_decomp_t> decomp_loc(new eigen_decomp_t(A));
        decomp_sp=decomp_loc;
      }
      A_inv=decomp_sp->inverse();
      return 0;
    }
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const eigen_mat_t &A, eigen_mat_t &A_inv) {
      eigen_mat_t A2=A;
      invert_dest(n,A2,A_inv);
      return 0;
    }
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const eigen_mat_t &A,
                           eigen_mat_t &A_inv, double &A_det) {
      if (&A!=Aptr) {
        eigen_mat_t A2=A;
        std::shared_ptr<eigen_decomp_t> decomp_loc(new eigen_decomp_t(A2));
        Aptr=0;
        decomp_sp=decomp_loc;
      }
      A_inv=decomp_sp->inverse();
      A_det=decomp_sp->determinant();
      return 0;
    }
    
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const eigen_mat_t &A) {
      if (&A!=Aptr) {
        eigen_mat_t A2=A;
        std::shared_ptr<eigen_decomp_t> decomp_loc(new eigen_decomp_t(A2));
        Aptr=0;
        decomp_sp=decomp_loc;
      }
      return decomp_sp->determinant();
    }
    
    /// Inver matrix \c A in place
    virtual int invert_inplace(size_t n, eigen_mat_t &A) {
      if (&A!=Aptr) {
        Aptr=&A;
        std::shared_ptr<eigen_decomp_t> decomp_loc(new eigen_decomp_t(A));
        decomp_sp=decomp_loc;
      }
      A=decomp_sp->inverse();
      return 0;
    }
    
  };
  
  // End of #if defined (O2SCL_SET_EIGEN) || defined (DOXYGEN)
#endif

#ifdef O2SCL_NEVER_DEFINED
  
  class matrix_invert_cholesky_fast :
    public matrix_invert_det<class mat_t=o2scl::tensor2<>,
                             class fp_t=double> {

  protected:

    /// Native O2scl Cholesky decomposition (slowest)
    matrix_invert_det_cholesky<o2scl::tensor2<>,double> o2;

#ifdef O2SCL_SET_ARMA    
    /// Cholesky decomposition from Armadillo
    matrix_invert_det_sympd_arma<> arma;
#endif

#ifdef O2SCL_SET_CUDA
    /// GPU-based Cholesky decomposition
    matrix_invert_det_cholesky_cuda cuda;
#endif

  public:

    /// The size over which to prefer Cuda over Armadillo
    int n_cuda_arma;
    
    /// The size over which to prefer Cuda over O2scl's native
    int n_cuda_o2;

    /// The mode
    //@{
    static const int fast=0;
    static const int force_o2=1;
    static const int force_arma=2;
    static const int force_cuda=3;
    int mode;
    //@}

    matrix_invert_choleksy_fast() {
      n_cuda_arma=400;
      n_cuda_o2=400;
      mode=fast;
    }
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const o2scl::tensor2<> &A,
                       o2scl::tensor2<> &A_inv) {

#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

      // Both cuda and Armadillo
      int ret;
      if (mode==force_cuda || n>n_cuda_arma) {
        vector<double> vd_inv(n*n);
        ret=cuda.invert(n,A.get_data(),vd_inv);
        A_inv.swap_data(vd_inv);
      } else if (mode==force_o2) {
        ret=o2.invert(n,A,A_inv);
      } else {

        double *Ap=&A.get(0,0);
        arma::mat am(Ap,n,n,false);
        double *Ap_inv=&A_inv.get(0,0);
        arma::mat am_inv(Ap_inv,n,n,false);
        
        ret=arma.invert(n,am,am_inv);
        
      }
      
#else

      // Cuda only
      int ret;
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda || n>n_cuda_o2) {
        vector<double> vd_inv(n*n);
        ret=cuda.invert(n,A,A_inv);
        A_inv.swap_data(vd_inv);
      } else {
        ret=o2.invert(n,A,A_inv);
      }
      
#endif
#else
#ifdef O2SCL_SET_ARMA

      // Armadillo only
      
      double *Ap=&A.get(0,0);
      arma::mat am(Ap,n,n,false);
      double *Ap_inv=&A_inv.get(0,0);
      arma::mat am_inv(Ap_inv,n,n,false);
      
      ret=arma.invert(n,am,am_inv);
        
#else

      // Neither cuda norm Armadillo
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda) {
        O2SCL_ERR("Mode is force_cuda but O2SCL_SET_CUDA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      ret=o2.invert(n,A,A_inv);
      
#endif
#endif
      
      return ret;
    }
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const o2scl::tensor2<> &A,
                           o2scl::tensor2<> &A_inv, double &A_det) {
      
#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

      // Both cuda and Armadillo
      int ret;
      if (mode==force_cuda || n>n_cuda_arma) {
        ret=cuda.invert_det(n,A,A_inv,A_det);
      } else if (mode==force_o2) {
        ret=o2.invert_det(n,A,A_inv,A_det);
      } else {
        ret=arma.invert_det(n,A,A_inv,A_det);
      }
      
#else

      // Cuda only
      int ret;
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda || n>n_cuda_o2) {
        ret=cuda.invert_det(n,A,A_inv,A_det);
      } else {
        ret=o2.invert_det(n,A,A_inv,A_det);
      }
      
#endif
#else
#ifdef O2SCL_SET_ARMA

      // Armadillo only
      ret=arma.invert_det(n,A,A_inv,A_det);
      
#else

      // Neither cuda norm Armadillo
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda) {
        O2SCL_ERR("Mode is force_cuda but O2SCL_SET_CUDA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      ret=o2.invert_det(n,A,A_inv,A_det);
      
#endif
#endif
      
      return ret;
    }

    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const o2scl::tensor2<> &A) {

#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

      // Both cuda and Armadillo
      int ret;
      if (mode==force_cuda || n>n_cuda_arma) {
        ret=cuda.det(n,A);
      } else if (mode==force_o2) {
        ret=o2.det(n,A);
      } else {
        ret=arma.det(n,A);
      }
      
#else

      // Cuda only
      int ret;
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda || n>n_cuda_o2) {
        ret=cuda.det(n,A);
      } else {
        ret=o2.det(n,A);
      }
      
#endif
#else
#ifdef O2SCL_SET_ARMA

      // Armadillo only
      ret=arma.det(n,A);
      
#else

      // Neither cuda norm Armadillo
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda) {
        O2SCL_ERR("Mode is force_cuda but O2SCL_SET_CUDA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      ret=o2.det(n,A);
      
#endif
#endif
      
      return ret;
    }      
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, o2scl::tensor2<> &A) {

#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

      // Both cuda and Armadillo
      int ret;
      if (mode==force_cuda || n>n_cuda_arma) {
        ret=cuda.invert_inplace(n,A);
      } else if (mode==force_o2) {
        ret=o2.invert_inplace(n,A);
      } else {
        ret=arma.invert_inplace(n,A);
      }
      
#else

      // Cuda only
      int ret;
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda || n>n_cuda_o2) {
        ret=cuda.invert_inplace(n,A);
      } else {
        ret=o2.invert_inplace(n,A);
      }
      
#endif
#else
#ifdef O2SCL_SET_ARMA

      // Armadillo only
      ret=arma.invert_inplace(n,A);
      
#else

      // Neither cuda norm Armadillo
      if (force_arma) {
        O2SCL_ERR("Mode is force_arma but O2SCL_SET_ARMA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      if (force_cuda) {
        O2SCL_ERR("Mode is force_cuda but O2SCL_SET_CUDA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      ret=o2.invert_inplace(n,A);
      
#endif
#endif
      
      return ret;
    }
    
  };
  
#endif
  
}
  
#else

#include <o2scl/invert_special.h>

// End of #if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)
#endif

// End of #ifndef O2SCL_MATRIX_INVERT_H
#endif
