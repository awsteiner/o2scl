/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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
/** \file invert_fast.h
    \brief File for matrix inversion classes
*/
#ifndef O2SCL_INVERT_FAST_H
#define O2SCL_INVERT_FAST_H

#include <o2scl/invert.h>
#include <o2scl/invert_cuda.h>
#include <o2scl/tensor.h>
#include <o2scl/set_cuda.h>
#include <o2scl/vector_special.h>

namespace o2scl_linalg {

  /** \brief Desc
   */
  class matrix_invert_cholesky_fast :
    public matrix_invert_det<o2scl::tensor2<>,double> {

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
    size_t n_cuda_arma;
    
    /// The size over which to prefer Cuda over native inversion
    size_t n_cuda_o2;

    /// The size over which to prefer Armadillo over native inversion
    size_t n_arma_o2;

    /// The mode
    //@{
    static const size_t fast=0;
    static const size_t force_o2=1;
    static const size_t force_arma=2;
    static const size_t force_cuda=3;
    int mode;
    //@}

    matrix_invert_cholesky_fast() {
      n_cuda_arma=400;
      n_cuda_o2=100;
      n_arma_o2=15;
      mode=fast;
    }
    
    /// Invert matrix \c A, returning the inverse in \c A_inv
    virtual int invert(size_t n, const o2scl::tensor2<> &A,
                       o2scl::tensor2<> &A_inv) {

#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

      // Both cuda and Armadillo
      int ret;
      if (mode==force_o2) {
        
        //std::cout << "1";
        ret=o2.invert(n,A,A_inv);
        
      } else if (mode==force_arma || (mode!=force_cuda && n<n_cuda_arma)) {
        
        //std::cout << "2";
        // We have to cast away constness :(
        double *Ap=(double *)(&A.get(0,0));
        arma::mat am(Ap,n,n,false);
        double *Ap_inv=&A_inv.get(0,0);
        arma::mat am_inv(Ap_inv,n,n,false);
        
        ret=arma.invert(n,am,am_inv);
        
      } else {
        
        std::vector<double> vd_inv(n*n);
        //std::cout << "3";
        ret=cuda.invert(n,A.get_data(),vd_inv);
        A_inv.swap_data(vd_inv);

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

      if (force_cuda) {
        O2SCL_ERR("Mode is force_cuda but O2SCL_SET_CUDA is false ",
                  "in matrix_invert_fast::invert().",o2scl::exc_eunimpl);
                  
      }
      // We have to cast away constness :(
      double *Ap=(double *)(&(A.get(0,0)));
      arma::mat am(Ap,n,n,false);
      double *Ap_inv=&A_inv.get(0,0);
      arma::mat am_inv(Ap_inv,n,n,false);
      
      int ret=arma.invert(n,am,am_inv);
        
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
      int ret=o2.invert(n,A,A_inv);
      
#endif
#endif
      
      return ret;
    }
    
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    virtual int invert_det(size_t n, const o2scl::tensor2<> &A,
                           o2scl::tensor2<> &A_inv, double &A_det) {
#ifdef O2SCL_NEVER_DEFINED
      
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
#endif
      return 0;
    }

    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    virtual double det(size_t n, const o2scl::tensor2<> &A) {
#ifdef O2SCL_NEVER_DEFINED

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
#endif
      return 0;
    }      
    
    /// Invert matrix \c A in place
    virtual int invert_inplace(size_t n, o2scl::tensor2<> &A) {

#ifdef O2SCL_NEVER_DEFINED
      
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
#endif
      return 0;
    }
    
  };
  
}
  
// End of #ifndef O2SCL_INVERT_FAST_H
#endif
