/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#ifndef O2SCL_INVERT_SPECIAL_H
#define O2SCL_INVERT_SPECIAL_H

#if !defined (O2SCL_COND_FLAG) || defined (O2SCL_ARMA)
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

#endif

#endif
