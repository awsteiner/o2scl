/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#ifndef O2SCL_LINEAR_SPECIAL_H
#define O2SCL_LINEAR_SPECIAL_H

#include <eigen3/Eigen/Dense>

namespace o2scl_linalg {

  /** \brief Eigen linear solver using QR decomposition with 
      column pivoting

      This class is only defined if Eigen support was enabled during
      installation.
  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_houseQR : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.householderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using QR decomposition with 
      column pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_colQR : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.colPivHouseholderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using QR decomposition with 
      full pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_fullQR : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.fullPivHouseholderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LU decomposition with 
      partial pivoting
      
      This requires the matrix \c A to be invertible.

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_partLU : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.partialPivLu().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LU decomposition with 
      full pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_fullLU : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.fullPivLu().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LLT decomposition with 
      full pivoting
      
      This requires the matrix \c A to be positive definite.

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_LLT : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.llt().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LDLT decomposition with 
      full pivoting
      
      This requires the matrix \c A to be positive or negative
      semidefinite.

      This class is only defined if Eigen support was enabled during
      installation.
  */
  template<class eigen_vec_t, class eigen_mat_t>
    class linear_solver_eigen_LDLT : 
  public linear_solver<eigen_vec_t,eigen_mat_t> {
    virtual void solve(size_t n, eigen_mat_t &A, eigen_vec_t &b,
		       eigen_vec_t &x) {
      x=A.ldlt().solve(b);
      return;
    }
  };

}

#endif
