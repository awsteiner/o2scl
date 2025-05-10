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
/** \file invert_cuda.h
    \brief File for CUDA solver
*/
#ifndef O2SCL_INVERT_CUDA_H
#define O2SCL_INVERT_CUDA_H

#include <iostream>
#include <vector>

namespace o2scl_linalg {

  /** \brief Use CUDA to invert a symmetric positive definite matrix
      stored as a <tt>std::vector</tt> on the GPU
  */
  class matrix_invert_det_cholesky_cuda {
  
  public:
  
    /// Invert matrix \c A, returning the inverse in \c A_inv
    int invert(size_t n, const std::vector<double> &A,
               std::vector<double> &A_inv);
  
    /** \brief Invert matrix \c A, returning the inverse in \c A_inv, 
        and the determinant in \c A_det
    */
    int invert_det(size_t n, const std::vector<double> &A,
                   std::vector<double> &A_inv, double &A_det);
  
    /** \brief Determine the determinant of the matrix \c A without
        inverting
    */
    double det(size_t n, const std::vector<double> &A);
  
    /// Invert matrix \c A in place
    int invert_inplace(size_t n, std::vector<double> &A);
  
  };

}

#endif

