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
/** \file solve_cuda.h
    \brief File for CUDA solver
*/
#ifndef O2SCL_SOLVE_CUDA_H
#define O2SCL_SOLVE_CUDA_H

#include <vector>

namespace o2scl_linalg {
  
  /** \brief Use CUDA to invert a symmetric positive definite matrix
      stored as a <tt>std::vector</tt> on the GPU
  */
  class linear_solver_LU_cuda_base {
    
  public:

    /** \brief Solve square linear system \f$ A x = b \f$ of size \c n

        \note We name this something other than \c solve() to ensure
        no confusion because of multiple inheritance.
     */
    int solve_base(int n, const std::vector<double> &A,
                   const std::vector<double> &B,
                   std::vector<double> &x);
  
  };

}

#endif

