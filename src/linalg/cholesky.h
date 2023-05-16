/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_CHOLESKY_H
#define O2SCL_CHOLESKY_H

/** \file cholesky.h
    \brief Header wrapper for \ref cholesky_base.h
*/

#include <o2scl/err_hnd.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/vector.h>

namespace o2scl_linalg {
  
  using namespace o2scl_cblas;

#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/cholesky_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

namespace o2scl_linalg_bracket {
  
  using namespace o2scl_cblas_bracket;
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/cholesky_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)

#if defined (O2SCL_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  
  /** \brief Eigen specialization of \ref cholesky_decomp()
   */
  template<>
    int cholesky_decomp<Eigen::MatrixXd>
    (const size_t M, Eigen::MatrixXd &A, bool err_on_fail);
     
}
#endif

#else
#include <o2scl/cholesky_special.h>
#endif

#endif
