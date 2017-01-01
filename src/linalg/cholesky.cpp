/*
  -------------------------------------------------------------------

  Copyright (C) 2015-2017, Andrew W. Steiner

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
#include <o2scl/cholesky.h>

void blank_func2() {
  return;
}

#ifdef O2SCL_EIGEN
  
// (Eigen specialization)
template<>
void o2scl_linalg::cholesky_decomp<Eigen::MatrixXd>
(const size_t M, Eigen::MatrixXd &A) {
  
  Eigen::LLT<Eigen::MatrixXd> llt(A);
  A=llt.matrixL();
  for(size_t i=0;i<M;i++) {
    for(size_t j=i;j<M;j++) {
      A(i,j)=A(j,i);
    }
  }
  return;
}
  
#endif

