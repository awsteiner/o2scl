/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_CHOLESKY_SPECIAL_H
#define O2SCL_CHOLESKY_SPECIAL_H

#include <eigen3/Eigen/Dense>

namespace o2scl_linalg {
  
  // (Eigen specialization)
  template<>
    int cholesky_decomp<Eigen::MatrixXd>
    (const size_t M, Eigen::MatrixXd &A, bool err_on_fail);
     
}

#endif
