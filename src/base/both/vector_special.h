/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#ifndef O2SCL_VECTOR_SPECIAL_H
#define O2SCL_VECTOR_SPECIAL_H

#include <armadillo>
#include <eigen3/Eigen/Dense>

#if !defined (O2SCL_COND_FLAG) || defined (O2SCL_ARMA)
#include <armadillo>
namespace o2scl_linalg {
  
  /// Armadillo version of \ref matrix_max()
  double matrix_max(const arma::mat &data);

  /// Armadillo version of \ref matrix_row()
  arma::subview_row<double> matrix_row(arma::mat &M, size_t row);

  /// Armadillo version of \ref matrix_column()
  arma::subview_col<double> matrix_column(arma::mat &M, size_t column);

}
#endif

#if !defined (O2SCL_COND_FLAG) || defined (O2SCL_EIGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {

  /// Eigen version of \ref matrix_max()
  double matrix_max(const Eigen::MatrixXd &data);
  
  /// Eigen version of \ref matrix_row()
  Eigen::MatrixXd::RowXpr matrix_row(Eigen::MatrixXd &M, size_t row);
     
  /// Eigen version of \ref matrix_column()
  Eigen::MatrixXd::ColXpr matrix_column(Eigen::MatrixXd &M, size_t column);
     
}
#endif

#endif
