/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#include <o2scl/vector.h>

using namespace std;
using namespace o2scl;

void blank_func1() {
  return;
}

#ifdef O2SCL_ARMA

template<> arma::subview_row<double> 
o2scl::matrix_row<arma::mat,arma::subview_row<double> >
(arma::mat &M, size_t row) {
  return M.row(row);
}

template<> arma::subview_col<double> 
o2scl::matrix_column<arma::mat,arma::subview_col<double> >
(arma::mat &M, size_t column) {
  return M.col(column);
}

double o2scl::matrix_max(const arma::mat &data) {
  size_t m=data.n_rows;
  size_t n=data.n_cols;
  if (n==0 || m==0) {
    O2SCL_ERR("Sent size=0 to matrix_max().",exc_efailed);
  }
  double max=data(0,0);
  for(size_t i=0;i<n;i++) {
    for(size_t j=0;j<m;j++) {
      if (data(i,j)>max) {
	max=data(i,j);
      }
    }
  }
  return max;
}

double o2scl::matrix_min(const arma::mat &data) {
  size_t m=data.n_rows;
  size_t n=data.n_cols;
  if (n==0 || m==0) {
    O2SCL_ERR("Sent size=0 to matrix_min().",exc_efailed);
  }
  double min=data(0,0);
  for(size_t i=0;i<n;i++) {
    for(size_t j=0;j<m;j++) {
      if (data(i,j)<min) {
	min=data(i,j);
      }
    }
  }
  return min;
}

#endif

#ifdef O2SCL_EIGEN

template<> Eigen::MatrixXd::RowXpr 
o2scl::matrix_row<Eigen::MatrixXd,Eigen::MatrixXd::RowXpr>
(Eigen::MatrixXd &M, size_t row) {
  return M.row(row);
}

template<> Eigen::MatrixXd::ColXpr 
o2scl::matrix_column<Eigen::MatrixXd,Eigen::MatrixXd::ColXpr>
(Eigen::MatrixXd &M, size_t column) {
  return M.col(column);
}

double o2scl::matrix_max(const Eigen::MatrixXd &data) {
  size_t m=data.rows();
  size_t n=data.cols();
  if (n==0 || m==0) {
    O2SCL_ERR("Sent size=0 to matrix_max().",exc_efailed);
  }
  double max=data(0,0);
  for(size_t i=0;i<n;i++) {
    for(size_t j=0;j<m;j++) {
      if (data(i,j)>max) {
	max=data(i,j);
      }
    }
  }
  return max;
}

double o2scl::matrix_min(const Eigen::MatrixXd &data) {
  size_t m=data.rows();
  size_t n=data.cols();
  if (n==0 || m==0) {
    O2SCL_ERR("Sent size=0 to matrix_min().",exc_efailed);
  }
  double min=data(0,0);
  for(size_t i=0;i<n;i++) {
    for(size_t j=0;j<m;j++) {
      if (data(i,j)<min) {
	min=data(i,j);
      }
    }
  }
  return min;
}

#endif
