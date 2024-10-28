/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2024, Andrew W. Steiner

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
#include <o2scl/qr.h>

void blank_func3() {
  return;
}

#ifdef O2SCL_SET_ARMA
  
// (Armadillo specialization)
template<>
void o2scl_linalg::QR_decomp_unpack<arma::mat,arma::mat,arma::mat,double>
(const size_t M, const size_t N,arma::mat &A, arma::mat &Q, arma::mat &R) {
  arma::qr_econ(Q,R,A);
  return;
}
  
#endif
#ifdef O2SCL_SET_EIGEN
  
// (Eigen specialization)
template<>
void o2scl_linalg::QR_decomp_unpack<Eigen::MatrixXd,Eigen::MatrixXd,
				    Eigen::MatrixXd,double>
(const size_t M, const size_t N, Eigen::MatrixXd &A, 
 Eigen::MatrixXd &Q, Eigen::MatrixXd &R) {
 
  Eigen::HouseholderQR<Eigen::MatrixXd> hqr(A);
  Q=hqr.householderQ();
  R=hqr.matrixQR().triangularView<Eigen::Upper>();
  return;
}
  
#endif

