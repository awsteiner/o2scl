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
#ifndef O2SCL_QR_SPECIAL_H
#define O2SCL_QR_SPECIAL_H

#if !defined (O2SCL_COND_FLAG) || defined (O2SCL_ARMA)
#include <armadillo>
namespace o2scl_linalg {
  
  // (Armadillo specialization)
  template<>
    void QR_decomp_unpack<arma::mat,arma::mat,arma::mat>
    (const size_t M, const size_t N, arma::mat &A, arma::mat &Q, 
     arma::mat &R);

}
#endif

#if !defined (O2SCL_COND_FLAG) || defined (O2SCL_EIGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  
  // (Eigen specialization)
  template<>
    void QR_decomp_unpack<Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd>
    (const size_t M, const size_t N, Eigen::MatrixXd &A, 
     Eigen::MatrixXd &Q, Eigen::MatrixXd &R);
     
}
#endif

#endif
