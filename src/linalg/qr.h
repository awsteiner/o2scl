/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file qr.h
    \brief Header wrapper for \ref qr_base.h
*/
#ifndef O2SCL_QR_H
#define O2SCL_QR_H

#include <o2scl/err_hnd.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/householder.h>
#include <o2scl/givens.h>

namespace o2scl_linalg {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#ifndef DOXYGEN
#include <o2scl/qr_base.h>  
#endif
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

namespace o2scl_linalg_bracket {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#ifndef DOXYGEN
#include <o2scl/qr_base.h>  
#endif
#undef O2SCL_IX
#undef O2SCL_IX2

}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)

#if defined (O2SCL_ARMA) || defined (DOXYGEN)
#include <armadillo>
namespace o2scl_linalg {
  
  /** \brief Armadillo specialization of \ref QR_decomp_unpack().
   */
  template<>
    void QR_decomp_unpack<arma::mat,arma::mat,arma::mat>
    (const size_t M, const size_t N, arma::mat &A, arma::mat &Q, 
     arma::mat &R);

}
#endif

#if defined (O2SCL_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  
  /** \brief Eigen specialization of \ref QR_decomp_unpack().
   */
  template<>
    void QR_decomp_unpack<Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd>
    (const size_t M, const size_t N, Eigen::MatrixXd &A, 
     Eigen::MatrixXd &Q, Eigen::MatrixXd &R);
     
}
#endif

#else
#include <o2scl/qr_special.h>
#endif

#endif
