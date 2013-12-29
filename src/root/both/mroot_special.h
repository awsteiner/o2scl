/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_MROOT_SPECIAL_H
#define O2SCL_MROOT_SPECIAL_H

#include <armadillo>
#include <Eigen/Dense>

namespace o2scl {
  /** \brief A version of \ref mroot_hybrids which uses 
      Armadillo for the QR decomposition
      
      \note This class template is only defined if Armadillo
      was enabled when \o2 was installed.
  */
  template<class func_t, class vec_t, class mat_t, class jfunc_t>
    class mroot_hybrids_arma_qr_econ :
  public mroot_hybrids<func_t,vec_t,mat_t,jfunc_t> {
    
    virtual void qr_decomp_unpack() {
      arma::qr_econ(this->q,this->r,this->J);
      return;
    }
    
  };
  /** \brief A version of \ref mroot_hybrids
      which uses Eigen for the QR decomposition

      \note This class template is only defined if Eigen 
      was enabled when \o2 was installed.
  */
  template<class func_t, class vec_t, class mat_t, class jfunc_t>
    class mroot_hybrids_eigen :
  public mroot_hybrids<func_t,vec_t,mat_t,jfunc_t> {
  
    virtual void qr_decomp_unpack() {
      Eigen::HouseholderQR<Eigen::MatrixXd> hqr(this->J);
      this->q=hqr.householderQ();
      this->r=hqr.matrixQR().triangularView<Eigen::Upper>();
      return;
    }
  
  };
}

#endif
