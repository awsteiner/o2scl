/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#ifndef O2SCL_LINEAR_SPECIAL_H
#define O2SCL_LINEAR_SPECIAL_H

#include <armadillo>

namespace o2scl_linalg {

  /** \brief Armadillo linear solver 

      This class is only defined if Armadillo support was enabled
      during installation
   */
  template<class arma_vec_t, class arma_mat_t> class linear_solver_arma : 
  public linear_solver<arma_vec_t,arma_mat_t> {
    virtual void solve(size_t n, arma_mat_t &A, arma_vec_t &b,
		       arma_vec_t &x) {
      x=arma::solve(A,b);
      return;
    }
  };

}

#endif
