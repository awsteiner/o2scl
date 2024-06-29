/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2024, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERPM_BASE_H
#define O2SCL_INTERPM_BASE_H

/** \file interpm_base.h
    \brief File defining \ref o2scl::interpm_base,
*/

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/table.h>

namespace o2scl {

  typedef boost::numeric::ublas::vector<double> ubvector;

  /** \brief Desc
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::matrix_view_table<>,
    class mat_y_t=o2scl::matrix_view_table_transpose<> >
    class interpm_base {
      
  public:
  
  /** \brief Desc
   */
  int set_data(size_t n_in, size_t n_out, size_t n_points,
               mat_x_t &user_x, mat_y_t &user_y, 
               bool rescale=false, bool err_on_fail=true)=0;
  
    };
  
}

  
