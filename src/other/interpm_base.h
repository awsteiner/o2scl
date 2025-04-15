/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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

  /** \brief Base class for multidimensional interpolation
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::const_matrix_view_table<>,
           class mat_y_t=o2scl::const_matrix_view_table<> >
  class interpm_base {

  protected:
    
    /// Number of parameters
    size_t n_params;
    
    /// Number of outputs
    size_t n_outputs;
    
    /// Number of points
    size_t n_points;
    
  public:
    
    /// If true, throw exceptions on convergence errors (default true)
    bool err_nonconv;
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    interpm_base() {
      verbose=0;
      n_params=0;
      n_outputs=0;
      n_points=0;
      err_nonconv=true;
    }

    /** \brief Set the data to be interpolated

        The variable \c user_x should be arranged so that parameters
        are indexed by the columns and data points are points are
        indexed by the rows. The variable \c user_y should be arranged
        so that output quantities are indexed by the columns and
        points are indexed by the rows. Children are allowed (but not
        required) to use \c swap to take over management of the input
        data. None of the values \c n_in, \c n_out, or \c n_pts are
        allowed to be zero. The matrix \c user_x should have \c n_pts
        rows and \c n_in columns, while the matrix \c user_y should
        have \c n_pts rows and \c n_out columns.
    */
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y)=0;

    /** \brief Add data to be interpolated

        The arguments in this function should be arranged similarly to
        the \ref set_data() function. Children may not implement this
        function (in which case it calls the error handler), or they
        may simply append the additional data set to the one
        previously specified in \ref set_data() call.
     */
    virtual int add_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y) {
      O2SCL_ERR2_RET("No add_data() function for this interpolator",
                     " in interpm_base().",o2scl::exc_eunimpl);
    }
    
    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int eval(const vec_t &x, vec_t &y) const=0;
    
    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int operator()(size_t nx, const vec_t &x,
                         size_t ny, vec_t &y) {
      if (nx!=n_params) {
        O2SCL_ERR2("Mismatch in number of parameters in ",
                   "interpm_base::operator().",o2scl::exc_einval);
      }
      if (ny!=n_outputs) {
        O2SCL_ERR2("Mismatch in number of outputs in ",
                   "interpm_base::operator().",o2scl::exc_einval);
      }
      return eval(x,y);
    }
    
    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int operator()(size_t nx, const vec_t &x,
                         vec_t &y) {
      if (nx!=n_params) {
        O2SCL_ERR2("Mismatch in number of parameters in ",
                   "interpm_base::operator().",o2scl::exc_einval);
      }
      return eval(x,y);
    }
    
    /** \brief Evaluate the interpolation at point \c x,
        returning a single value
    */
    virtual double operator()(size_t nx, const vec_t &x) {
      vec_t y;
      eval(x,y);
      return y[0];
    }
    
    /** \brief Evaluate the interpolation at point \c x,
        returning \c y and the uncertainties in \c y_unc
    */
    virtual int eval_unc(const vec_t &x, vec_t &y, vec_t &y_unc) const {
      for(size_t j=0;j<n_outputs;j++) {
        y_unc[j]=0.0;
      }
      return eval(x,y);
    }
    
  };
  
}

#endif
