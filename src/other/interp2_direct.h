/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Andrew W. Steiner
  
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
/** \file interp2_direct.h
    \brief File for definition of \ref o2scl::interp2_direct
*/
#ifndef O2SCL_INTERP2_DIRECT_H
#define O2SCL_INTERP2_DIRECT_H

#include <iostream>
#include <string>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/interp.h>
#include <o2scl/search_vec.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Bilinear or bicubic two-dimensional interpolation

      This class implements two-dimensional interpolation. First and
      second derivatives along both x- and y-directions can be
      computed. This class is likely a bit faster than \ref
      o2scl::interp2_seq but less flexible.

      The convention used by this class is that the first (row) index
      of the matrix enumerates the x coordinate and that the second
      (column) index enumerates the y coordinate. See the discussion
      in the User's guide in the section called \ref rowcol_subsect.

      The function set_data() does not copy the data, it stores
      pointers to the data. If the data is modified, then the function
      reset_interp() must be called to reset the interpolation
      information with the original pointer information. The storage
      for the data, including the arrays \c x_grid and \c y_grid are
      all managed by the user.

      By default, cubic spline interpolation with natural boundary
      conditions is used. This can be changed by calling set_interp()
      again with the same data and the new interpolation type.
      Only cubic spline and linear interpolation are supported.

      Based on D. Zaslavsky's routines at
      https://github.com/diazona/interp2d (licensed under GPLv3).
   */
  class interp2_direct {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

    interp2_direct();

    /** \brief Initialize the data for the 2-dimensional interpolation
     */
    void set_data(size_t n_x, size_t n_y, ubvector &x_grid,
		  ubvector &y_grid, ubmatrix &data, 
		  size_t interp_type=itp_cspline);
    
    /** \brief Perform the 2-d interpolation 
     */
    double eval(double x, double y);

    /** \brief Compute the partial derivative in the x-direction
     */
    double deriv_x(double x, double y);

    /** \brief Compute the partial second derivative in the x-direction
     */
    double deriv_xx(double x, double y);

    /** \brief Compute the partial derivative in the y-direction
     */
    double deriv_y(double x, double y);

    /** \brief Compute the partial second derivative in the y-direction
     */
    double deriv_yy(double x, double y);

    /** \brief Compute the mixed partial derivative 
	\f$ \frac{\partial^2 f}{\partial x \partial y} \f$
    */
    double deriv_xy(double x, double y);

#ifndef DOXYGEN_NO_O2NS
    
  protected:

    /// The number of x grid points
    size_t nx;

    /// The number of y grid points
    size_t ny;

    /// True if the data has been specified by the user
    bool data_set;

    /// The x grid
    ubvector *xfun;

    /// The y grid
    ubvector *yfun;

    /// The data
    ubmatrix *datap;

    /// Interpolation type
    size_t itype;

    /// Partial derivative with respect to x
    ubmatrix zx;

    /// Partial derivative with respect to y
    ubmatrix zy;

    /// Mixed partial derivative
    ubmatrix zxy;

    /// Searching object for x-direction
    search_vec<ubvector> svx;

    /// Searching object for y-direction
    search_vec<ubvector> svy;
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



