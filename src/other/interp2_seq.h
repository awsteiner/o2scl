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
/** \file interp2_seq.h
    \brief File for definition of \ref o2scl::interp2_seq
*/
#ifndef O2SCL_INTERP2_SEQ_H
#define O2SCL_INTERP2_SEQ_H

#include <iostream>
#include <string>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/interp.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Two-dimensional interpolation class by successive
      one-dimensional interpolation

      This class implements two-dimensional interpolation by iterating
      the \o2 one-dimensional interpolation routines. Derivatives and
      integrals along both x- and y-directions can be computed. This
      class is likely a bit slower than \ref interp2_seq but more
      flexible.

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

      There is an example for the usage of this class given
      in <tt>examples/ex_interp2_seq.cpp</tt>.

      \future Implement an improved caching system in case, for example
      \c xfirst is true and the last interpolation used the same
      value of \c x.

      \comment

      1. One could generalize the matrix and vector types. This would
      demand adding template types for the row and column objects.

      2. One could copy the user data instead of storing pointers to
      it.

      3. One could include get() and set() methods, but one would
      have to be careful about when to reset the interpolation.

      \endcomment

  */
  class interp2_seq {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

    interp2_seq();

    virtual ~interp2_seq();
    
    /** \brief Initialize the data for the 2-dimensional interpolation
	
	If \c x_first is true, then set_data() creates interpolation
	objects for each of the rows. Calls to interp() then uses
	these to create a column at the specified value of \c x. An
	interpolation object is created at this column to find the
	value of the function at the specified value \c y. If \c
	x_first is false, the opposite strategy is employed. These two
	options may give slightly different results.
    */
    void set_data(size_t n_x, size_t n_y, ubvector &x_grid,
		 ubvector &y_grid, ubmatrix &data, 
		 size_t interp_type=itp_cspline, bool x_first=true);
    
    /** \brief Reset the stored interpolation since the data has changed
	
        This will throw an exception if the set_data() has not been
        called.
    */
    void reset_interp();

    /** \brief Perform the 2-d interpolation 
     */
    double eval(double x, double y) const;

    /** \brief Perform the 2-d interpolation 
     */
    double operator()(double x, double y) const;

    /** \brief Compute the partial derivative in the x-direction
     */
    double deriv_x(double x, double y) const;

    /** \brief Compute the partial second derivative in the x-direction
     */
    double deriv_xx(double x, double y) const;

    /** \brief Compute the integral in the x-direction between x=x0
	and x=x1
    */
    double integ_x(double x0, double x1, double y) const;

    /** \brief Compute the partial derivative in the y-direction
     */
    double deriv_y(double x, double y) const;

    /** \brief Compute the partial second derivative in the y-direction
     */
    double deriv_yy(double x, double y) const;

    /** \brief Compute the integral in the y-direction between y=y0
	and y=y1
    */
    double integ_y(double x, double y0, double y1) const;

    /** \brief Compute the mixed partial derivative 
	\f$ \frac{\partial^2 f}{\partial x \partial y} \f$
    */
    double deriv_xy(double x, double y) const;

    /** \brief Compute a general interpolation result

	This computes
	\f[
	\frac{\partial^m}{\partial x^m}
	\frac{\partial^n}{\partial y^n} f(x,y) 
	\f]
	for \f$ m \in (-1,0,1,2) \f$ and \f$ n \in (-1,0,1,2) \f$ with
	the notation
	\f{eqnarray*}
	\frac{\partial^{-1}}{\partial x^{-1}} 
	&\equiv & \int_{x_0}^{x_1} f~dx \nonumber \\
	\frac{\partial^0}{\partial x^0} &\equiv & 
	\left.f\right|_{x=x_0} \nonumber \\
	\frac{\partial^1}{\partial x^1} &\equiv &
	\left(\frac{\partial f}{\partial x}\right)_{x=x_0} \nonumber \\
	\frac{\partial^2}{\partial x^2} &\equiv &
	\left(\frac{\partial^2 f}{\partial x^2}\right)_{x=x_0}
	\f}
	and the value of \f$ x_1 \f$ is ignored when \f$ m \geq 0 \f$
	and the value of \f$ y_1 \f$ is ignored when \f$ n \geq 0 \f$.

    */
    double eval_gen(int m, int n, double x0, double x1, 
		    double y0, double y1) const;

#ifndef DOXYGEN_NO_O2NS
    
  protected:

    /// The array of interpolation objects
    std::vector<interp_vec<ubvector> *> itps;

    /// An array of rows or columns
    std::vector<ubvector> vecs;

    /// The number of x grid points
    size_t nx;

    /// The number of y grid points
    size_t ny;

    /// True if the x interpolation should be done first
    bool xfirst;

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

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



