/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014-2022, Andrew W. Steiner
  
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
/** \file interp2.h
    \brief File defining \ref o2scl::interp2_base
*/
#ifndef O2SCL_INTERP2_H
#define O2SCL_INTERP2_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Two-dimensional interpolation base class [abstract]
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> >
    class interp2_base {

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
  public:
    
    interp2_base() {
    }
    
    virtual ~interp2_base() {
    }
    
    /** \brief Perform the 2-d interpolation 
     */
    virtual double eval(double x, double y) const=0;
    
    /** \brief Perform the 2-d interpolation 
     */
    virtual double operator()(double x, double y) const {
      return eval(x,y);
    }
    
    /** \brief Compute the partial derivative in the x-direction
     */
    virtual double deriv_x(double x, double y) const=0;

    /** \brief Compute the partial second derivative in the x-direction
     */
    virtual double deriv_xx(double x, double y) const=0;

    /** \brief Compute the integral in the x-direction between x=x0
	and x=x1
    */
    virtual double integ_x(double x0, double x1, double y) const=0;

    /** \brief Compute the partial derivative in the y-direction
     */
    virtual double deriv_y(double x, double y) const=0;

    /** \brief Compute the partial second derivative in the y-direction
     */
    virtual double deriv_yy(double x, double y) const=0;

    /** \brief Compute the integral in the y-direction between y=y0
	and y=y1
    */
    virtual double integ_y(double x, double y0, double y1) const=0;

    /** \brief Compute the mixed partial derivative 
	\f$ \frac{\partial^2 f}{\partial x \partial y} \f$
    */
    virtual double deriv_xy(double x, double y) const=0;

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
    virtual double eval_gen(int ix, int iy, double x0, double x1, 
		    double y0, double y1) const=0;

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The number of x grid points
    size_t nx;

    /// The number of y grid points
    size_t ny;

    /// The x grid
    vec_t *xfun;

    /// The y grid
    vec_t *yfun;

    /// The data
    mat_t *datap;
    
  private:
    
    interp2_base<vec_t,mat_t>(const interp2_base<vec_t,mat_t> &);
    interp2_base<vec_t,mat_t>& operator=(const interp2_base<vec_t,mat_t>&);

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



