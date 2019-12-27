/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
    \brief File defining \ref o2scl::interp2_seq
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
#include <o2scl/interp2.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Two-dimensional interpolation class by successive
      one-dimensional interpolation

      This class implements two-dimensional interpolation by iterating
      the \o2 one-dimensional interpolation routines. Derivatives and
      integrals along both x- and y-directions can be computed. This
      class is likely a bit slower than \ref interp2_direct but more
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

      Because of the way this class creates pointers to the
      data, copy construction is not currently allowed. 
      
      \future Implement an improved caching system in case, for example
      \c xfirst is true and the last interpolation used the same
      value of \c x.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double>,
    class mat_row_t=boost::numeric::ublas::matrix_row<mat_t> >
    class interp2_seq : public interp2_base<vec_t,mat_t> {

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
    interp2_seq() {
      data_set=false;
      itype=itp_cspline;
    }
    
    virtual ~interp2_seq() {
      for(size_t i=0;i<itps.size();i++) {
	delete itps[i];
	delete vecs[i];
      }
      itps.clear();
    }
    
    /** \brief Initialize the data for the 2-dimensional interpolation
	
	If \c x_first is true, then set_data() creates interpolation
	objects for each of the rows. Calls to interp() then uses
	these to create a column at the specified value of \c x. An
	interpolation object is created at this column to find the
	value of the function at the specified value \c y. If \c
	x_first is false, the opposite strategy is employed. These two
	options may give slightly different results.
    */
    void set_data(size_t n_x, size_t n_y, vec_t &x_grid,
		  vec_t &y_grid, mat_t &data, 
		  size_t interp_type=itp_cspline) {
      
      // Set new data
      itype=interp_type;
      nx=n_x;
      ny=n_y;
      xfun=&x_grid;
      yfun=&y_grid;
      datap=&data;
      
      // Set interpolation objects
      
      for(size_t i=0;i<itps.size();i++) {
	delete itps[i];
	delete vecs[i];
      }
      itps.clear();
      
      // If we interpolate along the x-axis first, then we want to fix the
      // first index, to get nx rows of size ny
      vecs.resize(nx);
      itps.resize(nx);
      for(size_t i=0;i<nx;i++) {
	vecs[i]=new mat_row_t
	  (o2scl::matrix_row<mat_t,mat_row_t>(*datap,i));
	itps[i]=new interp_vec<vec_t,mat_row_t>(ny,*yfun,*vecs[i],itype);
      }
      data_set=true;
      
      return;
    }
    
    /** \brief Reset the stored interpolation since the data has changed
	
        This will throw an exception if the set_data() has not been
        called.
    */
    void reset_interp() {
      if (data_set) {
	
	for(size_t i=0;i<itps.size();i++) {
	  delete itps[i];
	  delete vecs[i];
	}
	itps.clear();
	
	// Set interpolation objects
	vecs.resize(nx);
	itps.resize(nx);
	for(size_t i=0;i<nx;i++) {
	  vecs[i]=new mat_row_t
	    (o2scl::matrix_row<mat_t,mat_row_t>(*datap,i));
	  itps[i]=new interp_vec<vec_t,mat_row_t>(ny,*yfun,*vecs[i],itype);
	}
	
      } else {
	O2SCL_ERR("Data not set in interp2_seq::reset_interp().",exc_einval);
      }
      
      return;
    }
    
    /** \brief Perform the 2-d interpolation 
     */
    double eval(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->eval(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.eval(x);
      return result;
    }
    
    /** \brief Perform the 2-d interpolation 
     */
    double operator()(double x, double y) const {
      return eval(x,y);
    }
    
    /** \brief Compute the partial derivative in the x-direction
     */
    double deriv_x(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->eval(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.deriv(x);
      return result;
    }

    /** \brief Compute the partial second derivative in the x-direction
     */
    double deriv_xx(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::deriv_xx().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->eval(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.deriv2(x);
      return result;
    }

    /** \brief Compute the integral in the x-direction between x=x0
	and x=x1
    */
    double integ_x(double x0, double x1, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::integ_x().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->eval(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.integ(x0,x1);
      return result;

    }

    /** \brief Compute the partial derivative in the y-direction
     */
    double deriv_y(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::deriv_y().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->deriv(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.eval(x);
      return result;
    }

    /** \brief Compute the partial second derivative in the y-direction
     */
    double deriv_yy(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::deriv_yy().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->deriv2(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.eval(x);
      return result;
    }

    /** \brief Compute the integral in the y-direction between y=y0
	and y=y1
    */
    double integ_y(double x, double y0, double y1) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::integ_y().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->integ(y0,y1);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.eval(x);
      return result;
    }

    /** \brief Compute the mixed partial derivative 
	\f$ \frac{\partial^2 f}{\partial x \partial y} \f$
    */
    double deriv_xy(double x, double y) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
	return 0.0;
      }
      double result;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	icol[i]=itps[i]->deriv(y);
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      result=six.deriv(x);
      return result;
    }

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
    double eval_gen(int ix, int iy, double x0, double x1, 
		    double y0, double y1) const {
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_seq::interp_gen().",exc_efailed);
	return 0.0;
      }
      // Initialize to prevent un'inited var. warnings
      double result=0.0;
      ubvector icol(nx);
      for(size_t i=0;i<nx;i++) {
	if (iy==-1) {
	  icol[i]=itps[i]->integ(y0,y1);
	} else if (iy==0) {
	  icol[i]=itps[i]->eval(y0);
	} else if (iy==1) {
	  icol[i]=itps[i]->deriv(y0);
	} else if (iy==2) {
	  icol[i]=itps[i]->deriv2(y0);
	} else {
	  O2SCL_ERR2("Invalid value of 'iy' for interp2_seq::",
		     "interp_gen(). (xfirst=true)",exc_einval);
	}
      }
      interp_vec<vec_t,ubvector> six(nx,*xfun,icol,itype);
      if (ix==-1) {
	result=six.integ(x0,x1);
      } else if (ix==0) {
	result=six.eval(x0);
      } else if (ix==1) {
	result=six.deriv(x0);
      } else if (ix==2) {
	result=six.deriv2(x0);
      } else {
	O2SCL_ERR2("Invalid value of 'ix' for interp2_seq::",
		   "interp_gen(). (xfirst=true)",exc_einval);
      }
      return result;
    }

#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// The array of interpolation objects
    std::vector<interp_vec<vec_t,mat_row_t> *> itps;

    /// An array of rows
    std::vector<mat_row_t *> vecs;

    /// The number of x grid points
    size_t nx;

    /// The number of y grid points
    size_t ny;

    /// True if the data has been specified by the user
    bool data_set;

    /// The x grid
    vec_t *xfun;

    /// The y grid
    vec_t *yfun;

    /// The data
    mat_t *datap;
    
    /// Interpolation type
    size_t itype;

  private:

    interp2_seq<vec_t,mat_t,mat_row_t>
      (const interp2_seq<vec_t,mat_t,mat_row_t> &);
    interp2_seq<vec_t,mat_t,mat_row_t>& operator=
      (const interp2_seq<vec_t,mat_t,mat_row_t>&);

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



