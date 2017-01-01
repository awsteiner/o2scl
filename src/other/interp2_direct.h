/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2017, Andrew W. Steiner
  
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
/*
 * Some of the code in this file was originally part of interp2d, a
 * GSL-compatible two-dimensional interpolation library.
 * <http://www.ellipsix.net/interp2d.html>
 *
 * Copyright 2012 Thomas Beutlich, David Zaslavsky
 * Portions based on alglib 3.6 interpolation code,
 *  copyright Sergey Bochkanov
 * Portions based on GNU GSL interpolation code,
 *  copyright 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file interp2_direct.h
    \brief File defining \ref o2scl::interp2_direct
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
#include <o2scl/interp2.h>
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
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double>,
    class mat_row_t=boost::numeric::ublas::matrix_row<mat_t>,
    class mat_column_t=boost::numeric::ublas::matrix_column<mat_t> >
    class interp2_direct : public interp2_base<vec_t,mat_t> {
    
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
  public:
    
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_col;

    interp2_direct() {
      data_set=false;
      itype=itp_cspline;
    }

    /** \brief Initialize the data for the 2-dimensional interpolation
     */
    void set_data(size_t n_x, size_t n_y, vec_t &x_grid,
		  vec_t &y_grid, mat_t &data, 
		  size_t interp_type=itp_cspline) {
  
      if (interp_type!=itp_linear && interp_type!=itp_cspline &&
	  interp_type!=itp_cspline_peri) {
	O2SCL_ERR2("Unsupported interpolation type in ",
		   "interp2_direct::set_data().",exc_eunimpl);
      }
  
      this->nx=n_x;
      this->ny=n_y;
      this->xfun=&x_grid;
      this->yfun=&y_grid;
      this->datap=&data;
      itype=interp_type;

      svx.set_vec(n_x,x_grid);
      svy.set_vec(n_y,y_grid);
  
      if (interp_type==itp_cspline || interp_type==itp_cspline_peri) {
    
	zx.resize(n_x,n_y);
	zy.resize(n_x,n_y);
	zxy.resize(n_x,n_y);
  
	// Partial derivative with respect to x
	for(size_t j=0;j<n_y;j++) {
	  mat_column_t col=
	    o2scl::matrix_column<mat_t,mat_column_t>(data,j);
	  interp_vec<vec_t,mat_column_t> itp(n_x,x_grid,col,interp_type);
	  for(size_t i=0;i<n_x;i++) {
	    zx(i,j)=itp.deriv(x_grid[i]);
	  }
	}

	// Partial derivative with respect to y
	for(size_t i=0;i<n_x;i++) {
	  mat_row_t row=
	    o2scl::matrix_row<mat_t,mat_row_t>(data,i);
	  interp_vec<vec_t,mat_row_t> itp(n_y,y_grid,row,interp_type);
	  for(size_t j=0;j<n_y;j++) {
	    zy(i,j)=itp.deriv(y_grid[j]);
	  }
	}

	// Mixed partial derivative
	for(size_t j=0;j<n_y;j++) {
	  ubmatrix_col col=
	    o2scl::matrix_column<ubmatrix,ubmatrix_col>(zy,j);
	  interp_vec<vec_t,ubmatrix_col> itp(n_x,x_grid,col,interp_type);
	  for(size_t i=0;i<n_x;i++) {
	    zxy(i,j)=itp.deriv(x_grid[i]);
	  }
	}

      }

      data_set=true;

      return;
    }
    
    /** \brief Perform the 2-d interpolation 
     */
    virtual double eval(double x, double y) const {

      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::eval().",exc_einval);
      }
  
      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;

      if (itype==itp_linear) {
	return (1.0-t)*(1.0-u)*zminmin+t*(1.0-u)*zmaxmin+
	  (1.0-t)*u*zminmax+t*u*zmaxmax;
      }

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;

      double z=0.0;
      double v=zminmin;
      z+=v*t0*u0;
      v=zyminmin;
      z+=v*t0*u1;
      v=-3*zminmin+3*zminmax-2*zyminmin-zyminmax;
      z+=v*t0*u2;
      v=2*zminmin-2*zminmax+zyminmin+zyminmax;
      z+=v*t0*u3;
      v=zxminmin;
      z+=v*t1*u0;
      v=zxyminmin;
      z+=v*t1*u1;
      v=-3*zxminmin+3*zxminmax-2*zxyminmin-zxyminmax;
      z+=v*t1*u2;
      v=2*zxminmin-2*zxminmax+zxyminmin+zxyminmax;
      z+=v*t1*u3;
      v=-3*zminmin+3*zmaxmin-2*zxminmin-zxmaxmin;
      z+=v*t2*u0;
      v=-3*zyminmin+3*zymaxmin-2*zxyminmin-zxymaxmin;
      z+=v*t2*u1;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z+=v*t2*u2;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z+=v*t2*u3;
      v=2*zminmin-2*zmaxmin+zxminmin+zxmaxmin;
      z+=v*t3*u0;
      v=2*zyminmin-2*zymaxmin+zxyminmin+zxymaxmin;
      z+=v*t3*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z+=v*t3*u2;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z+=v*t3*u3;
 
      return z;
    }

    /** \brief Compute the partial derivative in the x-direction
     */
    virtual double deriv_x(double x, double y) const {
  
      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::deriv_x().",exc_einval);
      }

      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;
  
      if (itype==itp_linear) {
	return dt*(-(1.0-u)*zminmin+(1.0-u)*zmaxmin-u*zminmax+u*zmaxmax);
      }

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;

      double z_p=0.0;
      double v=zxminmin;
      z_p+=v*t0*u0;
      v=zxyminmin;
      z_p+=v*t0*u1;
      v=-3*zxminmin+3*zxminmax-2*zxyminmin-zxyminmax;
      z_p+=v*t0*u2;
      v=2*zxminmin-2*zxminmax+zxyminmin+zxyminmax;
      z_p+=v*t0*u3;
      v=-3*zminmin+3*zmaxmin-2*zxminmin-zxmaxmin;
      z_p+=2*v*t1*u0;
      v=-3*zyminmin+3*zymaxmin-2*zxyminmin-zxymaxmin;
      z_p+=2*v*t1*u1;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z_p+=2*v*t1*u2;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z_p+=2*v*t1*u3;
      v=2*zminmin-2*zmaxmin+zxminmin+zxmaxmin;
      z_p+=3*v*t2*u0;
      v=2*zyminmin-2*zymaxmin+zxyminmin+zxymaxmin;
      z_p+=3*v*t2*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z_p+=3*v*t2*u2;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z_p+=3*v*t2*u3;
      z_p*=dt;

      return z_p;
    }

    /** \brief Compute the partial second derivative in the x-direction
     */
    virtual double deriv_xx(double x, double y) const {

      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::deriv_xx().",exc_einval);
      }

      if (itype==itp_linear) {
	return 0.0;
      }
  
      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;

      double z_pp=0.0;
      double v=-3*zminmin+3*zmaxmin-2*zxminmin-zxmaxmin;
      z_pp+=2*v*t0*u0;
      v=-3*zyminmin+3*zymaxmin-2*zxyminmin-zxymaxmin;
      z_pp+=2*v*t0*u1;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z_pp+=2*v*t0*u2;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z_pp+=2*v*t0*u3;
      v=2*zminmin-2*zmaxmin+zxminmin+zxmaxmin;
      z_pp+=6*v*t1*u0;
      v=2*zyminmin-2*zymaxmin+zxyminmin+zxymaxmin;
      z_pp+=6*v*t1*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z_pp+=6*v*t1*u2;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z_pp+=6*v*t1*u3;
      z_pp*=dt*dt;

      return z_pp;
    }

    /** \brief Compute the partial derivative in the y-direction
     */
    virtual double deriv_y(double x, double y) const {
  
      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::deriv_y().",exc_einval);
      }

      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;

      if (itype==itp_linear) {
	return du*(-(1.0-t)*zminmin-t*zmaxmin+(1.0-t)*zminmax+t*zmaxmax);
      }

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;
  
      double z_p=0.0;
      double v=zyminmin;
      z_p+=v*t0*u0;
      v=-3*zminmin+3*zminmax-2*zyminmin-zyminmax;
      z_p+=2*v*t0*u1;
      v=2*zminmin-2*zminmax+zyminmin+zyminmax;
      z_p+=3*v*t0*u2;
      v=zxyminmin;
      z_p+=v*t1*u0;
      v=-3*zxminmin+3*zxminmax-2*zxyminmin-zxyminmax;
      z_p+=2*v*t1*u1;
      v=2*zxminmin-2*zxminmax+zxyminmin+zxyminmax;
      z_p+=3*v*t1*u2;
      v=-3*zyminmin+3*zymaxmin-2*zxyminmin-zxymaxmin;
      z_p+=v*t2*u0;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z_p+=2*v*t2*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z_p+=3*v*t2*u2;
      v=2*zyminmin-2*zymaxmin+zxyminmin+zxymaxmin;
      z_p+=v*t3*u0;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z_p+=2*v*t3*u1;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z_p+=3*v*t3*u2;
      z_p*=du;

      return z_p;
    }

    /** \brief Compute the partial second derivative in the y-direction
     */
    virtual double deriv_yy(double x, double y) const {
  
      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::deriv_yy().",exc_einval);
      }

      if (itype==itp_linear) {
	return 0.0;
      }

      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;

      double z_pp=0.0;
      double v=-3*zminmin+3*zminmax-2*zyminmin-zyminmax;
      z_pp+=2*v*t0*u0;
      v=2*zminmin-2*zminmax+zyminmin+zyminmax;
      z_pp+=6*v*t0*u1;
      v=-3*zxminmin+3*zxminmax-2*zxyminmin-zxyminmax;
      z_pp+=2*v*t1*u0;
      v=2*zxminmin-2*zxminmax+zxyminmin+zxyminmax;
      z_pp+=6*v*t1*u1;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z_pp+=2*v*t2*u0;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z_pp+=6*v*t2*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z_pp+=2*v*t3*u0;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z_pp+=6*v*t3*u1;
      z_pp*=du*du;

      return z_pp;
    }

    /** \brief Compute the mixed partial derivative 
	\f$ \frac{\partial^2 f}{\partial x \partial y} \f$
    */
    virtual double deriv_xy(double x, double y) const {
  
      if (!data_set) {
	O2SCL_ERR("Data not set in interp2_direct::deriv_xy().",exc_einval);
      }

      size_t xi=svx.find(x);
      size_t yi=svy.find(y);

      double xmin=(*this->xfun)[xi];
      double xmax=(*this->xfun)[xi+1];
      double ymin=(*this->yfun)[yi];
      double ymax=(*this->yfun)[yi+1];

      double zminmin=(*this->datap)(xi,yi);
      double zminmax=(*this->datap)(xi,yi+1);
      double zmaxmin=(*this->datap)(xi+1,yi);
      double zmaxmax=(*this->datap)(xi+1,yi+1);

      double dx=xmax-xmin;
      double dy=ymax-ymin;

      double t=(x-xmin)/dx;
      double u=(y-ymin)/dy;
      double dt=1.0/dx;
      double du=1.0/dy;

      if (itype==itp_linear) {
	return dt*du*(zminmin-zmaxmin-zminmax+zmaxmax);
      }

      double zxminmin=zx(xi,yi)/dt;
      double zxminmax=zx(xi,yi+1)/dt;
      double zxmaxmin=zx(xi+1,yi)/dt;
      double zxmaxmax=zx(xi+1,yi+1)/dt;

      double zyminmin=zy(xi,yi)/du;
      double zyminmax=zy(xi,yi+1)/du;
      double zymaxmin=zy(xi+1,yi)/du;
      double zymaxmax=zy(xi+1,yi+1)/du;

      double zxyminmin=zxy(xi,yi)/du/dt;
      double zxyminmax=zxy(xi,yi+1)/du/dt;
      double zxymaxmin=zxy(xi+1,yi)/du/dt;
      double zxymaxmax=zxy(xi+1,yi+1)/du/dt;

      double t0=1.0;
      double t1=t;
      double t2=t*t;
      double t3=t*t2;
      double u0=1.0;
      double u1=u;
      double u2=u*u;
      double u3=u*u2;

      double z_pp=0.0;
      double v=zxyminmin;
      z_pp+=v*t0*u0;
      v=-3*zxminmin+3*zxminmax-2*zxyminmin-zxyminmax;
      z_pp+=2*v*t0*u1;
      v=2*zxminmin-2*zxminmax+zxyminmin+zxyminmax;
      z_pp+=3*v*t0*u2;
      v=-3*zyminmin+3*zymaxmin-2*zxyminmin-zxymaxmin;
      z_pp+=2*v*t1*u0;
      v=9*zminmin-9*zmaxmin+9*zmaxmax-9*zminmax+6*zxminmin+3*zxmaxmin-
	3*zxmaxmax-6*zxminmax+6*zyminmin-6*zymaxmin-3*zymaxmax+
	3*zyminmax+4*zxyminmin+2*zxymaxmin+zxymaxmax+2*zxyminmax;
      z_pp+=4*v*t1*u1;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-4*zxminmin-2*zxmaxmin+
	2*zxmaxmax+4*zxminmax-3*zyminmin+3*zymaxmin+3*zymaxmax-
	3*zyminmax-2*zxyminmin-zxymaxmin-zxymaxmax-2*zxyminmax;
      z_pp+=6*v*t1*u2;
      v=2*zyminmin-2*zymaxmin+zxyminmin+zxymaxmin;
      z_pp+=3*v*t2*u0;
      v=-6*zminmin+6*zmaxmin-6*zmaxmax+6*zminmax-3*zxminmin-3*zxmaxmin+
	3*zxmaxmax+3*zxminmax-4*zyminmin+4*zymaxmin+2*zymaxmax-
	2*zyminmax-2*zxyminmin-2*zxymaxmin-zxymaxmax-zxyminmax;
      z_pp+=6*v*t2*u1;
      v=4*zminmin-4*zmaxmin+4*zmaxmax-4*zminmax+2*zxminmin+2*zxmaxmin-
	2*zxmaxmax-2*zxminmax+2*zyminmin-2*zymaxmin-2*zymaxmax+
	2*zyminmax+zxyminmin+zxymaxmin+zxymaxmax+zxyminmax;
      z_pp+=9*v*t2*u2;
      z_pp*=dt*du;

      return z_pp;
    }

    virtual double integ_x(double x0, double x1, double y) const {
      O2SCL_ERR("Integration unimplemented in interp2_direct.",
		exc_eunimpl);
      return 0.0;
    }

    virtual double integ_y(double x, double y0, double y1) const {
      O2SCL_ERR("Integration unimplemented in interp2_direct.",
		exc_eunimpl);
      return 0.0;
    }

    virtual double eval_gen(int ix, int iy, double x0, double x1, 
			    double y0, double y1) const {
      O2SCL_ERR("Function eval_gen() unimplemented in interp2_direct.",
		exc_eunimpl);
      return 0.0;
    }


#ifndef DOXYGEN_NO_O2NS
    
  protected:

    /// True if the data has been specified by the user
    bool data_set;

    /// Interpolation type
    size_t itype;

    /// Partial derivative with respect to x
    ubmatrix zx;

    /// Partial derivative with respect to y
    ubmatrix zy;

    /// Mixed partial derivative
    ubmatrix zxy;

    /// Searching object for x-direction
    search_vec<vec_t> svx;

    /// Searching object for y-direction
    search_vec<vec_t> svy;
    
  private:

    interp2_direct<vec_t,mat_t,mat_row_t,mat_column_t>
      (const interp2_direct<vec_t,mat_t,mat_row_t,mat_column_t> &);
    interp2_direct<vec_t,mat_t,mat_row_t,mat_column_t>& operator=
      (const interp2_direct<vec_t,mat_t,mat_row_t,mat_column_t>&);

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



