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
 * This file is part of interp2d, a GSL-compatible two-dimensional
 * interpolation library. <http://www.ellipsix.net/interp2d.html>
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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/interp2_direct.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

using namespace std;
using namespace o2scl;

interp2_direct::interp2_direct() {
  data_set=false;
  itype=itp_cspline;
}

void interp2_direct::set_data(size_t n_x, size_t n_y, ubvector &x_grid,
			   ubvector &y_grid, ubmatrix &data, 
			   size_t interp_type) {
  
  if (interp_type!=itp_linear && interp_type!=itp_cspline &&
      interp_type!=itp_cspline_peri) {
    O2SCL_ERR2("Unsupported interpolation type in ",
	       "interp2_direct::set_data().",exc_eunimpl);
  }
  
  nx=n_x;
  ny=n_y;
  xfun=&x_grid;
  yfun=&y_grid;
  datap=&data;
  itype=interp_type;

  svx.set_vec(n_x,x_grid);
  svy.set_vec(n_y,y_grid);
  
  if (interp_type==itp_cspline || interp_type==itp_cspline_peri) {
    
    zx.resize(n_x,n_y);
    zy.resize(n_x,n_y);
    zxy.resize(n_x,n_y);
  
    // Partial derivative with respect to x
    for(size_t j=0;j<n_y;j++) {
      ubmatrix_column col=
	o2scl::matrix_column<ubmatrix,ubmatrix_column>(data,j);
      interp_vec<ubvector,ubmatrix_column> itp(n_x,x_grid,col,interp_type);
      for(size_t i=0;i<n_x;i++) {
	zx(i,j)=itp.deriv(x_grid[i]);
      }
    }

    // Partial derivative with respect to y
    for(size_t i=0;i<n_x;i++) {
      ubmatrix_row row=
	o2scl::matrix_row<ubmatrix,ubmatrix_row>(data,i);
      interp_vec<ubvector,ubmatrix_row> itp(n_y,y_grid,row,interp_type);
      for(size_t j=0;j<n_y;j++) {
	zy(i,j)=itp.deriv(y_grid[j]);
      }
    }

    // Mixed partial derivative
    for(size_t j=0;j<n_y;j++) {
      ubmatrix_column col=
	o2scl::matrix_column<ubmatrix,ubmatrix_column>(zy,j);
      interp_vec<ubvector,ubmatrix_column> itp(n_x,x_grid,col,interp_type);
      for(size_t i=0;i<n_x;i++) {
	zxy(i,j)=itp.deriv(x_grid[i]);
      }
    }

  }

  data_set=true;

  return;
}

double interp2_direct::eval(double x, double y) {

  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::eval().",exc_einval);
  }
  
  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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

double interp2_direct::deriv_x(double x, double y) {
  
  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::deriv_x().",exc_einval);
  }

  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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

double interp2_direct::deriv_y(double x, double y) {
  
  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::deriv_y().",exc_einval);
  }

  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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

double interp2_direct::deriv_xx(double x, double y) {

  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::deriv_xx().",exc_einval);
  }

  if (itype==itp_linear) {
    return 0.0;
  }
  
  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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

double interp2_direct::deriv_xy(double x, double y) {
  
  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::deriv_xy().",exc_einval);
  }

  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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

double interp2_direct::deriv_yy(double x, double y) {
  
  if (!data_set) {
    O2SCL_ERR("Data not set in interp2_direct::deriv_yy().",exc_einval);
  }

  if (itype==itp_linear) {
    return 0.0;
  }

  size_t xi=svx.find(x);
  size_t yi=svy.find(y);

  double xmin=(*xfun)[xi];
  double xmax=(*xfun)[xi+1];
  double ymin=(*yfun)[yi];
  double ymax=(*yfun)[yi+1];

  double zminmin=(*datap)(xi,yi);
  double zminmax=(*datap)(xi,yi+1);
  double zmaxmin=(*datap)(xi+1,yi);
  double zmaxmax=(*datap)(xi+1,yi+1);

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
