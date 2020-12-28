/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/interp2_direct.h>
#include <o2scl/interp2_seq.h>
#include <o2scl/columnify.h>
#include <o2scl/test_mgr.h>
#include <o2scl/tensor_grid.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

double f(double x, double y) {
  return pow(sin(0.1*x+0.3*y),2.0);
}

double fx(double x, double y) {
  return 2.0*sin(0.1*x+0.3*y)*0.1*cos(0.1*x+0.3*y);
}

double fy(double x, double y) {
  return 2.0*sin(0.1*x+0.3*y)*0.3*cos(0.1*x+0.3*y);
}

double f2x(double x, double y) {
  return 2.0*0.1*0.1*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

double f2y(double x, double y) {
  return 2.0*0.3*0.3*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

double fxy(double x, double y) {
  return 2.0*0.1*0.3*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  interp2_seq<ubvector,ubmatrix,ubmatrix_row> it;
  interp2_direct<ubvector,ubmatrix,ubmatrix_row,ubmatrix_column> it2;
  double x0, x1, y0, y1;
  double tol=0.05;
  double tol2=0.4;

  {
    // N=M
    
    size_t M=40;
    size_t N=40;
    ubvector x2(M), y2(N);
    ubmatrix data2(M,N);
    for(size_t ii=0;ii<M;ii++) {
      x2[ii]=((double)ii)/10.0;
    }
    for(size_t jj=0;jj<N;jj++) {
      y2[jj]=((double)jj)/20.0;
    }
    for(size_t ii=0;ii<M;ii++) {
      for(size_t jj=0;jj<N;jj++) {
	data2(ii,jj)=f(x2[ii],y2[jj]);
      }
    }
  
    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;
    
    it2.set_data(M,N,x2,y2,data2);
    it.set_data(M,N,x2,y2,data2);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it2.eval(x0,y0),f(x0,y0),tol,"cspline M=N i");
      t.test_rel(it2.deriv_x(x0,y0),fx(x0,y0),tol,"cspline M=N dx");
      t.test_rel(it2.deriv_y(x0,y0),fy(x0,y0),tol,"cspline M=N dy");
      t.test_rel(it2.deriv_xy(x0,y0),fxy(x0,y0),tol,"cspline M=N dyx");
      t.test_rel(it2.deriv_xx(x0,y0),f2x(x0,y0),tol2,"cspline M=N d2x");
      t.test_rel(it2.deriv_yy(x0,y0),f2y(x0,y0),tol2,"cspline M=N d2y");
      t.test_rel(it2.eval(x0,y0),it.eval(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
      t.test_rel(it2.deriv_x(x0,y0),it.deriv_x(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
      t.test_rel(it2.deriv_y(x0,y0),it.deriv_y(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
      t.test_rel(it2.deriv_xy(x0,y0),it.deriv_xy(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
      t.test_rel(it2.deriv_xx(x0,y0),it.deriv_xx(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
      t.test_rel(it2.deriv_yy(x0,y0),it.deriv_yy(x0,y0),1.0e-9,
		 "cspline M=N it2 vs. it");
    }

  }

  {
    // M>N
    
    size_t M=50;
    size_t N=40;
    ubvector x2(M), y2(N);
    ubmatrix data2(M,N);
    for(size_t ii=0;ii<M;ii++) {
      x2[ii]=((double)ii)/10.0;
    }
    for(size_t jj=0;jj<N;jj++) {
      y2[jj]=((double)jj)/20.0;
    }
    for(size_t ii=0;ii<M;ii++) {
      for(size_t jj=0;jj<N;jj++) {
	data2(ii,jj)=f(x2[ii],y2[jj]);
      }
    }
  
    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

    it2.set_data(M,N,x2,y2,data2);
    it.set_data(M,N,x2,y2,data2);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it2.eval(x0,y0),f(x0,y0),tol,"cspline M>N dx");
      t.test_rel(it2.deriv_x(x0,y0),fx(x0,y0),tol,"cspline M>N dx");
      t.test_rel(it2.deriv_y(x0,y0),fy(x0,y0),tol,"cspline M>N dy");
      t.test_rel(it2.deriv_xy(x0,y0),fxy(x0,y0),tol,"cspline M>N dyx");
      t.test_rel(it2.deriv_xx(x0,y0),f2x(x0,y0),tol2,"cspline M>N d2x");
      t.test_rel(it2.deriv_yy(x0,y0),f2y(x0,y0),tol2,"cspline M>N d2y");
      t.test_rel(it2.eval(x0,y0),it.eval(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
      t.test_rel(it2.deriv_x(x0,y0),it.deriv_x(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
      t.test_rel(it2.deriv_y(x0,y0),it.deriv_y(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
      t.test_rel(it2.deriv_xy(x0,y0),it.deriv_xy(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
      t.test_rel(it2.deriv_xx(x0,y0),it.deriv_xx(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
      t.test_rel(it2.deriv_yy(x0,y0),it.deriv_yy(x0,y0),1.0e-9,
		 "cspline M>N it2 vs. it");
    }

  }

  {
    // N>M
    
    size_t M=40;
    size_t N=50;
    ubvector x2(M), y2(N);
    ubmatrix data2(M,N);
    for(size_t ii=0;ii<M;ii++) {
      x2[ii]=((double)ii)/10.0;
    }
    for(size_t jj=0;jj<N;jj++) {
      y2[jj]=((double)jj)/20.0;
    }
    for(size_t ii=0;ii<M;ii++) {
      for(size_t jj=0;jj<N;jj++) {
	data2(ii,jj)=f(x2[ii],y2[jj]);
      }
    }
  
    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

    it2.set_data(M,N,x2,y2,data2);
    it.set_data(M,N,x2,y2,data2);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it2.eval(x0,y0),f(x0,y0),tol,"cspline M<N dx");
      t.test_rel(it2.eval(x0,y0),f(x0,y0),tol,"cspline M<N dx");
      t.test_rel(it2.deriv_x(x0,y0),fx(x0,y0),tol,"cspline M<N dx");
      t.test_rel(it2.deriv_y(x0,y0),fy(x0,y0),tol,"cspline M<N dy");
      t.test_rel(it2.deriv_xy(x0,y0),fxy(x0,y0),tol,"cspline M<N dyx");
      t.test_rel(it2.deriv_xx(x0,y0),f2x(x0,y0),tol2,"cspline M<N d2x");
      t.test_rel(it2.deriv_yy(x0,y0),f2y(x0,y0),tol2,"cspline M<N d2y");
      t.test_rel(it2.eval(x0,y0),it.eval(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
      t.test_rel(it2.deriv_x(x0,y0),it.deriv_x(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
      t.test_rel(it2.deriv_y(x0,y0),it.deriv_y(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
      t.test_rel(it2.deriv_xy(x0,y0),it.deriv_xy(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
      t.test_rel(it2.deriv_xx(x0,y0),it.deriv_xx(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
      t.test_rel(it2.deriv_yy(x0,y0),it.deriv_yy(x0,y0),1.0e-9,
		 "cspline M<N it2 vs. it");
    }

  }

  {
    // N=M with linear interpolation
    
    size_t M=40;
    size_t N=40;
    ubvector x2(M), y2(N);
    ubmatrix data2(M,N);
    for(size_t ii=0;ii<M;ii++) {
      x2[ii]=((double)ii)/10.0;
    }
    for(size_t jj=0;jj<N;jj++) {
      y2[jj]=((double)jj)/20.0;
    }
    for(size_t ii=0;ii<M;ii++) {
      for(size_t jj=0;jj<N;jj++) {
	data2(ii,jj)=f(x2[ii],y2[jj]);
      }
    }
  
    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;
    
    it2.set_data(M,N,x2,y2,data2,itp_linear);
    it.set_data(M,N,x2,y2,data2,itp_linear);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it2.eval(x0,y0),f(x0,y0),tol,"linear i");
      t.test_rel(it2.deriv_x(x0,y0),fx(x0,y0),tol,"linear dx");
      t.test_rel(it2.deriv_y(x0,y0),fy(x0,y0),tol,"linear dy");
      t.test_rel(it2.deriv_xy(x0,y0),fxy(x0,y0),1.5,"linear dyx");
      t.test_rel(it2.eval(x0,y0),it.eval(x0,y0),1.0e-9,
		 "linear it2 vs. it");
      t.test_rel(it2.deriv_x(x0,y0),it.deriv_x(x0,y0),1.0e-9,
		 "linear it2 vs. it");
      t.test_rel(it2.deriv_y(x0,y0),it.deriv_y(x0,y0),1.0e-9,
		 "linear it2 vs. it");
      t.test_rel(it2.deriv_xy(x0,y0),it.deriv_xy(x0,y0),1.0e-9,
		 "linear it2 vs. it");
      t.test_rel(it2.deriv_xx(x0,y0),it.deriv_xx(x0,y0),1.0e-9,
		 "linear it2 vs. it");
      t.test_rel(it2.deriv_yy(x0,y0),it.deriv_yy(x0,y0),1.0e-9,
		 "linear it2 vs. it");
    }

  }

  {
    // Show how to slice a tensor
    tensor_grid3<> tg(3,2,1);
    double grid[6]={4,5,6,7,8,9};
    tg.set_grid_packed(grid);
    for(size_t j=0;j<3;j++) {
      for(size_t k=0;k<2;k++) {
	for(size_t ell=0;ell<1;ell++) {
	  tg.set(j,k,ell,((double)(j+k+ell)));
	}
      }
    }
    std::function<double &(size_t,size_t)> slice=
      std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
		(&tensor_grid3<>::get),
		&tg,std::placeholders::_1,0,std::placeholders::_2);
    
  }

  t.report();
  return 0;
}

