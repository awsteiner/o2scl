/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

#include <o2scl/interp2_seq.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

double f(double x, double y) {
  return pow(sin(0.1*x+0.3*y),2.0);
}

double fix(double x0, double x1, double y) {
  return (x1-x0)/2.0+(sin(2.0*(0.1*x0+0.3*y))-sin(2.0*(0.1*x1+0.3*y)))/0.4;
}

double fiy(double x, double y0, double y1) {
  return (y1-y0)/2.0+(sin(2.0*(0.1*x+0.3*y0))-sin(2.0*(0.1*x+0.3*y1)))/1.2;
}

double fx(double x, double y) {
  return 2.0*sin(0.1*x+0.3*y)*0.1*cos(0.1*x+0.3*y);
}

double fy(double x, double y) {
  return 2.0*sin(0.1*x+0.3*y)*0.3*cos(0.1*x+0.3*y);
}

double fxx(double x, double y) {
  return 2.0*0.1*0.1*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

double fyy(double x, double y) {
  return 2.0*0.3*0.3*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

double fxy(double x, double y) {
  return 2.0*0.1*0.3*(pow(cos(0.1*x+0.3*y),2.0)-
		      pow(sin(0.1*x+0.3*y),2.0));
}

double fxxy(double x, double y) {
  return -8.0*0.1*0.1*0.3*cos(0.1*x+0.3*y)*sin(0.1*x+0.3*y);
}

double fxyy(double x, double y) {
  return -8.0*0.1*0.3*0.3*cos(0.1*x+0.3*y)*sin(0.1*x+0.3*y);
}

double fxxyy(double x, double y) {
  return 8.0*0.1*0.1*0.3*0.3*(pow(sin(0.1*x+0.3*y),2.0)-
			      pow(cos(0.1*x+0.3*y),2.0));
}

template<class vec_t, class mat_t> 
int print_data(int nx, int ny, vec_t &x, vec_t &y, mat_t &data) {
  int j, k;
  cout.setf(ios::showpos);
  cout.precision(3);
  cout << "              ";
  for(k=0;k<ny;k++) {
    string stx="y";
    stx+=itos(k);
    cout.width(11);
    cout << stx;
  }
  cout << endl;
  cout << "               ";
  for(k=0;k<ny;k++) {
    cout << y[k] << " ";
  }
  cout << endl;
  for(j=0;j<nx;j++) {
    string stx="x";
    stx+=itos(j);
    cout.width(3);
    cout << stx << " ";
    cout << x[j] << " ";
    for(k=0;k<ny;k++) {
      cout << data(j,k) << " ";
    }
    cout << endl;
  }
  cout.unsetf(ios::showpos);
  cout.precision(6);
  return 0;
}

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  interp2_seq<ubvector,ubmatrix,ubmatrix_row> it;
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
  
    // x-first

    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

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
      t.test_rel(it.eval(x0,y0),f(x0,y0),1.0e-6,"2dintp x 1 i");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),1.0e-4,"2dintp x 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),5.0e-4,"2dintp x 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),7.0e-2,"2dintp x 1 dyx");
      t.test_abs(it.deriv_xx(x0,y0),fxx(x0,y0),1.0e-3,"2dintp x 1 d2x");
      t.test_abs(it.deriv_yy(x0,y0),fyy(x0,y0),5.0e-2,"2dintp x 1 d2y");
      t.test_rel(it.eval_gen(2,1,x0,x1,y0,y1),fxxy(x0,y0),
		 5.0e-1,"gen vs. spec.");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),5.0e-6,"2dintp x 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),5.0e-6,"2dintp x 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),1.0e-6,"2dintp x 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),1.0e-6,"2dintp x 1 diy");
      }
    }

    // y-first

    x0=2.42;
    y0=1.41;
    y0=1.41;
    y1=0.62;

    //it.set_data(M,N,x2,y2,data2,itp_cspline,false);
    it.set_data(M,N,x2,y2,data2,itp_cspline);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[N-1]) x0/=2.0;
      while (x1>x2[N-1]) x1/=2.0;
      while (y0>y2[M-1]) y0/=2.0;
      while (y1>y2[M-1]) y1/=2.0;
      t.test_rel(it.eval(x0,y0),f(x0,y0),1.0e-6,"2dintp y 1 i");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),1.0e-4,"2dintp y 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),5.0e-4,"2dintp y 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),7.0e-2,"2dintp y 1 dyx");
      t.test_abs(it.deriv_xx(x0,y0),fxx(x0,y0),1.0e-3,"2dintp y 1 d2x");
      t.test_abs(it.deriv_yy(x0,y0),fyy(x0,y0),5.0e-2,"2dintp y 1 d2y");
      t.test_rel(it.eval_gen(2,1,x0,x1,y0,y1),fxxy(x0,y0),
		 5.0e-1,"gen vs. spec.");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),5.0e-6,"2dintp y 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),5.0e-6,"2dintp y 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),1.0e-6,"2dintp y 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),1.0e-6,"2dintp y 1 diy");
      }
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
  
    // x-first

    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

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
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpM x 1 dx");
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpM x 1 dx");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintpM x 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintpM x 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),tol,"2dintpM x 1 dyx");
      t.test_rel(it.deriv_xx(x0,y0),fxx(x0,y0),tol2,"2dintpM x 1 d2x");
      t.test_rel(it.deriv_yy(x0,y0),fyy(x0,y0),tol2,"2dintpM x 1 d2y");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintpM x 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintpM x 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintpM x 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintpM x 1 diy");
      }
    }

    // y-first

    x0=2.42;
    y0=1.41;
    y0=1.41;
    y1=0.62;

    //it.set_data(M,N,x2,y2,data2,itp_cspline,false);
    it.set_data(M,N,x2,y2,data2,itp_cspline);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpM y 1 i");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintpM y 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintpM y 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),tol,"2dintpM y 1 dyx");
      t.test_rel(it.deriv_xx(x0,y0),fxx(x0,y0),tol2,"2dintpM y 1 d2x");
      t.test_rel(it.deriv_yy(x0,y0),fyy(x0,y0),tol2,"2dintpM y 1 d2y");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintpM y 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintpM y 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintpM y 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintpM y 1 diy");
      }
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
  
    // x-first

    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

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
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpN x 1 dx");
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpN x 1 dx");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintpN x 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintpN x 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),tol,"2dintpN x 1 dyx");
      t.test_rel(it.deriv_xx(x0,y0),fxx(x0,y0),tol2,"2dintpN x 1 d2x");
      t.test_rel(it.deriv_yy(x0,y0),fyy(x0,y0),tol2,"2dintpN x 1 d2y");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintpN x 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintpN x 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintpN x 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintpN x 1 diy");
      }
    }

    // y-first

    x0=2.42;
    y0=1.41;
    y0=1.41;
    y1=0.62;

    //it.set_data(M,N,x2,y2,data2,itp_cspline,false);
    it.set_data(M,N,x2,y2,data2,itp_cspline);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintpN y 1 i");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintpN y 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintpN y 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),tol,"2dintpN y 1 dyx");
      t.test_rel(it.deriv_xx(x0,y0),fxx(x0,y0),tol2,"2dintpN y 1 d2x");
      t.test_rel(it.deriv_yy(x0,y0),fyy(x0,y0),tol2,"2dintpN y 1 d2y");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintpN y 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintpN y 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintpN y 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintpN y 1 diy");
      }
    }

  }

  {
    // N=M with linear interpolaiton
    
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
  
    // x-first

    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;
    
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
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintp x 1 dx");
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintp x 1 dx");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintp x 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintp x 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),1.5,"2dintp x 1 dyx");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintp x 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintp x 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintp x 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintp x 1 diy");
      }
    }

    // y-first

    x0=2.42;
    y0=1.41;
    y0=1.41;
    y1=0.62;

    //it.set_data(M,N,x2,y2,data2,itp_linear,false);
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
      t.test_rel(it.eval(x0,y0),f(x0,y0),tol,"2dintp y 1 i");
      t.test_rel(it.deriv_x(x0,y0),fx(x0,y0),tol,"2dintp y 1 dx");
      t.test_rel(it.deriv_y(x0,y0),fy(x0,y0),tol,"2dintp y 1 dy");
      t.test_rel(it.deriv_xy(x0,y0),fxy(x0,y0),1.5,"2dintp y 1 dyx");
      if (x0<x1) {
	t.test_rel(it.integ_x(x0,x1,y0),fix(x0,x1,y0),tol2,"2dintp y 1 dix");
      } else {
	t.test_rel(it.integ_x(x1,x0,y0),fix(x1,x0,y0),tol2,"2dintp y 1 dix");
      }
      if (y0<y1) {
	t.test_rel(it.integ_y(x0,y0,y1),fiy(x0,y0,y1),tol2,"2dintp y 1 diy");
      } else {
	t.test_rel(it.integ_y(x0,y1,y0),fiy(x0,y1,y0),tol2,"2dintp y 1 diy");
      }
    }

  }

#ifdef O2SCL_ARMA

  {
    // N=M
    
    interp2_seq<arma::rowvec,arma::mat,arma::subview_row<double> > ait;
    
    size_t M=40;
    size_t N=40;
    arma::rowvec ax2(M), ay2(N);
    ubvector x2(M), y2(N);
    arma::mat adata2(M,N);
    ubmatrix data2(M,N);
    for(size_t ii=0;ii<M;ii++) {
      x2[ii]=((double)ii)/10.0;
      ax2[ii]=((double)ii)/10.0;
    }
    for(size_t jj=0;jj<N;jj++) {
      y2[jj]=((double)jj)/20.0;
      ay2[jj]=((double)jj)/20.0;
    }
    for(size_t ii=0;ii<M;ii++) {
      for(size_t jj=0;jj<N;jj++) {
	data2(ii,jj)=f(x2[ii],y2[jj]);
	adata2(ii,jj)=f(x2[ii],y2[jj]);
      }
    }
  
    // x-first

    x0=2.42;
    x1=3.22;
    y0=1.41;
    y1=0.62;

    it.set_data(M,N,x2,y2,data2);
    ait.set_data(M,N,ax2,ay2,adata2);

    for(size_t im=0;im<10;im++) {
      x0*=sqrt(5.0);
      x1*=sqrt(5.0);
      y0*=sqrt(7.0);
      y1*=sqrt(7.0);
      while (x0>x2[M-1]) x0/=2.0;
      while (x1>x2[M-1]) x1/=2.0;
      while (y0>y2[N-1]) y0/=2.0;
      while (y1>y2[N-1]) y1/=2.0;
      t.test_rel(it.eval(x0,y0),ait.eval(x0,y0),
		 1.0e-12,"arma i");
      t.test_rel(it.deriv_x(x0,y0),ait.deriv_x(x0,y0),
		 1.0e-12,"arma dx");
      t.test_rel(it.deriv_y(x0,y0),ait.deriv_y(x0,y0),
		 1.0e-12,"arma dy");
      t.test_rel(it.deriv_xy(x0,y0),ait.deriv_xy(x0,y0),
		 1.0e-12,"arma dyx");
      t.test_abs(it.deriv_xx(x0,y0),ait.deriv_xx(x0,y0),
		 1.0e-12,"arma d2x");
      t.test_abs(it.deriv_yy(x0,y0),ait.deriv_yy(x0,y0),
		 1.0e-12,"arma d2y");
      t.test_rel(it.eval_gen(2,1,x0,x1,y0,y1),ait.eval_gen(2,1,x0,x1,y0,y1),
		 1.0e-12,"arma gen");
    }
  }

#endif

  t.report();
  return 0;
}

