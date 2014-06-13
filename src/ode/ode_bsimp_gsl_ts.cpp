/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_rkck_gsl.h>
#include <o2scl/ode_bsimp_gsl.h>
#include <gsl/gsl_odeiv.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=-y[1]+x*y[0]/10.0;
  dydx[1]=y[0]+x;
  return 0;
}

int jac(double x, size_t nv, const ubvector &y, 
	ubmatrix &dfdy, ubvector &dfdx) {
  dfdy(0,0)=x/10.0;
  dfdy(0,1)=-1.0;
  dfdy(1,0)=1.0;
  dfdy(1,1)=0.0;
  dfdx[0]=y[0]/10.0;
  dfdx[1]=1.0;
  return 0;
}

int derivs_gsl(double x, const double y[], double dydx[], void *pa) {
  dydx[0]=-y[1]+x*y[0]/10.0;
  dydx[1]=y[0]+x;
  return 0;
}

int jac_gsl(double x, const double y[], double *dfdy,
	    double dfdt[], void *pa) {
  dfdy[0]=x/10.0;
  dfdy[1]=-1.0;
  dfdy[2]=1.0;
  dfdy[3]=0.0;
  dfdt[0]=y[0]/10.0;
  dfdt[1]=1.0;
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);
  cout.precision(4);
  
  test_mgr t;
  t.set_output_level(1);

  double x1, x2, dx=1.0e-1, x3;
  ubvector y1(2), dydx1(2), yout1(2), yerr1(2), dydx_out1(2);
  ubvector y2(2), dydx2(2), yout2(2), yerr2(2), dydx_out2(2);
  double y3[2], dydx3[2], yout3[2], yerr3[2], dydx_out3[2];
  
  ode_rkck_gsl<ubvector,ubvector,ubvector,ode_funct11> rk;
  ode_bsimp_gsl<ode_funct11,ode_jac_funct11> gb;

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  
  gsl_odeiv_step *s=gsl_odeiv_step_alloc(gsl_odeiv_step_bsimp,2);
  
  ode_funct11 od=derivs;
  ode_jac_funct11 oj=jac;
  gsl_odeiv_system sys={derivs_gsl,jac_gsl,2,0};
  
  x1=1.0;
  x2=1.0;
  x3=1.0;
  y1[0]=1.0;
  y1[1]=-1.0;
  y2[0]=1.0;
  y2[1]=-1.0;
  y3[0]=1.0;
  y3[1]=-1.0;

  derivs(x1,2,y1,dydx1);
  derivs(x2,2,y2,dydx2);
  derivs_gsl(x3,y3,dydx3,0);

  //gb.verbose=2;

  for(size_t i=1;i<=40;i++) {

    rk.step(x1,dx,2,y1,dydx1,y1,yerr1,dydx1,od);
    gb.step(x2,dx,2,y2,dydx2,y2,yerr2,dydx2,od,oj);
    gsl_odeiv_step_apply(s,x3,dx,y3,yerr3,dydx3,dydx3,&sys);

    if (fabs(y1[0])>4.0e-2) {
      t.test_rel(y1[0],y2[0],1.0e-1,"rk gb 1a");
    } else {
      t.test_rel(y1[0],y2[0],8.0e-1,"rk gb 1b");
    }
    if (fabs(y1[1])>1.5e-2) {
      t.test_rel(y1[1],y2[1],1.0e-1,"rk gb 2a");
    } else {
      t.test_rel(y1[1],y2[1],8.0e-1,"rk gb 2b");
    }

    t.test_rel(y3[0],y2[0],1.0e-10,"gsl gb 3");
    t.test_rel(y3[1],y2[1],1.0e-10,"gsl gb 4");

  }

  t.report();

  return 0;
}
