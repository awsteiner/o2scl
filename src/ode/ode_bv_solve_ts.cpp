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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/ode_bv_solve.h>
#include <o2scl/ode_funct.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector<int> ubvector_int;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

/* As an example, we solve y''[x]==-y[x] with the boundary
   conditions y[0]==1 and y'[1]==2. The solution is
   y=cos(x)+(2.0*sec(1)+tan(1))*sin(x)
*/
 
int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[1];
  dydx[1]=-y[0];
  return 0;
}

int exact_sol(double x, double &y, double &dydx) {
  y=cos(x)+(2.0/cos(1.0)+tan(1.0))*sin(x);
  dydx=-sin(x)+(2.0/cos(1.0)+tan(1.0))*cos(x);
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  ode_bv_shoot<> obs;
  ode_bv_shoot_grid<> obsg;
  
  ode_funct off=derivs;
  double x0=0.0, x1=1.0, h=(x1-x0)/20.0;
  ubvector ystart(2), yend(2), yerr(2), dydx_out(2);
  ubvector_int index(2);
  index[0]=ode_bv_solve::left;
  index[1]=ode_bv_solve::right;

  size_t n_sol=100;
  ubvector x_sol(n_sol);
  ubmatrix y_sol(n_sol,2), dydx_sol(n_sol,2), yerr_sol(n_sol,2);
  double tol=1.0e-7;
  double y, dydx;

  // ----------------------------------------------
  // Test solve_final_value()

  ystart[0]=1.0;
  ystart[1]=0.0;
  yend[0]=0.0;
  yend[1]=2.0;

  obs.solve_final_value(x0,x1,h,2,ystart,yend,index,yerr,dydx_out,off);
  
  cout << "Boundaries from shooting: " << endl;
  cout << ystart[0] << " " << ystart[1] << endl;
  cout << yend[0] << " " << yend[1] << endl;
  cout << endl;

  cout << "Boundaries from exact solution: " << endl;
  exact_sol(x0,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(ystart[1],dydx,tol,"ystart[1]");
  exact_sol(x1,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(yend[0],y,tol,"yend[0]");
  t.test_rel(yend[1],dydx,tol,"yend[1]");
  cout << endl;

  // ----------------------------------------------
  // Test solve_store()

  ystart[0]=1.0;
  ystart[1]=0.0;
  yend[0]=0.0;
  yend[1]=2.0;

  obs.solve_store<ubmatrix,ubmatrix_row>
    (x0,x1,h,2,ystart,yend,index,n_sol,x_sol,y_sol,yerr_sol,dydx_sol,off);
  
  cout << "Boundaries from shooting: " << endl;
  cout << ystart[0] << " " << ystart[1] << endl;
  cout << yend[0] << " " << yend[1] << endl;
  cout << endl;

  cout << "Boundaries from exact solution: " << endl;
  exact_sol(x0,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(ystart[1],dydx,tol,"ystart[1]");
  exact_sol(x1,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(yend[0],y,tol,"yend[0]");
  t.test_rel(yend[1],dydx,tol,"yend[1]");
  cout << endl;

  cout << "Stored solution" << endl;
  for(size_t i=0;i<n_sol;i++) {
    cout << x_sol[i] << " " << y_sol(i,0) << " " << y_sol(i,1) << " ";
    cout.setf(ios::showpos);
    cout << yerr_sol(i,0) << " " << yerr_sol(i,1) << endl;
    cout.unsetf(ios::showpos);
    exact_sol(x_sol[i],y,dydx);
    t.test_rel(y_sol(i,0),y,tol,"0");
    t.test_rel(y_sol(i,1),dydx,tol,"1");
    t.test_abs(yerr_sol(i,0),0.0,1.0e-6,"0e");
    t.test_abs(yerr_sol(i,1),0.0,1.0e-6,"1e");
  }
  cout << endl;

  // ----------------------------------------------
  // Test solve_grid()

  n_sol=20;
  for(size_t i=0;i<n_sol;i++) {
    x_sol[i]=x0+((double)i)/((double)(n_sol-1))*(x1-x0);
  }

  // Set boundaries and initial guess
  ystart[0]=1.0;
  ystart[1]=0.0;
  yend[0]=0.0;
  yend[1]=2.0;

  obsg.solve_grid(x0,x1,h,2,ystart,yend,index,n_sol,x_sol,y_sol,
		  yerr_sol,dydx_sol,off);

  ubmatrix_row y0(y_sol,0);
  ubmatrix_row dydx0(dydx_sol,0);

  cout << "Boundaries from exact solution: " << endl;
  exact_sol(x0,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(y_sol(0,0),y,tol,"ystart[1]");
  t.test_rel(dydx_sol(0,0),dydx,tol,"ystart[1]");
  exact_sol(x1,y,dydx);
  cout << y << " " << dydx << endl;
  t.test_rel(y_sol(n_sol-1,0),y,tol,"ystart[1]");
  t.test_rel(dydx_sol(n_sol-1,0),dydx,tol,"ystart[1]");
  cout << endl;

  cout << "Stored solution" << endl;
  for(size_t i=0;i<n_sol;i++) {
    cout << x_sol[i] << " " << y_sol(i,0) << " " << y_sol(i,1) << " ";
    cout.setf(ios::showpos);
    cout << yerr_sol(i,0) << " " << yerr_sol(i,1) << endl;
    cout.unsetf(ios::showpos);
    exact_sol(x_sol[i],y,dydx);
    t.test_rel(y_sol(i,0),y,tol,"0");
    t.test_rel(y_sol(i,1),dydx,tol,"1");
    t.test_abs(yerr_sol(i,0),0.0,1.0e-6,"0e");
    t.test_abs(yerr_sol(i,1),0.0,1.0e-6,"1e");
  }

  t.report();
  return 0;
}

