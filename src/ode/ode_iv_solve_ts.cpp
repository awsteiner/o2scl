/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#include <o2scl/vector.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/test_mgr.h>
#include <o2scl/ode_funct.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[1];
  dydx[1]=-y[0];
  return 0;
}

int exact_sol(double x, double &y, double &dydx, double &d2ydx2) {
  y=cos(x)+(2.0/cos(1.0)+tan(1.0))*sin(x);
  dydx=-sin(x)+(2.0/cos(1.0)+tan(1.0))*cos(x);
  d2ydx2=-cos(x)-(2.0/cos(1.0)+tan(1.0))*sin(x);
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  ode_iv_solve<> ivs;
  ode_iv_solve_grid<> ivsg;

  ode_funct od=derivs;

  ubvector y(2), dydx(2), yout(2), yerr(2), yend(2);

  double ey, edydx, ed2ydx2;

  cout << endl;
  
  // ------------------------------------------------

  cout << "Final value: " << endl;

  // Try enough values of 'h' so that we test both 
  // an even and an odd number of steps
  for(double h=0.1;h<=0.2;h*=1.1) {
    y[0]=1.0;
    y[1]=(2.0/cos(1.0)+tan(1.0));
    
    ivs.solve_final_value(0.0,1.0,h,2,y,yend,yerr,dydx,od);
    exact_sol(1.0,ey,edydx,ed2ydx2);
    
    t.test_rel(yend[0],ey,1.0e-7,"ey");
    t.test_rel(yend[1],edydx,1.0e-7,"edydx");
    t.test_rel(dydx[0],edydx,1.0e-7,"edydx");
    t.test_rel(dydx[1],ed2ydx2,1.0e-7,"ed2ydx2");
    cout << h << " " << ivs.nsteps << endl;
  }
  cout << endl;
  
  // ------------------------------------------------

  cout << "Final value, h<0: " << endl;

  for(double h=-0.1;h>=-0.2;h*=1.1) {
    exact_sol(1.0,ey,edydx,ed2ydx2);
    y[0]=ey;
    y[1]=edydx;
    
    ivs.solve_final_value(1.0,0.0,h,2,y,yend,yerr,dydx,od);
    exact_sol(0.0,ey,edydx,ed2ydx2);
    
    t.test_rel(yend[0],ey,1.0e-7,"ey");
    t.test_rel(yend[1],edydx,1.0e-7,"edydx");
    t.test_rel(dydx[0],edydx,1.0e-7,"edydx");
    t.test_rel(dydx[1],ed2ydx2,1.0e-7,"ed2ydx2");
    cout << h << " " << ivs.nsteps << endl;
  }
  cout << endl;
  
  // ------------------------------------------------
  // Two-equation test of solve_grid

  {
    ubvector xgrid;
    ubmatrix ygrid, dydxgrid, err_grid;
    size_t ngrid;

    cout << "Grid: " << endl;
  
    ngrid=11;
    xgrid.resize(ngrid);
    ygrid.resize(ngrid,2);
    dydxgrid.resize(ngrid,2);
    err_grid.resize(ngrid,2);
    vector_grid(uniform_grid_end<>(0.0,1.0,ngrid-1),xgrid);
    ygrid(0,0)=1.0;
    ygrid(0,1)=(2.0/cos(1.0)+tan(1.0));
  
    ivsg.solve_grid<ubmatrix>(1.0,2,ngrid,xgrid,ygrid,err_grid,dydxgrid,od);
  
    for(size_t i=0;i<ngrid;i++) {
      exact_sol(xgrid[i],ey,edydx,ed2ydx2);
      cout << xgrid[i] << " " << ygrid(i,0) << " " << ey << " "
	   << dydxgrid(i,0) << " " << edydx << endl;
      t.test_rel(ygrid(i,0),ey,1.0e-8,"y g2");
      t.test_rel(ygrid(i,1),edydx,1.0e-8,"y g2");
      t.test_rel(dydxgrid(i,0),edydx,1.0e-8,"y g2");
      t.test_rel(dydxgrid(i,1),ed2ydx2,1.0e-8,"y g2");
    }
  
    cout << endl;
  }

  // ------------------------------------------------
  // Test solve_store

  ubvector xgrid;
  ubmatrix ygrid, dydxgrid, err_grid;
  size_t ngrid;

  ngrid=15;

  cout << "Store:" << endl;
  xgrid.resize(ngrid);
  ygrid.resize(ngrid,2);
  dydxgrid.resize(ngrid,2);
  err_grid.resize(ngrid,2);

  ygrid(0,0)=1.0;
  ygrid(0,1)=(2.0/cos(1.0)+tan(1.0));
  
  ivs.solve_store<ubmatrix>
    (0.0,1.0,0.1,2,ngrid,xgrid,ygrid,err_grid,dydxgrid,od);
  
  for(size_t i=0;i<ngrid;i++) { 
    exact_sol(xgrid[i],ey,edydx,ed2ydx2);
    cout << xgrid[i] << " " << ygrid(i,0) << " " << ey << " "
	 << dydxgrid(i,0) << " " << edydx << endl;
    t.test_rel(ygrid(i,0),ey,2.0e-7,"y g2");
    t.test_rel(ygrid(i,1),edydx,2.0e-7,"y g2");
    t.test_rel(dydxgrid(i,0),edydx,2.0e-7,"y g2");
    t.test_rel(dydxgrid(i,1),ed2ydx2,2.0e-7,"y g2");
  }
  size_t ns1=ivs.nsteps;

  // Compare with solve_final_value

  y[0]=1.0;
  y[1]=(2.0/cos(1.0)+tan(1.0));
  ivs.solve_final_value(0.0,1.0,0.1,2,y,yend,yerr,dydx,od);

  t.test_gen(ns1==ivs.nsteps,"0");
  t.test_rel(yend[0],ygrid(ngrid-1,0),1.0e-14,"1");
  t.test_rel(yend[1],ygrid(ngrid-1,1),1.0e-14,"2");
  t.test_rel(dydx[0],dydxgrid(ngrid-1,0),1.0e-14,"3");
  t.test_rel(dydx[1],dydxgrid(ngrid-1,1),1.0e-14,"4");

  cout << endl;
  
  // ------------------------------------------------
  // Test solve_store with non-zero istart

  cout << "Store with istart: " << endl;

  // Variable ngrid was modified by solve_store() above, 
  // so we return it to 15.
  ngrid=15;
  size_t ix=5;

  ygrid(ix,0)=1.0;
  ygrid(ix,1)=(2.0/cos(1.0)+tan(1.0));
  
  ivs.solve_store<ubmatrix>
    (0.0,1.0,0.1,2,ngrid,xgrid,ygrid,err_grid,dydxgrid,od,ix);
  
  for(size_t i=ix;i<ngrid;i++) { 
    exact_sol(xgrid[i],ey,edydx,ed2ydx2);
    cout << xgrid[i] << " " << ygrid(i,0) << " " << ey << " "
	 << dydxgrid(i,0) << " " << edydx << endl;
    t.test_rel(ygrid(i,0),ey,2.0e-7,"y g2");
    t.test_rel(ygrid(i,1),edydx,2.0e-7,"y g2");
    t.test_rel(dydxgrid(i,0),edydx,2.0e-7,"y g2");
    t.test_rel(dydxgrid(i,1),ed2ydx2,2.0e-7,"y g2");
  }

  // Compare with solve_final_value

  t.test_gen(ns1==ivs.nsteps,"0");
  t.test_rel(yend[0],ygrid(ngrid-1,0),1.0e-14,"1");
  t.test_rel(yend[1],ygrid(ngrid-1,1),1.0e-14,"2");
  t.test_rel(dydx[0],dydxgrid(ngrid-1,0),1.0e-14,"3");
  t.test_rel(dydx[1],dydxgrid(ngrid-1,1),1.0e-14,"4");

  cout << endl;

  // ------------------------------------------------

  t.report();
  return 0;
}

