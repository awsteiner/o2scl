/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

/* Example: ex_ode_it.cpp
   -------------------------------------------------------------------
   Demonstrate the iterative method for solving ODEs
*/

#include <o2scl/ode_it_solve.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/linear_solver.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;

// The three-dimensional ODE system
class ode_system {
  
public:
  
  // The LHS boundary conditions
  double left(size_t ieq, double x, ovector_base &yleft) {
    if (ieq==0) return yleft[0]-1.0;
    return yleft[1]*yleft[1]+yleft[2]*yleft[2]-2.0;
  }
  
  // The RHS boundary conditions
  double right(size_t ieq, double x, ovector_base &yright) {
    return yright[1]-3.0;
  }
  
  // The differential equations
  double derivs(size_t ieq, double x, ovector_base &y) {
    if (ieq==1) return y[0]+y[1];
    else if (ieq==2) return y[0]+y[2];
    return y[1];
  }

  // This is the alternative specification for ode_iv_solve for
  // comparison
  int shoot_derivs(double x, size_t nv, const ovector_base &y,
		   ovector_base &dydx) {
    dydx[0]=y[1];
    dydx[1]=y[0]+y[1];
    dydx[2]=y[0]+y[2];
    return 0;
  }
  
};

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // The ODE solver
  ode_it_solve<ode_it_funct<ovector_base>,ovector_base,
    omatrix_base,omatrix_row,ovector_base,omatrix_base> oit;

  // The class which contains the functions to solve
  ode_system os;

  // Make function objects for the derivatives and boundary conditions
  ode_it_funct_mfptr<ode_system> f_d(&os,&ode_system::derivs);
  ode_it_funct_mfptr<ode_system> f_l(&os,&ode_system::left);
  ode_it_funct_mfptr<ode_system> f_r(&os,&ode_system::right);

  // Grid size
  size_t ngrid=40;

  // Number of ODEs
  size_t neq=3;

  // Number of LHS boundary conditions
  size_t nbleft=2;

  // Create space for the solution and make an initial guess
  ovector x(ngrid);
  omatrix y(ngrid,neq);
  for(size_t i=0;i<ngrid;i++) {
    x[i]=((double)i)/((double)(ngrid-1));
    y[i][0]=1.0+x[i]+1.0;
    y[i][1]=3.0*x[i];
    y[i][2]=-0.1*x[i]-1.4;
  }

  int pa=0;

  // Workspace objects
  omatrix A(ngrid*neq,ngrid*neq);
  ovector rhs(ngrid*neq), dy(ngrid*neq);

  // Perform the solution
  oit.verbose=1;
  oit.solve(ngrid,neq,nbleft,x,y,f_d,f_l,f_r,A,rhs,dy);
  
  // Compare with the initial value solver ode_iv_solve
  ode_iv_solve<> ois;
  ode_funct_mfptr<ode_system> f_sd(&os,&ode_system::shoot_derivs);
  ovector ystart(neq), yend(neq);
  for(size_t i=0;i<neq;i++) ystart[i]=y[0][i];
  ois.solve_final_value(0.0,1.0,0.01,neq,ystart,yend,f_sd);

  // Test the result
  t.test_rel(y[0][0],1.0,1.0e-3,"ya");
  t.test_rel(y[ngrid-1][0],yend[0],1.0e-3,"yb");
  t.test_rel(y[0][1],0.25951,1.0e-3,"yc");
  t.test_rel(y[ngrid-1][1],yend[1],1.0e-3,"yd");
  t.test_rel(y[0][2],-1.3902,1.0e-3,"ye");
  t.test_rel(y[ngrid-1][2],yend[2],1.0e-3,"yf");
  
  t.report();

  return 0;
}
// End of example

