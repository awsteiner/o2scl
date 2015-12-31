/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/ode_funct.h>
#include <o2scl/ode_bv_mshoot.h>

using namespace std;
using namespace o2scl;

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

#ifdef O2SCL_NEVER_DEFINED
  
  ode_bv_mshoot<> obs;
  
  ode_funct_fptr<> off(derivs);
  ubvector xbound(5);
  ubmatrix ybound(5,2);
  ubvector_int index(2);
  double h=0.05;
  xbound[0]=0.0;
  xbound[1]=0.25;
  xbound[2]=0.5;
  xbound[3]=0.75;
  xbound[4]=1.0;
  ybound[0][0]=1.0;
  ybound[0][1]=6.0;
  ybound[1][0]=2.0;
  ybound[1][1]=5.0;
  ybound[2][0]=3.0;
  ybound[2][1]=4.0;
  ybound[3][0]=4.0;
  ybound[3][1]=3.0;
  // The initial value of ybound[4][0] is ignored so we leave it at zero
  ybound[4][0]=0.0;
  ybound[4][1]=2.0;
  index[0]=ode_bv_solve<int>::left;
  index[1]=ode_bv_solve<int>::right;

  obs.solve_final_value(h,2,5,xbound,ybound,index,off);
  
  ybound[0][0]=1.0;
  ybound[0][1]=6.0;
  ybound[1][0]=2.0;
  ybound[1][1]=5.0;
  ybound[2][0]=3.0;
  ybound[2][1]=4.0;
  ybound[3][0]=4.0;
  ybound[3][1]=3.0;
  // The initial value of ybound[4][0] is ignored so we leave it at zero
  ybound[4][0]=0.0;
  ybound[4][1]=2.0;

  size_t n_sol=101;
  ubvector x_sol(n_sol);
  ubmatrix y_sol(n_sol,2), dydx_sol(n_sol,2), yerr_sol(n_sol,2);
  
  obs.solve_store<ubmatrix_row>(h,2,5,xbound,ybound,index,n_sol,x_sol,
			       y_sol,dydx_sol,yerr_sol,off);

  double tol=1.0e-7, y, dydx;

  cout << "Full solution" << endl;
  for(size_t i=0;i<n_sol;i++) {
    cout << x_sol[i] << " " << y_sol[i][0] << " " << y_sol[i][1] << " ";
    cout.setf(ios::showpos);
    cout << yerr_sol[i][0] << " " << yerr_sol[i][1] << endl;
    cout.unsetf(ios::showpos);
    exact_sol(x_sol[i],y,dydx);
    t.test_rel(y_sol[i][0],y,tol,"0");
    t.test_rel(y_sol[i][1],dydx,tol,"1");
    t.test_abs(yerr_sol[i][0],0.0,1.0e-6,"0e");
    t.test_abs(yerr_sol[i][1],0.0,1.0e-6,"1e");
  }

#endif
  
  t.report();
  return 0;
}

