/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/test_mgr.h>
#include <o2scl/ode_funct.h>
#include <o2scl/astep_boost.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  double x, dx;
  ubvector y(1), yerr(1), dydx(1), yout(1), dydx_out(1);
  int i;
  
  astep_boost<boost::numeric::odeint::runge_kutta_cash_karp54<ubvector>> na;

  ode_funct od=derivs;
  
  // Test astep()
  x=1.0;
  y[0]=1.0;
  dx=0.1;
  cout << "x            y(calc)      y(exact)     diff" << endl;
  cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
       << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  derivs(x,1,y,dydx);
  for(i=1;i<=40;i++) {
    na.astep(x,5.0,dx,1,y,dydx,yerr,od);
    t.test_rel(y[0],exp(x-1.0),1.0e-8,"y_calculated-y_exact");
    cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
	 << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  }

  ode_funct od11=derivs;
  astep_boost<boost::numeric::odeint::runge_kutta_cash_karp54<ubvector>> na11;

  // Test astep()
  x=1.0;
  y[0]=1.0;
  dx=0.1;
  cout << "x            y(calc)      y(exact)     diff" << endl;
  cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
       << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  derivs(x,1,y,dydx);
  for(i=1;i<=40;i++) {
    na11.astep(x,5.0,dx,1,y,dydx,yerr,od11);
    t.test_rel(y[0],exp(x-1.0),1.0e-8,"y_calculated-y_exact");
    cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
	 << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  }

  astep_boost<boost::numeric::odeint::runge_kutta_cash_karp54<ubvector>> na12;

  // Test astep_full()
  x=1.0;
  y[0]=1.0;
  dx=0.1;
  cout << "x            y(calc)      y(exact)     diff" << endl;
  cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
       << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  derivs(x,1,y,dydx);
  for(i=1;i<=40;i++) {
    na12.astep_full(x,5.0,x,dx,1,y,dydx,yout,yerr,dydx_out,od);
    y[0]=yout[0];
    dydx[0]=dydx_out[0];
    t.test_rel(y[0],exp(x-1.0),1.0e-8,"y_calculated-y_exact");
    cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
	 << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
  }

  t.report();

  return 0;
}
