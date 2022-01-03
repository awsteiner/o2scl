/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_rk8pd_gsl.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  return 0;
}

int derivs2(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  dydx[1]=y[1];
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  ubvector y(2), dydx(2), yerr(2);
  
  ode_rk8pd_gsl<> rk;

  ode_funct od=derivs;

  double x=1.0;
  y[0]=1.0;
  double dx=0.1;
  cout << "x            y(calc)      y(exact)     diff         yerr"
       << "          y'" << endl;
  cout << 1.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << -0.0 
       << " " << 1.0 << endl;
  for(size_t i=1;i<=40;i++) {
    derivs(x,1,y,dydx);
    rk.step(x,dx,1,y,dydx,y,yerr,dydx,od);
    t.test_rel(y[0],exp(x+dx)/exp(1.0),1.0e-13,"y_calculated-y_exact");
    t.test_abs(yerr[0],0.0,5.0e-13,"y_err");
    cout << x+dx << " " << y[0] << " " << exp(x+dx)/exp(1.0) << " "
	 << fabs(y[0]-exp(x+dx)/exp(1.0)) << " " << yerr[0] << " ";
    x+=dx;
    // Evaluate the derivative at the new point
    derivs(x,1,y,dydx);
    cout << dydx[0] << endl;
  }
  cout << endl;
  
  t.report();

  return 0;
}
