/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/gsl_astep_old.h>

using namespace std;
using namespace o2scl;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  gsl_astep_old<ode_funct_fptr<>,ubvector,ubvector,ubvector_alloc> ga;
  double x, dx;
  ubvector y(1),dydx(1),yerr(1),yout(1),dydx_out(1);
  int i;
  int vpx=0;
  
  ode_funct_fptr<ubvector> od(derivs);
  
  cout.setf(ios::scientific);

  // Test astep_derivs()
  x=1.0;
  y[0]=1.0;
  dx=0.4;
  derivs(x,1,y,dydx);
  for(i=1;i<=40;i++) {
    ga.astep(x,5.0,dx,1,y,dydx,yerr,od);
    t.test_rel(y[0],exp(x-1.0),6.0e-6,"y_calculated-y_exact");
    if (i%10==0 || x>5.0-1.0e-4) {
      cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
	   << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
    }
    if (x>5.0-1.0e-4) i=50;
  }

  // Test astep_full()
  x=1.0;
  y[0]=1.0;
  dx=0.4;
  derivs(x,1,y,dydx);
  for(i=1;i<=40;i++) {
    ga.astep_full(x,5.0,x,dx,1,y,dydx,yout,yerr,dydx_out,od);
    y[0]=yout[0];
    dydx[0]=dydx_out[0];
    t.test_rel(y[0],exp(x-1.0),6.0e-6,"y_calculated-y_exact");
    if (i%10==0 || x>5.0-1.0e-4) {
      cout << x << " " << y[0] << " " << exp(x-1.0) << " " 
	   << fabs(y[0]-exp(x-1.0))/exp(x-1.0) << endl;
    }
    if (x>5.0-1.0e-4) i=50;
  }

  t.report();

  return 0;
}
