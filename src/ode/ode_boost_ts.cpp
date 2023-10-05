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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_odeiv.h>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_rkck_gsl.h>
#include <o2scl/ode_boost.h>

typedef boost::numeric::ublas::vector<double> ubvector;

int expon(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  return 0;
}

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  {
    double x, dx=1.0e-1, x2;
    ode_funct of=expon;

    ubvector y(1), dydx(1), yout(1), yerr(1), dydx_out(1);
    ubvector y2(1), dydx2(1), yout2(1), yerr2(1), dydx_out2(1);

    ode_boost<boost::numeric::odeint::runge_kutta_cash_karp54<ubvector> > ob;
    ode_boost<boost::numeric::odeint::runge_kutta_dopri5<ubvector> > ob2;
    ode_boost<boost::numeric::odeint::runge_kutta_fehlberg78<ubvector> > ob3;
    //ode_boost<boost::numeric::odeint::rosenbrock4<ubvector> > ob4;
    ode_rkck_gsl<> org;

    for(size_t kk=0;kk<3;kk++) {
    
      x=1.0;
      y[0]=1.0;
      x2=1.0;
      y2[0]=1.0;

      cout << "x            y(calc)      y(exact)     diff         yerr"
           << "          y'" << endl;
      cout << 1.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << -0.0 
           << " " << 1.0 << endl;

      // Compute initial derivative
      expon(x,1,y,dydx);
      expon(x2,1,y2,dydx2);

      for(size_t i=1;i<=40;i++) {

        // Make a step
        if (kk==0) {
          ob.step(x,dx,1,y,dydx,y,yerr,dydx,of);
        } else if (kk==1) {
          ob2.step(x,dx,1,y,dydx,y,yerr,dydx,of);
        } else if (kk==2) {
          ob3.step(x,dx,1,y,dydx,y,yerr,dydx,of);
        } else {
          //ob4.step(x,dx,1,y,dydx,y,yerr,dydx,of);
        }

        // Output results
        cout << x+dx << " " << y[0] << " " << exp(x+dx)/exp(1.0) << " "
             << fabs(y[0]-exp(x+dx)/exp(1.0)) << " " << yerr[0] << " ";
        cout << dydx[0] << endl;

        // Make a step
        org.step(x2,dx,1,y2,dydx2,y2,yerr2,dydx2,of);

        // Output results
        cout << x2+dx << " " << y2[0] << " " << exp(x2+dx)/exp(1.0) << " "
             << fabs(y2[0]-exp(x2+dx)/exp(1.0)) << " " << yerr2[0] << " ";
        cout << dydx2[0] << endl;

        x+=dx;
        x2+=dx;

      }

      if (kk==0) {
        t.test_rel(y[0],y2[0],1.0e-15,"Boost vs. GSL");
      } else if (kk==1) {
        t.test_rel(y[0],y2[0],1.0e-7,"Boost vs. GSL");
      } else if (kk==2) {
        t.test_rel(y[0],y2[0],1.0e-7,"Boost vs. GSL");
      } else {
        t.test_rel(y[0],y2[0],1.0e-8,"Boost vs. GSL");
      }

    }

  }
  
  t.report();


  return 0;
}
