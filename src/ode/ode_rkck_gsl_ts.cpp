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

#include <gsl/gsl_odeiv.h>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_rkck_gsl.h>

typedef boost::numeric::ublas::vector<double> ubvector;

int expon(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  return 0;
}

int exponx(double x, size_t nv, 
	   const std::vector<double> &y,
	   boost::numeric::ublas::vector<double> &dydx) {
  dydx[0]=y[0];
  return 0;
}

int expon_gsl(double x, const double y[], double dydx[], 
	      void *params) {
  dydx[0]=y[0];
  return 0;
}

int expon2(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[0];
  dydx[1]=y[1];
  return 0;
}

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  {

    double x, dx=1.0e-1, x_gsl;

    ubvector y(2), dydx(2), yout(2), yerr(2), dydx_out(2);
    double y_gsl[1], yerr_gsl[1], dydx_gsl[1], dydx_out_gsl[1];

    // We have to keep the full type specification to specify
    // that we want ode_funct_fptr and not ode_funct11 independent
    // of whether or not O2SCL_CPP11 is defined
    ode_rkck_gsl<ubvector,ubvector,ubvector,
		 ode_funct<ubvector,ubvector> > rk;

    // GSL setup
    const gsl_odeiv_step_type *T=gsl_odeiv_step_rkck;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,1);
    gsl_odeiv_system sys={expon_gsl,0,1,0};

    // O2scl function setup
    ode_funct_fptr<ubvector> od(expon);

    x=1.0;
    y[0]=1.0;

    x_gsl=1.0;
    y_gsl[0]=1.0;

    cout << "x            y(calc)      y(exact)     diff         yerr"
	 << "          y'" << endl;
    cout << 1.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << -0.0 
	 << " " << 1.0 << endl;

    // Compute initial derivative
    expon(x,1,y,dydx);
    expon_gsl(x_gsl,y_gsl,dydx_gsl,0);

    for(size_t i=1;i<=40;i++) {

      // Make a step
      rk.step(x,dx,1,y,dydx,y,yerr,dydx,od);

      // Output results
      cout << x+dx << " " << y[0] << " " << exp(x+dx)/exp(1.0) << " "
	   << fabs(y[0]-exp(x+dx)/exp(1.0)) << " " << yerr[0] << " ";
      cout << dydx[0] << endl;

      // Test results
      t.test_rel(y[0],exp(x+dx)/exp(1.0),1.0e-6,"y_calculated-y_exact");
      t.test_abs(yerr[0],0.0,1.0e-6,"y_err");

      // Perform the same step with GSL to compare
      gsl_odeiv_step_apply(s,x_gsl,dx,y_gsl,yerr_gsl,
			   dydx_gsl,dydx_out_gsl,&sys);
      dydx_gsl[0]=dydx_out_gsl[0];

      // Test comparison
      t.test_rel(y_gsl[0],y[0],1.0e-10,"GSL vs. O2scl");
      t.test_rel(yerr_gsl[0],yerr[0],1.0e-9,"GSL vs. O2scl err");
      t.test_rel(dydx_gsl[0],dydx[0],1.0e-9,"GSL vs. O2scl deriv");

      x+=dx;
      x_gsl+=dx;

    }

  }

  // Show that we can use combinations of different vector
  // types
  {
    
    double x, dx=1.0e-1;

    typedef std::vector<double> vec1_t;
    typedef boost::numeric::ublas::vector<double> vec2_t;
    typedef boost::numeric::ublas::bounded_array<double,1> vec3_t;

    vec1_t y(1);
    vec2_t dydx(1);
    vec1_t yout(1);
    vec3_t yerr(1);
    vec2_t dydx_out(1);

    // We have to keep the full type specification to specify
    // that we want ode_funct and not ode_funct11 independent
    // of whether or not O2SCL_CPP11 is defined
    ode_rkck_gsl<vec1_t,vec2_t,vec3_t,ode_funct<vec1_t,vec2_t> > rk;

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    ode_funct_fptr<vec1_t,vec2_t> od(exponx);

    x=1.0;
    y[0]=1.0;
    cout << "x            y(calc)      y(exact)     diff         yerr"
	 << "          y'" << endl;
    cout << 1.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << -0.0 
	 << " " << 1.0 << endl;
    exponx(x,1,y,dydx);
    for(size_t i=1;i<=40;i++) {
      rk.step(x,dx,1,y,dydx,y,yerr,dydx,od);
      t.test_rel(y[0],exp(x+dx)/exp(1.0),1.0e-6,"y_calculated-y_exact 3 types");
      t.test_abs(yerr[0],0.0,1.0e-6,"y_err 3 types");
      cout << x+dx << " " << y[0] << " " << exp(x+dx)/exp(1.0) << " "
	   << fabs(y[0]-exp(x+dx)/exp(1.0)) << " " << yerr[0] << " ";
      x+=dx;
      cout << dydx[0] << endl;
    }

#ifndef O2SCL_NO_CPP11
    
    ode_rkck_gsl<vec1_t,vec2_t,vec3_t,ode_funct11<vec1_t,vec2_t> > rk11;
    ode_funct11<vec1_t,vec2_t> f11=exponx;

    x=1.0;
    y[0]=1.0;
    cout << "x            y(calc)      y(exact)     diff         yerr"
	 << "          y'" << endl;
    cout << 1.0 << " " << 1.0 << " " << 1.0 << " " << 0.0 << " " << -0.0 
	 << " " << 1.0 << endl;
    exponx(x,1,y,dydx);
    for(size_t i=1;i<=40;i++) {
      rk11.step(x,dx,1,y,dydx,y,yerr,dydx,f11);
      t.test_rel(y[0],exp(x+dx)/exp(1.0),1.0e-6,"y_calculated-y_exact c++11");
      t.test_abs(yerr[0],0.0,1.0e-6,"y_err c++11");
      cout << x+dx << " " << y[0] << " " << exp(x+dx)/exp(1.0) << " "
	   << fabs(y[0]-exp(x+dx)/exp(1.0)) << " " << yerr[0] << " ";
      x+=dx;
      cout << dydx[0] << endl;
    }
    
#endif

  }

  t.report();


  return 0;
}
