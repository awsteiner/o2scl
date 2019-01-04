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

/* Example: ex_chebapp.cpp
   -------------------------------------------------------------------
*/
  
#include <iostream>
#include <o2scl/constants.h>
#include <o2scl/test_mgr.h>
#include <o2scl/cheb_approx.h>
#include <o2scl/deriv_cern.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;

double func(double x) {
  return sin(1.0/(x+0.08));
}

double dfunc(double x) {
  return -cos(1.0/(x+0.08))/pow(x+0.08,2.0);
}

// Simple function to output information to file for plotting
void write_file(cheb_approx &gc);

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  funct tf=func;

  cheb_approx gc;
  deriv_cern<> cd;
  inte_qag_gsl<> gi;

  double res, err;
  double x0=0.55;

  // Initialize the Chebyshev approximation
  gc.init(func,100,0.0,2.0*o2scl_const::pi);

  // Evaluate the approximation and compare with the exact result
  cout << "f(0.55)" << endl;
  cout << "Exact         : " << func(x0) << endl;
  gc.eval_err(x0,res,err);
  cout << "Approx (n=100): " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << " Act. Error   : " << fabs(res-func(x0)) << endl;

  // Evaluate the approximation at lower order
  gc.eval_n_err(50,x0,res,err);
  cout << "Approx (n=50) : " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << " Act. Error   : " << fabs(res-func(x0)) << endl;
  gc.eval_n_err(25,x0,res,err);
  cout << "Approx (n=25) : " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << " Act. Error   : " << fabs(res-func(x0)) << endl;
  cout << endl;

  t.test_rel(gc.eval(x0),func(x0),1.0e-4,"eval");

  // Show how to use operator=() to create a new approximation
  cheb_approx gc2=gc;
  cout << "Using operator=(): " << gc2.eval(x0) << " " << func(x0) << endl;
  cout << endl;

  t.test_rel(gc2.eval(x0),gc.eval(x0),1.0e-10,"op=");

  // Show how to compute the derivative
  cheb_approx gc_deriv;
  gc.deriv(gc_deriv);

  cout << "f'(0.55)" << endl;
  cout << "Exact         : " << dfunc(x0) << endl;
  gc_deriv.eval_err(x0,res,err);
  cout << "Approx (n=100): " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << " Act. Error   : " << fabs(res-dfunc(x0)) << endl;
  cd.deriv_err(x0,tf,res,err);
  cout << "Direct deriv  : " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << " Act. Error   : " << fabs(res-dfunc(x0)) << endl;
  cout << endl;

  t.test_abs(res,dfunc(x0),1.0e-12,"deriv with deriv_cern");
  t.test_abs(gc_deriv.eval(x0),dfunc(x0),5.0e-3,"deriv with cheb");

  // Show how to compute the integral
  cheb_approx gc_integ;
  gc.integ(gc_integ);

  cout << "int(f,0,0.55)" << endl;
  gc_integ.eval_err(x0,res,err);
  cout << "Approx (n=100): " << res << endl;
  cout << " Est. Error   : " << err << endl;
  gi.integ_err(tf,0.0,x0,res,err);
  cout << "Direct integ  : " << res << endl;
  cout << " Est. Error   : " << err << endl;
  cout << "Rel. Error    : " << fabs(res-gc_integ.eval(x0)) << endl;
  cout << endl;

  t.test_abs(gc_integ.eval(x0),gi.integ(tf,0.0,x0),1.0e-6,"integral");

  write_file(gc);
  
  t.report();
  return 0;
}
// End of example

// Simple function to output information to file for plotting
void write_file(cheb_approx &gc) {

  ofstream fout;
  fout.open("ex_chebapp.out");
  fout.setf(ios::scientific);
  
  for(double x=0.0;x<1.0001;x+=0.01) {
    fout << x << " " << func(x) << " " << gc.eval(x) << " "
	 << gc.eval_n(50,x) << " " << gc.eval_n(25,x) << endl;
  }

  fout.close();

  return;
}
