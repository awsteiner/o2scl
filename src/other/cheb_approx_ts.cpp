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

#include <iostream>
#include <o2scl/constants.h>
#include <o2scl/test_mgr.h>
#include <o2scl/cheb_approx.h>

using namespace std;
using namespace o2scl;

double func(double x) {
  return sin(x);
}

double func2(double x) {
  double a=1.0/(exp(1.0)-1.0);
  double b=1.0/(1.0-exp(1.0));
  return exp(x)*a+b;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cheb_approx gc;

  // Initialize a Chebyshev approximation of the function y=sin(x)
  gc.init(func,10,0.0,2.0*o2scl_const::pi);
  t.test_rel(gc.eval(0.55),sin(0.55),5.0e-2,"cheb1");
  cout << gc.eval(0.55) << " " << sin(0.55) << endl;
  cout << gc.eval_n(5,0.55) << " " << sin(0.55) << endl;
  double res, err;
  gc.eval_err(0.55,res,err);
  cout << res << " " << err << endl;
  gc.eval_n_err(5,0.55,res,err);
  cout << res << " " << err << endl;

  // Try operator=()
  cheb_approx gc2=gc;
  t.test_rel(gc2.eval(0.55),sin(0.55),5.0e-2,"cheb1");
  cout << gc2.eval(0.55) << " " << sin(0.55) << endl;

  // Show how to get the derivative
  cheb_approx gc3;
  gc.deriv(gc3);
  t.test_rel(gc3.eval(0.55),cos(0.55),5.0e-2,"dcheb2");
  cout << gc3.eval(0.55) << " " << cos(0.55) << endl;

  // Show how to get the integral
  cheb_approx gc4;
  gc.integ(gc4);
  t.test_rel(gc4.eval(0.55)-gc4.eval(0.45),-cos(0.55)+cos(0.45),
	     5.0e-2,"icheb2");
  cout << gc4.eval(0.55) << " " << -cos(0.55)+cos(0.45) << endl;

  // Show how to use the eval() function to represent
  // the Chebyshev approximation as a funct object
  //funct_cmfptr<cheb_approx> fmn(&gc,&cheb_approx::eval);
  funct fmn=std::bind(std::mem_fn<double(double) const>
			(&cheb_approx::eval),&gc,
			std::placeholders::_1);
  
  double y;
  for(double x=0.0;x<1.01;x+=0.2) {
    y=fmn(x);
    cout << y << " " << sin(x) << endl;
    t.test_rel(y,sin(x),5.0e-2,"as func");
  }

  // Show that the endpoints are not exact
  gc.init(func2,10,0.0,1.0);
  cout << gc.eval(0.0) << " " << gc.eval(1.0) << endl;

  // Try to make them exact
  for(size_t i=0;i<5;i++) {
    double zval=gc.eval(0.0);
    double oval=gc.eval(1.0);
    cout << gc.get_coefficient(0) << " ";
    gc.set_coefficient(0,gc.get_coefficient(0)-zval);
    cout << gc.eval(0.0) << " " << gc.eval(1.0) << " "
	 << gc.get_coefficient(0) << endl;
  }
  cout << gc.eval(0.0) << " " << gc.eval(1.0) << endl;

  t.report();
  return 0;
}


