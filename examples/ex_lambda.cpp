/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2020, Andrew W. Steiner
  
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

// sphinx-example-start
/* Example: ex_lambda.cpp
   -------------------------------------------------------------------
   Demonstrate how to use standard library and lambda function objects
   with O2scl.
 
*/
#include <iostream>
#include <functional>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

// A global function
double gfn(double x) {
  return sin(x)-0.1;
}

// A global function with a parameter
double gfn_param(double x, double a) {
  return sin(x)-a;
}

class a_class {
public:
  // A member function
  double mfn(double x) {
    return sin(x)-0.1;
  }
  // A member function with a parameter
  double mfn_param(double x, double a) {
    return sin(x)-a;
  }
  // A member function with a const and non-const reference parameter
  double mfn_param_ref(double x, const double &a, double &b) {
    b=x*x;
    return sin(x)-a;
  }
};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  // The O2scl solver. Note that we use the same solver for 
  // all the examples below.
  root_brent_gsl<std::function<double(double)> > grb;

  // For the initial bracket 
  double a, b;

  // With a global function
  {
    a=-0.9, b=0.9;
    std::function<double(double)> f=gfn;
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function");
  }

  // O2scl defines 'funct' as std::function<double(double)>, so this
  // shorter notation also works. 
  {
    a=-0.9, b=0.9;
    funct f=gfn;
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function (take 2)");
  }

  // With a global function
  {
    a=-0.9, b=0.9;
    funct f=std::bind(gfn_param,std::placeholders::_1,0.1);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function with parameter");
  }

  // With a member function
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double)>(&a_class::mfn),
		      ac,std::placeholders::_1);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Member function");
  }

  // With a member function which has a fixed parameter
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double,double)>
		      (&a_class::mfn_param),
		      ac,std::placeholders::_1,0.1);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Member function with parameter");
  }

  // Inline specification of the function
  {
    a=-0.9, b=0.9;
    funct f=[](double x) -> double { double z=sin(x)-0.1; return z; };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Inline 1");
  }

  // A bit of a shorter notation 
  {
    a=-0.9, b=0.9;
    funct f=[](double x){ return sin(x)-0.1; };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Inline 2");
    // Show copy construction works
    a=-0.9, b=0.9;
    funct f2=f;
    grb.solve_bkt(a,b,f2);
    t.test_rel(a,asin(0.1),1.0e-12,"Inline 3");
  }

  // With a member function and reference parameters. Note that
  // we use std::cref for the const reference and std::ref
  // for the non-const reference
  {
    a=-0.9, b=0.9;
    a_class ac;
    double param1=0.1, param2;
    funct f=std::bind(std::mem_fn<double(double,const double &,double &)>
		      (&a_class::mfn_param_ref),
		      ac,std::placeholders::_1,std::cref(param1),
		      std::ref(param2));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Member function with references");
    // The last function call isn't always at the final root!
    t.test_rel(param2,0.01,1.0e-2,"Reference 2");
  }

  t.report();
  return 0;
}
