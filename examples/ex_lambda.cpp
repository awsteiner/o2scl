/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Andrew W. Steiner
  
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

// A global function with a value parameter
double gfn_param(double x, double p) {
  return sin(x)-p;
}

// A global function with a reference parameter
double gfn_rparam(double x, double &p) {
  return sin(x)-p;
}

// A global function with a const reference parameter
double gfn_crparam(double x, const double &p) {
  return sin(x)-p;
}

class a_class {
  
public:
  
  // A member function
  double mfn(double x) {
    return sin(x)-0.1;
  }

  // A const member function
  double cmfn(double x) const {
    return sin(x)-0.1;
  }
  
  // A member function with a value parameter
  double mfn_param(double x, double p) {
    return sin(x)-p;
  }
  
  // A member function with a reference parameter
  double mfn_rparam(double x, double &p) {
    return sin(x)-p;
  }
  
  // A const member function with a reference parameter
  double cmfn_rparam(double x, double &p) const {
    return sin(x)-p;
  }
  
  // A const member function with a const reference parameter
  double cmfn_crparam(double x, const double &p) const {
    return sin(x)-p;
  }
  
};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  // The O2scl solver. Note that we use the same solver for 
  // all the examples below.
  root_brent_gsl<> grb;

  // For the initial bracket 
  double a, b;

  // The parameter
  double p=0.1;
  
  // With a global function
  {
    a=-0.9, b=0.9;
    std::function<double(double)> f=gfn;
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function using std::function");
  }

  // O2scl defines 'funct' as std::function<double(double)>, so this
  // shorter notation also works. 
  {
    a=-0.9, b=0.9;
    funct f=gfn;
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function");
  }

  // A global function with a parameter
  {
    a=-0.9, b=0.9;
    funct f=std::bind(gfn_param,std::placeholders::_1,p);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Global function with parameter");
  }

  // A global function with a reference parameter, note that
  // std::ref() is used to indicate a reference parameter
  {
    a=-0.9, b=0.9;
    funct f=std::bind(gfn_rparam,std::placeholders::_1,std::ref(p));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Global function with reference parameter");
  }

  // A global function with a const reference parameter,
  // note that std::cref() is used for the const reference
  {
    a=-0.9, b=0.9;
    funct f=std::bind(gfn_crparam,std::placeholders::_1,std::cref(p));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Global function with const reference parameter");
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

  // With a const member function, note the extra const in the
  // template parameter for std::mem_fn
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double) const>(&a_class::cmfn),
		      ac,std::placeholders::_1);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Const member function");
  }

  // With a member function which has a value parameter.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double,double)>
		      (&a_class::mfn_param),
		      ac,std::placeholders::_1,0.1);
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Member function with parameter");
  }

  // With a member function which has a reference parameter.
  // Note the use of std::ref() for the reference parameter.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double,double &)>
		      (&a_class::mfn_rparam),
		      ac,std::placeholders::_1,std::ref(p));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Member function with reference parameter");
  }

  // With a const member function which has a reference parameter Note
  // the use of const in the template parameter for std::mem_fn.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double,double &) const>
		      (&a_class::cmfn_rparam),
		      ac,std::placeholders::_1,std::ref(p));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Const member function with reference parameter");
  }

  // With a const member function which has a const reference
  // parameter. Note the use of const in the template parameter
  // for std::mem_fn and the use of std::cref() for the reference
  // parameter.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=std::bind(std::mem_fn<double(double,const double &) const>
		      (&a_class::cmfn_crparam),
		      ac,std::placeholders::_1,std::cref(p));
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Const member function with const reference parameter");
  }

  // Lambda expression with inline specification
  {
    a=-0.9, b=0.9;
    funct f=[](double x) -> double { double z=sin(x)-0.1; return z; };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Lambda inline");
  }

  // Lambda expression for global function with parameter.
  {
    a=-0.9, b=0.9;
    funct f=[p](double x) -> double { return gfn_param(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for global function with parameter");
  }

  // Lambda expression for global function with reference parameter.
  // We require the 'mutable' keyword for the parameter 'p' because
  // default captures are const, and p is a non-const reference.
  {
    a=-0.9, b=0.9;
    funct f=[p](double x) mutable -> double { return gfn_rparam(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for global function with reference parameter");
  }

  // Lambda expression for global function with const reference parameter
  {
    a=-0.9, b=0.9;
    funct f=[p](double x)-> double { return gfn_crparam(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for global function with const reference parameter");
  }

  // Lambda expression for member function. We require the 'mutable'
  // keyword for the class instance 'ac' because default captures are
  // const, and mfn is not a const member function.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac](double x) mutable -> double { return ac.mfn(x); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for member function");
  }

  // Lambda expression for const member function 
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac](double x) -> double { return ac.cmfn(x); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,"Lambda for const member function");
  }

  // Lambda expression for member function with value parameter. We
  // require the 'mutable' keyword for the class instance 'ac' because
  // default captures are const, and mfn is not a const member
  // function.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac,p](double x) mutable -> double { return ac.mfn_param(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for member function with parameter");
  }

  // Lambda expression for member function with reference parameter.
  // We require the 'mutable' keyword for both 'ac' and 'p'.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac,p](double x) mutable -> double { return ac.mfn_rparam(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for member function with reference parameter");
  }

  // Lambda expression for const member function with reference
  // parameter. We require the 'mutable' keyword for the parameter 'p'
  // because default captures are const, and p is a non-const
  // reference.
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac,p](double x) mutable -> double { return ac.cmfn_rparam(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               "Lambda for const member function with reference parameter");
  }
  
  // Lambda expression for const member function with const reference
  // parameter. This is the case when 'mutable' is not required, because
  // both the reference and the member function are const. 
  {
    a=-0.9, b=0.9;
    a_class ac;
    funct f=[ac,p](double x)-> double { return ac.cmfn_crparam(x,p); };
    grb.solve_bkt(a,b,f);
    t.test_rel(a,asin(0.1),1.0e-12,
               ((std::string)"Lambda for const member function with ")+
               "const reference parameter");
  }

  t.report();
  return 0;
}
