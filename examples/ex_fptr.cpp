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

/* Example: ex_fptr.cpp
   -------------------------------------------------------------------
   This gives an example of the how member functions and external
   parameters are supplied to numerical routines. In this case, a
   member function with two parameters is passed to the root_brent_gsl
   class, which solves the equation. One of the parameters is member
   data, and the other is specified using the extra parameter argument
   to the function.
*/

#include <o2scl/funct.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

class my_class {

private:

  double parameter;

public:
  
  void set_parameter() { parameter=0.01; }
  
  // A function demonstrating the different ways of implementing
  // function parameters
  double function_to_solve(double x, double &p) {
    return atan((x-parameter)*4)*(1.0+sin((x-parameter)*50.0)/p);
  }
  
};

// Simple code to write the function to a file
int write_file(double sol);

int main(void) {
  
  cout.setf(ios::scientific);
  
  test_mgr t;
  // Only print something out if one of the tests fails
  t.set_output_level(1);

  // The solver, specifying the type of the parameter (double)
  // and the function type (funct<double>)
  root_brent_gsl<> solver;

  my_class c;
  c.set_parameter();
  
  double p=1.1;

  // This is the code that allows specification of class member
  // functions as functions to solve. This approach avoids the use of
  // static variables and functions and multiple inheritance at the
  // expense of a little overhead. We need to provide the address of
  // an instantiated object and the address of the member function.
  funct_mfptr_param<my_class,double> 
    function(&c,&my_class::function_to_solve,p);
  
  double x1=-1;
  double x2=2;
  
  // The value verbose=1 prints out iteration information
  // and verbose=2 requires a keypress between iterations.
  solver.verbose=1;
  solver.solve_bkt(x1,x2,function);

  // This is actually a somewhat difficult function to solve because
  // of the sinusoidal behavior.
  cout << "Solution: " << x1 
       << " Function value: " << c.function_to_solve(x1,p) << endl;

  // Write the function being solved to a file (see source code 
  // in examples directory for details)
  write_file(x1);

  // Obtain and summarize test results
  t.test_abs(c.function_to_solve(x1,p),0.0,1.0e-10,"ex_fptr");
  t.report();

  return 0;
}
// End of example

// Simple output function
int write_file(double x1) {
  my_class c;
  c.set_parameter();
  ofstream fout;
  double p=1.1;
  fout.open("ex_fptr.out");
  fout << "x y" << endl;
  for(double x=-1.0;x<=2.00001;x+=0.001) {
    fout << x << " " << c.function_to_solve(x,p) << endl;
  }
  fout.close();
  return 0;
}
