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

/* Example: ex_twod_intp.cpp
   -------------------------------------------------------------------
   A simple example for two-dimensional interpolation using
   the interp_twod class.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/interp2_seq.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

// A function for filling the data and comparing results
double f(double x, double y) {
  return pow(sin(0.1*x+0.3*y),2.0);
}

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  int i,j;

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  // Create the sample data

  ubvector x(3), y(3);
  ubmatrix data(3,3);

  // Set the grid
  x[0]=0.0;
  x[1]=1.0;
  x[2]=2.0;
  y[0]=3.0;
  y[1]=2.0;
  y[2]=1.0;

  // Set and print out the data
  cout << endl;
  cout << "Data: " << endl;
  cout << "          x  | ";
  for(i=0;i<3;i++) cout << x[i] << " ";
  cout << endl;
  cout << " y           |" << endl;
  cout << "-------------|-";
  for(i=0;i<3;i++) cout << "-------------";
  cout << endl;
  for(i=0;i<3;i++) {
    cout << y[i] << " | ";
    for(j=0;j<3;j++) {
      data(i,j)=f(x[j],y[i]);
      cout << data(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;

  // Perform the interpolation

  cout << "x            y            Calc.        Exact" << endl;

  interp2_seq ti;

  // Interpolation, x-first
  double tol=0.05;
  double tol2=0.4;

  ti.set_data(3,3,x,y,data,itp_cspline,true);

  double x0, y0, x1, y1;

  x0=0.5; y0=1.5;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  x0=0.99; y0=1.99;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  x0=1.0; y0=2.0;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  cout << endl;

  // Interpolation, y-first

  ti.set_data(3,3,x,y,data,itp_cspline,false);

  x0=0.5; y0=1.5;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  x0=0.99; y0=1.99;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  x0=1.0; y0=2.0;
  cout << x0 << " " << y0 << " "
       << ti.eval(x0,y0) << " " << f(x0,y0) << endl;

  cout << endl;

  t.report();
  return 0;
}
// End of example
