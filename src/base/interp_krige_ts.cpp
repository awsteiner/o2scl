/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017, Andrew W. Steiner
  
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
#include <o2scl/interp_krige.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double covar(double x, double y) {
  return exp(-2.0*(x-y)*(x-y));
}

double covar_noise(double x, double y) {
  if (x==y) return exp(-2.0*(x-y)*(x-y))+0.5;
  return exp(-2.0*(x-y)*(x-y));    
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // ---------------------------------------------------------------
  // Create test data

  ubvector x(4), y(4);
  for(size_t i=0;i<4;i++) {
    x[i]=((double)i)+1.0;
  }
  y[0]=5.0;
  y[1]=6.0;
  y[2]=2.0;
  y[3]=3.0;

  // ---------------------------------------------------------------
  //

  interp_krige<ubvector> ik;
  std::function<double(double,double)> f=covar;
  std::function<double(double,double)> f_noise=covar_noise;

  // ---------------------------------------------------------------
  // Test normal interpolation

  ik.set(4,x,y,f);
  cout << ik.eval(1.0) << endl;
  cout << ik.eval(1.5) << endl;
  cout << ik.eval(2.5) << endl;
  cout << ik.eval(3.5) << endl;
  cout << endl;

  // ---------------------------------------------------------------
  // Test interpolation with noise

  ik.set(4,x,y,f_noise);
  cout << ik.eval(1.0) << endl;
  cout << ik.eval(1.5) << endl;
  cout << ik.eval(2.5) << endl;
  cout << ik.eval(3.5) << endl;
  cout << endl;

  t.report();

  return 0;
}
