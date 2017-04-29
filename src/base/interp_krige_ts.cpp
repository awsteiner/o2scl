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
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double covar(double x, double y) {
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

  // ---------------------------------------------------------------
  // Test normal interpolation

  ik.set_covar(4,x,y,f);
  cout << ik.eval(1.0) << endl;
  cout << ik.eval(1.5) << endl;
  cout << ik.eval(2.5) << endl;
  cout << ik.eval(3.5) << endl;
  cout << endl;

  // ---------------------------------------------------------------
  // Test interpolation with noise
  
  ik.set_covar_noise(4,x,y,f,0.5);
  cout << ik.eval(1.0) << endl;
  cout << ik.eval(1.5) << endl;
  cout << ik.eval(2.5) << endl;
  cout << ik.eval(3.5) << endl;
  cout << endl;

  // ---------------------------------------------------------------
  // Second set of test data

  cout.setf(ios::showpos);
  
  ubvector x2(10), y2(10);
  for(size_t i=0;i<10;i++) {
    x2[i]=((double)i)/2.0;
    y2[i]=sin(x2[i]);
  }
    
  interp_vec<ubvector> io;
  io.set(10,x2,y2);
  cout << io.eval(1.01) << " " << sin(1.01) << " "
       << fabs((io.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
  cout << io.eval(1.0) << " " << sin(1.0) << " "
       << fabs((io.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
  cout << io.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
       << fabs(io.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
  cout << endl;
  
  interp_krige_optim<ubvector> iko;
  iko.set(10,x2,y2);
  cout << iko.eval(1.01) << " " << sin(1.01) << " "
       << fabs((iko.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
  cout << iko.eval(1.0) << " " << sin(1.0) << " "
       << fabs((iko.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
  cout << iko.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
       << fabs(iko.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
  cout << endl;

  iko.full_min=true;
  iko.set(10,x2,y2);
  cout << iko.eval(1.01) << " " << sin(1.01) << " "
       << fabs((iko.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
  cout << iko.eval(1.0) << " " << sin(1.0) << " "
       << fabs((iko.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
  cout << iko.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
       << fabs(iko.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
  cout << endl;

  // ---------------------------------------------------------------
  // Third set of test data

  cout.setf(ios::showpos);
  rng_gsl rg;

  double err=1.0e-2;
  ubvector x3(30), y3(30);
  for(size_t i=0;i<30;i++) {
    x3[i]=((double)i)/6.0;
    y3[i]=sin(x3[i])+err*(2.0*rg.random()-1.0);
  }
    
  io.set(30,x3,y3);
  cout << io.eval(1.01) << " " << sin(1.01) << " "
       << fabs((io.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
  cout << io.eval(1.0) << " " << sin(1.0) << " "
       << fabs((io.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
  cout << io.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
       << fabs(io.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
  cout << endl;

  iko.full_min=false;
  iko.set_noise(30,x3,y3,err*err);
  cout << iko.eval(1.01) << " " << sin(1.01) << " "
       << fabs((iko.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
  cout << iko.eval(1.0) << " " << sin(1.0) << " "
       << fabs((iko.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
  cout << iko.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
       << fabs(iko.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
  cout << endl;

  /*
    iko.full_min=true;
    iko.set_noise(30,x3,y3,err*err);
    cout << iko.eval(1.01) << " " << sin(1.01) << " "
    << fabs((iko.eval(1.01)-sin(1.01))/sin(1.01)) << endl;
    cout << iko.eval(1.0) << " " << sin(1.0) << " "
    << fabs((iko.eval(1.0)-sin(1.0))/sin(1.0)) << endl;
    cout << iko.eval(acos(-1.0)) << " " << sin(acos(-1.0)) << " "
    << fabs(iko.eval(acos(-1.0))-sin(acos(-1.0))) << endl;
    cout << endl;
  */
  
  t.report();

  return 0;
}
