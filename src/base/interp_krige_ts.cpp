/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2018, Andrew W. Steiner
  
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
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double covar(double x, double y) {
  return exp(-2.0*(x-y)*(x-y));
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

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
  t.test_rel(ik.eval(1.0),5.0,1.0e-6,"ik 1");
  t.test_rel(ik.eval(1.5),5.5,0.1,"ik 2");
  t.test_rel(ik.eval(2.5),4.0,0.1,"ik 3");
  t.test_rel(ik.eval(3.5),3.0,0.5,"ik 4");
  cout << endl;

  // ---------------------------------------------------------------
  // Test interpolation with noise
  
  ik.set_covar_noise(4,x,y,f,0.5);
  t.test_rel(ik.eval(1.0),5.0,0.7,"ik 1");
  t.test_rel(ik.eval(1.5),5.5,0.7,"ik 2");
  t.test_rel(ik.eval(2.5),4.0,0.7,"ik 3");
  t.test_rel(ik.eval(3.5),3.0,0.7,"ik 4");
  cout << endl;

  // ---------------------------------------------------------------
  // Second set of test data

  cout.setf(ios::showpos);
  
  ubvector x2(10), y2(10);
  for(size_t i=0;i<10;i++) {
    x2[i]=((double)i)/2.0;
    y2[i]=sin(x2[i]);
  }

  double exact, res;
  interp_vec<ubvector> io;
  io.set(10,x2,y2);
  
  exact=sin(1.01);
  res=io.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"io 1");
  exact=sin(1.0);
  res=io.eval(1.0);
  t.test_rel(exact,res,1.0e-8,"io 2");
  exact=sin(o2scl_const::pi);
  res=io.eval(o2scl_const::pi);
  t.test_rel(exact,res,2.0,"io 3");
  cout << endl;
  
  interp_krige_optim<ubvector> iko;
  iko.set(10,x2,y2);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-2,"iko 1");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-6,"iko 2");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_rel(exact,res,2.0,"iko 3");
  cout << endl;

  iko.full_min=true;
  iko.set(10,x2,y2);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-2,"iko 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-8,"iko 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_rel(exact,res,2.0,"iko 6");
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
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-2,"io 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-8,"io 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_rel(exact,res,2.0,"io 6");
  cout << endl;
  
  iko.full_min=false;
  iko.set_noise(30,x3,y3,err*err);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-2,"iko 7");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-2,"iko 8");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_rel(exact,res,2.0,"iko 9");
  cout << endl;
  
  t.report();

  return 0;
}
