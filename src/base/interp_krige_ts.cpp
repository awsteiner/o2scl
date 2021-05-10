/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2021, Andrew W. Steiner
  
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
  cout << "Data: " << endl;
  for(size_t i=0;i<4;i++) {
    cout.width(2);
    cout << x[i] << " " <<y[i] << endl;
  }
  cout << endl;

  // ---------------------------------------------------------------
  //

  interp_krige<ubvector> ik;
  std::function<double(double,double)> f=covar;

  // ---------------------------------------------------------------
  // Test normal interpolation

  cout << "Normal interpolation:" << endl;
  ik.set_covar(4,x,y,f);
  t.test_rel(ik.eval(1.0),5.0,1.0e-6,"ik 1");
  t.test_rel(ik.eval(1.5),5.5,0.1,"ik 2");
  t.test_rel(ik.eval(2.5),4.0,0.1,"ik 3");
  t.test_rel(ik.eval(3.5),3.0,0.5,"ik 4");
  cout << endl;
  
  // ---------------------------------------------------------------
  // Test normal interpolation with rescaling

  cout << "Normal interpolation with rescaling:" << endl;
  ik.set_covar(4,x,y,f,true);
  t.test_rel(ik.eval(1.0),5.0,1.0e-6,"ikr 1");
  t.test_rel(ik.eval(1.5),5.5,0.1,"ikr 2");
  t.test_rel(ik.eval(2.5),4.0,0.1,"ikr 3");
  t.test_rel(ik.eval(3.5),3.0,0.5,"ikr 4");
  cout << endl;
  
  // ---------------------------------------------------------------
  // Test interpolation with noise
  
  cout << "Noisy interpolation:" << endl;
  ik.set_covar_noise(4,x,y,f,0.5);
  t.test_rel(ik.eval(1.0),5.0,0.7,"ik 1");
  t.test_rel(ik.eval(1.5),5.5,0.7,"ik 2");
  t.test_rel(ik.eval(2.5),4.0,0.7,"ik 3");
  t.test_rel(ik.eval(3.5),3.0,0.7,"ik 4");
  cout << endl;

  // ---------------------------------------------------------------
  // Test interpolation with noise and rescaling
  
  cout << "Noisy interpolation with rescaling:" << endl;
  ik.set_covar_noise(4,x,y,f,0.5,true);
  t.test_rel(ik.eval(1.0),5.0,0.7,"ikr 1");
  t.test_rel(ik.eval(1.5),5.5,0.7,"ikr 2");
  t.test_rel(ik.eval(2.5),4.0,0.7,"ikr 3");
  t.test_rel(ik.eval(3.5),3.0,0.7,"ikr 4");
  cout << endl;

  // ---------------------------------------------------------------
  // Second set of test data

  cout << "Data: " << endl;
  ubvector x2(10), y2(10);
  for(size_t i=0;i<10;i++) {
    x2[i]=((double)i)/2.0;
    y2[i]=sin(x2[i]);
    cout.width(2);
    cout << i << " " << x2[i] << " " << y2[i] << endl;
  }
  cout << endl;

  // ---------------------------------------------------------------
  // Output default interpolation (cubic spline) results
  
  double exact, res;
  interp_vec<ubvector> io;
  io.set(10,x2,y2);

  cout << "Class interp_vec with cubic spline results: " << endl;
  exact=sin(1.01);
  res=io.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"io 1");
  exact=sin(1.0);
  res=io.eval(1.0);
  t.test_rel(exact,res,1.0e-8,"io 2");
  exact=sin(o2scl_const::pi);
  res=io.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-3,"io 3");
  cout << endl;
  
  // ---------------------------------------------------------------
  // Test interp_krige_optim interface
  
  interp_krige_optim<ubvector> iko;
  
  iko.set(10,x2,y2);

  cout << "Class interp_krige_optim with simple interface." << endl;
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"iko 1");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-4,"iko 2");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-4,"iko 3");
  cout << endl;

  cout << "Class interp_krige_optim with simple interface, "
       << "rescaled version." << endl;
  iko.set(10,x2,y2,true);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"ikor 1");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-4,"ikor 2");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-4,"ikor 3");
  cout << endl;

  // ---------------------------------------------------------------
  // Test interp_krige_optim interface with full minimization

  iko.full_min=true;
  
  iko.set(10,x2,y2);

  cout << "Class interp_krige_optim with full minimization" << endl;
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"iko 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-4,"iko 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-5,"iko 6");
  cout << endl;

  iko.set(10,x2,y2,true);

  cout << "Class interp_krige_optim with full minimization "
       << "and rescaling" << endl;
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-4,"iko 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-4,"iko 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-5,"iko 6");
  cout << endl;

  iko.full_min=false;
  
  // ---------------------------------------------------------------
  // Third set of test data

  cout.setf(ios::showpos);
  rng_gsl rg;

  double err=1.0e-2;
  ubvector x3(30), y3(30);
  cout << "Noisy data: " << endl;
  for(size_t i=0;i<30;i++) {
    x3[i]=((double)i)/6.0;
    y3[i]=sin(x3[i])+err*(2.0*rg.random()-1.0);
    cout.width(2);
    cout << i << " " << x3[i] << " " << y3[i] << endl;
  }
  cout << endl;

  cout << "Class interp_krige_optim without noise" << endl;
  
  iko.set(30,x3,y3);
  
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"io 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"io 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-2,"io 6");
  cout << endl;
  
  cout << "Class interp_krige_optim with noise" << endl;

  iko.set_noise(30,x3,y3,err*err);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"iko 7");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"iko 8");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-2,"iko 9");
  cout << endl;
  
  cout << "Class interp_krige_optim without noise but with rescaling" << endl;
  
  iko.set(30,x3,y3,true);
  
  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"io 4");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"io 5");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-2,"io 6");
  cout << endl;
  
  cout << "Class interp_krige_optim with noise and with rescaling" << endl;

  iko.set_noise(30,x3,y3,err*err,true);

  exact=sin(1.01);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"iko 7");
  exact=sin(1.0);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"iko 8");
  exact=sin(o2scl_const::pi);
  res=iko.eval(o2scl_const::pi);
  t.test_abs(exact,res,1.0e-2,"iko 9");
  cout << endl;
  
  t.report();

  return 0;
}
