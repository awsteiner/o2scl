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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_gen_comp.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int nfunc;

double lower_limit(size_t nv, const ubvector &x) {
  return 0.0;
}

double upper_limit(size_t nv, const ubvector &x) {
  if (nv==0) return 1.0;
  else if (nv==1) return sqrt(1.0-x[0]*x[0]);
  else if (nv==2) return sqrt(1.0-x[0]*x[0]-x[1]*x[1]);
  else if (nv==3) return sqrt(1.0-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
  return 0.0;
}

double test_fun(size_t nv, const ubvector &x) {
  nfunc++;
  return 1.0;
}

double test_fun2(size_t nv, const ubvector &x) {
  nfunc++;
  if (nv==1) return x[0]*x[0];
  else if (nv==2) return x[1]*x[1];
  return x[2]*x[2];
}

int main(void) {
  inte_gen_comp<multi_funct<> > ci;

  /// The integrator
  typedef inte_gen_comp<multi_funct<>,
    multi_funct<>,multi_funct<>,ubvector> ior_type;

  inte_qag_gsl<funct11 > gl1;
  inte_qag_gsl<funct11 > gl2;
  inte_qag_gsl<funct11 > gl3;
  inte_qag_gsl<funct11 > gl4;
  
  int vp=0;
  inte<funct11 > **ip;
  double res;
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
 
  multi_funct_fptr<> l(lower_limit);
  multi_funct_fptr<> u(upper_limit);
  multi_funct_fptr<> tf(test_fun);
  multi_funct_fptr<> tf2(test_fun2);

  ci.set_oned_inte(gl1,0);
  ci.set_oned_inte(gl2,1);
  ci.set_oned_inte(gl3,2);
  ci.set_oned_inte(gl4,3);

  // Calculate the volume of one octant of a 3d sphere of unit radius
  // = (4/3 pi)/8 = pi/6
  
  nfunc=0;
  res=ci.ginteg(tf,3,l,u);
  cout << acos(-1.0)/6.0 << " " << res << " " << nfunc << endl;
  t.test_rel(acos(-1.0)/6.0,res,5.0e-10,"value");

  // Calculate the volume of one 16th of a 4d sphere of unit radius
  // = pi^2/2/16 = pi/32

  // This doesn't work as well. Why?
  
  nfunc=0;
  res=ci.ginteg(tf,4,l,u);
  cout << pow(acos(-1.0),2.0)/32.0 << " " << res << " " << nfunc << endl;
  t.test_rel(pow(acos(-1.0),2.0)/32.0,res,5.0e-2,"value");
  
  // Calculate the volume of one octant of a sphere of unit radius
  // weighted by x^2*y^2*z^2 = pi/1890

  nfunc=0;
  res=ci.ginteg(tf2,3,l,u);
  cout << acos(-1.0)/1890.0 << " " << res << " " << nfunc << endl;
  t.test_rel(acos(-1.0)/1890.0,res,5.0e-10,"value");

  t.report();
  return 0;
}

