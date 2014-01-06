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
#include <o2scl/mm_funct.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int func(size_t nv, const ubvector &x, ubvector &y, double &pa) {
  y[0]=pa+x[0]-x[1];
  y[1]=pa+x[0]+x[1];
  return 0;
}

int funca(size_t nv, const double x[2], double y[2], double &pa) {
  y[0]=pa+x[0]-x[1];
  y[1]=pa+x[0]+x[1];
  return 0;
}

typedef double arr_t[2];

class ac {
public:
  int mfunc(size_t nv, const ubvector &x, ubvector &y, double &pa) {
    y[0]=pa+x[0]-x[1];
    y[1]=pa+x[0]+x[1];
    return 0;
  }
  int mfunca(size_t nv, const arr_t &x, arr_t &y, double &pa) {
    y[0]=pa+x[0]-x[1];
    y[1]=pa+x[0]+x[1];
    return 0;
  }

};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  double a=2.0;
  mm_funct<ubvector> *basep;
  mm_funct_fptr_param<double,ubvector> f1(func,a);
  
  ubvector ap(1);
  ap[0]=2.0;
  ubvector x(2), y(2);
  x[0]=1.2;
  x[1]=3.5;
  
  double xa[2]={1.2,3.5}, ya[2];

  // ---------------------------------------------------------
  // ubvector section
  // ---------------------------------------------------------
  
  // Normal function pointer
  y[0]=0.0; y[1]=0.0;
  f1(2,x,y);
  t.test_rel(y[0],-0.3,1.0e-6,"fptr");
  t.test_rel(y[1],6.7,1.0e-6,"fptr");

  // Check that a base class pointer also works
  y[0]=0.0; y[1]=0.0;
  basep=&f1;
  (*basep)(2,x,y);
  t.test_rel(y[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(y[1],6.7,1.0e-6,"mfptr");

  // Member function pointer
  y[0]=0.0; y[1]=0.0;
  ac c1;
  mm_funct_mfptr_param<ac,double,ubvector> f3(&c1,&ac::mfunc,a);
  f3(2,x,y);
  t.test_rel(y[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(y[1],6.7,1.0e-6,"mfptr");

  // Check that a base class pointer also works
  y[0]=0.0; y[1]=0.0;
  basep=&f3;
  (*basep)(2,x,y);
  t.test_rel(y[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(y[1],6.7,1.0e-6,"mfptr");
  
  // The following line will not compile under gcc:
  mm_funct_mfptr_param<ac,double,arr_t> f3a(&c1,&ac::mfunca,a);
  // This is the reason we need to have a separate set of classes
  // for double[nv].
  f3a(2,xa,ya);
  t.test_rel(ya[0],-0.3,1.0e-6,"mfptra");
  t.test_rel(ya[1],6.7,1.0e-6,"mfptra");

#ifdef O2SCL_NEVER_DEFINED

  // ---------------------------------------------------------
  // double[] section
  // ---------------------------------------------------------
  
  mm_vfunct<double,2> *basepa;
  mm_vfunct_fptr<double,2> f1a(funca);
  //mm_vfunct_strings<double,2> f2a(2,fs,"x,x2",1,"a");

  // Normal function pointer
  ya[0]=0.0; ya[1]=0.0;
  f1a(2,xa,ya);
  t.test_rel(ya[0],-0.3,1.0e-6,"fptr");
  t.test_rel(ya[1],6.7,1.0e-6,"fptr");

  // Check that a base class pointer also works
  ya[0]=0.0; ya[1]=0.0;
  basepa=&f1a;
  (*basepa)(2,xa,ya);
  t.test_rel(ya[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(ya[1],6.7,1.0e-6,"mfptr");

  // Member function pointer
  ya[0]=0.0; ya[1]=0.0;
  mm_vfunct_mfptr<ac,double,2> f3a(&c1,&ac::mfunca);
  f3a(2,xa,ya);
  t.test_rel(ya[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(ya[1],6.7,1.0e-6,"mfptr");

  // Check that a base class pointer also works
  ya[0]=0.0; ya[1]=0.0;
  basepa=&f3a;
  (*basepa)(2,xa,ya);
  t.test_rel(ya[0],-0.3,1.0e-6,"mfptr");
  t.test_rel(ya[1],6.7,1.0e-6,"mfptr");

#endif
  
  t.report();
  return 0;
}

