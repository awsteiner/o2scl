/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/multi_funct.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double func2(size_t nv, const ubvector &x, double &pa);

double func2(size_t nv, const ubvector &x, double &pa) {
  return pa+x[0]-x[1];
}

class ac {
public:
  double mfunc2(size_t nv, const ubvector &x, double &pa) {
    return pa+x[0]-x[1];
  }
};

typedef double arr_t[2];

double vfunc2(size_t nv, const arr_t &x, double &pa);

double vfunc2(size_t nv, const arr_t &x, double &pa) {
  return pa+x[0]-x[1];
}

class vac {
public:
  double mvfunc2(size_t nv, const arr_t &x, double &pa) {
    return pa+x[0]-x[1];
  }
};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

#ifdef O2SCL_NEVER_DEFINED

  // --------------------------------------------------------------------
  // test ubvector versions
  // --------------------------------------------------------------------

  double a=2.0, y;
  ac c1;
  vac c12;

  multi_funct_fptr_param<double,ubvector> f4(func2,a);
  multi_funct_mfptr_param<ac,double,ubvector> f5(&c1,&ac::mfunc2,a);

  ubvector x(2);
  x[0]=1.2;
  x[1]=3.5;
  void *vp=0;
  
  y=f4(2,x);
  t.test_rel(y,-0.3,1.0e-6,"fptr");
  
  y=f5(2,x);
  t.test_rel(y,-0.3,1.0e-6,"mfptr");

  // --------------------------------------------------------------------
  // test array versions
  // --------------------------------------------------------------------

  multi_funct_fptr_param<double,arr_t> f42(vfunc2,a);
  multi_funct_mfptr_param<vac,double,arr_t> f52(&c12,&vac::mvfunc2,a);

  double x2[2];
  x2[0]=1.2;
  x2[1]=3.5;
  
  y=f42(2,x2);
  t.test_rel(y,-0.3,1.0e-6,"fptr");
  
  y=f52(2,x2);
  t.test_rel(y,-0.3,1.0e-6,"mfptr");

#endif

  t.report();
  return 0;
}
