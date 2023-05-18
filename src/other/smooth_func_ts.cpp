/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2016-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/smooth_func.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/rng.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class ns_test {

public:

  double eps;

  ns_test() {
    eps=1.0e-5;
  }

  int f(size_t nv, const ubvector &x, ubvector &y) {
    // Generate deterministic noise from a high-frequency trig.
    // function (if it was truly random the test result would
    // be unpredictable).
    y[0]=x[0]*x[0]*x[0]-3.0+eps*sin(1.0e9*(x[0]+x[1]));
    y[1]=x[0]*x[1]-4.0+eps*cos(1.0e9*(x[0]+x[1]));
    return 0;
  }

};

double od_test(double x) {
  return 2.0+sin(x*1.0e20);
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  ns_test nst;

  smooth_func<ubvector,mm_funct> sf;

  mm_funct mf=std::bind
    (std::mem_fn<int(size_t,const ubvector &, ubvector &)>
     (&ns_test::f),&nst,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  sf.set_func(mf);

  vector<double> step={1.0e-4};
  sf.set_step(step);

  mm_funct sff=std::bind
    (std::mem_fn<int(size_t,const ubvector &, ubvector &)>
     (&smooth_func<ubvector,mm_funct>::operator()),&sf,
     std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  mroot_hybrids<> mh;
  mh.err_nonconv=false;
  mh.verbose=1;
  ubvector x(2);
  int mret;

  x[0]=1.0;
  x[1]=1.0;
  mret=mh.msolve(2,x,mf);
  cout << "ret: " << mret << endl;

  cout << endl;
  
  x[0]=1.0;
  x[1]=1.0;
  mret=mh.msolve(2,x,sff);
  cout << "ret: " << mret << endl;
  cout << x[0]*x[0]*x[0] << " " << x[0]*x[1] << endl;

  t.test_rel(x[0]*x[0]*x[0],3.0,1.0e-4,"x[0]");
  t.test_rel(x[0]*x[1],4.0,1.0e-4,"x[1]");

  funct fx=od_test;
  gauss_filter<> gf;
  gf.set_func(fx);
  gf.set_K(11);
  cout << gf(2.0) << endl;
  cout << fx(2.0) << endl;
  
  t.report();
  
  return 0;
}


