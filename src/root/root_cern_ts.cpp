/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/root_cern.h>

using namespace std;
using namespace o2scl;

double gfn(double x) {
  return sin(x-0.2);
}

double gfn2(double x) {
  return atan((x-0.2)*4)*(1.0+sin((x-0.2)*50.0)/1.1);
}

double gfn3(double x) {
  return tanh(100.0*(x-0.2));
}

class cl {
public:
  double mfn(double x) {
    return gfn(x);
  }
};

template<class fp_t> fp_t cbrt_fun(fp_t x) {
  return x*x*x-5;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);
  
  cl acl;
  double a, b;
  
  typedef double (*gfnt)(double);
  
  // 1 - Non-templated access through a funct object 
  funct fmf=std::bind(std::mem_fn<double(double)>
                      (&cl::mfn),&acl,std::placeholders::_1);
  root_cern<> cr1;
  
  a=1.0e-5;
  cr1.solve(a,fmf);
  t.test_rel(a,0.2,1.0e-6,"1");
  
  // 4 - Templated access through a global function pointer
  root_cern<gfnt> cr4;
  gfnt gfnv=&gfn;
  a=1.0;
  cr4.solve(a,gfnv);

  t.test_rel(a,0.2,1.0e-6,"4");

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  cout << "Using adaptive multiprecision with a simple function and a\n"
       << "  lambda expression:" << endl;
  double am=1.0, valm, errm;
  cr1.verbose=1;
  int amret=cr1.solve_multip(am,[](auto &&tx) mutable
  { return cbrt_fun(tx); },errm);
  cout << dtos(am,0) << " " << dtos(cbrt(5.0),0) << endl;
  t.test_rel(am,cbrt(5.0),2.0e-15,"multiprecision");
  
#endif
  
  t.report();
  return 0;
}

