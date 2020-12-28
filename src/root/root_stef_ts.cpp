/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/root_stef.h>

double gfn(double x) {
  return sin(x-0.2);
}

double gfnd(double x) {
  return cos(x-0.2);
}

class cl {
public:
  
  double mfn(double x) {
    return sin(x-0.2);
  }
  
  double mfnd(double x) {
    return cos(x-0.2);
  }
};

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;

  cl acl;
  double a;
  
  t.set_output_level(2);

  // 1 - Non-templated access through a funct object 
  funct fmf=std::bind(std::mem_fn<double(double)>
			(&cl::mfn),&acl,std::placeholders::_1);
  funct fmfd=std::bind(std::mem_fn<double(double)>
			 (&cl::mfnd),&acl,std::placeholders::_1);
  root_stef<funct,funct> cr1;
  a=1.0e-5;
  cr1.verbose=1;
  cr1.solve_de(a,fmf,fmfd);
  t.test_rel(a,0.2,1.0e-6,"1");
    
  // 4 - Templated access through a global function pointer
  typedef double (*gfnt)(double);
  root_stef<gfnt,gfnt> cr4;
  gfnt gfnv=&gfn;
  gfnt gfnvd=&gfnd;
  a=1.0e-5;
  cr4.verbose=1;
  cr4.solve_de(a,gfnv,gfnvd);
  t.test_rel(a,0.2,1.0e-6,"6");

  t.report();
  return 0;
}
