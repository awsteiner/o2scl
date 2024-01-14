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

class cl {
public:
  double mfn(double x) {
    return sin(x-0.2);
  }
};

int main(void) {
  cl acl;
  double a, b;
  int i;
  int vp=0;
  size_t tmp;
  int N=1;
  int t1=0, t2=0, t3=0, t4=0;
  test_mgr t;

  t.set_output_level(2);
  
  typedef double (*gfnt)(double);

  for(int kk=0;kk<1;kk++) {
    
    // 1 - Non-templated access through a funct object 
    funct fmf=std::bind(std::mem_fn<double(double)>
			  (&cl::mfn),&acl,std::placeholders::_1);
    root_cern<> cr1;

    tmp=clock();
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
	a=1.0e-5;
	cr1.solve(a,fmf);
      }
    }
    t1+=(clock()-tmp)/10000;
    cout << (clock()-tmp)/10000 << " " << a << endl;
    t.test_rel(a,0.2,1.0e-6,"1");
    
    // 4 - Templated access through a global function pointer
    root_cern<gfnt> cr4;
    gfnt gfnv=&gfn;
    tmp=clock();
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
	a=1.0;
	cr4.solve(a,gfnv);
      }
    }
    t4+=(clock()-tmp)/10000;
    cout << (clock()-tmp)/10000 << " " << a << endl;
    t.test_rel(a,0.2,1.0e-6,"4");
    cout << endl;
  }

  cout << t1 << endl;
  cout << t4 << endl;

  t.report();
  return 0;
}

