 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_cern.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int gfn(size_t nv, const ubvector &x, ubvector &y) {
  y[0]=sin(x[1]-0.2);
  y[1]=sin(x[0]-0.25);
  return 0;
}

class cl {
public:

  int mfn(size_t nv, const ubvector &x, ubvector &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }
  
  // Placeholder jacobian (ignored by mroot_cern)
  int mfnd(size_t nv, ubvector &x, size_t ny, 
	   ubvector &y, ubmatrix&j) {
    return 0;
  }
  
};

using namespace std;
using namespace o2scl;

int main(void) {
  cl acl;
  ubvector x(2);
  double xa[2];
  int i;
  int vp=0;
  size_t tmp;
  //  int N=700;
  int N=1;
  int t1=0, t2=0, t3=0, t4=0, t5=0;
  test_mgr t;

  t.set_output_level(1);

  for(int kk=0;kk<1;kk++) {

    // 1 - Non-templated access through a funct object 
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&cl::mfn),&acl,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    mroot_cern<mm_funct,ubvector,jac_funct> cr1;
    tmp=clock();
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
	x[0]=0.5;
	x[1]=0.5;
	cr1.msolve(2,x,fmf);
      }
    }
    t1+=(clock()-tmp)/10000;
    cout << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
    t.test_rel(x[0],0.25,1.0e-6,"1a");
    t.test_rel(x[1],0.2,1.0e-6,"1b");

    // Show that msolve_de() works
    jac_funct fmfd=std::bind
      (std::mem_fn<int(size_t,ubvector &,size_t,ubvector &,ubmatrix &)>
       (&cl::mfnd),&acl,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,std::placeholders::_4,std::placeholders::_5);

    x[0]=0.5;
    x[1]=0.5;
    cr1.msolve_de(2,x,fmf,fmfd);
    t.test_rel(x[0],0.25,1.0e-6,"1c");
    t.test_rel(x[1],0.2,1.0e-6,"1d");

    // 4 - Templated access through a global function pointer
    typedef int (*gfnt)(size_t, const ubvector &, ubvector &);
    mroot_cern<gfnt,ubvector> cr4;
    gfnt gfnv=&gfn;
    tmp=clock();
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
	x[0]=0.5;
	x[1]=0.5;
	cr4.msolve(2,x,gfnv);
      }
    }
    t4+=(clock()-tmp)/10000;
    cout << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
    t.test_rel(x[0],0.25,1.0e-6,"5a");
    t.test_rel(x[1],0.2,1.0e-6,"5b");

  }

  cout << t1 << endl;
  cout << t2 << endl;
  cout << t4 << endl;
  cout << t5 << endl;

  t.report();
  return 0;
}

