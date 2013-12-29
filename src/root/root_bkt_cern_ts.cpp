/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/funct.h>

double gfn(double x) {
  return sin(x-0.2);
}

class cl {
public:
  
  virtual ~cl() {}

  virtual double mfn(double x) {
    return sin(x-0.2);
  }

};

using namespace std;
using namespace o2scl;

int main(void) {
  cl acl;
  double a, b;
  int i;
  size_t tmp;
  test_mgr t;

  t.set_output_level(1);

  // 1 - Non-templated access through a funct object 
  funct_mfptr<cl> fmf(&acl,&cl::mfn);
  root_bkt_cern<funct_mfptr<cl> > cr1;
  tmp=clock();
  a=1.0e-5;
  b=1.0;
  cr1.solve_bkt(a,b,fmf);
  t.test_rel(a,0.2,1.0e-6,"1");
    
  // 4 - Templated access through a global function pointer
  typedef double (*gfnt)(double);
  root_bkt_cern<gfnt> cr4b;
  gfnt gfnv=&gfn;
  tmp=clock();
  a=1.0e-5;
  b=1.0;
  cr4b.solve_bkt(a,b,gfnv);
  t.test_rel(a,0.2,1.0e-6,"4");

  t.report();
  return 0;
}

