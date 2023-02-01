/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/funct.h>
#include <o2scl/root_robbins_monro.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double gfn(double x) {
  return sin((x-0.2)*4)*(1.0+sin((x-0.2)*50.0)/1.1);
}

double gfnp(double x) {
  return tanh(100.0*(x-0.2));
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  cl acl;
  double a, b, r;

  funct fmg1=gfn;
  funct fmg2=gfn2;
  funct_gsl fg1(fmg1);
  funct_gsl fg2(fmg2);

  for(size_t k=0;k<2;k++) {

    // GSL 
    {
      int status;
      int iter=0, max_iter=100;

      a=-1.0;
      b=1.0;
    
      const gsl_root_fsolver_type *T=gsl_root_fsolver_brent;
      gsl_root_fsolver *s=gsl_root_fsolver_alloc (T);
      
      if (k==0) gsl_root_fsolver_set(s,&fg1,a,b);
      else gsl_root_fsolver_set(s,&fg2,a,b);
     
      do {

	iter++;
	status=gsl_root_fsolver_iterate(s);
	r=gsl_root_fsolver_root(s);
	a=gsl_root_fsolver_x_lower(s);
	b=gsl_root_fsolver_x_upper(s);
	status=gsl_root_test_interval(a,b,0,1.0e-6);

	cout << iter << " " << a << " " << r << " " << b << endl;

      } while (status==gsl_continue && iter<max_iter);
     
      gsl_root_fsolver_free(s);

      if (k==0) {
	cout << r << " " << gfn(r) << endl;
      } else {
	cout << r << " " << gfn2(r) << endl;
      }
      t.test_rel(r,0.2,1.0e-10,"1");
    }

    // Using set() and iterate() with a function pointer
    typedef double (*gfnt)(double);
    root_robbins_monro<gfnt> grb1;
    gfnt gfnv=&gfn;
    if (k==1) gfnv=&gfn2;
    a=-1.0;
    b=1.0;
    grb1.set(gfnv,a,b);
    int status=gsl_continue;
    for(size_t iter=0;iter<100 && status==gsl_continue;iter++) {

      grb1.iterate(gfnv);

      r=grb1.get_root();
      a=grb1.get_lower();
      b=grb1.get_upper();
      status=gsl_root_test_interval(a,b,0,1.0e-6);

      cout << iter << " " << a << " " << r << " " << b << endl;
    }
    
    if (k==0) {
      cout << r << " " << gfn(r) << endl;
    } else {
      cout << r << " " << gfn2(r) << endl;
    }
    t.test_rel(grb1.get_root(),0.2,1.0e-10,"1");

  }

  // Using a funct object and the solve_bkt() interface
  funct fmf=std::bind(std::mem_fn<double(double)>
			(&cl::mfn),&acl,std::placeholders::_1);
  root_robbins_monro<> grb2;
  a=-1.0;
  b=1.0;
  grb2.solve_bkt(a,b,fmf);
  t.test_rel(a,0.2,1.0e-10,"1");
  cout << a << " " << gfn(a) << endl;

  t.report();
  return 0;
}

