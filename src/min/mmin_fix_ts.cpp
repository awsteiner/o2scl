/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
#include <o2scl/mmin_fix.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mmin_conf.h>
#include <o2scl/mmin_bfgs2.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double minfun(size_t n, const ubvector &x) {
  double ret;
  ret=(x[0]+x[1]/10.0)*(x[0]+x[1]/10.0)+
    (x[1]-2.0+x[0])*(x[1]-2.0+x[0])+3.0;
  return ret;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);

  double min, min2;
  ubvector x(3);

  mmin_fix_params<> g;
  g.tol_abs/=1.0e2;

  multi_funct mf=minfun;
  
  // Test mmin() which doesn't fix any parameters
  x[0]=1.0;
  x[1]=1.0;
  g.verbose=1;
  int ret=g.mmin(2,x,min,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],-2.0/9.0,1.0e-4,"mmin_fix 1");
  t.test_rel(x[1],20.0/9.0,2.0e-4,"mmin_fix 2");
  g.verbose=0;
  
  bool fix[2]={false,false};
  bool *fixp=&(fix[0]);

  // Test mmin_fix() with all entries to false
  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin_fix(2,x,min,fixp,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],-2.0/9.0,1.0e-4,"mmin_fix 3");
  t.test_rel(x[1],20.0/9.0,2.0e-4,"mmin_fix 4");
  
  // Test mmin_fix() while fixing 1st parameter
  x[0]=0.0;
  x[1]=2.0;
  fix[0]=true;
  ret=g.mmin_fix(2,x,min,fixp,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_rel(x[0],0.0,1.0e-6,"mmin_fix 5");
  t.test_rel(x[1],200.0/101.0,1.0e-6,"mmin_fix 6");
  fix[0]=false;

  // Test mmin_fix() while fixing 2nd parameter
  x[0]=0.0;
  x[1]=2.0;
  fix[1]=true;
  ret=g.mmin_fix(2,x,min,fixp,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_rel(x[0],-0.1,1.0e-4,"mmin_fix 7");
  t.test_rel(x[1],2.0,1.0e-6,"mmin_fix 8");
  fix[1]=false;
  
#ifdef O2SCL_NEVER_DEFINED

  // AWS 3/30/2020: this doesn't work anymore because mmin_fix_params
  // currently requires dfunc_t and func_t to be the same. I think the
  // way to fix this is to add another template parameter to
  // mmin_fix_params which specifies the new dfunc_t type.
  
  // Test using a different minimizer
  mmin_conf<mmin_fix_params<>,ubvector,grad_funct,
    gradient<mmin_fix_params<>,ubvector>,
    gradient_gsl<mmin_fix_params<>,ubvector> > gmc;
  g.set_mmin(gmc);

  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin_fix(2,x,min,fixp,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],-2.0/9.0,1.0e-4,"mmin_fix 9");
  t.test_rel(x[1],20.0/9.0,2.0e-4,"mmin_fix 10");

  x[0]=0.0;
  x[1]=2.0;
  fix[1]=true;
  ret=g.mmin_fix(2,x,min,fixp,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_rel(x[0],-0.1,1.0e-4,"mmin_fix 11");
  t.test_rel(x[1],2.0,1.0e-6,"mmin_fix 12");
  fix[1]=false;

#endif
  
  t.report();
  return 0;
}
 
