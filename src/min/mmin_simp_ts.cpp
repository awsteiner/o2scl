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
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/multi_funct.h>
#include <o2scl/mmin_simp.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double minfun(size_t n, const ubvector &x) {
  return x[0]*x[0]+(x[1]-2.0)*(x[1]-2.0)+3.0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  double min=0.0, min2;
  ubvector x(2), x2(2);
  double xa[2];
  int vp=0;
  mmin_simp<multi_funct<> > g;
  
  cout.setf(ios::scientific);
  
  multi_funct_fptr<> mf(minfun);
  
  int ret;

  // Standard function

  x[0]=1.1;
  x[1]=0.9;
  ret=g.mmin(2,x,min,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],0.0,1.0e-4,"gsl_mmin_nmsimplex 1");
  t.test_rel(x[1],2.0,2.0e-4,"gsl_mmin_nmsimplex 2");
  
  // GSL-like interface

  x[0]=1.0;
  x[1]=1.0;
  ubvector step_size(2);
  step_size[0]=0.1;
  step_size[1]=0.1;
  g.set(mf,2,x,step_size);
  for(size_t it=0;it<100 && g.size>g.tol_abs;it++) {
    g.iterate();
  }
  x[0]=g.x[0];
  x[1]=g.x[1];
  min=g.fval;
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],0.0,1.0e-4,"gsl_mmin_nmsimplex 1");
  t.test_rel(x[1],2.0,2.0e-4,"gsl_mmin_nmsimplex 2");

  // Specify full simplex

  ubmatrix simp(3,2);
  simp(0,0)=1.0;
  simp(0,1)=1.0;
  simp(1,0)=1.1;
  simp(1,1)=1.1;
  simp(2,0)=2.0;
  simp(2,1)=1.0;
  ret=g.mmin_simplex(2,simp,min,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(simp(0,0),0.0,1.0e-4,"gsl_mmin_nmsimplex 1");
  t.test_rel(simp(0,1),2.0,2.0e-4,"gsl_mmin_nmsimplex 2");

  t.report();
  return 0;
}
 
