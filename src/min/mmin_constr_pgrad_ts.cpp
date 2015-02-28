/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#include <gsl/gsl_blas.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mmin_constr_pgrad.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double quadratic(size_t nv, const ubvector &x) {
  double y=0.0;
  for(size_t i=0;i<nv;i++) {
    double t=x[i]-((double)i)/10.0;
    y+=(i+1.0)*t*t;
  }
  return y;
}

int quadratic_df(size_t nv, ubvector &x, ubvector &g) {

  for(size_t i=0;i<nv;i++) {
    double t=x[i]-((double)i)/10.0;
    g[i]=2.0*(i+1.0)*t;
  }
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  std::cout.setf(ios::scientific);

  size_t nv=100;

  multi_funct11 mff=quadratic;
  grad_funct11 gff=quadratic_df;
  
  mmin_constr_pgrad<multi_funct11,grad_funct11,ubvector> omp;

  ubvector c1(nv), c2(nv), x(nv);

  double fmin;
  int vp=0;
    
  for(size_t i=0;i<nv;i++) {
    c1[i]=-3.0;
    c2[i]=3.0;
  }
    
  for(size_t i=0;i<nv;i++) {
    x[i]=1.0+((double)i);
  }

  omp.ntrial=1000;
  omp.set_constraints(nv,c1,c2);
  omp.mmin_de(nv,x,fmin,mff,gff);
  t.test_gen(omp.last_ntrial==135,"Num iters.");
  t.test_rel(x[0],4.8555e-5,1.0e-4,"result");
  t.test_rel(fmin,9.300970e+4,1.0e-6,"min");

  t.report();

  return 0;
}
