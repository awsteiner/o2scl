/*
   ───────────────────────────────────────────────────────────────────

   Copyright (C) 2026, Andrew W. Steiner

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
#include <iostream>
#include <cmath>

#include <o2scl/multi_funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/cma_es.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

// A simple, separable, unimodal test function
double sphere(size_t nvar, const ubvector &x) {
  double s=0.0;
  for(size_t i=0;i<nvar;i++) s+=x[i]*x[i];
  return s;
}

// The classic non-separable, narrow-valley Rosenbrock function,
// with global minimum 0 at x_i=1 for all i
double rosenbrock(size_t nvar, const ubvector &x) {
  double s=0.0;
  for(size_t i=0;i+1<nvar;i++) {
    double t1=x[i+1]-x[i]*x[i];
    double t2=1.0-x[i];
    s+=100.0*t1*t1+t2*t2;
  }
  return s;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // Plain (mu/mu_w,lambda)-CMA-ES, no restarts, on the sphere
  // function
  {
    cma_es<> ce;
    ce.rg.set_seed(1);

    ubvector x(5);
    for(size_t i=0;i<5;i++) x[i]=1.0;
    double fmin;
    multi_funct mf=sphere;

    ce.mmin(5,x,fmin,mf);

    cout << "cma_es, sphere, no restarts: fmin=" << fmin << endl;
    t.test_rel(fmin,0.0,1.0e-6,"sphere fmin");
    for(size_t i=0;i<5;i++) {
      t.test_abs(x[i],0.0,1.0e-3,"sphere x");
    }
  }

  // Rosenbrock without restarts
  {
    cma_es<> ce;
    ce.rg.set_seed(2);
    ce.ntrial=2000;

    ubvector x(4);
    for(size_t i=0;i<4;i++) x[i]=0.0;
    double fmin;
    multi_funct mf=rosenbrock;

    ce.mmin(4,x,fmin,mf);

    cout << "cma_es, rosenbrock, no restarts: fmin=" << fmin << endl;
    t.test_gen(fmin<1.0,"rosenbrock fmin (no restart)");
  }

  // Rosenbrock with IPOP restarts, using a small initial
  // population and step size to make the restart logic more
  // likely to be exercised
  {
    cma_es<> ce;
    ce.rg.set_seed(3);
    ce.restart_mode=cma_es<>::cma_es_restart_ipop;
    ce.lambda=6;
    ce.max_evals=40000;

    ubvector x(4);
    for(size_t i=0;i<4;i++) x[i]=0.0;
    double fmin;
    multi_funct mf=rosenbrock;

    ce.mmin(4,x,fmin,mf);

    cout << "cma_es, rosenbrock, IPOP: fmin=" << fmin
         << ", n_restarts=" << ce.n_restarts << endl;
    t.test_gen(fmin<1.0e-3,"rosenbrock fmin (IPOP)");
    t.test_gen(ce.n_restarts>0,"IPOP performed at least one restart");
  }

  // Rosenbrock with BIPOP restarts
  {
    cma_es<> ce;
    ce.rg.set_seed(4);
    ce.restart_mode=cma_es<>::cma_es_restart_bipop;
    ce.lambda=6;
    ce.max_evals=40000;

    ubvector x(4);
    for(size_t i=0;i<4;i++) x[i]=0.0;
    double fmin;
    multi_funct mf=rosenbrock;

    ce.mmin(4,x,fmin,mf);

    cout << "cma_es, rosenbrock, BIPOP: fmin=" << fmin
         << ", n_restarts=" << ce.n_restarts << endl;
    t.test_gen(fmin<1.0e-2,"rosenbrock fmin (BIPOP)");
    t.test_gen(ce.n_restarts>0,"BIPOP performed at least one restart");
  }

#ifdef O2SCL_SET_EIGEN

  // Same sphere test, but with the Eigen-backed eigendecomposition
  {
    cma_es_eigen<> ce;
    ce.rg.set_seed(5);

    ubvector x(5);
    for(size_t i=0;i<5;i++) x[i]=1.0;
    double fmin;
    multi_funct mf=sphere;

    ce.mmin(5,x,fmin,mf);

    cout << "cma_es_eigen, sphere, no restarts: fmin=" << fmin << endl;
    t.test_rel(fmin,0.0,1.0e-6,"sphere fmin (eigen)");
    t.test_gen(((std::string)ce.type())==((std::string)"cma_es_eigen"),
              "cma_es_eigen::type()");
  }

#endif

  t.report();

  return 0;
}
