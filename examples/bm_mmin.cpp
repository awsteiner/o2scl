/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2026, Andrew W. Steiner

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
/*
  Benchmark and compare the multidimensional minimizers in src/min ,
  following the same style as bm_min.cpp (timing statistics via
  expval_scalar) and using the current minimizer API and test
  function, as in examples/ex_mmin.cpp .
*/

#include <cmath>
#include <ctime>
#include <functional>

#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/constants.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/mmin_conf.h>
#include <o2scl/mmin_conp.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/diff_evo.h>
#include <o2scl/diff_evo_adapt.h>
#include <o2scl/cma_es.h>
#include <o2scl/rng.h>
#include <o2scl/expval.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

/** \brief A difficult "spring" function to minimize

    The global minimum is at (1,0,0), with f(1,0,0)=1. This is the
    same test function used in the original version of this file;
    see the <tt>spring_two()</tt> function in examples/ex_mmin.cpp
    for a related, harder variant with an adjustable parameter.
*/
double spring(size_t nv, const ubvector &x) {
  double theta=atan2(x[1],x[0]);
  double r=hypot(x[0],x[1]);
  double z=x[2];
  while (z>pi) z-=2.0*pi;
  while (z<-pi) z+=2.0*pi;
  double tmz=theta-z;
  double rm1=r-1.0;
  return exp(tmz*tmz+rm1*rm1)+fabs(x[2]/10.0);
}

/** \brief Run \c n_reps randomized trials of minimizer \c mn on the
    \ref spring() function, reporting average wall-clock time,
    function evaluations, and iterations (with block-averaged
    error estimates from \ref o2scl::expval_scalar), and the
    fraction of trials which converged close to the true minimum

    If \c max_evals_budget is nonzero, it is used to set \c
    mn.max_evals before each trial, following the max_evals/
    last_n_evals convention described in the documentation of
    \ref o2scl::mmin_base -- this makes it possible to compare
    very different minimizer classes (e.g. \ref o2scl::mmin_simp2
    and \ref o2scl::cma_es) on a common, algorithm-neutral
    evaluation budget, rather than only on their (very differently
    priced) iteration counts.
*/
template<class min_t>
void bench_alg(std::string name, min_t &mn, multi_funct &mf,
                rng<> &gr, size_t n_reps, size_t max_evals_budget) {

  expval_scalar time_ev, evals_ev, iters_ev;
  time_ev.set_blocks(10,1);
  evals_ev.set_blocks(10,1);
  iters_ev.set_blocks(10,1);

  mn.max_evals=max_evals_budget;
  mn.err_nonconv=false;

  size_t n_success=0;

  for(size_t k=0;k<n_reps;k++) {

    ubvector x(3);
    x[0]=1.0+(gr.random()-0.5)*0.1;
    x[1]=0.0+(gr.random()-0.5)*0.1;
    x[2]=2.0*pi+(gr.random()-0.5)*0.4;

    double fmin;
    clock_t clk1=clock();
    mn.mmin(3,x,fmin,mf);
    clock_t clk2=clock();

    time_ev.add(((double)(clk2-clk1))/CLOCKS_PER_SEC);
    evals_ev.add((double)mn.last_n_evals);
    iters_ev.add((double)mn.last_ntrial);

    double dist=sqrt(pow(x[0]-1.0,2.0)+x[1]*x[1]+x[2]*x[2]);
    if (dist<0.05) n_success++;
  }

  double avg, sd, avge;

  cout.width(16);
  cout << std::left << name << std::right;
  cout.setf(ios::scientific);
  cout.precision(3);

  time_ev.current_avg(avg,sd,avge);
  cout << " " << avg << " +/- " << avge << " ";

  evals_ev.current_avg(avg,sd,avge);
  cout << avg << " +/- " << avge << " ";

  iters_ev.current_avg(avg,sd,avge);
  cout << avg << " +/- " << avge << " ";

  cout.unsetf(ios::scientific);
  cout << n_success << "/" << n_reps << endl;

  return;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  rng<> gr;
  gr.clock_seed();

  multi_funct mf=spring;

  mmin_simp2<> g1;
  mmin_conf<> g2;
  mmin_conp<> g3;
  mmin_bfgs2<> g4;
  diff_evo<> g5;
  diff_evo_adapt<> g6;
  cma_es<> g7;

  // The spring function is difficult for the local minimizers, so
  // give them more trials than the default
  g1.ntrial*=4;
  g2.ntrial*=4;
  g3.ntrial*=4;
  g4.ntrial*=4;

  static const size_t n_reps=40;

  cout << "Unconstrained comparison (default ntrial, "
       << "max_evals=0, i.e. unbounded):" << endl;
  cout << "Algorithm        Avg. time (s)        "
       << "Avg. evals           Avg. iters           Success"
       << endl;

  bench_alg("mmin_simp2",g1,mf,gr,n_reps,0);
  bench_alg("mmin_conf",g2,mf,gr,n_reps,0);
  bench_alg("mmin_conp",g3,mf,gr,n_reps,0);
  bench_alg("mmin_bfgs2",g4,mf,gr,n_reps,0);
  bench_alg("diff_evo",g5,mf,gr,n_reps,0);
  bench_alg("diff_evo_adapt",g6,mf,gr,n_reps,0);
  bench_alg("cma_es",g7,mf,gr,n_reps,0);

  cout << endl;

  // Now use max_evals to put every minimizer on the same,
  // algorithm-neutral evaluation budget (see the documentation of
  // o2scl::mmin_base for why this is the right way to compare
  // minimizers -- an "iteration" costs a very different number of
  // function evaluations from one algorithm to the next). Since
  // max_evals now bounds the number of function evaluations
  // directly, ntrial is set generously so that it doesn't cut a
  // minimizer off before its evaluation budget is spent.
  size_t budget=300;

  g1.ntrial=100000;
  g2.ntrial=100000;
  g3.ntrial=100000;
  g4.ntrial=100000;
  g5.ntrial=100000;
  g6.ntrial=100000;
  g7.ntrial=100000;

  cout << "Comparison on a common evaluation budget "
       << "(max_evals=" << budget << "):" << endl;
  cout << "Algorithm        Avg. time (s)        "
       << "Avg. evals           Avg. iters           Success"
       << endl;

  bench_alg("mmin_simp2",g1,mf,gr,n_reps,budget);
  bench_alg("mmin_conf",g2,mf,gr,n_reps,budget);
  bench_alg("mmin_conp",g3,mf,gr,n_reps,budget);
  bench_alg("mmin_bfgs2",g4,mf,gr,n_reps,budget);
  bench_alg("diff_evo",g5,mf,gr,n_reps,budget);
  bench_alg("diff_evo_adapt",g6,mf,gr,n_reps,budget);
  bench_alg("cma_es",g7,mf,gr,n_reps,budget);

  // A basic sanity check that the local minimizers can still find
  // the minimum of this well-worn test function without a budget
  // constraint
  {
    ubvector x(3);
    x[0]=1.0; x[1]=0.0; x[2]=2.0*pi;
    double fmin;
    mmin_simp2<> gs;
    gs.ntrial*=4;
    gs.mmin(3,x,fmin,mf);
    t.test_rel(x[0],1.0,1.0e-3,"mmin_simp2 x[0]");
    t.test_rel(fmin,1.0,1.0e-3,"mmin_simp2 fmin");
  }

  t.report();

  return 0;
}
