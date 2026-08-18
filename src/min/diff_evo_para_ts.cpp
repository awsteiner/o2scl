/*
   ───────────────────────────────────────────────────────────────────

   Copyright (C) 2006-2026, Andrew W. Steiner and Edwin van Leeuwen

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

#include <gsl/gsl_sf_bessel.h>

#include <o2scl/multi_funct.h>
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/diff_evo_adapt.h>
#include <o2scl/diff_evo_para.h>
#include <o2scl/rng.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

// A simple function with many local minima. A "greedy" minimizer
// would likely fail to find the correct minimum.
double func(size_t nvar, const ubvector &x) {
  double a, b;
  a=(x[0]-2.0);
  b=(x[1]+3.0);
  return -gsl_sf_bessel_J0(a)*gsl_sf_bessel_J0(b);
}

rng<> gr;

int init_function(size_t dim, const ubvector &x, ubvector &y) {
  for (size_t i=0;i<dim;++i) {
    y[i]=20*gr.random()-10;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  multi_funct fx=func;
  mm_funct init_f=init_function;

  // ───────────────────────────────────────────────────────────
  // First, verify that diff_evo_adapt's new opt-in generational
  // (synchronous) update mode still finds the correct minimum,
  // independently of the OpenMP-parallel class below. This mode
  // defaults to off, so this also exercises a code path that
  // isn't touched by the pre-existing diff_evo_adapt_ts.cpp test.
  {
    diff_evo_adapt<multi_funct> de;
    double result;
    ubvector init(2);

    de.set_init_function(init_f);
    de.verbose=0;
    de.ntrial=1000;
    de.generational=true;

    de.mmin(2,init,result,fx);

    t.test_rel(init[0],2.0,1.0e-2,"generational - value 0");
    t.test_rel(init[1],-3.0,1.0e-2,"generational - value 1");
    t.test_rel(result,-1.0,1.0e-2,"generational - min");
  }

  // ───────────────────────────────────────────────────────────
  // Now test diff_evo_para with a range of thread counts,
  // including 1 (which should behave like the serial generational
  // update above) and several values greater than 1, to check
  // that parallelizing the per-generation trial evaluation doesn't
  // change the fact that the global minimum is found.
  for (size_t n_threads : {1, 2, 4}) {

    diff_evo_para<multi_funct> de;
    double result;
    ubvector init(2);

    de.set_init_function(init_f);
    de.verbose=0;
    de.ntrial=1000;
    de.n_threads=n_threads;

    if (n_threads==2) de.verbose=1;
    else de.verbose=0;
    de.mmin(2,init,result,fx);

    string desc=((string)"n_threads=")+itos(n_threads);
    t.test_rel(init[0],2.0,1.0e-2,(desc+" - value 0").c_str());
    t.test_rel(init[1],-3.0,1.0e-2,(desc+" - value 1").c_str());
    t.test_rel(result,-1.0,1.0e-2,(desc+" - min").c_str());

    // n_threads should never be reduced below what was requested
    // when built without OpenMP support it is silently forced to 1,
    // which is also a valid outcome here since n_threads<=4 in every
    // case tested
    t.test_gen(de.n_threads>=1,(desc+" - n_threads positive").c_str());
  }

  // ───────────────────────────────────────────────────────────
  // Check that leaving pop_size at its default (0) still produces
  // a population size which is both at least the usual 10*nvar
  // dimensionality-based value and a multiple of the number of
  // requested threads, by re-running with a larger, deliberately
  // awkward thread count and confirming convergence is still
  // obtained (an uneven, non-multiple pop_size would leave some
  // threads idle every generation but should not otherwise affect
  // correctness; this is primarily a smoke test that the
  // select_pop_size() override doesn't crash or infinite-loop for
  // odd thread counts).
  {
    diff_evo_para<multi_funct> de;
    double result;
    ubvector init(2);

    de.set_init_function(init_f);
    de.verbose=0;
    de.ntrial=1000;
    de.n_threads=3;

    de.mmin(2,init,result,fx);

    t.test_rel(init[0],2.0,1.0e-2,"n_threads=3 - value 0");
    t.test_rel(init[1],-3.0,1.0e-2,"n_threads=3 - value 1");
    t.test_rel(result,-1.0,1.0e-2,"n_threads=3 - min");
  }

  t.report();

  return 0;
}
