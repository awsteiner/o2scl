/*
  ───────────────────────────────────────────────────────────────────

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

  ───────────────────────────────────────────────────────────────────
*/
#include <iostream>
#include <cmath>

#include <o2scl/nucmass_fit_iso.h>
#include <o2scl/nucmass_two_interp.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

// A fast, non-LibTorch chain formula (base + interpm_idw residual)
// -- this deliberately avoids LibTorch so this test can run in any
// O2scl build, not just one configured with -DO2SCL_SET_LIBTORCH.
// The LibTorch-specific piece of nucmass_fit_iso's parallelism (the
// fix to interpm_libtorch's global RNG seeding) is instead verified
// by hand against the nucML project, since interpm_libtorch.h isn't
// compiled as part of the normal O2scl library or test suite.
typedef nucmass_two_interp<nucmass_semi_empirical,interpm_idw<> > chain_t;

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  nucmass_ame ame;
  ame.load("20");

  int minZ=20, maxZ=50;

  // Serial reference run
  chain_t proto1;
  nucmass_fit_iso nfi1;
  nfi1.n_threads=1;
  size_t n_fit1=nfi1.fit(proto1,ame,minZ,maxZ);
  size_t n_interp1=nfi1.fit_interp_pass(ame);

  // Same fit, with 4 OpenMP threads (or serially, if O2scl wasn't
  // built with OpenMP support -- nucmass_fit_iso::n_threads is
  // silently forced back to 1 in that case)
  chain_t proto2;
  nucmass_fit_iso nfi2;
  nfi2.n_threads=4;
  size_t n_fit2=nfi2.fit(proto2,ame,minZ,maxZ);
  size_t n_interp2=nfi2.fit_interp_pass(ame);

  t.test_gen(n_fit1==n_fit2,"same number of chains fit");
  t.test_gen(n_interp1==n_interp2,
            "same number of interpolators trained");

  // nucmass_fit's simplex minimizer is deterministic (no RNG
  // involved for this chain_t, since interpm_idw has no random
  // training step either), so a correctly race-free parallel fit
  // should give numerically identical results to the serial one,
  // regardless of how many threads were used or how chains happened
  // to be distributed across them.
  int n_checked=0, n_mismatch=0;
  double max_diff=0.0;
  for(int Z=minZ;Z<=maxZ;Z++) {
    for(int N=0;N<150;N++) {
      bool inc1=nfi1.is_included(Z,N);
      bool inc2=nfi2.is_included(Z,N);
      if (inc1!=inc2) {
        n_mismatch++;
        continue;
      }
      if (inc1) {
        double diff=fabs(nfi1.mass_excess(Z,N)-nfi2.mass_excess(Z,N));
        if (diff>max_diff) max_diff=diff;
        n_checked++;
      }
    }
  }

  t.test_gen(n_checked>100,"a reasonable number of nuclei were checked");
  t.test_gen(n_mismatch==0,"no is_included() mismatches");
  t.test_abs(max_diff,0.0,1.0e-8,
            "serial and n_threads=4 predictions match exactly");

  t.report();

  return 0;
}
