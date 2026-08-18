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
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_two_interp.h>
#include <o2scl/nucmass_fit_iso.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  nucmass_ame ame;
  ame.load("20");

  // ------------------------------------------------------------
  // Check that clone() produces genuinely independent copies:
  // fitting one clone's base formula must not change another
  // clone's (or the prototype's) base-formula parameters. This is
  // the behavior that a naive shallow/pointer-based design would
  // get wrong. (The interpolator side isn't exercised here, since
  // interpm_idw::eval() requires set_data() to have been called
  // first -- that's covered below, via fit_interp_pass().)

  nucmass_two_interp<nucmass_dz_fit,interpm_idw<> > proto;

  nucmass_fit_base *c1=proto.clone();
  nucmass_fit_base *c2=proto.clone();

  nucmass_fit_base::ubvector x2_before(c2->nfit), x2_after(c2->nfit);
  c2->guess_fun(c2->nfit,x2_before);

  std::vector<nucleus> dist;
  nucdist_set(dist,ame,"Z>=46 && Z<=54");
  nucmass_fit nf;
  nf.dist=dist;
  // Whether or not this particular fit converges is irrelevant to
  // what's being tested here (clone() independence), so don't let a
  // non-convergence abort the test.
  //nf.def_mmin.err_nonconv=false;
  double res;
  nf.fit(*c1,res);

  c2->guess_fun(c2->nfit,x2_after);

  bool unchanged=true;
  for(size_t i=0;i<c2->nfit;i++) {
    if (x2_before[i]!=x2_after[i]) unchanged=false;
  }
  t.test_gen(unchanged,"clone() produces independent copies");

  delete c1;
  delete c2;

  // ------------------------------------------------------------
  // Exercise the full two-stage (base formula, then interpolator
  // on the residual) fit, driven through nucmass_fit_iso over a
  // small Z range for speed.

  nucmass_fit_iso nfi;
  nfi.x=3;
  nfi.mf.def_mmin.ntrial=60000;

  nucmass_two_interp<nucmass_dz_fit,interpm_idw<> > ti_proto;
  size_t n_fit=nfi.fit(ti_proto,ame,48,52);
  t.test_gen(n_fit>0,"at least one chain fit");

  size_t n_interp=nfi.fit_interp_pass(ame);
  t.test_gen(n_interp==n_fit,
             "fit_interp_pass() trains every fitted chain");

  // Every checked (Z,N) below is one of fit_interp_pass()'s own
  // training points, and interpm_idw is exact at its training
  // points, so the two-stage fit should reproduce them to within
  // numerical noise, regardless of how well the base formula's
  // simplex fit itself converged.
  size_t n_checked=0;
  double sum_sq=0.0;
  for(int Z=49;Z<=51;Z++) {
    for(int N=Z;N<Z+10;N++) {
      if (nfi.is_included(Z,N) && ame.is_included(Z,N)) {
        double me_fit=nfi.mass_excess(Z,N);
        double me_exp=ame.mass_excess(Z,N);
        t.test_gen(std::isfinite(me_fit),"two-interp mass is finite");
        sum_sq+=(me_fit-me_exp)*(me_fit-me_exp);
        n_checked++;
      }
    }
  }
  t.test_gen(n_checked>0,"at least one (Z,N) checked");
  if (n_checked>0) {
    double rms=sqrt(sum_sq/((double)n_checked));
    t.test_gen(rms<1.0e-6,
               "two-interp chains reproduce training data (RMS)");
  }

  t.report();
  return 0;
}
