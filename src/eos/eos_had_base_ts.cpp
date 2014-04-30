/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/mroot_hybrids.h>

#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/eos_had_rmf.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  eos_had_rmf he, *he2;
  test_mgr t;
  t.set_output_level(2);
  int vp=0;

  fermion n(939.0/hc_mev_fm,2.0), p(939.0/hc_mev_fm,2.0);

  he.eoa=1.0;
  he.n0=2.0;
  he.comp=3.0;
  he.esym=4.0;
  he.msom=5.0;
  he.kprime=6.0;

  thermo th;
  th.ed=1.0;
  th.pr=2.0;
  th.en=3.0;

  he.set_n_and_p(n,p);
  he.set_thermo(th);

  t.report();
  return 0;

}

