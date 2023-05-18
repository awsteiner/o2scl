/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/fermion.h>
#include <o2scl/eos_had_ddc.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  fermion n(939.0/hc_mev_fm,2.0);
  fermion p(939.0/hc_mev_fm,2.0);
  n.non_interacting=false;
  p.non_interacting=false;

  eos_had_ddc ddc;

  n.n=0.153/2.0;
  p.n=0.153/2.0;

  double sig=0.197375;
  double ome=0.129147;
  double rho=0.0;
  double f1, f2, f3;
  thermo th;
  
  ddc.calc_eq_e(n,p,sig,ome,rho,f1,f2,f3,th);
  
  t.test_rel(n.ms/n.m,0.555,5.0e-5,"n.ms");
  t.test_rel(p.ms/p.m,0.555,5.0e-5,"p.ms");
  t.test_rel((n.mu-n.m)*hc_mev_fm,-16.2,5.0e-3,"mun");
  t.test_rel((p.mu-p.m)*hc_mev_fm,-16.2,5.0e-3,"mup");
  t.test_rel(f1,0.0,1.0e-3,"f1");
  t.test_rel(f2,0.0,1.0e-3,"f2");
  t.test_rel(f3,0.0,1.0e-3,"f3");
  t.test_rel(th.pr,0.0,1.0e-4,"pr");
  t.test_rel((th.ed/0.153-n.m)*hc_mev_fm,-16.247,1.0e-4,"ed");
  
  t.report();
  return 0;
}
