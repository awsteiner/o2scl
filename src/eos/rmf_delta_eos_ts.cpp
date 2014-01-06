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
#include <o2scl/constants.h>
#include <o2scl/fermion.h>
#include <o2scl/rmf_delta_eos.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  rmf_delta_eos red;
  double gs;
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  fermion n(939.0/hc_mev_fm,2.0), p(939.0/hc_mev_fm,2.0);
  n.non_interacting=false;
  p.non_interacting=false;
  red.set_n_and_p(n,p);
  
  cout << "Test Set I and Set II from PRC 65, 045201:" << endl;

  red.mnuc=939.0/hc_mev_fm;
  red.ms=550.0/hc_mev_fm;
  red.mr=770.0/hc_mev_fm;
  red.mw=783.0/hc_mev_fm;
  red.md=980.0/hc_mev_fm;

  red.n0=0.16;
  red.comp=240.0/hc_mev_fm;
  red.eoa=-16.0/hc_mev_fm;
  red.esym=30.5/hc_mev_fm;
  red.msom=0.75;
  n.mu=4.8;
  p.mu=4.8;
  red.set_fields(0.1,0.07,-0.001,0.0);
  red.fix_saturation();
  t.test_rel(red.cs*red.cs,10.33,1.0e-3,"1");
  t.test_rel(red.cw*red.cw,5.42,1.0e-3,"2");
  t.test_rel(red.cr*red.cr/4.0,0.95,4.0e-2,"3");
  t.test_rel(red.b*red.mnuc,0.033,1.0e-3,"4");
  t.test_rel(red.c,-0.0048,1.0e-2,"5");
  
  red.cr=sqrt(3.15*4.0);
  red.cd=sqrt(2.50);
   
  n.mu=4.8;
  p.mu=4.8;
  red.set_n_and_p(n,p);
  red.set_fields(0.1,0.07,-0.001,0.0);
  red.saturation();
  t.test_rel(red.fesym(red.n0)*hc_mev_fm,30.5,1.0e-2,"6");

  t.report();
  cout << endl;
  
  cout << "Test \"stiff\" model from PLB 388 191." << endl;

  red.cs=sqrt(11.25);
  red.cw=sqrt(6.483);
  red.cr=sqrt(10.0);
  gs=red.cs*red.ms;
  red.cd=sqrt(1.3);
  red.b=0.003825;
  red.c=3.5e-6;

  n.mu=4.8;
  p.mu=4.8;
  red.set_n_and_p(n,p);
  red.set_fields(0.1,0.07,-0.001,0.0);
  red.saturation();

  t.test_rel(red.n0,0.145,4.0e-3,"7");
  t.test_rel(red.fesym(red.n0)*hc_mev_fm,34.0,1.0e-1,"8");

  t.report();

  return 0;
}

