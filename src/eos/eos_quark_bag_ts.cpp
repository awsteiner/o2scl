/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_quark_bag.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  eos_quark_bag ebg;
  thermo th;
  for(double nB=1.0e-5;nB<2.0;nB*=2.0) {
    for(double T=0.1;T<100.1;T*=10.0) {
      for(double nQ=1.0e-5;nB<2.0;nB*=2.0) {
        for(double nS=0.0;nS<1.0/3.0+1.0e-4;nS+=1.0/12.0) {
          int ret=ebg.calc_temp_f_gen(nB,nQ,nS,T,th);
          cout << nB << " " << nQ << " " << nS << " " << T << " "
               << th.ed << " " << ret << endl;
        }
      }
    }
  }

  t.report();
  return 0;
}

