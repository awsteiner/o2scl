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
#include <o2scl/boson_rel.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  double T, t1, t2, t3, t4, t5;

  boson b3(1.0,2.0);
  boson b(1.0,2.0);

  t.test_rel(b3.m,1.0,1.0e-6,"mass inheritance");

  boson_rel rb;

  b.non_interacting=true;
  b3.non_interacting=true;

  cout << "----------------------------------------------------" << endl;
  cout << "boson_rel tests:" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
    
  cout << "Degenerate: " << endl;
  b.m=1.1;
  b.mu=1.0;
  T=0.3;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  b3.m=1.1;
  b3.mu=1.0;
  T=0.3;
  cout << "boson_rel: calc_mu(T) vs. calc_density(T)" << endl;
  rb.calc_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  rb.calc_density(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl;
  
  cout << "Non-degenerate (large mass): " << endl;
  b.m=1.0;
  b.mu=0.11;
  T=1.0;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  b3.m=1.0;
  b3.mu=0.11;
  T=1.0;
  cout << "boson_rel: calc_mu(T) vs. calc_density(T)" << endl;
  rb.calc_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  rb.calc_density(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl;
  cout << endl;

  cout << "Non-degenerate (small mass): " << endl;
  b.m=0.11;
  b.mu=0.1;
  T=1.0;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  b.m=0.11;
  b.mu=0.1;
  T=1.0;
  cout << "boson_rel: pair_mu(T) vs. pair_density(T)" << endl;
  rb.pair_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  rb.pair_density(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl;
  cout << endl;
  /*
    rb.def_dit.tol_rel*=1.0e2;
    rb.def_dit.tol_abs*=1.0e2;
  */

  if (0) {
    part_calibrate_class pcc;
    double v1=pcc.part_calibrate<boson,boson_rel>
      (b,rb,0,"../../data/o2scl/boson_cal.o2",false,true,2,true);
  }
  
  t.report();
  return 0;

}
