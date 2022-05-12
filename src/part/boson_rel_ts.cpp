/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/boson_rel.h>
#include <o2scl/boson_eff.h>
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
  boson_eff eb;

  b.non_interacting=true;
  b3.non_interacting=true;

  cout << "----------------------------------------------------" << endl;
  cout << "boson_rel tests:" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
    
  cout << "Degenerate: " << endl;
  b.m=1.1;
  b.mu=1.0;
  T=0.1;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  cout << "boson_eff: calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl; 
  t1=b.n; t2=b.mu; t3=b.ed; t4=b.pr; t5=b.en;
  eb.calc_density(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  t.test_rel(b.n,t1,1.0e-10,"deg eb calc_mu vs calc_density density");
  t.test_rel(b.mu,t2,1.0e-10,"deg eb calc_mu vs calc_density chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"deg eb calc_mu vs calc_density energy");
  t.test_rel(b.pr,t4,1.0e-10,"deg eb calc_mu vs calc_density pressure");
  t.test_rel(b.en,t5,1.0e-10,"deg eb calc_mu vs calc_density entropy");
  cout << endl;
  
  b3.m=1.1;
  b3.mu=1.0;
  T=0.1;
  cout << "boson_rel: calc_mu(T) vs. calc_density(T)" << endl;
  rb.calc_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  t.test_rel(b.n,t1,1.0e-10,"deg rb calc_mu vs calc_density density");
  t.test_rel(b.mu,t2,2.0e-10,"deg rb calc_mu vs calc_density chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"deg rb calc_mu vs calc_density energy");
  t.test_rel(b.pr,t4,1.0e-10,"deg rb calc_mu vs calc_density pressure");
  t.test_rel(b.en,t5,1.0e-10,"deg rb calc_mu vs calc_density entropy");
  rb.calc_density(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl;
  t.test_rel(b.n,t1,1.0e-10,"deg eb vs. rb density");
  t.test_rel(b.mu,t2,2.0e-10,"deg eb vs. rb chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"deg eb vs. rb energy");
  t.test_rel(b.pr,t4,1.0e-10,"deg eb vs. rb pressure");
  t.test_rel(b.en,t5,1.0e-10,"deg eb vs. rb entropy");
  cout << endl;

  cout << "Non-degenerate (large mass): " << endl;
  b.m=1.0;
  b.mu=0.11;
  T=1.0;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  cout << "boson_eff: calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl; 
  t1=b.n; t2=b.mu; t3=b.ed; t4=b.pr; t5=b.en;
  eb.calc_density(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  t.test_rel(b.n,t1,1.0e-10,"ndeg lm eb calc_mu vs calc_density density");
  t.test_rel(b.mu,t2,2.0e-10,"ndeg lm eb calc_mu vs calc_density chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"ndeg lm eb calc_mu vs calc_density energy");
  t.test_rel(b.pr,t4,1.0e-10,"ndeg lm eb calc_mu vs calc_density pressure");
  t.test_rel(b.en,t5,1.0e-10,"ndeg lm eb calc_mu vs calc_density entropy");
  cout << endl;
  
  b3.m=1.0;
  b3.mu=0.11;
  T=1.0;
  cout << "boson_rel: calc_mu(T) vs. calc_density(T)" << endl;
  rb.calc_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  t.test_rel(b.n,t1,1.0e-10,"ndeg lm rb calc_mu vs calc_density density");
  t.test_rel(b.mu,t2,2.0e-10,"ndeg lm rb calc_mu vs calc_density chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"ndeg lm rb calc_mu vs calc_density energy");
  t.test_rel(b.pr,t4,1.0e-10,"ndeg lm rb calc_mu vs calc_density pressure");
  t.test_rel(b.en,t5,1.0e-10,"ndeg lm rb calc_mu vs calc_density entropy");
  rb.calc_density(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl;
  t.test_rel(b.n,t1,1.0e-10,"ndeg lm rb vs. eb density");
  t.test_rel(b.mu,t2,2.0e-10,"ndeg lm rb vs. eb chem. pot.");
  t.test_rel(b.ed,t3,1.0e-10,"ndeg lm rb vs. eb energy");
  t.test_rel(b.pr,t4,1.0e-10,"ndeg lm rb vs. eb pressure");
  t.test_rel(b.en,t5,1.0e-10,"ndeg lm rb vs. eb entropy");
  cout << endl;

  cout << "Non-degenerate (small mass): " << endl;
  b.m=0.11;
  b.mu=0.1;
  T=1.0;
  cout << "(m=" << b.m << ", mu=" << b.mu << ", T=" << T << ")" << endl;
  cout << endl;
  
  cout << "boson_eff: pair_mu(T) vs. pair_density(T)" << endl;
  eb.pair_mu(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl; 
  t1=b.n; t2=b.mu; t3=b.ed; t4=b.pr; t5=b.en;
  eb.pair_density(b,T);
  cout << b.n << " " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  //  t.test_rel(b.n,t1,1.0e-10,"calc_mu vs calc_density density");
  //  t.test_rel(b.mu,t2,1.0e-10,"calc_mu vs calc_density chem. pot.");
  //  t.test_rel(b.ed,t3,1.0e-10,"calc_mu vs calc_density energy");
  //  t.test_rel(b.pr,t4,1.0e-10,"calc_mu vs calc_density pressure");
  //  t.test_rel(b.en,t5,1.0e-10,"calc_mu vs calc_density entropy");
  b.m=0.11;
  b.mu=0.1;
  T=1.0;
  /*
  cout << "boson_rel: pair_mu(T) vs. pair_density(T)" << endl;
  ret=rb.pair_mu(b3,T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
       << b3.pr << " " << b3.en << endl; 
  t1=b.n; t2=b.mu; t3=b.ed; t4=b.pr; t5=b.en;
  ret=b3.pair_density(T);
  cout << b3.n << " " << b3.mu << " " << b3.ed << " " 
  << b3.pr << " " << b3.en << endl;
  */
  /*
    t.test_rel(b.n,t1,1.0e-10,"density");
    t.test_rel(b.mu,t2,1.0e-10,"chem. pot.");
    t.test_rel(b.ed,t3,1.0e-10,"energy");
    t.test_rel(b.pr,t4,1.0e-10,"pressure");
    t.test_rel(b.en,t5,1.0e-10,"entropy");
    cout << endl;
  */

  /*
    rb.def_dit.tol_rel*=1.0e2;
    rb.def_dit.tol_abs*=1.0e2;
    part_calibrate(b,rb,0,"../../data/o2scl/boson_cal.o2",false,1,1);
  */
  
  t.report();
  return 0;

}
