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
#include <o2scl/test_mgr.h>

#include <o2scl/nonrel_fermion.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  test_mgr t;
  t.set_output_level(1);
  
  double t1, t2, t3, t4, T;

  nonrel_fermion nrf;
  fermion e2(5.0,2.0);
  t.test_rel(e2.m,5.0,1.0e-6,"mass inheritance");
  e2.non_interacting=true;

  cout.setf(ios::scientific);

  // -----------------------------------------------------------------
  // nonrel_fermion tests

  cout << "----------------------------------------------------" << endl;
  cout << "nonrel_fermion tests:" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
    
  cout << "calc_mu_zerot() vs. calc_mu(T) vs. calc_density(T) vs. "
       << "Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=5.1;
  T=0.01;
  // math2c den1
  t1=0.034194453306012236;
  // math2c end
  // math2c ed1
  t2=0.17312331286395766;
  // math2c end
  // math2c en1
  t3=0.016563189255283155;
  // math2c end
  // math2c pr1
  t4=0.0014340308892575593;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu_zerot(e2);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.pr << " " 
       << e2.en << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.pr << " " 
       << e2.en << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.pr << " " 
       << e2.en << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t4 << " " 
       << t3 << endl;
  cout << endl;

  cout << "calc_mu(T) vs. calc_density(T) vs. Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=5.1;
  T=0.1;
  // math2c den2
  t1=0.07074119812209832;
  // math2c end
  // math2c ed2
  t2=0.3671902323113951;
  // math2c end
  // math2c en2
  t3=0.1539961635579667;
  // math2c end
  // math2c pr2
  t4=0.00898949446710301;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t3 << " " 
       << t4 << endl;
  cout << endl;

  cout << "calc_mu(T) vs. calc_density(T) vs. Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=5.1;
  T=1.0;
  // math2c den3
  t1=1.1749327822805162;
  // math2c end
  // math2c ed3
  t2=7.891014509994989;
  // math2c end
  // math2c en3
  t3=3.2430910527220838;
  // math2c end
  // math2c pr3
  t4=1.3442337323577283;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t3 << " " 
       << t4 << endl;
  cout << endl;

  cout << "calc_mu(T) vs. calc_density(T) vs. Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=6.0;
  T=0.05;
  // math2c den4
  t1=1.0713203764363604;
  // math2c end
  // math2c ed4
  t2=6.007288246759999;
  // math2c end
  // math2c en4
  t3=0.263137956873924;
  // math2c end
  // math2c pr4
  t4=0.43379090970185885;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t3 << " " 
       << t4 << endl;
  cout << endl;

  cout << "calc_mu(T) vs. calc_density(T) vs. Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=6.0;
  T=0.1;
  // math2c den5
  t1=1.0813235579115335;
  // math2c end
  // math2c ed5
  t2=6.08683836632324;
  // math2c end
  // math2c en5
  t3=0.5237740336317366;
  // math2c end
  // math2c pr5
  t4=0.45348038450913464;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t3 << " " 
       << t4 << endl;
  cout << endl;

  cout << "calc_mu(T) vs. calc_density(T) vs. Mathematica: " << endl;
  e2.m=5.0;
  e2.mu=6.0;
  T=1.0;
  // math2c den6
  t1=2.2370331047788405;
  // math2c end
  // math2c ed6
  t2=15.449257153402318;
  // math2c end
  // math2c en6
  t3=4.8697862776430725;
  // math2c end
  // math2c pr6
  t4=2.842727752913799;
  // math2c end
  cout << "(m=" << e2.m << ", mu=" << e2.mu << ", T=" << T << ")" << endl;
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(t1,e2.n,5.0e-9,"density");
  t.test_rel(t2,e2.ed,5.0e-9,"energy density");
  t.test_rel(t3,e2.en,5.0e-9,"entropy");
  t.test_rel(t4,e2.pr,5.0e-9,"pressure");
  cout << t1 << " " << e2.mu << " " << t2 << " " << t3 << " " 
       << t4 << endl;
  cout << endl;

  cout << "Compare with inc_rest_mass=false: " << endl;
  cout << e2.n << " " << e2.mu-e2.m << " " << e2.ed-e2.n*e2.m 
       << " " << e2.en << " " << e2.pr << endl;
  t1=e2.n;
  t2=e2.mu-e2.m;
  t3=e2.ed-e2.n*e2.m;
  t4=e2.en;
  e2.inc_rest_mass=false;
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(e2.n,t1,5.0e-9,"density");
  t.test_rel(e2.mu,t2,5.0e-9,"chemical potential");
  t.test_rel(e2.ed,t3,5.0e-9,"energy density");
  t.test_rel(e2.en,t4,5.0e-9,"entropy");
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.mu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(e2.n,t1,5.0e-9,"density");
  t.test_rel(e2.mu,t2,5.0e-9,"chemical potential");
  t.test_rel(e2.ed,t3,5.0e-9,"energy density");
  t.test_rel(e2.en,t4,5.0e-9,"entropy");
  cout << endl;

  cout << "Compare with inc_rest_mass=false when non_interacting is false: " 
       << endl;
  e2.non_interacting=false;
  e2.inc_rest_mass=false;
  nrf.calc_density(e2,T);
  cout << e2.n << " " << e2.nu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(e2.n,t1,5.0e-9,"density");
  t.test_rel(e2.mu,t2,5.0e-9,"chemical potential");
  t.test_rel(e2.ed,t3,5.0e-9,"energy density");
  t.test_rel(e2.en,t4,5.0e-9,"entropy");
  nrf.calc_mu(e2,T);
  cout << e2.n << " " << e2.nu << " " << e2.ed << " " << e2.en << " " 
       << e2.pr << endl;
  t.test_rel(e2.n,t1,5.0e-9,"density");
  t.test_rel(e2.mu,t2,5.0e-9,"chemical potential");
  t.test_rel(e2.ed,t3,5.0e-9,"energy density");
  t.test_rel(e2.en,t4,5.0e-9,"entropy");
  cout << endl;

  t.set_output_level(2);
  t.report();

  return 0;
}
