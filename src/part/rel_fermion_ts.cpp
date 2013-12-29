/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

#include <o2scl/rel_fermion.h>
#include <o2scl/eff_fermion.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  fermion e(1.0,2.0);
  fermion e3(1.0,2.0);
  eff_fermion ef;
  rel_fermion rf;
  t.test_rel(e3.m,1.0,1.0e-6,"mass inheritance");
  e.non_interacting=true;
  e3.non_interacting=true;

  double T;
  double t1, t2, t3, t4, t5;

  cout.setf(ios::scientific);

  // -----------------------------------------------------------------
  // rel_fermion tests
  
  cout << "----------------------------------------------------" << endl;
  cout << "rel_fermion tests:" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  cout << "Non-degenerate with M/T large: " << endl;
  e.m=4.5;
  e.n=1.0e-9;
  e3.m=4.5;
  e3.n=1.0e-9;
  T=1.0e-3;
  ef.calc_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  ef.calc_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  
  e3.mu=0.0;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  cout << "0 " << rf.unc.n << " " << rf.unc.mu << " " << rf.unc.ed << " " 
       << rf.unc.pr << " " << rf.unc.en << endl;

  cout << "Degenerate: " << endl;
  e.m=0.1;
  e.mu=1.0;
  T=0.1;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  ef.calc_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.calc_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-10,"density 1");
  t.test_rel(t2,e.mu,1.0e-10,"chem. pot. 1");
  t.test_rel(t3,e.ed,1.0e-10,"energy 1");
  t.test_rel(t4,e.pr,1.0e-10,"pressure 1");
  t.test_rel(t5,e.en,1.0e-10,"entropy 1");
  
  e3.m=0.1;
  e3.mu=1.0;
  T=0.1;

  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-9,"density 2");
  t.test_rel(t2,e3.mu,1.0e-9,"chem. pot. 2");
  t.test_rel(t3,e3.ed,1.0e-9,"energy 2");
  t.test_rel(t4,e3.pr,1.0e-9,"pressure 2");
  t.test_rel(t5,e3.en,1.0e-9,"entropy 2");
  cout << endl;

  cout << "Non-degenerate: " << endl;
  e.m=0.1;
  e.mu=0.11;
  T=1.0;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  ef.calc_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.calc_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-9,"density 3");
  t.test_rel(t2,e.mu,1.0e-9,"chem. pot. 3");
  t.test_rel(t3,e.ed,1.0e-9,"energy 3");
  t.test_rel(t4,e.pr,1.0e-9,"pressure 3");
  t.test_rel(t5,e.en,1.0e-9,"entropy 3");

  e3.m=0.1;
  e3.mu=0.11;
  T=1.0;
  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 4");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 4");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 4");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 4");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 4");
  cout << endl;

  cout << "Zero-temperature: " << endl;
  e3.m=0.1;
  e3.mu=0.11;
  T=0.001;
  cout << "(m=" << e3.m << ", mu=" << e3.mu << ", T=" << T << ")" << endl;
  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=e3.m;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 5");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 5");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 5");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 5");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 5");

  cout << "rel_fermion: calc_mu(T=0) vs. calc_density(T=0)" << endl;
  ef.calc_mu_zerot(e3);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  ef.calc_density_zerot(e3);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 6");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 6");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 6");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 6");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 6");
  cout << endl;

  cout << "Non-degenerate: " << endl;
  e.m=0.1;
  e.mu=0.11;
  T=1.0;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  ef.pair_mu(e,T);
  cout << endl;
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.pair_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  cout << endl;
  t.test_rel(t1,e.n,1.0e-10,"density 7");
  t.test_rel(t2,e.mu,1.0e-10,"chem. pot. 7");
  t.test_rel(t3,e.ed,1.0e-10,"energy 7");
  t.test_rel(t4,e.pr,1.0e-10,"pressure 7");
  t.test_rel(t5,e.en,1.0e-10,"entropy 7");

  e3.m=0.1;
  e3.mu=0.11;
  T=1.0;
  cout << "rel_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  rf.pair_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.pair_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 8");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 8");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 8");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 8");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 8");
  cout << endl;

  cout << "Comparing answers near deg_limit: " << endl;
  e.m=0.1;
  T=1.0;
  e.mu=(rf.deg_limit+1.0e-4)*T+e.m;
  e3.m=0.1;
  T=1.0;
  e3.mu=(rf.deg_limit+1.0e-4)*T+e3.m;
  cout << "(m=" << e3.m << ", mu=" << e3.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  ef.pair_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.pair_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-9,"density 9");
  t.test_rel(t2,e.mu,1.0e-9,"chem. pot. 9");
  t.test_rel(t3,e.ed,1.0e-9,"energy 9");
  t.test_rel(t4,e.pr,1.0e-9,"pressure 9");
  t.test_rel(t5,e.en,1.0e-9,"entropy 9");

  cout << "rel_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  rf.pair_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.pair_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-8,"density 10");
  t.test_rel(t2,e3.mu,1.0e-8,"chem. pot. 10");
  t.test_rel(t3,e3.ed,1.0e-8,"energy 10");
  t.test_rel(t4,e3.pr,1.0e-8,"pressure 10");
  t.test_rel(t5,e3.en,1.0e-8,"entropy 10");

  e.m=0.1;
  T=1.0;
  e.mu=(rf.deg_limit-1.0e-4)*T+e.m;
  e3.m=0.1;
  T=1.0;
  e3.mu=(rf.deg_limit-1.0e-4)*T+e3.m;
  cout << "(m=" << e3.m << ", mu=" << e3.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  ef.pair_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.pair_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-9,"density 11");
  t.test_rel(t2,e.mu,1.0e-9,"chem. pot. 11");
  t.test_rel(t3,e.ed,1.0e-9,"energy 11");
  t.test_rel(t4,e.pr,1.0e-9,"pressure 11");
  t.test_rel(t5,e.en,1.0e-9,"entropy 11");

  cout << "rel_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  rf.pair_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.pair_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-9,"density 12");
  t.test_rel(t2,e3.mu,1.0e-9,"chem. pot. 12");
  t.test_rel(t3,e3.ed,1.0e-9,"energy 12");
  t.test_rel(t4,e3.pr,1.0e-9,"pressure 12");
  t.test_rel(t5,e3.en,1.0e-9,"entropy 12");
  cout << endl;


  // -----------------------------------------------------------------
  // rel_fermion tests (with errors)

  cout << "----------------------------------------------------" << endl;
  cout << "rel_fermion tests (with errors):" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
    
  cout << "Degenerate: " << endl;
  e.m=0.1;
  e.mu=1.0;
  T=0.1;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  ef.calc_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.calc_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-10,"density 13");
  t.test_rel(t2,e.mu,1.0e-10,"chem. pot. 13");
  t.test_rel(t3,e.ed,1.0e-10,"energy 13");
  t.test_rel(t4,e.pr,1.0e-10,"pressure 13");
  t.test_rel(t5,e.en,1.0e-10,"entropy 13");

  e3.m=0.1;
  e3.mu=1.0;
  T=0.1;
  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t.test_rel(t1,e3.n,1.0e-9,"density 14");
  t.test_rel(t2,e3.mu,1.0e-9,"chem. pot. 14");
  t.test_rel(t3,e3.ed,1.0e-9,"energy 14");
  t.test_rel(t4,e3.pr,1.0e-9,"pressure 14");
  t.test_rel(t5,e3.en,1.0e-9,"entropy 14");
  cout << endl;

  cout << "Non-degenerate: " << endl;
  e.m=0.1;
  e.mu=0.11;
  T=1.0;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  ef.calc_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.calc_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-9,"density 15");
  t.test_rel(t2,e.mu,1.0e-9,"chem. pot. 15");
  t.test_rel(t3,e.ed,1.0e-9,"energy 15");
  t.test_rel(t4,e.pr,1.0e-9,"pressure 15");
  t.test_rel(t5,e.en,1.0e-9,"entropy 15");

  e3.m=0.1;
  e3.mu=0.11;
  T=1.0;
  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 16");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 16");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 16");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 16");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 16");
  cout << endl;

  cout << "Zero-temperature: " << endl;
  e3.m=0.1;
  e3.mu=0.11;
  T=0.001;
  cout << "(m=" << e3.m << ", mu=" << e3.mu << ", T=" << T << ")" << endl;
  cout << "rel_fermion: calc_mu(T) vs. calc_density(T)" << endl;
  rf.calc_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=e3.m;
  rf.calc_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 17");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 17");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 17");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 17");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 17");

  cout << "rel_fermion: calc_mu(T=0) vs. calc_density(T=0)" << endl;
  ef.calc_mu_zerot(e3);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  ef.calc_density_zerot(e3);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 18");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 18");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 18");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 18");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 18");
  cout << endl;

  cout << "Non-degenerate: " << endl;
  e.m=0.1;
  e.mu=0.11;
  T=1.0;
  cout << "(m=" << e.m << ", mu=" << e.mu << ", T=" << T << ")" << endl;
  cout << "eff_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  ef.pair_mu(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl; 
  t1=e.n;
  t2=e.mu;
  t3=e.ed;
  t4=e.pr;
  t5=e.en;
  e.mu=0.0;
  ef.pair_density(e,T);
  cout << e.n << " " << e.mu << " " << e.ed << " " 
       << e.pr << " " << e.en << endl;
  t.test_rel(t1,e.n,1.0e-10,"density 19");
  t.test_rel(t2,e.mu,1.0e-10,"chem. pot. 19");
  t.test_rel(t3,e.ed,1.0e-10,"energy 19");
  t.test_rel(t4,e.pr,1.0e-10,"pressure 19");
  t.test_rel(t5,e.en,1.0e-10,"entropy 19");

  e3.m=0.1;
  e3.mu=0.11;
  T=1.0;
  cout << "rel_fermion: pair_mu(T) vs. pair_density(T)" << endl;
  rf.pair_mu(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl; 
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t1=e3.n;
  t2=e3.mu;
  t3=e3.ed;
  t4=e3.pr;
  t5=e3.en;
  e3.mu=0.0;
  rf.pair_density(e3,T);
  cout << e3.n << " " << e3.mu << " " << e3.ed << " " 
       << e3.pr << " " << e3.en << endl;
  cout << "U " << rf.unc.n << " " << 0.0 << " " << rf.unc.ed << " "
       << rf.unc.pr << " " << rf.unc.en << endl;
  t.test_rel(t1,e3.n,1.0e-10,"density 20");
  t.test_rel(t2,e3.mu,1.0e-10,"chem. pot. 20");
  t.test_rel(t3,e3.ed,1.0e-10,"energy 20");
  t.test_rel(t4,e3.pr,1.0e-10,"pressure 20");
  t.test_rel(t5,e3.en,1.0e-10,"entropy 20");
  cout << endl;

  cout << "----------------------------------------------------" << endl;
  cout << "Four combinations of inc_rest_mass & non_interacting" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  double den;
  for(den=1.0e-6;den<=1.01e2;den*=1.0e2) {

    cout << "density: " << den << endl;

    // -----------------------------------------------------------------
    // Test calc_density() and pair_density() with inc_rest_mass is 
    // both true and false, and non_interacting is both true and false
    // -----------------------------------------------------------------
  
    fermion xr(0.511/hc_mev_fm,2.0);
    double temper=10.0/hc_mev_fm, d1, d2, mu;

    // -----------------------------------------------------------------

    xr.non_interacting=false;
    cout << "Interacting: " << endl;
    cout << endl;
    
    xr.inc_rest_mass=true;
    cout << "Include rest mass." << endl;
    
    xr.n=den;
    xr.nu=xr.m;
    rf.calc_density(xr,temper);
    cout << xr.nu << " ";
    t1=xr.nu;
    rf.calc_mu(xr,temper);
    cout << xr.n << " " << den << endl;
    t.test_rel(xr.n,den,1.0e-6,"den1");

    xr.n=den;
    xr.nu=xr.m;
    rf.pair_density(xr,temper);
    mu=xr.nu;
    cout << xr.nu << " ";
    t2=xr.nu;
    rf.pair_mu(xr,temper);
    cout << xr.n << " " << den << endl;
    t.test_rel(xr.n,den,1.0e-6,"den2");
    
    xr.nu=mu;
    rf.calc_mu(xr,temper);
    d1=xr.n;
    if (xr.inc_rest_mass) {
      xr.nu=-mu;
    } else {
      xr.nu=-mu-2.0*xr.m;
    }
    rf.calc_mu(xr,temper);
    d2=xr.n;
    cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
    t3=d1;
    t4=d2;
    t.test_rel(den,d1-d2,1.0e-6,"den3");
    cout << endl;

    // -----------------------------------------------------------------

    xr.non_interacting=true;
    cout << "Non-interacting: " << endl;
    cout << endl;

    xr.inc_rest_mass=true;
    cout << "Include rest mass." << endl;
    
    xr.n=den;
    xr.mu=0.0;
    rf.calc_density(xr,temper);
    cout << xr.mu << " ";
    t.test_rel(xr.nu,t1,1.0e-6,"nu1");
    rf.calc_mu(xr,temper);
    cout << xr.n << " " << den << endl;
    t.test_rel(xr.n,den,1.0e-6,"den1");
    
    xr.n=den;
    xr.mu=0.0;
    rf.pair_density(xr,temper);
    mu=xr.mu;
    cout << xr.mu << " ";
    t.test_rel(xr.nu,t2,1.0e-6,"nu2");
    rf.pair_mu(xr,temper);
    cout << xr.n << " " << den << endl;
    t.test_rel(xr.n,den,1.0e-6,"den2");
    
    xr.mu=mu;
    rf.calc_mu(xr,temper);
    d1=xr.n;
    xr.mu=-mu;
    rf.calc_mu(xr,temper);
    d2=xr.n;
    cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
    t.test_rel(den,d1-d2,1.0e-6,"den3");
    t.test_rel(t3,d1,1.0e-6,"den4");
    t.test_rel(t4,d2,1.0e-6,"den5");
    cout << endl;

  }

  // -----------------------------------------------------------------

  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate()." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
  
  double v1=rf.calibrate(e3,1);
  t.test_rel(v1,0.0,1.0e-6,"calibrate");

  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate() with better limits." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  // These seem to improve the accuracy. It's not clear that even
  // smaller tolerances will improve results.
  rf.upper_limit_fac=40.0;
  rf.dit->tol_abs=1.0e-13;
  rf.dit->tol_rel=1.0e-13;
  rf.nit->tol_abs=1.0e-13;
  rf.nit->tol_rel=1.0e-13;
  rf.density_root->tol_rel=1.0e-10;

  double v2=rf.calibrate(e3,1);
  t.test_rel(v2,0.0,1.0e-7,"calibrate 2");
  cout << v1 << " " << v2 << endl;

  // -----------------------------------------------------------------
  // Downcast the shared_ptr to the default integration type 

  inte_qag_gsl<> *qag=dynamic_cast<inte_qag_gsl<> *>(rf.dit.get());
  cout << qag->type() << " " << qag->get_rule() << endl;
  
  // -----------------------------------------------------------------
  
  t.report();

  return 0;
}

