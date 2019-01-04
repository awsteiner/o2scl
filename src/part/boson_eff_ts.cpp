/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#include <o2scl/boson_eff.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/boson_rel.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  boson_eff eb;
  boson b(1.0,2.0);
  t.test_rel(b.m,1.0,1.0e-6,"mass inheritance");
  b.non_interacting=true;

  double T;

  cout << "----------------------------------------------------" << endl;
  cout << "boson_eff tests:" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  cout << "massless, non-degenerate (mu=m=0) "
       << "(m=0.011, mu=0.01, T=0.1)" << endl;
  b.m=0.011;
  b.mu=0.01;
  T=0.1;
  cout << "calc_mu(T) vs. calc_density(T) vs. massless_calc(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl; 
  double tn=b.n, ted=b.ed, tpr=b.pr, ten=b.en;
  eb.calc_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  t.test_rel(b.n,tn,3.0e-9,"density 1");
  t.test_rel(b.ed,ted,3.0e-9,"energy 1");
  t.test_rel(b.pr,tpr,3.0e-9,"pressure 1");
  t.test_rel(b.en,ten,3.0e-9,"entropy 1");
  b.m=0.011;
  b.mu=0.01;
  T=0.1;
  b.massless_calc(T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  t.test_rel(b.n,tn,2.0e0,"density 2");
  t.test_rel(b.ed,ted,2.0e-1,"energy 2");
  t.test_rel(b.pr,tpr,2.0e-1,"pressure 2");
  t.test_rel(b.en,ten,2.0e-1,"entropy 2");
  cout << endl;

  cout << "degenerate relativistic (mu=m, mu>>T) "
       << "(m=0.11, mu=0.1, T=0.01)" << endl;
  b.m=0.11;
  b.mu=0.1;
  T=0.01;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  tn=b.n, ted=b.ed, tpr=b.pr, ten=b.en;
  eb.calc_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << endl;
  t.test_rel(b.n,tn,1.0e-10,"density 3");
  t.test_rel(b.ed,ted,1.0e-10,"energy 3");
  t.test_rel(b.pr,tpr,1.0e-10,"pressure 3");
  t.test_rel(b.en,ten,1.0e-10,"entropy 3");

  cout << "Testing mu=m limit "
       << "(m=0.101, mu=0.1, T=0.01)" << endl;
  b.m=0.101;
  b.mu=0.1;
  T=0.01;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << "(m=0.1, mu=0.1, T=0.01)" << endl;
  b.m=0.1;
  b.mu=0.1;
  T=0.01;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << endl;

  cout << "degenerate non-relativistic (m>>mu>>T) "
       << "(m=0.5, mu=0.1, T=0.05)" << endl;
  b.m=0.5;
  b.mu=0.1;
  T=0.05;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  eb.calc_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << endl;

  cout << "non-degenerate relativistic (mu=m, mu<<T) "
       << "(m=0.011, mu=0.01, T=0.1)" << endl;
  b.m=0.011;
  b.mu=0.01;
  T=0.1;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  eb.calc_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << "pair_mu(T) vs. pair_density(T)" << endl;
  b.m=0.011;
  b.mu=0.01;
  T=0.1;
  eb.pair_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  tn=b.n, ted=b.ed, tpr=b.pr, ten=b.en;
  eb.pair_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << endl;
  //  t.test_rel(b.n,tn,2.0e-6,"density");
  //  t.test_rel(b.ed,ted,2.0e-6,"energy");
  //  t.test_rel(b.pr,tpr,2.0e-6,"pressure");
  //  t.test_rel(b.en,ten,2.0e-6,"entropy");

  cout << "non-degenerate non-relativistic (mu<<m, mu<<T) "
       << "(m=0.1, mu=0.01, T=1.0)" << endl;
  b.m=0.1;
  b.mu=0.01;
  T=1.0;
  cout << "calc_mu(T) vs. calc_density(T)" << endl;
  eb.calc_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl; 
  tn=b.n, ted=b.ed, tpr=b.pr, ten=b.en;
  eb.calc_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  t.test_rel(b.n,tn,2.0e-6,"density");
  t.test_rel(b.ed,ted,2.0e-6,"energy");
  t.test_rel(b.pr,tpr,2.0e-6,"pressure");
  t.test_rel(b.en,ten,2.0e-6,"entropy");
  cout << "pair_mu(T) vs. pair_density(T)" << endl;
  b.m=0.1;
  b.mu=0.01;
  T=1.0;
  eb.pair_mu(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  tn=b.n, ted=b.ed, tpr=b.pr, ten=b.en;
  eb.pair_density(b,T);
  cout << "mu,ed,pr,en: " << b.mu << " " << b.ed << " " 
       << b.pr << " " << b.en << endl;
  cout << endl;
  //  t.test_rel(b.n,tn,2.0e-6,"density");
  //  t.test_rel(b.ed,ted,2.0e-6,"energy");
  //  t.test_rel(b.pr,tpr,2.0e-6,"pressure");
  //  t.test_rel(b.en,ten,2.0e-6,"entropy");

  t.report();
  return 0;
}

