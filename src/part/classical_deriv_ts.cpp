/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

#include <o2scl/classical_deriv.h>
#include <o2scl/classical.h>
#include <o2scl/fermion_nonrel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);
  
  classical_deriv snc;
  classical cla;

  part_deriv sf(5.0,2.0);
  part ef(5.0,2.0);

  t.test_rel(sf.m,5.0,1.0e-5,"mass_inheritance");
  t.test_rel(ef.m,5.0,1.0e-5,"mass_inheritance");
  
  sf.inc_rest_mass=false;
  ef.inc_rest_mass=true;

  double T=1.0;
  sf.non_interacting=false;
  sf.m=2.5;
  sf.nu=1.99*T-sf.m+sf.ms;
  ef.non_interacting=false;
  ef.m=2.5;
  ef.nu=1.99*T+sf.ms;

  cout << ef.inc_rest_mass << " " << sf.inc_rest_mass << endl;
  
  cout << "\"Non-degenerate\": " << endl;
  snc.calc_mu(sf,T);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << " " << sf.nu << endl;
  cout << sf.ed+sf.pr-sf.n*sf.nu << endl;
  cla.calc_mu(ef,T);
  cout << "ef: " << ef.n << " " << ef.ed-ef.n*ef.m << " " << ef.en << " " 
       << ef.pr << " " << ef.nu-ef.m << endl;
  cout << ef.ed+ef.pr-ef.n*ef.nu << endl;
  t.test_rel(sf.n,ef.n,1.0e-4,"ndeg density");
  t.test_rel(sf.ed,ef.ed-ef.n*ef.m,1.0e-4,"ndeg energy density");
  t.test_rel(sf.pr,ef.pr,1.0e-4,"ndeg pressure");
  t.test_rel(sf.en,ef.en,1.0e-4,"ndeg entropy");
  cout << endl;
  
  cout << "\"Degenerate\": " << endl;
  sf.nu=2.01*T-sf.m+sf.ms;
  ef.nu=2.01*T+sf.ms;
  snc.calc_mu(sf,T);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << endl;
  cla.calc_mu(ef,T);
  cout << "ef: " << ef.n << " " << ef.ed-ef.n*ef.m << " " << ef.en << " " 
       << ef.pr << endl;
  t.test_rel(sf.n,ef.n,1.0e-4,"deg density");
  t.test_rel(sf.ed,ef.ed-ef.n*ef.m,1.0e-4,"deg energy density");
  t.test_rel(sf.pr,ef.pr,1.0e-4,"deg pressure");
  t.test_rel(sf.en,ef.en,1.0e-4,"deg entropy");
  cout << endl;

  cout << "Test derivatives (\"non-degenerate\", direct): " << endl;
  double d1, d2, eps=1.0e-4;
  double dndmu, dndT, dsdT, dndm;
  
  sf.nu=1.0*T-sf.m+sf.ms+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.nu=1.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndmu=(d1-d2)/eps;

  sf.nu=1.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T+eps);
  d1=sf.n;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndT=(d1-d2)/eps;

  snc.calc_mu(sf,T+eps);
  d1=sf.en;
  snc.calc_mu(sf,T);
  d2=sf.en;
  dsdT=(d1-d2)/eps;

  sf.ms=5.0+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.ms=5.0;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndm=(d1-d2)/eps;

  snc.calc_mu(sf,T);
  cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
  cout << "nu: " << dndmu << " " << dndT << " " << dsdT << " "
       << dndm << endl;
  t.test_rel(dndmu,sf.dndmu,5.0e-4,"ndeg dir dndmu");
  t.test_rel(dndT,sf.dndT,5.0e-4,"ndeg dir dndT");
  t.test_rel(dsdT,sf.dsdT,5.0e-4,"ndeg dir dsdT");
  cout << endl;

  cout << "Test derivatives (\"degenerate\", direct): " << endl;

  sf.nu=4.0*T-sf.m+sf.ms+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.nu=4.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndmu=(d1-d2)/eps;

  sf.nu=4.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T+eps);
  d1=sf.n;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndT=(d1-d2)/eps;

  snc.calc_mu(sf,T+eps);
  d1=sf.en;
  snc.calc_mu(sf,T);
  d2=sf.en;
  dsdT=(d1-d2)/eps;
  
  sf.ms=5.0+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.ms=5.0;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndm=(d1-d2)/eps;

  snc.calc_mu(sf,T);
  cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
  cout << "nu: " << dndmu << " " << dndT << " " << dsdT << " "
       << dndm << endl;
  t.test_rel(dndmu,sf.dndmu,5.0e-4,"deg dir dndmu");
  t.test_rel(dndT,sf.dndT,5.0e-4,"deg dir dndT");
  t.test_rel(dsdT,sf.dsdT,5.0e-4,"deg dir dsdT");
  cout << endl;

  cout << "Test derivatives (\"non-degenerate\", byparts): " << endl;

  sf.nu=1.0*T-sf.m+sf.ms+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.nu=1.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndmu=(d1-d2)/eps;

  sf.nu=1.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T+eps);
  d1=sf.n;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndT=(d1-d2)/eps;

  snc.calc_mu(sf,T+eps);
  d1=sf.en;
  snc.calc_mu(sf,T);
  d2=sf.en;
  dsdT=(d1-d2)/eps;
  
  sf.ms=5.0+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.ms=5.0;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndm=(d1-d2)/eps;

  snc.calc_mu(sf,T);
  cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
  cout << "nu: " << dndmu << " " << dndT << " " << dsdT << " "
       << dndm << endl;
  t.test_rel(dndmu,sf.dndmu,5.0e-4,"ndeg byp dndmu");
  t.test_rel(dndT,sf.dndT,5.0e-4,"ndeg byp dndT");
  t.test_rel(dsdT,sf.dsdT,5.0e-4,"ndeg byp dsdT");
  cout << endl;

  cout << "Test derivatives (\"degenerate\", byparts): " << endl;

  sf.nu=4.0*T-sf.m+sf.ms+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.nu=4.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndmu=(d1-d2)/eps;

  sf.nu=4.0*T-sf.m+sf.ms;
  snc.calc_mu(sf,T+eps);
  d1=sf.n;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndT=(d1-d2)/eps;

  snc.calc_mu(sf,T+eps);
  d1=sf.en;
  snc.calc_mu(sf,T);
  d2=sf.en;
  dsdT=(d1-d2)/eps;
  
  sf.ms=5.0+eps;
  snc.calc_mu(sf,T);
  d1=sf.n;
  sf.ms=5.0;
  snc.calc_mu(sf,T);
  d2=sf.n;
  dndm=(d1-d2)/eps;

  snc.calc_mu(sf,T);
  cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
  cout << "nu: " << dndmu << " " << dndT << " " << dsdT << " "
       << dndm << endl;
  t.test_rel(dndmu,sf.dndmu,5.0e-4,"deg byp dndmu");
  t.test_rel(dndT,sf.dndT,5.0e-4,"deg byp dndT");
  t.test_rel(dsdT,sf.dsdT,5.0e-4,"deg byp dsdT");
  cout << endl;

  cout << "Check calc_density()." << endl;
  snc.calc_mu(sf,1.5);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << " " << sf.nu << endl;
  snc.calc_density(sf,1.5);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << " " << sf.nu << endl;
  snc.calc_mu(sf,0.7);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << " " << sf.nu << endl;
  snc.calc_density(sf,0.7);
  cout << "sf: " << sf.n << " " << sf.ed << " " << sf.en << " " 
       << sf.pr << " " << sf.nu << endl;

  t.report();
  return 0;
}
