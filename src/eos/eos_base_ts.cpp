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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_base.h>
#include <o2scl/test_mgr.h>
#include <o2scl/deriv_gsl.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  thermo th;

  eos_base eo;
  
  th.ed=1.0;
  th.pr=2.0;
  th.en=3.0;

  eo.set_thermo(th);

  eos_leptons elep;
  elep.include_muons=false;
  elep.include_deriv=true;

  // First, do a test with just pair_mu():

  cout << "pair_mu(): " << endl;
  
  elep.e.mu=1.0;
  
  elep.default_acc();
  elep.pair_mu(0.1);
  cout << "n: " << dtos(elep.e.n,0) << " en: ";
  cout << dtos(elep.e.en,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;
  cout << "  dndT: " << dtos(elep.ed.dndT) << " dsdT: "
       << dtos(elep.ed.dsdT) << endl;

  // Check temperature derivatives
  fermion et=elep.e;
  part_deriv_press edt=elep.ed;
  elep.pair_mu(0.1+1.0e-4);
  t.test_rel(elep.ed.dndT,(elep.e.n-et.n)/1.0e-4,1.0e-3,"dndT 1");
  t.test_rel(elep.ed.dsdT,(elep.e.en-et.en)/1.0e-4,1.0e-3,"dsdT 1");
  
  elep.improved_acc();
  elep.pair_mu(0.1);
  cout << "n: " << dtos(elep.e.n,0) << " en: ";
  cout << dtos(elep.e.en,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;
  cout << "  dndT: " << dtos(elep.ed.dndT) << " dsdT: "
       << dtos(elep.ed.dsdT) << endl;

  elep.ld_acc();
  elep.pair_mu(0.1);
  cout << "n: " << dtos(elep.e.n,0) << " en: ";
  cout << dtos(elep.e.en,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;
  cout << "  dndT: " << dtos(elep.ed.dndT) << " dsdT: "
       << dtos(elep.ed.dsdT) << endl;

  elep.fp_25_acc();
  elep.pair_mu(0.1);
  cout << "n: " << dtos(elep.e.n,0) << " en: ";
  cout << dtos(elep.e.en,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;
  cout << "  dndT: " << dtos(elep.ed.dndT) << " dsdT: "
       << dtos(elep.ed.dsdT) << endl;

  cout << endl;

  // Now, with pair_density()

  cout << "pair_density(): " << endl;
  
  elep.e.n=1.0e-6;
  elep.pair_density(0.1);
  elep.pair_mu(0.1);
  cout << "mu: " << dtos(elep.e.mu,0) << " n: ";
  cout << dtos(elep.e.n,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;
  
  elep.improved_acc();
  elep.e.n=1.0e-6;
  elep.pair_density(0.1);
  elep.pair_mu(0.1);
  cout << "mu: " << dtos(elep.e.mu,0) << " n: ";
  cout << dtos(elep.e.n,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;

  elep.ld_acc();
  elep.e.n=1.0e-6;
  elep.pair_density(0.1);
  elep.pair_mu(0.1);
  cout << "mu: " << dtos(elep.e.mu,0) << " n: ";
  cout << dtos(elep.e.n,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;

  elep.fp_25_acc();
  elep.e.n=1.0e-6;
  elep.pair_density(0.1);
  elep.pair_mu(0.1);
  cout << "mu: " << dtos(elep.e.mu,0) << " n: ";
  cout << dtos(elep.e.n,0) << " dndmu: ";
  cout << dtos(elep.ed.dndmu,0) << endl;

  cout << endl;
  
  t.report();
  return 0;
}

