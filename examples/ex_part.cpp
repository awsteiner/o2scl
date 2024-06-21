/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
// sphinx-example-start
/* Example: ex_part.cpp
   -------------------------------------------------------------------
   See "License Information" section of the documentation for license
   information.
*/

#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/classical.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  // We work in units of inverse Fermis, so that energy density is
  // fm^{-4}. We also use a classical particle, to compare to the
  // nondegenerate approximation.
  fermion_rel relf;
  classical_thermo cla;
  
  fermion e2(o2scl_settings.get_convert_units().convert
             ("kg","1/fm",o2scl_const::mass_electron_f<double>()),2.0);
  fermion e3(e2.m,2.0);
  
  // We provide an initial guess to the chemical potential. This
  // is not a great guess for nondegenerate matter, but O2scl
  // aims to be successful even with bad guesses.
  e2.mu=e2.m;
  e3.mu=e3.m;

  // Compute the pressure at a density of 0.0001 fm^{-3} and a
  // temperature of 10 MeV. At these temperatures, the electrons are
  // non-degenerate, and Boltzmann statistics nearly applies.
  e2.n=0.0001;
  relf.calc_density(e2,10.0/hc_mev_fm);
  e3.n=0.0001;
  cla.calc_density(e3,10.0/hc_mev_fm);

  cout << e2.pr << " " << e3.pr << endl;

  // Test
  t.test_rel(e2.pr,e3.pr,4.0e-1,"classical vs. exact");

  // Compute the pressure at a density of 0.1 fm^{-3} and a
  // temperature of 1 MeV. At these temperatures, the electrons are
  // strongly degenerate
  e2.n=0.0001;
  relf.calc_density(e2,10.0/hc_mev_fm);
  cout << e2.pr << endl;

  // Now add the contribution to the pressure from positrons using the
  // implementation of part::pair_density()
  e2.n=0.0001;
  relf.pair_density(e2,10.0/hc_mev_fm);
  cout << e2.pr << endl;

  t.report();
  return 0;
}
// End of example
