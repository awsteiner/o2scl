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
  eos_base eo, eo2, eo3;
  eos_base *eo4, *eo5, *eo6;
  test_mgr t;
  t.set_output_level(2);

  thermo th, *th2;
  int vp;
  string stype;

  cout.setf(ios::scientific);

  th.ed=1.0;
  th.pr=2.0;
  th.en=3.0;

  eo.set_thermo(th);

  t.report();
  return 0;
}

