/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#include <o2scl/classical.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  part n(5.0,2.0);
  classical cl;
  t.test_rel(n.m,5.0,1.0e-6,"mass inheritance");
  double temper=0.1;
  n.n=0.1;

  cl.calc_density(n,temper);
  t.test_rel(n.pr,temper*n.n,1.0e-8,"ideal gas");
  t.test_rel(n.ed,n.m*n.n+1.5*temper*n.n,1.0e-8,"energy equipartition");
  t.test_rel(n.pr+n.ed-temper*n.en-n.n*n.mu,0.0,1.0e-8,
	     "thermodynamic identity");
  cl.calc_mu(n,temper);
  t.test_rel(n.n,0.1,1.0e-8,"calc_mu(calc_density)");

  t.report();
  return 0;
}
