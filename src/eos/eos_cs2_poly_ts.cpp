/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2019, Andrew W. Steiner
  
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
#include <iostream>

#include <o2scl/test_mgr.h>
#include <o2scl/eos_cs2_poly.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  
  eos_cs2_poly ecp;

  ecp.fix_params(1.2,1.0,2.0,0.5,1.0,0.1);
  ecp.fix_integ_consts(1.2,5.0,1.2,6.0);

  double eps=1.0e-5;
  for(double nb=1.2;nb<2.01;nb+=0.1) {
    double mu0=ecp.mu_from_nb(nb);
    double mu1=ecp.mu_from_nb(nb+eps);
    double dmudnb=(mu1-mu0)/eps;
    double ed0=ecp.ed_from_nb(nb);
    double ed1=ecp.ed_from_nb(nb+eps);
    double mu_test=(ed1-ed0)/eps;
    cout.precision(3);
    cout << nb << " ";
    cout.precision(6);
    cout << ecp.cs2_from_nb(nb) << " "
	 << nb/ecp.mu_from_nb(nb)*dmudnb << " "
	 << ecp.mu_from_nb(nb) << " " << mu_test << " " 
	 << ecp.ed_from_nb(nb) << endl;
    t.test_rel(nb/ecp.mu_from_nb(nb)*dmudnb,
	       ecp.cs2_from_nb(nb),eps*4.0,"cs2 vs. numerical.");
    t.test_rel(ecp.mu_from_nb(nb),
	       mu_test,eps*4.0,"mu vs. numerical.");
  }

  t.report();
  return 0;
}
