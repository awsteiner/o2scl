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
#include <o2scl/test_mgr.h>
#include <o2scl/fermion.h>
#include <o2scl/fermion_eff.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  
  fermion e(1.0,2.0);
  fermion_eff fet;
  e.non_interacting=true;

  double alpha, two13, alpha16, cbt, alpha2, temper;
  
  // This section is discussed in the documentation for
  // the massless_pair_density() function
  cout << "Testing the limits of the expansions in"
       << " massless_pair_density()." << endl;
  temper=0.5;
  cout << "alpha        Large alpha  Exact        Small alpha" << endl;
  t.set_output_level(1);

  for(alpha=1.0e10;alpha>=9.9e-11;alpha/=2.0) {

    // Use the large alpha expansion to determine the chemical potential
    e.mu=(2.0/3.0/sqrt(alpha)-8.0/81.0/pow(alpha,1.5)+
	  32.0/729.0/pow(alpha,2.5))*pi*temper/sqrt(3.0);
    // Then compute the density 
    fet.massless_pair_mu(e,temper);
    // And from the density the new value of alpha
    alpha2=4.0*pi2*pow(temper,6.0)/243.0/e.n/e.n;
    cout << alpha << " ";
    cout << fabs((alpha-alpha2)/alpha) << " ";
    if (alpha>1.0e4) t.test_rel(fabs((alpha-alpha2)/alpha),0.0,
				5.0e-13,"large alpha");

    // Use the exact expression for the chemical potential
    cbt=pow(-1.0+sqrt(1.0+alpha),1.0/3.0)/pow(alpha,1.0/6.0);
    e.mu=pi*temper/sqrt(3.0)*(1.0/cbt-cbt);
    // Then compute the density 
    fet.massless_pair_mu(e,temper);
    // And from the density the new value of alpha
    alpha2=4.0*pi2*pow(temper,6.0)/243.0/e.n/e.n;
    cout << fabs((alpha-alpha2)/alpha) << " ";
    if (alpha>=3.0e-4 && alpha<=1.0e4) 
      t.test_rel(fabs((alpha-alpha2)/alpha),0.0,
		 5.0e-13,"moderate alpha");

    // Use the small alpha expansion to determine the chemical potential
    two13=cbrt(2.0);
    alpha16=pow(alpha,1.0/6.0);
    e.mu=(two13/alpha16-alpha16/two13+alpha/alpha16/6.0/two13/two13
	  +alpha*alpha16/12.0/two13-alpha*alpha/alpha16/18.0/two13/two13-
	  5.0*alpha*alpha*alpha16/144.0/two13+
	  77.0/2592.0*alpha*alpha*alpha/alpha16/two13/two13)*
      pi*temper/sqrt(3.0);
    // Then compute the density 
    fet.massless_pair_mu(e,temper);
    // And from the density the new value of alpha
    alpha2=4.0*pi2*pow(temper,6.0)/243.0/e.n/e.n;
    cout << fabs((alpha-alpha2)/alpha) << endl;
    if (alpha<3.0e-4) t.test_rel(fabs((alpha-alpha2)/alpha),0.0,
				 5.0e-13,"small alpha");
  }
  t.set_output_level(2);
  
  // Testing ndnr expansion

  /*
    e.m=5.0;
    e.mu=4.95;
    e.exp_precision=1.0e-8;
    ret=e.ndnr_expansion(1.0e-2);
    t.test_gen(ret==0,"ndnr_expansion(0)");
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ret=e.calc_mu(1.0e-2);
    cout << ret << " " << err_hnd->get_str() << endl;
    t.test_gen(ret==0,"calc_mu(1)");
    t.test_rel(e.n,tmp1,2.0e-3,"ndnr_expansion(1)");
    t.test_rel(e.mu,tmp2,1.0e-3,"ndnr_expansion(2)");
    t.test_rel(e.ed,tmp3,2.0e-3,"ndnr_expansion(3)");
    t.test_rel(e.pr,tmp4,2.0e-3,"ndnr_expansion(4)");
    t.test_rel(e.en,tmp5,2.0e-3,"ndnr_expansion(5)");
  */

  t.report();

  return 0;
}
