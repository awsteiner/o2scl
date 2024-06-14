/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/test_mgr.h>
#include <o2scl/fermion.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/fermion_rel.h>

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

  e.mu=2.0;
  e.n=0.01;
  fet.massless_pair_density(e,1.5);
  
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
  cout << endl;

  // -----------------------------------------------------------------
  // Test the ndeg_terms function in the non-degenerate limit 
  
  if (true) {
    
    cout << "Testing ndeg_terms():" << endl;
    cout.setf(ios::showpos);

    e.init(0.511/197.33,2.0);
    fermion_rel ft;
    double mux=1.0e-7;

    size_t j=1;
    double pterm, nterm, enterm, edterm;
    double pterm1, nterm1, enterm1, edterm1;
    double pterm2, nterm2, enterm2, edterm2;
    double T=20.0/197.33;

    double tt=T/e.m;
    double psi;
    bool inc_antip;

    inc_antip=false;

    e.mu=-mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm1,nterm1,enterm1,edterm1);
    cout << pterm1 << " " << nterm1 << " " << enterm1 << " "
         << edterm1 << endl;
    e.mu=mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm2,nterm2,enterm2,edterm2);
    cout << pterm2 << " " << nterm2 << " " << enterm2 << " "
         << edterm2 << endl;
    cout << pterm1+pterm2 << " " << nterm1-nterm2 << " "
         << enterm1+enterm2 << " " << edterm1+edterm2 << endl;

    inc_antip=true;

    e.mu=-mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm,nterm,enterm,edterm);
    cout << pterm << " " << nterm << " " << enterm << " "
         << edterm << endl;

    t.test_rel(pterm1+pterm2,pterm,1.0e-12,"pterm j=1");
    t.test_rel(nterm1-nterm2,nterm,1.0e-9,"nterm j=1");
    t.test_rel(enterm1+enterm2,enterm,1.0e-10,"enterm j=1");
    t.test_rel(enterm1+enterm2,enterm,1.0e-10,"enterm j=1");
    
    e.mu=mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm,nterm,enterm,edterm);
    cout << pterm << " " << nterm << " " << enterm << " "
         << edterm << endl;

    t.test_rel(pterm1+pterm2,pterm,1.0e-12,"pterm antip. j=1");
    t.test_rel(nterm2-nterm1,nterm,1.0e-9,"nterm antip. j=1");
    t.test_rel(enterm1+enterm2,enterm,1.0e-10,"enterm antip. j=1");
    t.test_rel(enterm1+enterm2,enterm,1.0e-10,"enterm antip. j=1");
    
    cout << endl;

    j=2;

    inc_antip=false;

    e.mu=-mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm1,nterm1,enterm1,edterm1);
    cout << pterm1 << " " << nterm1 << " " << enterm1 << " "
         << edterm1 << endl;
    e.mu=mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm2,nterm2,enterm2,edterm2);
    cout << pterm2 << " " << nterm2 << " " << enterm2 << " "
         << edterm2 << endl;
    cout << pterm1+pterm2 << " " << nterm1-nterm2 << " "
         << enterm1+enterm2 << " " << edterm1+edterm2 << endl;

    inc_antip=true;

    e.mu=-mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm,nterm,enterm,edterm);
    cout << pterm << " " << nterm << " " << enterm << " "
         << edterm << endl;

    t.test_rel(pterm1+pterm2,pterm,1.0e-11,"pterm j=2");
    t.test_rel(nterm1-nterm2,nterm,1.0e-9,"nterm j=2");
    t.test_rel(enterm1+enterm2,enterm,1.0e-9,"enterm j=2");
    t.test_rel(enterm1+enterm2,enterm,1.0e-9,"enterm j=2");
    
    e.mu=mux;
    psi=(e.mu-e.m)/T;
    ft.ndeg_terms(j,tt,psi*tt,e.m,e.inc_rest_mass,inc_antip,
               pterm,nterm,enterm,edterm);
    cout << pterm << " " << nterm << " " << enterm << " "
         << edterm << endl;

    t.test_rel(pterm1+pterm2,pterm,1.0e-11,"pterm antip. j=2");
    t.test_rel(nterm2-nterm1,nterm,1.0e-9,"nterm antip. j=2");
    t.test_rel(enterm1+enterm2,enterm,1.0e-9,"enterm antip. j=2");
    t.test_rel(enterm1+enterm2,enterm,1.0e-9,"enterm antip. j=2");
    
    cout << endl;

  }

  if (true) {
    cout << "An alternative way of testing ndeg_terms(), by\n  "
         << "comparing with numerical values in fermion.ipynb." << endl;
      
    double tt, psi;
    e.ms=0.511/197.33;
    e.inc_rest_mass=true;
    tt=0.1/0.511;
    psi=(8.0e-12-e.ms)/0.1*197.33;
    for(size_t j=1;j<6;j++) {
      double pterm, nterm, enterm, edterm;
      fet.ndeg_terms(j,tt,psi*tt,e.ms,e.inc_rest_mass,false,
                 pterm,nterm,enterm,edterm);
      std::cout.setf(std::ios::showpos);
      std::cout.precision(8);
      std::cout << j << " " << pterm << " " << nterm << " "
                << enterm << " " << edterm << std::endl;
    }
    std::cout << endl;
    
    for(size_t j=1;j<6;j++) {
      double pterm, nterm, enterm, edterm;
      fet.ndeg_terms(j,tt,psi*tt,e.ms,e.inc_rest_mass,true,
                     pterm,nterm,enterm,edterm);
      std::cout.setf(std::ios::showpos);
      std::cout.precision(8);
      std::cout << j << " " << pterm << " " << nterm << " "
                << enterm << " " << edterm << std::endl;
    }
  }

  if (false) {
    // Test pair_den_ndeg()
    e.mu=1.0e-14;
    fermion_rel ft;
    double Tx=1.0;
    ft.pair_mu(e,Tx);
    cout << "e.n: " << e.n << endl;
    ft.pair_den_ndeg(e,Tx);
    ft.calc_mu_ndeg(e,Tx,1.0e-18,true,1);
  }
  
  t.report();

  return 0;
}
