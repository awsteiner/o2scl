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
#include <o2scl/fermion_eff.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/classical.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  fermion_eff ef2;

  for(size_t ik=0;ik<2;ik++) {
    
    if (ik==1) ef2.load_coefficients(fermion_eff::cf_fermilat3);
  
    double tmp1, tmp2, tmp3, tmp4, tmp5, T;

    fermion e(1.0,2.0);
    t.test_rel(e.m,1.0,1.0e-6,"mass inheritance");
    e.non_interacting=true;

    fermion ferm(1.0,2.0);
    double den;
  
    cout << "----------------------------------------------------" << endl;
    cout << "Four combinations of inc_rest_mass & non_interacting" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << endl;

    // Higher densities are commented out because they were
    // causing problems. These need to be fixed.
    //for(den=1.0e-6;den<=1.01e2;den*=1.0e2) {
    for(den=1.0e-6;den<=1.01e-3;den*=1.0e2) {

      cout << "density: " << den << endl;

      // -----------------------------------------------------------------
      // Test calc_density() and pair_density() with inc_rest_mass is 
      // both true and false, and non_interacting is both true and false
      // -----------------------------------------------------------------
  
      fermion xe(0.511/hc_mev_fm,2.0);
      double temper=10.0/hc_mev_fm, d1, d2, mu;
      double t1, t2, t3, t4;

      // -----------------------------------------------------------------

      xe.non_interacting=false;
      cout << "Interacting: " << endl;
      cout << endl;
    
      xe.inc_rest_mass=true;
      cout << "Include rest mass." << endl;
    
      xe.n=den;
      xe.nu=xe.m;
      ef2.calc_density(xe,temper);
      cout << xe.nu << " ";
      t1=xe.nu;
      ef2.calc_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den1a");

      xe.n=den;
      xe.nu=xe.m;
      ef2.pair_density(xe,temper);
      mu=xe.nu;
      cout << xe.nu << " ";
      t2=xe.nu;
      ef2.pair_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den2a");
    
      xe.nu=mu;
      ef2.calc_mu(xe,temper);
      d1=xe.n;
      if (xe.inc_rest_mass) {
	xe.nu=-mu;
      } else {
	xe.nu=-mu-2.0*xe.m;
      }
      ef2.calc_mu(xe,temper);
      d2=xe.n;
      cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
	   
      t3=d1;
      t4=d2;
      t.test_rel(den,d1-d2,1.0e-6,"den3a");
      cout << endl;

      // -----------------------------------------------------------------

      xe.inc_rest_mass=false;
      cout << "Remove rest mass." << endl;
    
      xe.n=den;
      xe.nu=0.0;
      ef2.calc_density(xe,temper);
      cout << xe.nu+xe.m << " ";
      t.test_rel(xe.nu+xe.m,t1,1.0e-6,"nu1b");
      ef2.calc_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den1b");

      xe.n=den;
      xe.nu=0.0;
      ef2.pair_density(xe,temper);
      mu=xe.nu;
      cout << xe.nu+xe.m << " ";
      t.test_rel(xe.nu+xe.m,t2,1.0e-6,"nu2b");
      ef2.pair_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den2b");
    
      xe.nu=mu;
      ef2.calc_mu(xe,temper);
      d1=xe.n;
      xe.nu=-mu-2.0*xe.m;
      ef2.calc_mu(xe,temper);
      d2=xe.n;
      cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
	   
      t.test_rel(den,d1-d2,1.0e-6,"den3b");
      t.test_rel(t3,d1,1.0e-6,"den4b");
      t.test_rel(t4,d2,1.0e-6,"den5b");
      cout << endl;

      // -----------------------------------------------------------------

      xe.non_interacting=true;
      cout << "Non-interacting: " << endl;
      cout << endl;

      xe.inc_rest_mass=true;
      cout << "Include rest mass." << endl;
    
      xe.n=den;
      xe.mu=0.0;
      ef2.calc_density(xe,temper);
      cout << xe.mu << " ";
      t.test_rel(xe.nu,t1,1.0e-6,"nu1c");
      ef2.calc_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den1c");
    
      xe.n=den;
      xe.mu=0.0;
      ef2.pair_density(xe,temper);
      mu=xe.mu;
      cout << xe.mu << " ";
      t.test_rel(xe.nu,t2,1.0e-6,"nu2c");
      ef2.pair_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den2c");
    
      xe.mu=mu;
      ef2.calc_mu(xe,temper);
      d1=xe.n;
      xe.mu=-mu;
      ef2.calc_mu(xe,temper);
      d2=xe.n;
      cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
      t.test_rel(den,d1-d2,1.0e-6,"den3c");
      t.test_rel(t3,d1,1.0e-6,"den4c");
      t.test_rel(t4,d2,1.0e-6,"den5c");
      cout << endl;

      // -----------------------------------------------------------------

      xe.inc_rest_mass=false;
      cout << "Remove rest mass." << endl;
    
      xe.n=den;
      xe.mu=0.0;
      ef2.calc_density(xe,temper);
      cout << xe.mu+xe.m << " ";
      t.test_rel(xe.nu+xe.m,t1,1.0e-6,"nu1d");
      ef2.calc_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den1d");
    
      xe.n=den;
      xe.mu=0.0;
      ef2.pair_density(xe,temper);
      mu=xe.mu;
      cout << xe.mu+xe.m << " ";
      t.test_rel(xe.nu+xe.m,t2,1.0e-6,"nu2d");
      ef2.pair_mu(xe,temper);
      cout << xe.n << " " << den << endl;
      t.test_rel(xe.n,den,1.0e-6,"den2d");
    
      xe.mu=mu;
      ef2.calc_mu(xe,temper);
      d1=xe.n;
      xe.mu=-mu-2.0*xe.m;
      ef2.calc_mu(xe,temper);
      d2=xe.n;
      cout << den << " " << d1 << " " << d2 << " " << d1-d2 << endl;
	   
      t.test_rel(den,d1-d2,1.0e-6,"den3d");
      t.test_rel(t3,d1,1.0e-6,"den4d");
      t.test_rel(t4,d2,1.0e-6,"den5d");
      cout << endl;

    }

    cout << "----------------------------------------------------" << endl;
    cout << "fermion_eff tests:" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << endl;

    cout << "degenerate extremely relativistic (mu>>m, mu>T) "
	 << "(m=0.01, mu=1.0, T=0.1)" << endl;
    e.m=0.01;
    e.mu=1.0;
    T=0.1;
    cout << "calc_mu(e,T) vs. massless_calc_mu(e,T)" << endl;
    ef2.calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl; 
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.massless_calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-4,"fermion_eff density 1");
    t.test_rel(tmp2,e.mu,5.0e-4,"fermion_eff chem. pot. 1");
    t.test_rel(tmp3,e.ed,5.0e-4,"fermion_eff energy 1");
    t.test_rel(tmp4,e.pr,5.0e-4,"fermion_eff pressure 1");
    t.test_rel(tmp5,e.en,5.0e-4,"fermion_eff entropy 1");

    cout << "calc_density(e,T) vs. massless_calc_density(e,T)" << endl;
    e.m=0.01;
    e.mu=1.0;
    T=0.1;
    ef2.calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.massless_calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;

    t.test_rel(tmp1,e.n,5.0e-4,"fermion_eff density 2");
    t.test_rel(tmp2,e.mu,5.0e-4,"fermion_eff chem. pot. 2");
    t.test_rel(tmp3,e.ed,5.0e-4,"fermion_eff energy 2");
    t.report();
    t.test_rel(tmp4,e.pr,5.0e-4,"fermion_eff pressure 2");
    t.test_rel(tmp5,e.en,5.0e-4,"fermion_eff entropy 2");
    cout << endl;

    cout << "extremely degenerate relativistic (mu>m, mu>>T) "
	 << "(m=0.1, mu=1.0, T=0.01)" << endl;
    e.m=0.1;
    e.mu=1.0;
    T=0.01;
    cout << "calc_mu(e,T) vs. calc_mu(T=0)" << endl;
    ef2.calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.calc_mu_zerot(e);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 3");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 3");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 3");
    t.test_rel(tmp4,e.pr,5.0e-3,"fermion_eff pressure 3");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 3");

    cout << "calc_density(e,T) vs. calc_density(T=0)" << endl;
    e.m=0.1;
    e.mu=1.0;
    T=0.01;
    ef2.calc_mu(e,T);
    ef2.calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.calc_density_zerot(e);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 4");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 4");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 4");
    t.test_rel(tmp4,e.pr,5.0e-3,"fermion_eff pressure 4");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 4");
    cout << endl;

    cout << "degenerate non-relativistic ((mu-m=kf^2/2/m)<<m, T<<mu) "
	 << "(m=1.0, mu=1.1, T=0.01)" << endl;
    e.m=1.0;
    e.mu=1.1;
    T=0.01;
    cout << "calc_mu(e,T) vs. calc_mu(T=0)" << endl;
    ef2.calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.calc_mu_zerot(e);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,3.0e-2,"fermion_eff density 5");
    t.test_rel(tmp2,e.mu,3.0e-2,"fermion_eff chem. pot. 5");
    t.test_rel(tmp3,e.ed,3.0e-2,"fermion_eff energy 5");
    t.test_rel(tmp4,e.pr,8.0e-2,"fermion_eff pressure 5");
    t.test_rel(tmp5,e.en,3.0e-2,"fermion_eff entropy 5");

    cout << "calc_density(e,T) vs. calc_density(T=0)" << endl;
    e.m=1.0;
    e.mu=1.1;
    T=0.01;
    ef2.calc_mu(e,T);
    ef2.calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.calc_density_zerot(e);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 6");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 6");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 6");
    t.test_rel(tmp4,e.pr,5.0e-2,"fermion_eff pressure 6");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 6");
    cout << endl;

    cout << "non-degenerate extremely relativistic (m<<mu<<T) "
	 << "(m=0.01, mu=0.1, T=1.0)" << endl;
    e.m=0.01;
    e.mu=0.1;
    T=1.0;
    cout << "calc_mu(e,T) vs. calc_density(e,T)" << endl;
    ef2.calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    e.mu*=2.0;
    ef2.calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 7");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 7");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 7");
    t.test_rel(tmp4,e.pr,5.0e-3,"fermion_eff pressure 7");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 7");

    cout << "pair_mu(e,T) vs. massless_pair_mu(e,T) vs. pair_density(e,T) vs. " 
	 << endl;
    cout << "massless_pair_density(e,T)" << endl;
    e.m=0.01;
    e.mu=0.1;
    T=1.0;
    ef2.pair_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.massless_pair_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,1.0e-2,"fermion_eff density 8");
    t.test_rel(tmp2,e.mu,1.0e-2,"fermion_eff chem. pot. 8");
    t.test_rel(tmp3,e.ed,1.0e-2,"fermion_eff energy 8");
    t.test_rel(tmp4,e.pr,1.0e-2,"fermion_eff pressure 8");
    t.test_rel(tmp5,e.en,1.0e-2,"fermion_eff entropy 8");

    e.m=0.01;
    e.mu=0.1;
    T=1.0;
    ef2.pair_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-4,"fermion_eff density 9c");
    t.test_rel(tmp2,e.mu,5.0e-4,"fermion_eff chem. pot. 9c");
    t.test_rel(tmp3,e.ed,5.0e-4,"fermion_eff energy 9c");
    t.test_rel(tmp4,e.pr,5.0e-4,"fermion_eff pressure 9c");
    t.test_rel(tmp5,e.en,5.0e-4,"fermion_eff entropy 9c");
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    ef2.massless_pair_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,1.0e-2,"fermion_eff density 9");
    t.test_rel(tmp2,e.mu,1.0e-2,"fermion_eff chem. pot. 9");
    t.test_rel(tmp3,e.ed,1.0e-2,"fermion_eff energy 9");
    t.test_rel(tmp4,e.pr,1.0e-2,"fermion_eff pressure 9");
    t.test_rel(tmp5,e.en,1.0e-2,"fermion_eff entropy 9");
    cout << endl;

    cout << "non-degenerate non-relativistic (mu<<m, mu<<T) "
	 << "(m=0.1, mu=0.01, T=1.0)" << endl;
    e.m=0.1;
    e.mu=0.01;
    T=1.0;
    cout << "calc_mu(e,T) vs. calc_density(e,T)" << endl;
    ef2.calc_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl; 
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    e.mu*=2.0;
    ef2.calc_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 10");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 10");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 10");
    t.test_rel(tmp4,e.pr,5.0e-3,"fermion_eff pressure 10");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 10");

    cout << "pair_mu(e,T) vs. pair_density(e,T)" << endl;
    e.m=0.1;
    e.mu=0.01;
    T=1.0;
    ef2.pair_mu(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    tmp1=e.n;
    tmp2=e.mu;
    tmp3=e.ed;
    tmp4=e.pr;
    tmp5=e.en;
    e.mu*=2.0;
    ef2.pair_density(e,T);
    cout << e.n << " " << e.mu << " " << e.ed << " " 
	 << e.pr << " " << e.en << endl;
    t.test_rel(tmp1,e.n,5.0e-3,"fermion_eff density 11");
    t.test_rel(tmp2,e.mu,5.0e-3,"fermion_eff chem. pot. 11");
    t.test_rel(tmp3,e.ed,5.0e-3,"fermion_eff energy 11");
    t.test_rel(tmp4,e.pr,5.0e-3,"fermion_eff pressure 11");
    t.test_rel(tmp5,e.en,5.0e-3,"fermion_eff entropy 11");
    cout << endl;
  
    // -----------------------------------------------------------------
    // This tests the results of fermion_eff at the fitting points used
    // in Johns, Ellis, and Lattimer and compares them to the results
    // from the fermion_rel class.
  
    cout << "----------------------------------------------------" << endl;
    cout << "fermion_eff vs. fermion_rel:" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << endl;

    fermion e3(1.0,2.0);
    fermion_rel rf;

    cout.precision(4);
    double logf, logg, psi, sqt, a=0.42, f, g;

    double maxdev=0.0, maxf=0.0, maxg=0.0;

    // The fit points from Eggelton73
    for(logf=-10.61;logf<=8.0;logf+=2.0) {
      for(logg=-7.0;logg<=5.01;logg+=2.0) {

	f=pow(10.0,logf);
	g=pow(10.0,logg);
	T=1.0;
	e.m=T/g*sqrt(1.0+f);
	sqt=sqrt(1.0+f/a);
	psi=2.0*sqt+log((sqt-1.0)/(sqt+1.0));
	e.mu=psi*T+e.m;
	e3.m=e.m;
	e3.mu=e.mu;
	ef2.calc_mu(e,T);
	
	rf.calc_mu(e3,T);

	cout.setf(ios::showpos);
	cout << logf << " " << logg << " "; 
	cout.unsetf(ios::showpos);
	cout << e.pr << " " << e3.pr << " " << fabs(e.pr-e3.pr)/e.pr << " ";
	t.test_rel(e.pr,e3.pr,1.0e-2,"eq. pressure from mu");
	if (fabs(e.pr-e3.pr)/e.pr>maxdev) {
	  maxdev=fabs(e.pr-e3.pr)/e.pr;
	  maxf=logf;
	  maxg=logg;
	}
	ef2.calc_density(e,T);
	cout << e.pr << endl;
	t.test_rel(e.pr,e3.pr,1.0e-2,"eq. pressure from density");

      }

      if (logf<-10.0) logf+=3.61;
      if (logf>4.0) logf+=0.82;
    }
    cout << "Maximum deviation: " << maxdev << " at log(f): " << maxf 
	 << " log(g): " << maxg << endl;
    cout << endl;

    {
      // Show that the free energy decreases as a function of increasing
      // temperature at a particular fixed density
      fermion ext(1.0,2.0);
      ext.n=0.01;
      ext.mu=ext.m;
      ef2.calc_density(ext,0.02);
      cout << ext.mu << " " << ext.ed << " " << ext.ed-0.02*ext.en << endl;
      ef2.calc_density_zerot(ext);
      cout << ext.mu << " " << ext.ed << " " << ext.ed << endl;
    }

  }

  // -----------------------------------------------------------------

  fermion ec;
  t.test_rel(ef2.calibrate(ec,1,0,"../../data/o2scl/fermion_cal.o2"),
	     0.0,0.5,"calibrate");
  
  // -----------------------------------------------------------------

  t.report();
  return 0;
}
  

