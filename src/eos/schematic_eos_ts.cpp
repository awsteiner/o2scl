/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/schematic_eos.h>
#include <o2scl/skyrme_eos.h>
#include <o2scl/cold_nstar.h>
#include <o2scl/apr_eos.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  schematic_eos se, *se2;
  double dtemp, n0, nb;
  test_mgr t;
  t.set_output_level(1);

  thermo th;
  fermion n(939.0/hc_mev_fm,2.0), p(939.0/hc_mev_fm,2.0);

  n.n=0.159;
  p.n=0.001;
  se.calc_e(n,p,th);
  cout << (th.ed/0.16-n.m)*hc_mev_fm << endl;

  n.n=0.16;
  p.n=0.00;
  se.calc_e(n,p,th);
  cout << (th.ed/0.16-n.m)*hc_mev_fm << endl;

  n.n=0.08;
  p.n=0.08;
  se.calc_e(n,p,th);
  cout << (th.ed/0.16-n.m)*hc_mev_fm << endl;

  se.eoa=-16.0/hc_mev_fm;
  se.comp=200.0/hc_mev_fm;
  se.a=18.0/hc_mev_fm;
  se.b=12.0/hc_mev_fm;
  se.kprime=-2000.0/hc_mev_fm;
  se.n0=0.16;

  cout.setf(ios::scientific);

  n0=se.fn0(0.0,dtemp);
  t.test_rel(n0,0.16,1.0e-8,"n0");
  t.test_rel(se.feoa(n0)*hc_mev_fm,-16.0,1.0e-8,"eoa");
  t.test_rel(se.fcomp(n0)*hc_mev_fm,200.0,2.0e-8,"comp");
  t.test_rel(se.fesym(n0)*hc_mev_fm,30.0,4.0e-8,"esym");
  t.test_rel(se.fkprime(n0)*hc_mev_fm,-2.0e3,1.0e-6,"kprime");

  n.n=n0/2.0;
  p.n=n0/2.0;
  se.calc_e(n,p,th);
  t.test_rel((th.ed/n0-(n.m+p.m)/2.0)*hc_mev_fm,
	     se.eoa*hc_mev_fm,1.0e-8,"eoa");
  t.test_rel(n.mu*hc_mev_fm,p.mu*hc_mev_fm,1.0e-8,"mus");

  // Double check to ensure that the chemical potentials are correct

  double h=1.0e-5;

  {
    n.n=0.1;
    p.n=0.12;
    se.calc_e(n,p,th);
    double der1=th.ed;
    double mun=n.mu;
    double mup=p.mu;
    n.n+=h;
    se.calc_e(n,p,th);
    double der2=th.ed;
    n.n-=h;
    p.n+=h;
    se.calc_e(n,p,th);
    double der3=th.ed;
    p.n-=h;
    t.test_rel(mun,(der2-der1)/h,1.0e-5,"mun");
    t.test_rel(mup,(der3-der1)/h,1.0e-5,"mup");
  }

  se.gamma=1.5;
  {
    n.n=0.1;
    p.n=0.12;
    se.calc_e(n,p,th);
    double der1=th.ed;
    double mun=n.mu;
    double mup=p.mu;
    n.n+=h;
    se.calc_e(n,p,th);
    double der2=th.ed;
    n.n-=h;
    p.n+=h;
    se.calc_e(n,p,th);
    double der3=th.ed;
    p.n-=h;
    t.test_rel(mun,(der2-der1)/h,1.0e-5,"mun");
    t.test_rel(mup,(der3-der1)/h,1.0e-5,"mup");
  }

  // See how calc_p() works with an object of type schematic_eos
  n.n=0.16;
  p.n=0.16;
  for(double mu=4.52;mu<=4.68;mu+=0.01) {
    n.mu=mu;
    p.mu=mu;
    se.calc_p(n,p,th);
    se.calc_e(n,p,th);
    t.test_rel(n.mu,mu,1.0e-6,"neutron mu");
    t.test_rel(p.mu,mu,1.0e-6,"proton mu");
  }

  se.set_a_from_mstar(0.5,939.0/197.33);
  cout << se.a*hc_mev_fm << endl;
  se.set_a_from_mstar(1.0,939.0/197.33);
  cout << se.a*hc_mev_fm << endl;

  t.report();
  return 0;
}

