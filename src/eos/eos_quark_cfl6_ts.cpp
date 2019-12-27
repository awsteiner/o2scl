/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

#include <o2scl/eos_quark_bag.h>
#include <o2scl/eos_quark_cfl6.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mroot_cern.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  {

    double ss1,ss2,ss3,gap1,gap2,gap3,n3,n8;
    double ss12,ss22,ss32,gap12,gap22,gap32,n32,n82;
    eos_quark_cfl6 cfl2;
    thermo th, th2;
    quark u2(cfl2.up_default_mass,6.0);
    quark d2(cfl2.down_default_mass,6.0);
    quark s2(cfl2.strange_default_mass,6.0);
    cfl2.set_quarks(u2,d2,s2);
    cfl2.set_thermo(th2);

    // Parameters
    cfl2.set_parameters();
    cfl2.KD=cfl2.K;

    // Eigenvalues and derivatives
    u2.mu=2.0;
    d2.mu=2.1;
    s2.mu=2.2;
    u2.qq=-1.0;
    d2.qq=-1.1;
    s2.qq=-1.2;
    u2.del=0.5;
    d2.del=0.6;
    s2.del=0.7;
    cfl2.test_derivatives(1.5,0.0,0.0,t);

    // Compare new vs. old when KD=0

    cfl2.KD=0.0;

    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,4.0/hc_mev_fm);

    cfl2.KD=cfl2.K;

    u2.mu=2.5;
    d2.mu=2.6;
    s2.mu=2.7;
    u2.qq=-1.0;
    d2.qq=-1.1;
    s2.qq=-1.2;
    u2.del=0.5;
    d2.del=0.6;
    s2.del=0.7;

    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,4.0/hc_mev_fm);

    /*  
	We can't test everything here, because the densities in the old
	version of eos_quark_cfl6 only correspond to the deriative of the
	pressure with respect to the chemical potential when the gap
	equations have been properly solved.
      
	t.test_rel(ss1,ss12,1.0e-6,"ss1");
	t.test_rel(ss2,ss22,1.0e-6,"ss2");
	t.test_rel(ss3,ss32,1.0e-6,"ss3");
	t.test_rel(gap1,gap12,1.0e-6,"gap1");
	t.test_rel(gap2,gap22,1.0e-6,"gap2");
	t.test_rel(gap3,gap32,1.0e-6,"gap3");
	t.test_rel(u.n,u2.n,1.0e-6,"nu");
	t.test_rel(d.n,d2.n,1.0e-6,"nd");
	t.test_rel(s.n,s2.n,1.0e-6,"ns");
	t.test_rel(u.mu,u2.mu,1.0e-6,"muu");
	t.test_rel(d.mu,d2.mu,1.0e-6,"mud");
	t.test_rel(s.mu,s2.mu,1.0e-6,"mus");
	t.test_rel(u.ed,u2.ed,1.0e-6,"edu");
	t.test_rel(d.ed,d2.ed,1.0e-6,"edd");
	t.test_rel(s.ed,s2.ed,1.0e-6,"eds");
	t.test_rel(u.pr,u2.pr,1.0e-6,"pru");
	t.test_rel(d.pr,d2.pr,1.0e-6,"prd");
	t.test_rel(s.pr,s2.pr,1.0e-6,"prs");
	t.test_rel(th.ed,th2.ed,1.0e-6,"thed");
	t.test_rel(th.en,th2.en,1.0e-2,"then");
    */

    // ----------------------------------------------
    // In a case where KD=0, ensure that the quark 
    // condensate gap equation corresponds to the 
    // derivative of the pressure wrt the quark condensate

    cfl2.G=0.2166586;
    cfl2.GD=0.1575699;
    cfl2.K=0.0;
    cfl2.KD=0.0;
    cfl2.B0=20.07754;
    u2.m=0.0;
    d2.m=0.0;
    s2.m=0.0;
    cfl2.test_derivatives(1.6,0.0,0.0,t);
  
    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,2.0/hc_mev_fm);
    double dpdqq=ss12;
    double pr1=th2.pr;
    u2.qq+=1.0e-3;
    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,2.0/hc_mev_fm);
    u2.qq-=1.0e-3;
    double pr2=th2.pr;
    t.test_rel(-dpdqq,(pr2-pr1)/1.0e-3,1.0e-3,"dpdqq");
  
    // -----------------------------
    // The same test when KD != 0

    cfl2.KD=-cfl2.K;

    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,2.0/hc_mev_fm);
    dpdqq=ss12;
    pr1=th2.pr;
    u2.qq+=1.0e-3;
    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,2.0/hc_mev_fm);
    u2.qq-=1.0e-3;
    pr2=th2.pr;
    t.test_rel(-dpdqq,(pr2-pr1)/1.0e-3,1.0e-3,"dpdqq");

    // -----------------------------
    //0 8.0000e-02 -5.9692e-01 1.1750e-01 9.8591e-01 -2.1482e-01 4.5142e-01

    u2.m=0.0;
    d2.m=0.0;
    s2.m=0.0;

    cfl2.G=0.2166586*0.92;
    cfl2.GD=0.1575699*1.2;
    cfl2.K=0.007878494;
    cfl2.KD=-0.007878494;
    cfl2.B0=20.07754;

    u2.mu=0.98591;
    u2.del=0.11750;
    u2.qq=-0.56962;
    d2.mu=0.98591;
    d2.del=0.11750;
    d2.qq=-0.56962;
    s2.mu=0.98591;
    s2.del=0.11750;
    s2.qq=-0.56962;

    cfl2.calc_eq_temp_p(u2,d2,s2,ss12,ss22,ss32,gap12,gap22,gap32,
			0.0,0.0,n32,n82,th2,2.0/hc_mev_fm);
  }

  t.report();

  return 0;
}
