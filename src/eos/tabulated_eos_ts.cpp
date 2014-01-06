/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <o2scl/nonrel_fermion.h>
#include <o2scl/cold_nstar.h>
#include <o2scl/tabulated_eos.h>
#include <o2scl/schematic_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  tabulated_eos te;
  schematic_eos se;

  se.eoa=-16.0/hc_mev_fm;
  se.comp=200.0/hc_mev_fm;
  se.a=18.0/hc_mev_fm;
  se.b=12.0/hc_mev_fm;
  se.kprime=0.0;
  se.n0=0.16;

  fermion ne(939.0/hc_mev_fm,2.0), pr(939.0/hc_mev_fm,2.0);
  ne.non_interacting=false;
  pr.non_interacting=false;
  thermo hb;

  double tol=8.0e-4;

  ubvector rho(16), nuc(16), neut(16);
  for(size_t i=0;i<16;i++) {
    double barn=((double)i+1)/100.0;
    rho[i]=barn;
    ne.n=barn/2.0;
    pr.n=barn/2.0;
    se.calc_e(ne,pr,hb);
    nuc[i]=(hb.ed/barn-ne.m)*hc_mev_fm;
    ne.n=barn;
    pr.n=0.0;
    se.calc_e(ne,pr,hb);
    neut[i]=(hb.ed/barn-ne.m)*hc_mev_fm;
  }
  te.set_eos(16,rho,nuc,neut);

  for(double nb=0.01;nb<=0.1601;nb+=0.01) {
    ne.n=nb/2.0;
    pr.n=nb/2.0;
    se.calc_e(ne,pr,hb);
    double t1=ne.mu;
    double t2=pr.mu;
    double t3=hb.ed;
    te.calc_e(ne,pr,hb);
    t.test_rel(t1,ne.mu,tol,"mun");
    t.test_rel(t2,pr.mu,tol,"mup");
    t.test_rel(t3,hb.ed,tol,"ed");
  }
  for(double nb=0.01;nb<=0.1601;nb+=0.01) {
    ne.n=nb;
    pr.n=0.0;
    se.calc_e(ne,pr,hb);
    double t1=ne.mu;
    double t2=pr.mu;
    double t3=hb.ed;
    te.calc_e(ne,pr,hb);
    t.test_rel(t1,ne.mu,tol,"mun");
    t.test_rel(t2,pr.mu,tol,"mup");
    t.test_rel(t3,hb.ed,tol,"ed");
  }
  for(double nb=0.01;nb<=0.1601;nb+=0.01) {
    pr.n=nb/4.0;
    ne.n=nb-pr.n;
    se.calc_e(ne,pr,hb);
    double t1=ne.mu;
    double t2=pr.mu;
    double t3=hb.ed;
    te.calc_e(ne,pr,hb);
    t.test_rel(t1,ne.mu,tol,"mun");
    t.test_rel(t2,pr.mu,tol,"mup");
    t.test_rel(t3,hb.ed,tol,"ed");
  }

  t.report();
  return 0;
}
