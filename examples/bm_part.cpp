/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/sn_fermion.h>
#include <o2scl/eff_fermion.h>
#include <o2scl/rel_fermion.h>
#include <o2scl/sn_classical.h>
#include <o2scl/classical.h>
#include <o2scl/nonrel_fermion.h>

using namespace std;
using namespace o2scl;

// Things to check. Ensure that
// - calc_mu() and calc_density() both work with a bad guess
// - The fractional error in sn_fermion is sufficiently small
// - sn_fermion matches eff_fermion and rel_fermion
// - zero temperature works
// - classical matches sn_classical
//
// Should also compare with mathematica (now partially done
// in bm_part2.cpp)
//
// Remember we're testing performance not correctness (6/14/09 -
// Correctness testing ought to be done in *_ts.cpp, not here.)

class part_test {
public:

  classical cl;
  sn_classical sc;
  sn_fermion sf;
  eff_fermion ef;
  rel_fermion rf;
  nonrel_fermion nf;

  part_test() {
    cl.init(0.0,2.0);
    sc.init(0.0,2.0);
    sf.init(0.0,2.0);
    ef.init(0.0,2.0);
    rf.init(0.0,2.0);
    nf.init(0.0,2.0);
  }

  int copy(part &src, part &dest) {
    dest.n=src.n;
    dest.mu=src.mu;
    dest.ed=src.ed;
    dest.pr=src.pr;
    dest.en=src.en;
    return 0;
  }

  // Compare two particles individually setting val[0]...val[4]
  int compare(part &one, part &two, double *val) {
    if (fabs(one.n)>1.0e-99) {
      val[0]=fabs(one.n-two.n)/fabs(one.n);
    } else val[0]=-1.0;
    if (fabs(one.mu)>1.0e-99) {
      val[1]=fabs(one.mu-two.mu)/fabs(one.mu);
    } else val[1]=-1.0;
    if (fabs(one.ed)>1.0e-99) {
      val[2]=fabs(one.ed-two.ed)/fabs(one.ed);
    } else val[2]=-1.0;
    if (fabs(one.pr)>1.0e-99) {
      val[3]=fabs(one.pr-two.pr)/fabs(one.pr);
    } else val[3]=-1.0;
    if (fabs(one.en)>1.0e-99) {
      val[4]=fabs(one.en-two.en)/fabs(one.en);
    } else val[4]=-1.0;
    return 0;
  }
  
  // Run a test at the specified T, mu, and m, returning val[]
  int test(double T, double psi, double mot, double val[]) {
    int ret, ix=0;
    fermion temp;

    // Compare eff_fermion::calc_mu() to eff_fermion::calc_density()
    ef.m=T*mot;
    ef.mu=psi*T+sf.m;
    ret=ef.calc_mu(T);
    copy(ef,temp);
    val[ix++]=ret;
    ret=ef.calc_density(T);
    compare(temp,ef,&val[ix]); 
    ix+=5;
    val[ix++]=ret;
    // Try bad guess #1 for mu
    ef.nu=-1.0e6;
    ret=ef.calc_density(T);
    val[ix++]=ret;
    // Try bad guess #2 for mu
    ef.nu=1.0e6;
    ret=ef.calc_density(T);
    val[ix++]=ret;

    // Compare rel_fermion::calc_mu() to rel_fermion::calc_density()
    rf.m=T*mot;
    rf.mu=psi*T+sf.m;
    ret=rf.calc_mu(T);
    compare(temp,rf,&val[ix]);
    ix+=5;
    val[ix++]=ret;
    ret=rf.calc_density(T);
    compare(temp,rf,&val[ix]);
    ix+=5;
    val[ix++]=ret;
    // Try bad guess #1 for mu
    rf.nu=-1.0e6;
    ret=rf.calc_density(T);
    val[ix++]=ret;
    // Try bad guess #2 for mu
    rf.nu=1.0e6;
    ret=rf.calc_density(T);
    val[ix++]=ret;

    // Compare sn_fermion::calc_mu() to sn_fermion::calc_density()
    sf.m=T*mot;
    sf.mu=psi*T+sf.m;
    ret=sf.calc_mu(T);
    compare(temp,sf,&val[ix]);
    val[ix++]=ret;
    // Check uncertainties in sn_fermion
    val[ix++]=fabs(sf.unc.n/sf.n);
    val[ix++]=fabs(sf.unc.mu/sf.mu);
    val[ix++]=fabs(sf.unc.ed/sf.ed);
    val[ix++]=fabs(sf.unc.pr/sf.pr);
    val[ix++]=fabs(sf.unc.en/sf.en);
    ret=sf.calc_density(T);
    compare(temp,sf,&val[ix]);
    val[ix++]=ret;
    // Try bad guess #1 for mu
    sf.nu=-1.0e6;
    ret=sf.calc_density(T);
    val[ix++]=ret;
    // Try bad guess #2 for mu
    sf.nu=1.0e6;
    ret=sf.calc_density(T);
    val[ix++]=ret;

    return ix;
  }
};

int main(void) {

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);

  part_test pt;
  double vals[100];
  
  for(double T=1.0e-3;T<=1.01e3;T*=1.0e3) {
    double psi=-20.0;
    while (psi<800.0) {
      for(double mot=1.0e-2;mot<=1.01e1;mot*=10.0) {
	int ntot=pt.test(T,psi,mot,vals);
	cout.precision(3);
	cout << T << " " << psi << " " << mot << " " << ntot << " ";
	for(size_t k=0;k<9;k++) {
	  cout << vals[k] << " ";
	}
	cout << endl;
	cout << "  ";
	for(size_t k=9;k<20;k++) {
	  cout << vals[k] << " ";
	}
	cout << endl;
	cout << "  ";
	for(size_t k=20;k<31;k++) {
	  cout << vals[k] << " ";
	}
	cout << endl;
      }
      if (psi>0.0) psi*=4.0;
      if (psi<0.0) psi/=4.0;
      if (psi>-1.0 && psi<0.0) psi*=-1.0;
    }
  }

  return 0;
}

