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
#include <fstream>
#include <vector>
#include <o2scl/sn_fermion.h>
#include <o2scl/eff_fermion.h>
#include <o2scl/rel_fermion.h>

using namespace std;
using namespace o2scl;

typedef struct {
  double lmot;
  double lpsi;
  double n;
  double ed;
  double pr;
  double en1;
  double en2;
  double err;
} pm_entry;

class part_test {
public:

  sn_fermion sf;
  eff_fermion ef;
  rel_fermion rf;
  size_t n;

  std::vector<pm_entry> tab;

  part_test() {

    sf.init(0.0,2.0);
    ef.init(0.0,2.0);
    rf.init(0.0,2.0);

    if (0) {
      // This gives results to within 1 part in 10^{10}
      rf.def_dit.tolx/=1.0e2;
      rf.def_dit.tolf/=1.0e2;
      rf.def_nit.tolx/=1.0e2;
      rf.def_nit.tolf/=1.0e2;
      rf.upper_limit_fac=30.0;
    }
    
  }

  int read_file() {
    ifstream fin("../src/internal/partmath.dat");
    for(size_t i=0;i<35;i++) {
      pm_entry p;
      fin >> p.lmot >> p.lpsi >> p.n >> p.ed >> p.pr 
	  >> p.en1 >> p.en2 >> p.err;
      tab.push_back(p);
    }
    fin.close();
    n=tab.size();
    return 0;
  }

  int test() {

    size_t cnt=0;

    part sf_dev, ef_dev, rf_dev;
    part sf_bad, ef_bad, rf_bad;
    part exact;

    int verbose=0;

    // k=0 is with rest mass, k=1 is without
    for(size_t k=0;k<2;k++) {
      
      // Reset statistics

      sf_dev.n=0.0; sf_dev.ed=0.0; sf_dev.pr=0.0; sf_dev.en=0.0;
      ef_dev.n=0.0; ef_dev.ed=0.0; ef_dev.pr=0.0; ef_dev.en=0.0;
      rf_dev.n=0.0; rf_dev.ed=0.0; rf_dev.pr=0.0; rf_dev.en=0.0;
      
      sf_bad.n=0.0; sf_bad.ed=0.0; sf_bad.pr=0.0; sf_bad.en=0.0;
      ef_bad.n=0.0; ef_bad.ed=0.0; ef_bad.pr=0.0; ef_bad.en=0.0;
      rf_bad.n=0.0; rf_bad.ed=0.0; rf_bad.pr=0.0; rf_bad.en=0.0;
      
      // Temperature loop
      for(double T=1.0e-2;T<=1.001e2;T*=1.0e2) {

	// Loop over each point in the data file
	for(size_t i=0;i<n;i++) {

	  double m, mu;

	  if (k==0) {

	    sf.inc_rest_mass=true;
	    rf.inc_rest_mass=true;
	    ef.inc_rest_mass=true;
	    
	    m=pow(10.0,tab[i].lmot)*T;
	    mu=m+T*pow(10.0,tab[i].lpsi);

	  } else {

	    sf.inc_rest_mass=false;
	    rf.inc_rest_mass=false;
	    ef.inc_rest_mass=false;
	    
	    m=pow(10.0,tab[i].lmot)*T;
	    mu=T*pow(10.0,tab[i].lpsi);

	  }
      
	  if (verbose>0) {
	    cout << "T,m,mu: " << T << " " << m << " " << mu << endl;
	  }

	  sf.m=m;
	  sf.mu=mu;
	  sf.calc_mu(T);

	  ef.m=m;
	  ef.mu=mu;
	  ef.calc_mu(T);

	  rf.m=m;
	  rf.mu=mu;
	  rf.calc_mu(T);

	  exact.n=tab[i].n*pow(T,3.0);
	  if (k==0) {
	    exact.ed=tab[i].ed*pow(T,4.0);
	  } else {
	    exact.ed=tab[i].ed*pow(T,4.0)-exact.n*m;
	  }
	  exact.pr=tab[i].pr*pow(T,4.0);
	  exact.en=tab[i].en1*pow(T,3.0);

	  if (verbose>0) {
	    cout << "sf: " << sf.n << " " << sf.ed << " " << sf.pr << " " 
		 << sf.en << endl;
	    cout << "ef: " << ef.n << " " << ef.ed << " " << ef.pr << " " 
		 << ef.en << endl;
	    cout << "rf: " << rf.n << " " << rf.ed << " " << rf.pr << " " 
		 << rf.en << endl;
	    cout << "ex: " << exact.n << " " << exact.ed << " " 
		 << exact.pr << " " << exact.en << endl;
	    cout << endl;
	  }
      
	  sf_dev.n+=fabs((sf.n-exact.n)/exact.n);
	  sf_dev.ed+=fabs((sf.ed-exact.ed)/exact.ed);
	  sf_dev.pr+=fabs((sf.pr-exact.pr)/exact.pr);
	  sf_dev.en+=fabs((sf.en-exact.en)/exact.en);
	  ef_dev.n+=fabs((ef.n-exact.n)/exact.n);
	  ef_dev.ed+=fabs((ef.ed-exact.ed)/exact.ed);
	  ef_dev.pr+=fabs((ef.pr-exact.pr)/exact.pr);
	  ef_dev.en+=fabs((ef.en-exact.en)/exact.en);
	  rf_dev.n+=fabs((rf.n-exact.n)/exact.n);
	  rf_dev.ed+=fabs((rf.ed-exact.ed)/exact.ed);
	  rf_dev.pr+=fabs((rf.pr-exact.pr)/exact.pr);
	  rf_dev.en+=fabs((rf.en-exact.en)/exact.en);
	  cnt++;

	  if (fabs((sf.n-exact.n)/exact.n)>sf_bad.n) {
	    sf_bad.n=fabs((sf.n-exact.n)/exact.n);
	  }
	  if (fabs((sf.ed-exact.ed)/exact.ed)>sf_bad.ed) {
	    sf_bad.ed=fabs((sf.ed-exact.ed)/exact.ed);
	  }
	  if (fabs((sf.pr-exact.pr)/exact.pr)>sf_bad.pr) {
	    sf_bad.pr=fabs((sf.pr-exact.pr)/exact.pr);
	  }
	  if (fabs((sf.en-exact.en)/exact.en)>sf_bad.en) {
	    sf_bad.en=fabs((sf.en-exact.en)/exact.en);
	  }
	  if (fabs((ef.n-exact.n)/exact.n)>ef_bad.n) {
	    ef_bad.n=fabs((ef.n-exact.n)/exact.n);
	  }
	  if (fabs((ef.ed-exact.ed)/exact.ed)>ef_bad.ed) {
	    ef_bad.ed=fabs((ef.ed-exact.ed)/exact.ed);
	  }
	  if (fabs((ef.pr-exact.pr)/exact.pr)>ef_bad.pr) {
	    ef_bad.pr=fabs((ef.pr-exact.pr)/exact.pr);
	  }
	  if (fabs((ef.en-exact.en)/exact.en)>ef_bad.en) {
	    ef_bad.en=fabs((ef.en-exact.en)/exact.en);
	  }
	  if (fabs((rf.n-exact.n)/exact.n)>rf_bad.n) {
	    rf_bad.n=fabs((rf.n-exact.n)/exact.n);
	  }
	  if (fabs((rf.ed-exact.ed)/exact.ed)>rf_bad.ed) {
	    rf_bad.ed=fabs((rf.ed-exact.ed)/exact.ed);
	  }
	  if (fabs((rf.pr-exact.pr)/exact.pr)>rf_bad.pr) {
	    rf_bad.pr=fabs((rf.pr-exact.pr)/exact.pr);
	  }
	  if (fabs((rf.en-exact.en)/exact.en)>rf_bad.en) {
	    rf_bad.en=fabs((rf.en-exact.en)/exact.en);
	  }

	  if (false) {

	    // Now use calc_density()

	    sf.mu=0.0;
	    sf.n=exact.n;
	    sf.calc_density(T);
	  
	    ef.mu=0.0;
	    ef.n=exact.n;
	    ef.calc_density(T);
	  
	    rf.mu=0.0;
	    rf.n=exact.n;
	    rf.calc_density(T);

	    if (verbose>0) {
	      cout << "rf: " << rf.mu << " " << rf.ed << " " << rf.pr << " " 
		   << rf.en << endl;
	      cout << "ex: " << mu << " " << exact.ed << " " 
		   << exact.pr << " " << exact.en << endl;
	    }
	  
	    sf_dev.mu+=fabs((sf.mu-mu)/mu);
	    sf_dev.ed+=fabs((sf.ed-exact.ed)/exact.ed);
	    sf_dev.pr+=fabs((sf.pr-exact.pr)/exact.pr);
	    sf_dev.en+=fabs((sf.en-exact.en)/exact.en);
	    ef_dev.mu+=fabs((ef.mu-mu)/mu);
	    ef_dev.ed+=fabs((ef.ed-exact.ed)/exact.ed);
	    ef_dev.pr+=fabs((ef.pr-exact.pr)/exact.pr);
	    ef_dev.en+=fabs((ef.en-exact.en)/exact.en);
	    rf_dev.mu+=fabs((rf.mu-mu)/mu);
	    rf_dev.ed+=fabs((rf.ed-exact.ed)/exact.ed);
	    rf_dev.pr+=fabs((rf.pr-exact.pr)/exact.pr);
	    rf_dev.en+=fabs((rf.en-exact.en)/exact.en);
	    cnt++;

	    exit(-1);

	  }

	  // End of loop over points in data file
	}
	// End of temperature loop
      }

      sf_dev.n/=cnt;
      sf_dev.ed/=cnt;
      sf_dev.pr/=cnt;
      sf_dev.en/=cnt;
      ef_dev.n/=cnt;
      ef_dev.ed/=cnt;
      ef_dev.pr/=cnt;
      ef_dev.en/=cnt;
      rf_dev.n/=cnt;
      rf_dev.ed/=cnt;
      rf_dev.pr/=cnt;
      rf_dev.en/=cnt;

      if (k==0) {
	cout << "Include rest mass:\n" << endl;
      } else {
	cout << "Without rest mass:\n" << endl;
      }
    
      cout << "Average: " << endl;
      cout << "sf: " << sf_dev.n << " " << sf_dev.ed << " " 
	   << sf_dev.pr << " " << sf_dev.en << endl;
      cout << "ef: " << ef_dev.n << " " << ef_dev.ed << " " 
	   << ef_dev.pr << " " << ef_dev.en << endl;
      cout << "rf: " << rf_dev.n << " " << rf_dev.ed << " " 
	   << rf_dev.pr << " " << rf_dev.en << endl;
      cout << endl;

      cout << "Worst case: " << endl;
      cout << "sf: " << sf_bad.n << " " << sf_bad.ed << " " 
	   << sf_bad.pr << " " << sf_bad.en << endl;
      cout << "ef: " << ef_bad.n << " " << ef_bad.ed << " " 
	   << ef_bad.pr << " " << ef_bad.en << endl;
      cout << "rf: " << rf_bad.n << " " << rf_bad.ed << " " 
	   << rf_bad.pr << " " << rf_bad.en << endl;
      cout << endl;

    }

    return 0;
  }

};

int main(void) {

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);

  part_test pt;
  pt.read_file();
  pt.test();

  return 0;
}

