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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass.h>
#include <o2scl/ldrop_mass.h>
#include <o2scl/ldrop_shell.h>
#include <o2scl/apr_eos.h>
#include <o2scl/skyrme_eos.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  double d1, d2, qual;
  nucmass_ame au;
  ame_load(au,"03");
  ldrop_mass_pair ld;
  ldrop_shell ld2;
  apr_eos apr;

  fermion nrn(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  fermion nrp(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  
  nrn.non_interacting=false;
  nrp.non_interacting=false;
  nrn.inc_rest_mass=true;
  nrp.inc_rest_mass=true;

  skyrme_eos sk;
  skyrme_load(sk,"SLy4");

  cout << "-------------------------------------------------\n" << endl;

  // APR

  ld.set_hadronic_eos_temp(sk);
  ld.set_n_and_p(nrn,nrp);
  ld.n0=0.16;
  ld.n1=0.0;
  cout << "Lead from APR with ldrop_mass_pair: " << endl;
  cout << "Mass excess:\t\t " <<  ld.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
  cout << "Total mass:\t\t " << ld.total_mass(82,126) << endl;
  cout << endl;

  t.test_rel(au.total_mass(82,126),ld.total_mass(82,126),1.0e-3,
	     "ame vs ld2");

  ld2.set_hadronic_eos_temp(sk);
  ld2.set_n_and_p(nrn,nrp);
  ld2.n0=0.16;
  ld2.n1=0.0;
  cout << "Lead from APR with ldrop_shell: " << endl;
  cout << "Mass excess:\t\t " <<  ld2.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld2.binding_energy(82,126)/208.0 << endl;
  cout << "Total mass:\t\t " << ld2.total_mass(82,126) << endl;
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  // fit APR

  cout << "ldrop_mass_pair fit:" << endl;
  nucmass_fit mf;
  mf.set_exp_mass(au);
  mf.even_even=false;
  mf.def_mmin.ntrial*=10;

  ld.n0=0.16;
  ld.n1=0.0;
  ld.surften=1.1;
  ld.coul_coeff=1.0;
  cout << "Before: " << endl;
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << " " << ld.doi << "\n "
       << ld.ss << " " << ld.Epair << endl;
  
  mf.fit(ld,qual);
  mf.eval(ld,qual);
  cout << "After: " << endl;
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << " " << ld.doi << "\n "
       << ld.ss << " " << ld.Epair << endl;
  cout << "Quality: " << qual << endl;

  t.test_rel(qual,2.4629,1.0e-3,"pair qual");
  
  cout << endl;
  
  cout << "-------------------------------------------------\n" << endl;
  
  ld2.n1=ld.n1;
  ld2.n0=ld.n0;
  ld2.surften=ld.surften;
  ld2.coul_coeff=ld.coul_coeff;
  ld2.doi=ld.doi;
  ld2.ss=ld.ss;
  ld2.Epair=ld.Epair;

  cout << "ldrop_mass_shell2 fit:" << endl;

  cout << "Before: " << endl;
  cout << ld2.n1 << " " << ld2.n0 << " " << ld2.surften << " " 
       << ld2.coul_coeff << " " << ld2.doi << "\n "
       << ld2.ss << " " << ld2.Epair << endl;
  cout << ld2.s_a1 << " " << ld2.s_a2 << " " << ld2.s_a3 << " "
       << ld2.s_anp << endl;
  
  mf.def_mmin.ntrial*=2;

  for(size_t k=0;k<1;k++) {

    mf.fit(ld2,qual);
    mf.eval(ld2,qual);
    
    cout << "After: " << endl;
    cout << ld2.n1 << " " << ld2.n0 << " " << ld2.surften << " " 
	 << ld2.coul_coeff << " " << ld2.doi << "\n "
	 << ld2.ss << " " << ld2.Epair << endl;
    cout << ld2.s_a1 << " " << ld2.s_a2 << " " << ld2.s_a3 << " "
	 << ld2.s_anp << endl;
    cout << "Quality: " << qual << endl;

  }
  t.test_rel(qual,1.1189,1.0e-3,"ldrop_shell qual");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  frdm_shell fs;

  mf.eval(fs,qual);
  t.test_rel(qual,10.671,1.0e-3,"fs pre-fit");
  cout << "Quality: " << qual << endl;
  for(size_t i=0;i<3;i++) {
    mf.fit(fs,qual);
    cout << "Quality: " << qual << endl;
  }
  t.test_rel(qual,0.96559,1.0e-2,"fs post-fit");
  cout << endl;

  if (false) {
    ubvector xdm(10);
    dvi_mass dm;
    dm.guess_fun(10,xdm);
    
    nucmass_ame_exp ame13;
    ame_load(ame13);
    mf.set_exp_mass(ame13);
    mf.eval(dm,qual);
    vector_out(cout,10,xdm);
    cout << endl;
    //t.test_rel(qual,10.671,1.0e-3,"dm pre-fit");
    cout << "Quality: " << qual << endl;
    for(size_t i=0;i<3;i++) {
      mf.fit(dm,qual);
      cout << "Quality: " << qual << endl;
    }
    dm.guess_fun(10,xdm);
    vector_out(cout,10,xdm);
    cout << endl;
    //t.test_rel(qual,0.96559,1.0e-2,"dm post-fit");
  }

  /*
    int max_iso=30;
    ubvector vqual;
    ubvector_int n_qual;
    
    mf.eval_isospin(fs,n_qual,vqual,max_iso);
    for(int i=0;i<2*max_iso+1;i++) {
    cout << i-max_iso << " " << n_qual[i] << " ";
    if (n_qual[i]>0) cout << vqual[i] << endl;
    else cout << endl;
    }
    cout << endl;
  */

  t.report();
  return 0;
}

