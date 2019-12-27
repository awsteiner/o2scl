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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ldrop.h>
#include <o2scl/nucmass_ldrop_shell.h>
#include <o2scl/eos_had_skyrme.h>
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

#ifndef O2SCL_FAST_TEST

  double d1, d2, qual;
  nucmass_ame au;
  o2scl_hdf::ame_load_ext(au,"../../data/o2scl/nucmass/ame03.o2","ame03.o2");
  nucmass_ldrop_pair ld;
  nucmass_ldrop_shell ld2;

  fermion nrn(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  fermion nrp(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  
  nrn.non_interacting=false;
  nrp.non_interacting=false;
  nrn.inc_rest_mass=true;
  nrp.inc_rest_mass=true;

  eos_had_skyrme sk;
  o2scl_hdf::skyrme_load(sk,"../../data/o2scl/skdata/SLy4.o2",1);

  cout << "-------------------------------------------------\n" << endl;

  ld.set_eos_had_temp_base(sk);
  ld.set_n_and_p(nrn,nrp);
  ld.n0=0.16;
  ld.n1=0.0;
  cout << "Lead from SLy4 with nucmass_ldrop_pair: " << endl;
  cout << "Mass excess:\t\t " <<  ld.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
  cout << "Total mass:\t\t " << ld.total_mass(82,126) << endl;
  cout << endl;

  t.test_rel(au.total_mass(82,126),ld.total_mass(82,126),1.0e-3,
	     "ame vs ld2");

  ld2.set_eos_had_temp_base(sk);
  ld2.set_n_and_p(nrn,nrp);
  ld2.n0=0.16;
  ld2.n1=0.0;
  cout << "Lead from SLy4 with nucmass_ldrop_shell: " << endl;
  cout << "Mass excess:\t\t " <<  ld2.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld2.binding_energy(82,126)/208.0 << endl;
  cout << "Total mass:\t\t " << ld2.total_mass(82,126) << endl;
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  cout << "nucmass_ldrop_pair fit:" << endl;
  nucmass_fit mf;
  nucdist_set(mf.dist,au);
  mf.even_even=false;
  mf.def_mmin.ntrial*=10;

  ld.n0=0.16;
  ld.n1=0.0;
  ld.surften=1.1;
  ld.coul_coeff=1.0;
  cout << "Parameters before: ";
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << " " << ld.doi << "\n "
       << ld.ss << " " << ld.Epair << endl;
  
  mf.fit(ld,qual);
  mf.eval(ld,qual);
  cout << "Paramters after: ";
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << " " << ld.doi << "\n "
       << ld.ss << " " << ld.Epair << endl;
  cout << "Quality: " << qual << endl;

  t.test_rel(qual,2.4675,1.0e-2,"pair qual");
  
  cout << endl;
  
  cout << "-------------------------------------------------\n" << endl;
  
  ld2.n1=ld.n1;
  ld2.n0=ld.n0;
  ld2.surften=ld.surften;
  ld2.coul_coeff=ld.coul_coeff;
  ld2.doi=ld.doi;
  ld2.ss=ld.ss;
  ld2.Epair=ld.Epair;

  cout << "nucmass_nucmass_ldrop_shell fit:" << endl;

  cout << "Paramters before: ";
  cout << ld2.n1 << " " << ld2.n0 << " " << ld2.surften << " " 
       << ld2.coul_coeff << " " << ld2.doi << "\n "
       << ld2.ss << " " << ld2.Epair << endl;
  cout << ld2.s_a1 << " " << ld2.s_a2 << " " << ld2.s_a3 << " "
       << ld2.s_anp << endl;
  
  mf.def_mmin.ntrial*=2;

  for(size_t k=0;k<1;k++) {

    mf.fit(ld2,qual);
    mf.eval(ld2,qual);
    
    cout << "Paramters after: ";
    cout << ld2.n1 << " " << ld2.n0 << " " << ld2.surften << " " 
	 << ld2.coul_coeff << " " << ld2.doi << "\n "
	 << ld2.ss << " " << ld2.Epair << endl;
    cout << ld2.s_a1 << " " << ld2.s_a2 << " " << ld2.s_a3 << " "
	 << ld2.s_anp << endl;
    cout << "Quality: " << qual << endl;

  }
  t.test_rel(qual,1.1353,2.0e-2,"nucmass_ldrop_shell qual");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  cout << "nucmass_frdm_shell fit:" << endl;
  nucmass_frdm_shell fs;

  mf.eval(fs,qual);
  t.test_rel(qual,10.671,1.0e-3,"fs pre-fit");
  cout << "Quality: " << qual << endl;
  for(size_t i=0;i<3;i++) {
    mf.fit(fs,qual);
    cout << "Quality: " << qual << endl;
  }
  t.test_rel(qual,0.96559,1.0e-2,"fs post-fit");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  cout << "nucmass_dvi fit:" << endl;
  ubvector xdvi(10);
  nucmass_dvi dvi;
  dvi.guess_fun(10,xdvi);
    
  nucmass_ame ame13;
  o2scl_hdf::ame_load_ext(ame13,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2",true);
  nucdist_set(mf.dist,ame13);
  mf.eval(dvi,qual);
  cout << "Parameters before: ";
  vector_out(cout,10,xdvi);
  cout << endl;
  t.test_rel(qual,24.99316,1.0e-3,"dvi pre-fit");
  cout << "Quality: " << qual << endl;
  for(size_t i=0;i<3;i++) {
    mf.fit(dvi,qual);
    cout << "Quality: " << qual << endl;
  }
  dvi.guess_fun(10,xdvi);
  cout << "Parameters after: ";
  vector_out(cout,10,xdvi);
  cout << endl;
  // This is a little better fit than FRDM above, but maybe
  // because only the experimental masses are used
  t.test_rel(qual,0.8850699,1.0e-2,"dvi post-fit");

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

#endif
  
  t.report();
  return 0;
}

