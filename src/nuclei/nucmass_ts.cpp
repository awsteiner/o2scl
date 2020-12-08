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
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_ktuy.h>
#include <o2scl/nucmass_dglg.h>
#include <o2scl/nucmass_wlw.h>
#include <o2scl/nucmass_sdnp.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  // Test parse_elstring()
  {
    nucmass_info nmi;
    int tN, tZ, tA; 
    nmi.parse_elstring("Pb208",tZ,tN,tA);
    t.test_gen(tZ==82,"parse3");
    t.test_gen(tA==208,"parse4");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("Pb-208",tZ,tN,tA);
    t.test_gen(tZ==82,"parse5");
    t.test_gen(tA==208,"parse6");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("Pb 208",tZ,tN,tA);
    t.test_gen(tZ==82,"parse7");
    t.test_gen(tA==208,"parse8");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("pb-208",tZ,tN,tA);
    t.test_gen(tZ==82,"parse9");
    t.test_gen(tA==208,"parse10");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("pb 208",tZ,tN,tA);
    t.test_gen(tZ==82,"parse11");
    t.test_gen(tA==208,"parse12");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("H4",tZ,tN,tA);
    t.test_gen(tZ==1,"parse13");
    t.test_gen(tA==4,"parse14");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("208Pb",tZ,tN,tA);
    t.test_gen(tZ==82,"parse15");
    t.test_gen(tA==208,"parse16");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("208-Pb",tZ,tN,tA);
    t.test_gen(tZ==82,"parse17");
    t.test_gen(tA==208,"parse18");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("208 Pb",tZ,tN,tA);
    t.test_gen(tZ==82,"parse19");
    t.test_gen(tA==208,"parse20");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("208-pb",tZ,tN,tA);
    t.test_gen(tZ==82,"parse21");
    t.test_gen(tA==208,"parse22");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("208 pb",tZ,tN,tA);
    t.test_gen(tZ==82,"parse23");
    t.test_gen(tA==208,"parse24");
    tN=0;
    tZ=0;
    tA=0;
    nmi.parse_elstring("4H",tZ,tN,tA);
    t.test_gen(tZ==1,"parse25");
    t.test_gen(tA==4,"parse26");
  }

  // Create an instance of all the various mass formulae

  nucmass_ame ame;
  o2scl_hdf::ame_load_ext(ame,"../../data/o2scl/nucmass/ame03.o2","ame03.o2");
  nucmass_ame ame95rmd;
  o2scl_hdf::ame_load_ext(ame95rmd,"../../data/o2scl/nucmass/ame95rmd.o2",
		      "ame95rmd.o2");
  nucmass_ame ame03round;
  o2scl_hdf::ame_load_ext(ame03round,"../../data/o2scl/nucmass/ame03round.o2",
		      "ame03round.o2");
  nucmass_ame ame95exp;
  o2scl_hdf::ame_load_ext(ame95exp,"../../data/o2scl/nucmass/ame95exp.o2",
		      "ame95exp.o2");
  nucmass_ame ame12;
  o2scl_hdf::ame_load_ext(ame12,"../../data/o2scl/nucmass/ame12.o2","ame12.o2");

  // Output the references
  cout << "References: " << endl;
  cout << ame.get_reference() << endl;
  cout << ame03round.get_reference() << endl;
  cout << ame95rmd.get_reference() << endl;
  cout << ame95exp.get_reference() << endl;
  cout << ame12.get_reference() << endl;

  nucmass_semi_empirical sm;
  nucmass_mnmsk m95;
  o2scl_hdf::mnmsk_load(m95,"mnmsk97","../../data/o2scl/nucmass/mnmsk.o2");
  
  nucmass_ktuy kt;
  kt.load("../../data/o2scl/nucmass/ktuy04.o2",1);
  nucmass_ktuy kt2;
  kt2.load("../../data/o2scl/nucmass/ktuy05.o2",1);

  nucmass_wlw wlw1;
  wlw1.load("../../data/o2scl/nucmass/wlw10.o2",1);
  nucmass_wlw wlw2;
  wlw2.load("../../data/o2scl/nucmass/wllw10.o2",1);
  nucmass_wlw wlw3;
  wlw3.load("../../data/o2scl/nucmass/lwdw11.o2",1);
  //cout << "Here2b4." << endl;
  //nucmass_wlw wlw4("../../data/o2scl/nucmass/wl11.o2",1);
  //cout << "Here2b5." << endl;
  //  nucmass_wlw wlw5("../../data/o2scl/nucmass/wlwm14.o2",1);
  
  nucmass_sdnp sdnp1;
  sdnp1.load("../../data/o2scl/nucmass/sdnp03.o2",1);
  nucmass_sdnp sdnp2;
  sdnp2.load("../../data/o2scl/nucmass/sd_skp_04.o2",1);
  nucmass_sdnp sdnp3;
  sdnp3.load("../../data/o2scl/nucmass/sd_sly4_04.o2",1);
  
  nucmass_hfb hfb2;
  o2scl_hdf::hfb_load(hfb2,2,"../../data/o2scl/nucmass/");
  nucmass_hfb hfb8;
  o2scl_hdf::hfb_load(hfb8,8,"../../data/o2scl/nucmass/");
  nucmass_hfb hfb14;
  o2scl_hdf::hfb_load(hfb14,14,"../../data/o2scl/nucmass/");
  nucmass_hfb hfb14_v0;
  o2scl_hdf::hfb_load(hfb14_v0,15,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb17;
  o2scl_hdf::hfb_sp_load(hfb17,17,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb21;
  o2scl_hdf::hfb_sp_load(hfb21,21,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb22;
  o2scl_hdf::hfb_sp_load(hfb22,22,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb23;
  o2scl_hdf::hfb_sp_load(hfb23,23,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb24;
  o2scl_hdf::hfb_sp_load(hfb24,24,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb25;
  o2scl_hdf::hfb_sp_load(hfb25,25,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb26;
  o2scl_hdf::hfb_sp_load(hfb26,26,"../../data/o2scl/nucmass/");
  nucmass_hfb_sp hfb27;
  o2scl_hdf::hfb_sp_load(hfb27,27,"../../data/o2scl/nucmass/");
  //nucmass_dglg dglg;

  // Set up generic pointers for testing
  nucmass_table *nmd[27]={&ame,&ame95rmd,&ame03round,&ame95exp,
			  &m95,&kt,&kt2,&hfb2,&hfb8,
			  &hfb14,&hfb17,&hfb21,&hfb22,&hfb23,&hfb24,&hfb25,
			  &hfb26,&hfb27,&wlw1,&wlw1,&wlw2,&wlw3,&wlw3,
			  &wlw3,&sdnp1,&sdnp2,&sdnp3};

  // Test the spins obtained from mnmsk

  nucmass_mnmsk::entry mmk;
  mmk=m95.get_ZN(83,125);
  t.test_gen(mmk.spinp==((string)"9/2-"),"spinp");
  t.test_gen(mmk.spinn==((string)"1/2-"),"spinn");

  // Test the various formulae for the binding energy of lead

  t.test_rel(ame12.binding_energy(82,126)/208.0,-7.867,1.0e-4,"ame12 be");
  t.test_rel(ame.binding_energy(82,126)/208.0,-7.867,1.0e-4,"ame be");
  t.test_rel(sm.binding_energy(82,126)/208.0,-7.867,5.0e-2,"sm be");
  t.test_rel(m95.binding_energy(82,126)/208.0,-7.867,5.0e-4,"m95 be");
  t.test_rel(hfb2.binding_energy(82,126)/208.0,-7.867,2.0e-3,"hfb be");
  t.test_rel(hfb8.binding_energy(82,126)/208.0,-7.867,1.0e-3,"hfb2 be");
  t.test_rel(hfb14.binding_energy(82,126)/208.0,-7.867,1.0e-3,"hfb3 be");
  t.test_rel(kt.binding_energy(82,126)/208.0,-7.867,1.0e-3,"kt be");
  t.test_rel(kt2.binding_energy(82,126)/208.0,-7.867,1.0e-3,"kt2 be");

  // Test the binding energy and mass excess from get_nucleus()
  double mass_neutron=o2scl_mks::mass_neutron*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_proton=o2scl_mks::mass_proton*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_electron=o2scl_mks::mass_electron*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);
  double mass_amu=o2scl_mks::unified_atomic_mass*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);

  for(size_t i=0;i<27;i++) {
    nucleus n;
    nmd[i]->get_nucleus(82,126,n);
    t.test_rel(n.be*o2scl_const::hc_mev_fm/208.0,-7.867,4.0e-3,"ptr be");
    t.test_rel(n.m-126.0*mass_neutron-82.0*mass_proton,n.be,1.0e-12,"ptr be2");
    t.test_rel(n.m+82.0*mass_electron-208.0*mass_amu,n.mex,1.0e-11,
	       "ptr mex");
  }

  // Test size of ame95rmd
  
  t.test_gen(ame95rmd.get_nentries()==2931,"ame.n");

  // Test nucmass_radius
  nucmass_radius nr;
  double rho0, N, N_err;
  nr.eval_N_err(6.0,0.5,0.08,N,N_err);
  cout << N << " " << N_err << endl;
  t.test_rel(N,77.3433,1.0e-8,"nucmass_radius");

  // Test spinp_to_int()
  nucmass_info nmi;
  cout << nmi.spinp_to_int("9") << endl;
  cout << nmi.spinp_to_int("99") << endl;
  cout << nmi.spinp_to_int("9+") << endl;
  cout << nmi.spinp_to_int("+9") << endl;
  cout << nmi.spinp_to_int("99+") << endl;
  cout << nmi.spinp_to_int("+99") << endl;
  cout << nmi.spinp_to_int("9/2+") << endl;
  cout << nmi.spinp_to_int("+9/2") << endl;
  cout << nmi.spinp_to_int("99/2+") << endl;
  cout << nmi.spinp_to_int("+99/2") << endl;
  cout << nmi.spinp_to_int("9-") << endl;
  cout << nmi.spinp_to_int("-9") << endl;
  cout << nmi.spinp_to_int("99-") << endl;
  cout << nmi.spinp_to_int("-99") << endl;
  cout << nmi.spinp_to_int("9/2-") << endl;
  cout << nmi.spinp_to_int("-9/2") << endl;
  cout << nmi.spinp_to_int("99/2-") << endl;
  cout << nmi.spinp_to_int("-99/2") << endl;
  
  t.report();
  return 0;
}

