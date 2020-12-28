/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <o2scl/eos_had_apr.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/nucmass_ldrop.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  mmin_simp2<> gms;
  gms.err_nonconv=false;
  
  double d1, d2, qual;
  nucmass_ame au;
  o2scl_hdf::ame_load_ext(au,"../../data/o2scl/nucmass/ame12.o2","ame12.o2");

  nucmass_semi_empirical sm;
  nucmass_ldrop ld;
  eos_had_apr apr;

  cout << "\nTest rms_radius: " << endl;
  nucmass_radius rr;
  double r1, r2, r3;
  // First compute the radii assuming the central density is 0.092
  rr.eval_rms_rho(0.092,126,0.3,r1,r2,r3);
  // Now go back and recompute the central density from the radii
  rr.eval_rms_rsq(r2,126,0.3,r1,r2,r3);
  t.test_rel(r1,0.092,1.0e-6,"rms radius");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

#ifdef O2SCL_NEVER_DEFINED
  if (false) {
    nucmass_ldrop_skin tst;
    tst.new_skin_mode=true;
    tst.set_eos_had_temp_base(apr);
    tst.set_n_and_p(nrn,nrp);
    tst.n0=0.17;
    tst.n1=-0.05;
    double be;
    
    double chi=1.0e-2, Rn, Rws, nN, t1, t2, f1, f2;
    
    thermo th;
    nrn.n=1.0e-3;
    nrp.n=0.0;
    apr.calc_temp_e(nrn,nrp,1.0e-3,th);
    double fdrip=th.ed-1.0e-3*th.en-nrn.n*nrn.m;

    be=tst.drip_binding_energy_d(40,80,0.0,1.0e-3,chi,1.0e-3)/hc_mev_fm;
    cout << "be: " << be << endl;
    Rn=tst.Rn;
    Rws=Rn/cbrt(chi);
    nN=3.0/4.0/pi/pow(Rws,3.0);
    cout << nN << " " << Rws << endl;
    t1=be*nN;
    t2=(1.0-chi)*fdrip;
    f1=t1+t2;
    cout << t1 << " " << t2 << " " << f1 << endl;
    
    tst.rel_vacuum=false;
    be=tst.drip_binding_energy_d(40,80,0.0,1.0e-3,chi,1.0e-3)/hc_mev_fm;
    Rn=tst.Rn;
    Rws=Rn/cbrt(chi);
    nN=3.0/4.0/pi/pow(Rws,3.0);
    t1=be*nN;
    t2=fdrip;
    f2=t1+t2;
    cout << t1 << " " << t2 << " " << f2 << endl;
    cout << f1-f2 << endl;

    exit(-1);
  }
#endif

  cout << "-------------------------------------------------\n" << endl;

  cout << "Lead from AME: " << endl;
  cout <<  "Mass excess:\t\t " <<au.mass_excess(82,126) << endl;
  cout <<  "Binding energy:\t\t " <<au.binding_energy(82,126)/208.0 << endl;
  cout <<  "Total mass:\t\t " << au.total_mass(82,126) << endl;
  cout << endl;
  
  cout << "-------------------------------------------------\n" << endl;

  // RMF

  ld.def_had_eos.def_mroot.def_jac.set_epsmin(1.0e-15);

  double eoa;
  ld.n0=ld.def_had_eos.fn0(0.0,eoa);
  cout << "Lead from RMF: " << endl;
  cout << "Mass excess:\t\t " << ld.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
  cout <<  "Total mass:\t\t " << ld.total_mass(82,126) << endl;

  t.test_rel(au.total_mass(82,126),ld.total_mass(82,126),1.0e-3,
	     "ame vs ld");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  // fit RMF

  nucmass_fit mf3;
  mf3.even_even=false;
  nucmass_ame ame;
  o2scl_hdf::ame_load_ext(ame,"../../data/o2scl/nucmass/ame12.o2","ame12.o2");
  nucdist_set(mf3.dist,ame);
  //mf3.set_exp_mass(ame);

  cout << "RMF: " << endl;
  ld.n1=0.02;
  ld.n0=0.15;
  ld.surften=1.1;
  ld.coul_coeff=1.0;
  
  mf3.eval(ld,qual);
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << endl;
  cout << qual << endl;

  if (false) {
    
    cout << endl;
    cout << "Fit RMF: " << endl;
    mf3.fit(ld,qual);
    cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
	 << ld.coul_coeff << endl;
    cout << qual << endl;
    cout << endl;

    cout << "-------------------------------------------------\n" << endl;
      
    // RMF after fit
      
    cout << "RMF (after fit): " << endl;
    cout << "Saturation density:\t " << ld.n0 << endl;
      
    cout << "Mass excess:\t\t " << ld.mass_excess(82,126) << endl;
    cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
    cout <<  "Total mass:\t\t " << ld.total_mass(82,126) << endl;
      
  }

  cout << endl;
  cout << "-------------------------------------------------\n" << endl;

  // APR

  ld.set_eos_had_temp_base(apr);
  ld.n0=0.16;
  ld.n1=0.0;
  cout << "Lead from APR: " << endl;
  cout << "Mass excess:\t\t " <<  ld.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
  cout << "Total mass:\t\t " << ld.total_mass(82,126) << endl;

  t.test_rel(au.total_mass(82,126),ld.total_mass(82,126),1.0e-3,
	     "ame vs ld2");

  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  // fit APR

#ifdef O2SCL_NEVER_DEFINED

  nucmass_fit mf2;
  mf2.even_even=false;

  ld.n0=0.16;
  ld.n1=-0.05;
  ld.surften=1.1;
  ld.coul_coeff=1.0;
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << endl;
  
  mf2.eval(ld,qual);
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << endl;
  cout << qual << endl;

  cout << "Fit APR: " << endl;
  
  mf2.fit(ld,qual);
  cout << ld.n1 << " " << ld.n0 << " " << ld.surften << " " 
       << ld.coul_coeff << endl;
  cout << qual << endl;
  t.test_gen(qual<4.0,"APR Fit");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;
  
  // APR after fit
  
  cout << "APR (after fit): " << endl;
  cout << "Saturation density:\t " << ld.n0 << endl;
  
  cout << "Mass excess:\t\t " << ld.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ld.binding_energy(82,126)/208.0 << endl;
  cout <<  "Total mass:\t\t " << ld.total_mass(82,126) << endl;
  
  cout << "Densities:\t\t " << ld.nn << " " << ld.np << endl;
  cout << "Implied radii: " << cbrt(3*126/4.0/pi/ld.nn) << " "
       << cbrt(3*82/4.0/pi/ld.np) << endl;

  cout << endl;

#endif

  cout << "-------------------------------------------------\n" << endl;

  cout << "Compare nucmass_ldrop vs. nucmass_ldrop_skin:\n" << endl;

  nucmass_ldrop_skin ldf;

  ld.n1=-0.11;
  ld.n0=0.16;
  ldf.n1=-0.11;
  ldf.n0=0.16;
  ld.set_eos_had_temp_base(apr);
  ldf.set_eos_had_temp_base(apr);

  ldf.doi=1.0;
  ld.surften=1.1;
  ldf.surften=1.1;

  ld.coul_coeff=1.0;
  ldf.coul_coeff=1.0;

  ldf.ss=0.0;
  ldf.full_surface=false;

  cout << ld.mass_excess(82,126) << endl;
  cout << ld.binding_energy(82,126)/208.0 << endl;
  cout << ld.bulk << " " << ld.surf << " " << ld.coul << endl;
  cout << (ld.bulk+ld.surf+ld.coul) << endl;
  cout << ld.nn << " " << ld.np << endl;
  cout << endl;

  cout << ldf.mass_excess(82,126) << endl;
  cout << ldf.binding_energy(82,126)/208.0 << endl;
  cout << ldf.bulk << " " << ldf.surf << " " << ldf.coul << endl;
  cout << (ldf.bulk+ldf.surf+ldf.coul) << endl;
  cout << ldf.nn << " " << ldf.np << endl;
  cout << endl;

  t.test_rel(ld.mass_excess(82,126),ldf.mass_excess(82,126),1.0e-6,
	     "mass excesses");
  t.test_rel(ld.binding_energy(82,126),ldf.binding_energy(82,126),1.0e-6,
	     "binding energies");
  t.test_rel(ld.bulk,ldf.bulk,1.0e-6,
	     "bulk energies");
  t.test_rel(ld.surf,ldf.surf,1.0e-6,
	     "surface energies");
  t.test_rel(ld.nn,ldf.nn,1.0e-6,"neutron density");
  t.test_rel(ld.np,ldf.np,1.0e-6,"proton density");

  ldf.full_surface=true;
  ldf.ss=0.5;

  cout << "-------------------------------------------------\n" << endl;

  // With neutron skin
  
  ldf.set_eos_had_temp_base(apr);

  cout << "APR (with skin): " << endl;
  cout  << "Saturation density:\t " << ldf.n0 << endl;
  
  cout << "Mass excess:\t\t " <<  ldf.mass_excess(82,126) << endl;
  cout << "Binding energy:\t\t " << ldf.binding_energy(82,126)/208.0 << endl;
  d1=ldf.binding_energy(82,126)/208.0;
  cout << "Total mass:\t\t " <<  ldf.total_mass(82,126) << endl;
  
  t.test_rel(au.total_mass(82,126),ldf.total_mass(82,126),2.0e-2,
	     "ame vs ldf");
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;
  
  cout << "Fit with skin: " << endl;

  nucmass_fit mf;

#ifdef O2SCL_NEVER_DEFINED

  ldf.n0=0.1812;
  ldf.n1=-0.1249;
  ldf.doi=0.88284;
  ldf.surften=1.189;
  ldf.ss=1.66703;
  ldf.coul_coeff=0.88263;

  mf.eval(ldf,qual);

  cout << ldf.n0 << " " << ldf.n1 << " " << ldf.doi << " "
       << ldf.surften << " " << ldf.ss << " " << ldf.coul_coeff << endl;
  cout << qual << endl;
  t.test_gen(qual<5.0,"w/ Skin fit quality");

#endif

  if (false) {
    
    //ldf.new_skin_mode=true;
    mf.eval(ldf,qual);

    //mf.def_mmin.verbose=1;
    mf.def_mmin.err_nonconv=false;
    //mf.fit(ldf,qual);

    cout << ldf.n0 << " " << ldf.n1 << " " << ldf.doi << " "
	 << ldf.surften << " " << ldf.ss << " " << ldf.coul_coeff << endl;
    cout << qual << endl;

    mf.set_mmin(gms);

    mf.fit(ldf,qual);
    
    cout << ldf.n0 << " " << ldf.n1 << " " << ldf.doi << " "
	 << ldf.surften << " " << ldf.ss << " " << ldf.coul_coeff << endl;
    cout << qual << endl;
  }

#ifdef O2SCL_NEVER_DEFINED
  if (false) {

    mf.def_mmin.verbose=1;
    mf.fit(ldf,qual);
    
    cout << ldf.n0 << " " << ldf.n1 << " " << ldf.doi << " "
	 << ldf.surften << " " << ldf.ss << " " << ldf.coul_coeff << endl;
    cout << qual << endl;
    exit(-1);
  }
#endif

  cout << "-------------------------------------------------\n" << endl;

  cout << "Finite temperature: " << endl;

  ldf.set_eos_had_temp_base(apr);
  ldf.n0=0.184;
  ldf.n1=-0.05;
  ldf.doi=0.9;
  ldf.surften=1.17;
  ldf.ss=1.85;
  ldf.coul_coeff=0.89;

  cout << "ldf: " << ldf.binding_energy(82,126)/208.0 << endl;
  cout << endl;

  cout << "ldf: " << ldf.binding_energy(82,126)/208.0 << endl;
  cout << endl;

  cout << "ldf,drip: " 
       << ldf.drip_binding_energy_d(82,126,0.02,0.02,0.0,0.0)/208.0 << endl;
  cout << endl;

  cout << "ldf,drip,dim=3: " 
       << ldf.drip_binding_energy_d(82,126,0.02,0.02,0,0.0)/208.0 
       << endl;
  cout << "nn,np: " << ldf.nn << " " << ldf.np << endl;
  cout << "bulk,surf,coul: " 
       << ldf.bulk << " " << ldf.surf << " " << ldf.coul << endl;
  cout << endl;

  cout << "ldf,drip,dim,T=0.01: " 
       << ldf.drip_binding_energy_d(82,126,0.02,0.02,0,0.01)/208.0 
       << endl;
  cout << "nn,np: " << ldf.nn << " " << ldf.np << endl;
  cout << "bulk,surf,coul: " 
       << ldf.bulk << " " << ldf.surf << " " << ldf.coul << endl;
  cout << endl;

  cout << "-------------------------------------------------\n" << endl;

  t.report();
  return 0;
}

