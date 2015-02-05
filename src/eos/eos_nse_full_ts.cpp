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
#include <o2scl/eos_nse_full.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  eos_nse_full nse;
  nse.verbose=0;

  dense_matter dm;
  int ret;

  // Nuclear masses
  nucmass_ame_exp ame;
  o2scl_hdf::ame_load(ame,"12");

  // Load Skyrme EOS
  eos_had_skyrme sk;
  o2scl_hdf::skyrme_load(sk,"SLy4");
  nse.set_eos(sk);

  // Test nucmass_densmat

  nucmass_densmat &nd=nse.nuc_dens;
  nd.set_mass(ame);
  double t1, t2, t3, t4;
  nd.test_derivatives(1.0e-6,t1,t2,t3,t4);
  t.test_rel(t1,0.0,1.0e-6,"dEdnp");
  t.test_rel(t3,0.0,3.0e-6,"dEdnneg");
  
  // Create a distribution of three nuclei

  o2scl::nucleus nuc;

  dm.n.n=0.01;
  dm.p.n=0.01;
  dm.dist.clear();

  ame.get_nucleus(26,26,nuc);
  nuc.g=2.0;
  nuc.n=0.01/50.0;
  dm.dist.push_back(nuc);

  ame.get_nucleus(26,27,nuc);
  nuc.g=2.0;
  nuc.n=0.01/50.0;
  dm.dist.push_back(nuc);

  ame.get_nucleus(26,28,nuc);
  nuc.g=2.0;
  nuc.n=0.01/50.0;
  dm.dist.push_back(nuc);

  // Set the temperature
  dm.T=1.0/hc_mev_fm;

  // Compare analytical and numerical values from calc_density_noneq()
  // for eta_i with inc_prot_coul = false

  double eps=1.0e-6, fr1, fr2, eta_n, eta_p, eta_nuc[3];

  nse.inc_prot_coul=false;

  nse.calc_density_noneq(dm);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.n.n*=(1.0+eps);
  nse.calc_density_noneq(dm);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.n.n/=(1.0+eps);
  eta_n=(fr2-fr1)/(eps*dm.n.n);

  nse.calc_density_noneq(dm);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.p.n*=(1.0+eps);
  nse.calc_density_noneq(dm);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.p.n/=(1.0+eps);
  eta_p=(fr2-fr1)/(eps*dm.p.n);

  for(size_t i=0;i<3;i++) {
    nse.calc_density_noneq(dm);
    fr1=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n*=(1.0+eps);
    nse.calc_density_noneq(dm);
    fr2=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n/=(1.0+eps);
    eta_nuc[i]=(fr2-fr1)/(eps*dm.dist[i].n);
  }
  
  ret=nse.calc_density_noneq(dm);
  t.test_gen(ret==0,"ret 1");
  t.test_rel(dm.eta_n,eta_n,1.0e-6,"eta_n");
  t.test_rel(dm.eta_p,eta_p,1.0e-6,"eta_p");
  t.test_rel(dm.eta_nuc[0],eta_nuc[0],1.0e-6,"eta_nuc[0]");
  t.test_rel(dm.eta_nuc[1],eta_nuc[1],1.0e-6,"eta_nuc[1]");
  t.test_rel(dm.eta_nuc[2],eta_nuc[2],1.0e-6,"eta_nuc[2]");

  // Compare analytical and numerical values from calc_density_noneq()
  // for eta_i with inc_prot_coul = true (the default)

  nse.inc_prot_coul=true;

  nse.calc_density_noneq(dm);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.n.n*=(1.0+eps);
  nse.calc_density_noneq(dm);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.n.n/=(1.0+eps);
  eta_n=(fr2-fr1)/(eps*dm.n.n);

  nse.calc_density_noneq(dm);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.p.n*=(1.0+eps);
  nse.calc_density_noneq(dm);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.p.n/=(1.0+eps);
  eta_p=(fr2-fr1)/(eps*dm.p.n);

  for(size_t i=0;i<3;i++) {
    nse.calc_density_noneq(dm);
    fr1=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n*=(1.0+eps);
    nse.calc_density_noneq(dm);
    fr2=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n/=(1.0+eps);
    eta_nuc[i]=(fr2-fr1)/(eps*dm.dist[i].n);
  }
  
  ret=nse.calc_density_noneq(dm);
  t.test_gen(ret==0,"ret 1");
  t.test_rel(dm.eta_n,eta_n,1.0e-6,"eta_n");
  t.test_rel(dm.eta_p,eta_p,1.0e-6,"eta_p");
  t.test_rel(dm.eta_nuc[0],eta_nuc[0],1.0e-6,"eta_nuc[0]");
  t.test_rel(dm.eta_nuc[1],eta_nuc[1],1.0e-6,"eta_nuc[1]");
  t.test_rel(dm.eta_nuc[2],eta_nuc[2],1.0e-6,"eta_nuc[2]");
  
  // Now try to minimize the free energy for this composition
  // using calc_density_by_min()

  dm.nB=5.2e-2;
  dm.Ye=0.4943;

  ret=nse.calc_density_by_min(dm);
  t.test_gen(ret==0,"ret 2");
  t.test_rel(dm.baryon_density(),dm.nB,1.0e-6,"baryon density");
  t.test_rel(dm.electron_fraction(),dm.Ye,1.0e-6,"electron fraction");

  fr1=dm.th.ed-dm.th.en*dm.T;

  // Double check the free energy, baryon density, and electron
  // fraction with that returned by calc_density_noneq() 

  ret=nse.calc_density_noneq(dm);
  t.test_gen(ret==0,"ret 3");
  fr2=dm.th.ed-dm.th.en*dm.T;
  t.test_rel(fr1,fr2,1.0e-10,"free energies");

  // Minimize the free energy when inc_prot_coul is false

  nse.inc_prot_coul=false;
  // (This requires very accurate minimization to work correctly)
  nse.def_mmin.tol_rel/=1.0e4;
  nse.def_mmin.tol_abs/=1.0e4;

  ret=nse.calc_density_by_min(dm);
  t.test_gen(ret==0,"ret 4");

  for(size_t i=0;i<3;i++) {
    t.test_rel(dm.dist[i].be+dm.dist[i].mu,
	       dm.dist[i].Z*dm.p.mu+dm.dist[i].N*dm.n.mu,1.0e-2,"NSE 1");
    t.test_rel(dm.eta_n*dm.dist[i].N+dm.eta_p*dm.dist[i].Z,
	       dm.eta_nuc[i],1.0e-3,"NSE 2");
  }

#ifdef O2SCL_NEVER_DEFINED
  
  nse.calc_density_noneq(dm);
  nse.calc_density_fixnp(dm);

  // But even with this accurate minimization, it doesn't quite
  // get the right densities
  cout << dm.nB << " " << dm.baryon_density() << endl;
  cout << dm.Ye << " " << dm.electron_fraction() << endl;

  if (true) {
    // Add a full distribution and test mup_for_Ye
    dm.dist.clear();
    nucdist_set(dm.dist,ame,"(N+Z)>1");
    dm.nB=1.0e-10;
    dm.Ye=0.45;
    dm.T=0.1/hc_mev_fm;
    double mup_high=-5.0/hc_mev_fm;
    double mup_low=-5.5/hc_mev_fm;
    double mun_low=-12.0/hc_mev_fm;
    double mun_high=-11.5/hc_mev_fm;
    nse.verbose=1;
    nse.mup_for_Ye(-6.0/hc_mev_fm,mun_low,mun_high,dm);
    cout << dm.n.mu*hc_mev_fm << endl;
    cout << dm.baryon_density() << endl;
    cout << dm.electron_fraction() << endl;
    exit(-1);

    nse.bracket_mu_solve(mun_low,mun_high,mup_low,mup_high,dm);
  }
  
#endif
  
  t.report();
  return 0;
}


