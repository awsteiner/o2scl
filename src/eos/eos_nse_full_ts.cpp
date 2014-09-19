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

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  eos_nse_full nse;
  dense_matter dm;
  int ret;

  // Test nucmass_densmat

  nucmass_densmat &nd=nse.nuc_dens;
  double t1, t2, t3, t4;
  nd.test_derivatives(1.0e-6,t1,t2,t3,t4);
  t.test_rel(t1,0.0,1.0e-6,"dEdnp");
  t.test_rel(t3,0.0,3.0e-6,"dEdnneg");
  
  // Create a distribution of three nuclei
  o2scl::nucleus nuc;
  o2scl::nucmass_ame_exp &ame=nse.nuc_dens.ame;

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

  // Compare analytical and numerical values for eta_i with
  // inc_prot_coul = false

  double eps=1.0e-6, fr1, fr2, eta_n, eta_p, eta_nuc[3];

  nse.inc_prot_coul=false;

  nse.calc_density_noneq(dm,0);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.n.n*=(1.0+eps);
  nse.calc_density_noneq(dm,0);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.n.n/=(1.0+eps);
  eta_n=(fr2-fr1)/(eps*dm.n.n);

  nse.calc_density_noneq(dm,0);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.p.n*=(1.0+eps);
  nse.calc_density_noneq(dm,0);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.p.n/=(1.0+eps);
  eta_p=(fr2-fr1)/(eps*dm.p.n);

  for(size_t i=0;i<3;i++) {
    nse.calc_density_noneq(dm,0);
    fr1=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n*=(1.0+eps);
    nse.calc_density_noneq(dm,0);
    fr2=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n/=(1.0+eps);
    eta_nuc[i]=(fr2-fr1)/(eps*dm.dist[i].n);
  }
  
  ret=nse.calc_density_noneq(dm,0);
  t.test_gen(ret==0,"ret 1");
  t.test_rel(dm.eta_n,eta_n,1.0e-6,"eta_n");
  t.test_rel(dm.eta_p,eta_p,1.0e-6,"eta_p");
  t.test_rel(dm.eta_nuc[0],eta_nuc[0],1.0e-6,"eta_nuc[0]");
  t.test_rel(dm.eta_nuc[1],eta_nuc[1],1.0e-6,"eta_nuc[1]");
  t.test_rel(dm.eta_nuc[2],eta_nuc[2],1.0e-6,"eta_nuc[2]");

  // Compare analytical and numerical values for eta_i with
  // inc_prot_coul = true (the default)

  nse.inc_prot_coul=true;

  nse.calc_density_noneq(dm,0);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.n.n*=(1.0+eps);
  nse.calc_density_noneq(dm,0);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.n.n/=(1.0+eps);
  eta_n=(fr2-fr1)/(eps*dm.n.n);

  nse.calc_density_noneq(dm,0);
  fr1=dm.th.ed-dm.T*dm.th.en;
  dm.p.n*=(1.0+eps);
  nse.calc_density_noneq(dm,0);
  fr2=dm.th.ed-dm.T*dm.th.en;
  dm.p.n/=(1.0+eps);
  eta_p=(fr2-fr1)/(eps*dm.p.n);

  for(size_t i=0;i<3;i++) {
    nse.calc_density_noneq(dm,0);
    fr1=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n*=(1.0+eps);
    nse.calc_density_noneq(dm,0);
    fr2=dm.th.ed-dm.T*dm.th.en;
    dm.dist[i].n/=(1.0+eps);
    eta_nuc[i]=(fr2-fr1)/(eps*dm.dist[i].n);
  }
  
  ret=nse.calc_density_noneq(dm,0);
  t.test_gen(ret==0,"ret 1");
  t.test_rel(dm.eta_n,eta_n,1.0e-6,"eta_n");
  t.test_rel(dm.eta_p,eta_p,1.0e-6,"eta_p");
  t.test_rel(dm.eta_nuc[0],eta_nuc[0],1.0e-6,"eta_nuc[0]");
  t.test_rel(dm.eta_nuc[1],eta_nuc[1],1.0e-6,"eta_nuc[1]");
  t.test_rel(dm.eta_nuc[2],eta_nuc[2],1.0e-6,"eta_nuc[2]");

  // Now try to minimize the free energy for this composition

  dm.nB=5.2e-2;
  dm.Ye=0.4943;

  nse.calc_density_fixcomp(dm,0);

  fr1=dm.th.ed-dm.th.en*dm.T;

  // Double check the free energy

  ret=nse.calc_density_noneq(dm,1);
  t.test_gen(ret==0,"ret 2");
  fr2=dm.th.ed-dm.th.en*dm.T;
  t.test_rel(fr1,fr2,1.0e-6,"free energies");

  //dm.T=4.0/hc_mev_fm;
  //nse.calc_density(dm);
  
  t.report();
  return 0;
}


