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
#include <o2scl/eos_nse.h>
#include <o2scl/nucmass_frdm.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  double nb, Ye, mun, mup;
  
  double T=5.0/hc_mev_fm;
  thermo th;

  // Load nuclear masses
  nucmass_mnmsk mm;
  o2scl_hdf::mnmsk_load(mm);

  // Set the distribution of nuclei to use
  vector<nucleus> ad;
  nucdist_set(ad,mm,"Z>=24 & Z<=32 & N>=55 & N<=58");
  cout << ad.size() << endl;
  for(size_t i=0;i<ad.size();i++) {
    ad[i].g=1.0;
  }
  t.test_gen(ad.size()==36,"distribution size");

  eos_nse ne;

  nb=0.03;
  Ye=0.36;

  // Initial guess
  mun=-0.04;
  mup=-0.04;
  
  int ret=ne.calc_density(nb,Ye,T,mun,mup,th,ad);
  cout << ret << " " << mun << " " << mup << " " << nb << " " << Ye << endl;
  cout << endl;

  // Double check that the density and electron fraction are properly
  // reproduced
  double nbnew=0.0;
  double Yenew=0.0;
  for(size_t i=0;i<ad.size();i++) {
    cout << ad[i].Z << " " << mm.Ztoel(ad[i].Z) << " " 
	 << ad[i].A << " " << ad[i].n << " " 
	 << ((double)ad[i].Z)/((double)ad[i].A) << endl;
    nbnew+=ad[i].n*ad[i].A;
    Yenew+=ad[i].n*ad[i].Z;
  }
  Yenew/=nbnew;
  cout << endl;

  t.test_rel(nbnew,0.03,1.0e-6,"nb match.");
  t.test_rel(Yenew,0.36,1.0e-6,"Ye match.");

  for(;T>0.01/hc_mev_fm;T/=1.1) {
    cout << "Here." << nb << " " << Ye << endl;
    ret=ne.calc_density(nb,Ye,T,mun,mup,th,ad);
    cout << ret << " " << T << " "
	 << mun << " " << mup << " ";

    for(size_t i=0;i<ad.size();i++) {
      nbnew+=ad[i].n*ad[i].A;
      Yenew+=ad[i].n*ad[i].Z;
    }
    Yenew/=nbnew;

    cout << nbnew << " " << Yenew << endl;
    exit(-1);
  }
  cout << endl;

  t.test_rel(nbnew,0.03,1.0e-6,"nb match.");
  t.test_rel(Yenew,0.36,1.0e-6,"Ye match.");
  
  t.report();
  return 0;
}


