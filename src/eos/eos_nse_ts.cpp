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
#include <o2scl/nucdist_arb.h>
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

  mnmsk_mass mm;
  o2scl_hdf::mnmsk_load(mm);
  nucdist_arb ad;
  ad.set_dist(mm,"Z>=24 & Z<=32 & N>=55 & N<=58");
  eos_nse ne;

  nb=0.028;
  Ye=0.364;

  // Initial guess
  mun=-0.064;
  mup=-0.012;
  
  int ret=ne.calc_density(nb,Ye,T,mun,mup,th,ad);
  cout << ret << " " << mun << " " << mup << " " << nb << " " << Ye << endl;

  // Double check that the density and electron fraction are properly
  // reproduced
  double nbnew=0.0;
  double Yenew=0.0;
  for(nucdist::iterator ndi=ad.begin();ndi!=ad.end();ndi++) {
    cout << ndi->Z << " " << mm.Ztoel(ndi->Z) << " " 
	 << ndi->A << " " << ndi->n << endl;
    nbnew+=ndi->n*ndi->A;
    Yenew+=ndi->n*ndi->Z;
  }
  Yenew/=nbnew;
  t.test_rel(nbnew,0.028,1.0e-3,"nb match.");
  t.test_rel(Yenew,0.364,1.0e-3,"Ye match.");
  
  t.report();
  return 0;
}


