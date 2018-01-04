/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

  double nB, Ye, mun, mup, nn, np;
  eos_nse en;
  
  double T=5.0/hc_mev_fm;
  thermo th;

  // Load nuclear masses
  nucmass_mnmsk mm;
  o2scl_hdf::mnmsk_load(mm,"../../data/o2scl/nucmass/");

  // Set the distribution of nuclei to use
  vector<nucleus> ad;
  nucdist_set(ad,mm,"Z>=24 && Z<=32 && N>=55 && N<=58");
  for(size_t i=0;i<ad.size();i++) {
    ad[i].g=1.0;
  }
  t.test_gen(ad.size()==36,"distribution size");

  // ---------------------------------------------------------
  // Test make_guess(), show that it can solve for the
  // densities within a factor of two

  en.verbose=1;
  en.make_guess_iters=60;

  double fac=2.0;

  // Test make_guess() with very large chemical potentials
  mun=1.0;
  mup=1.0;
  en.make_guess(mun,mup,T,th,ad,1.0/fac,fac,1.0/fac,fac);

  // Test make_guess() with very small chemical potentials
  mun=-1.0;
  mup=-1.0;
  en.make_guess(mun,mup,T,th,ad,1.0/fac,fac,1.0/fac,fac);

  // Test make_guess() with large mun and small mup
  mun=1.0;
  mup=-1.0;
  en.make_guess(mun,mup,T,th,ad,1.0/fac,fac,1.0/fac,fac);

  // Test make_guess() with large mup and small mun
  mun=-1.0;
  mup=1.0;
  en.make_guess(mun,mup,T,th,ad,1.0/fac,fac,1.0/fac,fac);

  // ---------------------------------------------------------
  // Now test calc_density()
  
  nn=0.64*0.03;
  np=0.36*0.03;

  mun=1.0;
  mup=1.0;
  
  int ret=en.calc_density(nn,np,T,mun,mup,th,ad);

  en.verbose=0;

  // Double check that the density and electron fraction are properly
  // reproduced
  double nBnew=0.0;
  double Yenew=0.0;
  for(size_t i=0;i<ad.size();i++) {
    cout << ad[i].Z << " " << mm.Ztoel(ad[i].Z) << " " 
	 << ad[i].A << " " << ad[i].n << " " 
	 << ((double)ad[i].Z)/((double)ad[i].A) << endl;
    nBnew+=ad[i].n*ad[i].A;
    Yenew+=ad[i].n*ad[i].Z;
  }
  Yenew/=nBnew;
  cout << endl;

  t.test_rel(nBnew,0.03,1.0e-6,"nB match.");
  t.test_rel(Yenew,0.36,1.0e-6,"Ye match.");

  // Now proceed to lower temperatures
  for(;T>0.01/hc_mev_fm;T/=1.1) {
    ret=en.calc_density(nn,np,T,mun,mup,th,ad);
    cout << ret << " " << T << " "
	 << mun << " " << mup << " ";

    nBnew=0.0;
    Yenew=0.0;
    for(size_t i=0;i<ad.size();i++) {
      nBnew+=ad[i].n*ad[i].A;
      Yenew+=ad[i].n*ad[i].Z;
    }
    Yenew/=nBnew;

    cout << nBnew << " " << Yenew << endl;
  }
  cout << endl;

  // Output the distribution at the lowest temperature
  for(size_t i=0;i<ad.size();i++) {
    cout << ad[i].Z << " " << mm.Ztoel(ad[i].Z) << " " 
	 << ad[i].A << " " << ad[i].n << endl;
  }
  cout << endl;

  // Double check that the density and electron fraction are properly
  // reproduced
  t.test_rel(nBnew,0.03,1.0e-6,"nB match.");
  t.test_rel(Yenew,0.36,1.0e-6,"Ye match.");

  // ---------------------------------------------------------
  // Test with a more complete distribution at large and
  // small Ye

  o2scl::nucdist_set(ad,mm,"1");
  double min=1.0, max=0.0;
  for(size_t i=0;i<ad.size();i++) {
    if (((double)ad[i].Z)/((double)ad[i].A)<min) {
      min=((double)ad[i].Z)/((double)ad[i].A);
    }
    if (((double)ad[i].Z)/((double)ad[i].A)>max) {
      max=((double)ad[i].Z)/((double)ad[i].A);
    }
  }
  cout << "min Ye, max Ye: " << min << " " << max << endl;

  mun=1.0;
  mup=1.0;
  Ye=0.24;
  nn=1.0e-12*(1.0-Ye);
  np=1.0e-12*Ye;
  T=100.0/hc_mev_fm;
  double fact=2.0;
  int ret1=0, ret2, ret3;
  while (T>1.0/hc_mev_fm) {
    ret3=en.calc_density(nn,np,T,mun,mup,th,ad);
    cout << ret1 << " " << ret3 << " " << T << " " << mun << " "
	 << mup << endl;
    T/=1.5;
  }
  cout << endl;
  
  mun=1.0;
  mup=1.0;
  Ye=0.67;
  nn=1.0e-12*(1.0-Ye);
  np=1.0e-12*Ye;
  T=100.0/hc_mev_fm;
  while (T>1.0/hc_mev_fm) {
    ret3=en.calc_density(nn,np,T,mun,mup,th,ad);
    cout << ret1 << " " << ret3 << " " << T << " " << mun << " "
	 << mup << endl;
    T/=1.1;
  }
  cout << endl;
  
  t.report();
  return 0;
}


