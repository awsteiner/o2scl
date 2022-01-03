/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/nucmass_frdm.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  nucmass_mnmsk mo12;
  nucmass_mnmsk mo;
  nucmass_frdm frdm;
  nucmass_ame au;
  nucmass_semi_empirical sm;

  // Load experimental and theoretical masses
  o2scl_hdf::ame_load_ext(au,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2");
  o2scl_hdf::mnmsk_load(mo,"mnmsk97",
			"../../data/o2scl/nucmass/mnmsk.o2");
  o2scl_hdf::mnmsk_load(mo12,"msis16",
			"../../data/o2scl/nucmass/msis16.o2");

  // Show that the nucmass_frdm gives reasonable (but not great)
  // values for the binding energy of Lead 208
  cout << "AME2003   : ";
  cout << au.mass_excess(82,126) << " ";
  cout << au.binding_energy(82,126) << " ";
  cout << au.total_mass(82,126) << endl;
  cout << "Simple    : ";
  cout << sm.mass_excess(82,126) << " ";
  cout << sm.binding_energy(82,126) << " ";
  cout << sm.total_mass(82,126) << endl;
  cout << "FRDM      : ";
  cout << frdm.mass_excess(82,126) << " ";
  cout << frdm.binding_energy(82,126) << " ";
  cout << frdm.total_mass(82,126) << endl;
  cout << "MNMSK '97 : ";
  cout << mo.mass_excess(82,126) << " ";
  cout << mo.binding_energy(82,126) << " ";
  cout << mo.total_mass(82,126) << endl;
  cout << "MSIS '16  : ";
  cout << mo12.mass_excess(82,126) << " ";
  cout << mo12.binding_energy(82,126) << " ";
  cout << mo12.total_mass(82,126) << endl;

  // Explicitly add the microscopic part of the binding energy for 
  // Lead 208 from the table in Moller et al. (1995) 
  // and compare with the theoretical value quoted in that table
  t.test_rel(frdm.mass_excess(82,126)-12.84,-21.15,5.0e-4,"Lead 208");

  // ----------------------------------------------------------------
  // Compare nucmass_frdm with the macroscopic parts from
  // nucmass_mnmsk table and show that they're almost the same
  
  // Quality of fit
  double qual=0.0;
  
  // Total number of nuclei
  size_t nnuc=0;

  // Create the distribution
  vector<nucleus> fd;
  nucdist_set(fd,mo);

  // Old code to make mnmks97 table
  //ofstream fout;
  //fout.open("mnmsk97.dat");
  for(vector<nucleus>::iterator ndi=fd.begin();ndi!=fd.end();ndi++) {
    nucmass_mnmsk::entry mme=mo.get_ZN(ndi->Z,ndi->N);
    /*
    if (true) {
      fout << ((int)(mme.N+1.0e-12)) << " ";
      fout << ((int)(mme.Z+1.0e-12)) << " ";
      fout << ((int)(mme.A+1.0e-12)) << " ";
      fout << mme.eps2 << " ";
      fout << mme.eps3 << " ";
      fout << mme.eps4 << " ";
      fout << mme.eps6 << " ";
      fout << mme.eps6sym << " ";
      fout << mme.beta2 << " ";
      fout << mme.beta3 << " ";
      fout << mme.beta4 << " ";
      fout << mme.beta6 << " ";
      fout << mme.Emic << " ";
      fout << mme.Mth << " ";
      fout << mme.Mexp << " ";
      fout << mme.sigmaexp << " ";
      fout << mme.EmicFL << " ";
      fout << mme.MthFL << " ";
      string stmp(&mme.spinp[0]);
      fout << spin_to_double(stmp) << " ";
      string stmp2(&mme.spinn[0]);
      fout << spin_to_double(stmp2) << " ";
      fout << mme.gapp << " ";
      fout << mme.gapn << " ";
      fout << mme.be << " ";
      fout << mme.S1n << " ";
      fout << mme.S2n << " ";
      fout << mme.PA << " ";
      fout << mme.PAm1 << " ";
      fout << mme.PAm2 << " ";
      fout << mme.Qbeta << " ";
      fout << mme.Tbeta << " ";
      fout << mme.S1p << " ";
      fout << mme.S2p << " ";
      fout << mme.Qalpha << " ";
      fout << mme.Talpha << endl;
    }
    */
    if (ndi->N>=8 && ndi->Z>=8) {
      qual+=pow(frdm.mass_excess(ndi->Z,ndi->N)-(mme.Mth-mme.Emic),2.0);
      nnuc++;
    }
  }
  //fout.close();
  qual=sqrt(qual/nnuc);
  t.test_abs(qual,0.0,5.0e-3,"Compareare with Mth-Emic from full FRDM");
  
  cout << "Compared " << nnuc << " nuclei from nucmass_frdm and "
       << "nucmass_mnmsk. Fit quality factor: " << qual << " ." << endl;

  // ----------------------------------------------------------------
  // Fit the macroscopic part of FRDM in the nucmass_frdm class
  // to experiment
  
  nucmass_fit mf;
  nucdist_set(mf.dist,au);
  mf.eval(frdm,qual);
  cout << "Before fit: " << qual << endl;
  t.test_rel(qual,3.1860,1.0e-2,"FRDM pre-fit.");
  mf.def_mmin.ntrial*=10;
  mf.fit(frdm,qual);
  cout << "After fit: " << qual << endl;
  t.test_rel(qual,2.3836,1.0e-2,"FRDM post-fit.");

  nucmass_patch np;
  np.load();
  for(int Z=10;Z<100;Z++) {
    int N=110-Z;
    cout << Z << " " << N << " " << np.mass_excess(Z,N) << endl;
  }
  
  t.report();

  return 0;
}
