/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2026, Andrew W. Steiner

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

  ───────────────────────────────────────────────────────────────────
*/
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_fit_iso.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  nucmass_ame ame;
  ame.load("20");

  // A plain, single, global fit of the Duflo-Zuker
  // formula, for comparison.
  nucmass_dz_fit dz_global;
  nucmass_fit mf_global;
  nucdist_set(mf_global.dist,ame);
  //mf_global.def_mmin.err_nonconv=false;
  mf_global.def_mmin.ntrial*=6;
  double res_global;
  mf_global.fit(dz_global,res_global);
  cout << "Global DZ-33 fit, rms=" << res_global << " MeV." << endl;
  t.test_gen(res_global<2.0,"Global fit is reasonably good.");

  // Now the per-chain ensemble, restricted to a small Z range so
  // this test runs quickly
  nucmass_fit_iso nfi;

  nfi.mf.def_mmin.ntrial=60000;
  size_t n_fit=nfi.fit(dz_global,ame,48,52);
  cout << "Number of chains fit: " << n_fit << endl;
  t.test_gen(n_fit==5,"All 5 requested chains were fit.");
  t.test_gen(nfi.get_nentries()==5,"get_nentries() matches.");

  // Sn (Z=50) should be included for at least one of its
  // well-known stable isotopes, and the fitted mass excess should
  // be close to the AME value.
  t.test_gen(nfi.is_included(50,70),"Sn-120 is included.");
  double me_ame=ame.mass_excess(50,70);
  double me_nfi=nfi.mass_excess(50,70);
  cout << "Sn-120: AME=" << me_ame << " nfi=" << me_nfi << endl;
  t.test_abs(me_nfi,me_ame,1.0,"Sn-120 mass excess close to AME.");

  // A Z with no fitted chain should not be included.
  t.test_gen(nfi.is_included(60,70)==false,
             "Z=60 (not fit) correctly excluded.");

  // Two independently cloned chains shouldn't be the same object
  // and (since Sn and Cd are fit over different, only partially
  // overlapping, neighborhoods) shouldn't give identical
  // predictions far from their shared calibration region.
  double me_sn=nfi.mass_excess(50,90);
  double me_cd=nfi.mass_excess(48,90);
  cout << "N=90: Sn(Z=50)=" << me_sn << " Cd(Z=48)=" << me_cd << endl;
  t.test_gen(me_sn!=me_cd,"Independently-fit chains differ.");

  // HDF5 round-trip: write nfi to a temporary file, read it back
  // into a fresh nucmass_fit_iso (with a fresh, unfit prototype,
  // exactly as a real user would after starting a new program), and
  // check that the reloaded object reproduces the same chain count
  // and the same predictions.
  {
    hdf_file hf;
    hf.open_or_create("nucmass_fit_iso_ts.o2");
    nfi.hdf_output(hf,"nfi_test");
    hf.close();

    nucmass_dz_fit dz_proto;
    nucmass_fit_iso nfi2;
    hdf_file hf2;
    hf2.open("nucmass_fit_iso_ts.o2");
    nfi2.hdf_input(hf2,dz_proto,"nfi_test");
    hf2.close();

    t.test_gen(nfi2.get_nentries()==nfi.get_nentries(),
               "HDF5 round-trip: chain count matches.");
    t.test_gen(nfi2.x==nfi.x,"HDF5 round-trip: x matches.");
    t.test_gen(nfi2.is_included(50,70),
               "HDF5 round-trip: Sn-120 still included.");
    t.test_rel(nfi2.mass_excess(50,70),nfi.mass_excess(50,70),1.0e-10,
               "HDF5 round-trip: Sn-120 mass excess matches exactly.");
    t.test_rel(nfi2.mass_excess(48,90),nfi.mass_excess(48,90),1.0e-10,
               "HDF5 round-trip: Cd chain mass excess matches exactly.");
    t.test_gen(nfi2.is_included(60,70)==false,
               "HDF5 round-trip: unfit Z=60 still excluded.");

    // Auto-discovery: an empty name should find the (only) object
    // of type nucmass_fit_iso in the file.
    nucmass_fit_iso nfi3;
    hdf_file hf3;
    hf3.open("nucmass_fit_iso_ts.o2");
    // hdf_input()'s name parameter is by value (matching the
    // convention in e.g. interpm_libtorch), so the discovered name
    // isn't propagated back to the caller; hdf_input_n() must be
    // used directly to recover it.
    std::string found_name;
    nfi3.hdf_input_n(hf3,dz_proto,found_name);
    hf3.close();
    t.test_gen(found_name=="nfi_test",
               "HDF5 auto-discovery found the right group name.");
    t.test_gen(nfi3.get_nentries()==nfi.get_nentries(),
               "HDF5 auto-discovery: chain count matches.");
  }

  // clear() removes everything.
  nfi.clear();
  t.test_gen(nfi.get_nentries()==0,"clear() empties the table.");
  t.test_gen(nfi.is_included(50,70)==false,
             "Nothing included after clear().");

  t.report();

  return 0;
}
