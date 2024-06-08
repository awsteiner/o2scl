/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2024, Andrew W. Steiner
  
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
/* Example: ex_nstar_rot.cpp
   -------------------------------------------------------------------
*/

#include <o2scl/test_mgr.h>
#include <o2scl/nstar_rot.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  nstar_rot nst;
  
  eos_had_skyrme sk;
  // Load the SLy4 EOS
  o2scl_hdf::skyrme_load(sk,"SLy4");

  nstar_cold nco;
  nco.def_tov.verbose=0;
  nco.set_eos(sk);

  // Compute the Skyrme EOS in beta-equilibrium. The nstar_rot class
  // cannot handle fine grids, so we increase the grid spacing.
  nco.dnb*=2.0;
  nco.calc_eos();
  std::shared_ptr<table_units<> > eos=nco.get_eos_results();

  // The set_eos_fm() function adds its own crust, so we make
  // sure to remove the low-density part of the SLy4 EOS
  eos->delete_rows_func("nb<0.08");

  // Send the EOS to the nstar_rot object
  eos_nstar_rot_interp enri;
  enri.set_eos_fm(eos->get_nlines(),(*eos)["ed"],(*eos)["pr"],(*eos)["nb"]);
  nst.set_eos(enri);
  
  // First compute non-rotating configurations
  cout << "Non-rotating stars: " << endl;
  cout << "ed_cent      M (Msun)     R (km)" << endl;
  for (double ed_cent=3.0e14;ed_cent<3.3e15;ed_cent*=1.5) {
    nst.fix_cent_eden_non_rot(ed_cent);
    cout << ed_cent << " " << nst.Mass/nst.MSUN << " "
	 << nst.R_e/1.0e5 << endl;
  }
  cout << endl;

  // Compute configurations rotating at the Keplerian rotation rate
  cout << "Keplerian rotation rate:" << endl;
  cout << "ed_cent      M (Msun)     R (km)       "
       << "Omega        Omega_K" << endl;
  for (double ed_cent=3.0e14;ed_cent<3.3e15;ed_cent*=1.5) {
    if (ed_cent>3.01e14) {
      // After the first point, use the previous solution as an
      // initial guess
      nst.fix_cent_eden_with_kepler_alt(ed_cent,true);
    } else {
      nst.fix_cent_eden_with_kepler_alt(ed_cent);
    }
    cout << ed_cent << " " << nst.Mass/nst.MSUN << " "
	 << nst.R_e/1.0e5 << " " << nst.Omega << " "
	 << nst.Omega_K << endl;
  }
  cout << endl;

  // Compute configurations rotating at a frequency of 200 Hz
  double f=200.0, ang_vel=2.0*o2scl_const::pi*f;
  cout << "Rotating at 200 Hz: " << endl;
  cout << "ed_cent      M (Msun)     R (km)" << endl;
  for (double ed_cent=3.0e14;ed_cent<1.02e15;ed_cent*=1.5) {
    if (ed_cent>3.01e14) {
      nst.fix_cent_eden_ang_vel_alt(ed_cent,ang_vel,true);
    } else {
      nst.fix_cent_eden_ang_vel_alt(ed_cent,ang_vel);
    }
    cout << ed_cent << " " << nst.Mass/nst.MSUN << " "
	 << nst.R_e/1.0e5 << endl;
  }
  cout << endl;

  t.report();

  return 0;
}
// End of example
