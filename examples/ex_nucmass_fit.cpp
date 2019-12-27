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

/* Example: ex_nucmass_fit.cpp
   -------------------------------------------------------------------
*/

#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_fit.h>
#ifdef O2SCL_HDF
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_nucmass_io.h>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

#ifdef O2SCL_HDF

  // The RMS deviation of the fit
  double res;
  // The mass formula to be fitted
  nucmass_semi_empirical sem;
  // The fitting class
  nucmass_fit mf;

  // Load the experimental data
  nucmass_ame ame;
  o2scl_hdf::ame_load(ame,"12");
  nucdist_set(mf.dist,ame);

  // Perform the fit
  mf.fit(sem,res);

  // Output the results
  cout << sem.B << " " << sem.Sv << " " << sem.Ss << " " 
       << sem.Ec << " " << sem.Epair << endl;
  cout << res << endl;
  t.test_gen(res<4.0,"Successful fit.");

#else
  cout << "No fitting was performed because O2scl appears not to have been "
       << "compiled with HDF support." << endl;
#endif
  
  t.report();
  return 0;
}
// End of example
