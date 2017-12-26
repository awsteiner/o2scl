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
#include <o2scl/constants.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
   
  cout.setf(ios::scientific);

  // Load several mass formulae to make distributions out of

  nucmass_ame ame12;
  o2scl_hdf::ame_load(ame12,"../../data/o2scl/nucmass/ame12.o2","ame12.o2");

  nucmass_ame_exp amex12;
  o2scl_hdf::ame_load(amex12,"../../data/o2scl/nucmass/ame12.o2","ame12.o2");

  nucmass_mnmsk mth;
  o2scl_hdf::mnmsk_load(mth,"../../data/o2scl/nucmass/");

  nucmass_mnmsk_exp mexp;
  o2scl_hdf::mnmsk_load(mexp,"../../data/o2scl/nucmass/");

  nucmass_ame ame03;
  o2scl_hdf::ame_load(ame03,"../../data/o2scl/nucmass/ame03.o2","ame03.o2");
  
  nucmass_ame_exp amex03;
  o2scl_hdf::ame_load(amex03,"../../data/o2scl/nucmass/ame03.o2","ame03.o2");

  vector<nucleus> fdx;
  nucdist_set(fdx,mth,"1",400,true);

  t.report();
  return 0;
}
