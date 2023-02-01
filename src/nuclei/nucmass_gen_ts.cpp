/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/nucmass_gen.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  nucmass_gen ng;
  nucmass_ame ame;

  o2scl_hdf::ame_load_ext(ame,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2");

  string fnames[10]={"../../data/o2scl/nucmass/frib_mex/ddme2.o2",
		    "../../data/o2scl/nucmass/frib_mex/ddmed.o2",
		    "../../data/o2scl/nucmass/frib_mex/ddpc1.o2",
		    "../../data/o2scl/nucmass/frib_mex/nl3s.o2",
		    "../../data/o2scl/nucmass/frib_mex/sly4_all.o2",
		    "../../data/o2scl/nucmass/frib_mex/skms_all.o2",
		    "../../data/o2scl/nucmass/frib_mex/skp_all.o2",
		    "../../data/o2scl/nucmass/frib_mex/sv_min_all.o2",
		    "../../data/o2scl/nucmass/frib_mex/unedf0_all.o2",
		    "../../data/o2scl/nucmass/frib_mex/unedf1_all.o2"};
		    
  // Show that the nucmass_gen gives reasonable (but not great)
  // values for the binding energy of Lead 208
  cout << "AME2003 : ";
  cout << ame.mass_excess(82,126) << " ";
  cout << ame.binding_energy(82,126) << " ";
  cout << ame.total_mass(82,126) << endl;
  for(size_t i=0;i<10;i++) {
    if (i<4) {
      ng.load_be(fnames[i],"E",1.0,true);
    } else {
      ng.load_be(fnames[i],"Binding_Energy__MeV_",1.0,true);
    }
    cout << "GEN    : ";
    cout << ng.mass_excess(82,126) << " ";
    cout << ng.binding_energy(82,126) << " ";
    cout << ng.total_mass(82,126) << endl;
    t.test_abs(ng.mass_excess(82,126),ame.mass_excess(82,126),
	       6.0,"gen lead test");
  }
  t.report();

  return 0;
}
