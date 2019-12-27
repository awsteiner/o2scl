/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2020, Andrew W. Steiner
  
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
#include <o2scl/nucmass_densmat.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  dense_matter dm;
  int ret;

  // Nuclear masses
  nucmass_ame ame;
  o2scl_hdf::ame_load_ext(ame,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2",true);

  // Test dense_matter

  nucdist_set(dm.dist,ame);
  for(size_t i=0;i<dm.dist.size();i++) {
    if (i<10) dm.dist[i].n=((double)i)/100.0;
    else dm.dist[i].n=0.0;
  }
  dm.eta_nuc.resize(dm.dist.size());

  dm.output(cout,0);
  dm.output(cout,1);
  
  // Test nucmass_densmat

  /*
  nucmass_densmat &nd=nse.nuc_dens;
  nd.set_mass(ame);
  double t1, t2, t3, t4;
  nd.test_derivatives(1.0e-6,t1,t2,t3,t4);
  t.test_rel(t1,0.0,1.0e-6,"dEdnp");
  t.test_rel(t3,0.0,3.0e-6,"dEdnneg");
  */
  
  t.report();
  return 0;
}


