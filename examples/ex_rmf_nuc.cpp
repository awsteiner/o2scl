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

/* Example: ex_rmf_nuc.cpp
   -------------------------------------------------------------------

   This example uses the NL3 interaction to compute the structure
   of Lead-208
*/
#include <o2scl/rmf_nucleus.h>
#include <o2scl/test_mgr.h>
#ifdef O2SCL_HDF
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
#ifdef O2SCL_HDF
using namespace o2scl_hdf;
#endif

int main(int argv, char *argc[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  rmf_nucleus rn;
  rn.set_verbose(0);

  rmf_load(rn.def_rmf,"NL3");
  rn.run_nucleus(82,126,3,3);
  
  bool neutron=false;
  cout << "Proton orbitals: " << endl;
  cout << "Index State      Energy       Nodes Deg." << endl;
  for(int i=0;i<rn.nlevels;i++) {
    if (rn.levels[i].isospin<0.0 && neutron==false) {
      cout << endl;
      cout << "Neutron orbitals: " << endl;
      cout << "Index State      Energy       Nodes Deg." << endl;
      neutron=true;
    }
    cout.width(2);
    cout << i << " ";
    cout.width(12);
    cout.setf(ios::left);
    cout << rn.levels[i].state << " ";
    cout.unsetf(ios::left);
    cout << rn.levels[i].eigen << " "
	 << rn.levels[i].nodes << "     "
	 << rn.levels[i].twojp1 << endl;
  }
  cout << endl;

  cout << "Proton orbitals: " << endl;
  cout << "Index State      Energy       Nodes Deg." << endl;
  neutron=false;
  for(int i=0;i<rn.nuolevels;i++) {
    if (rn.unocc_levels[i].isospin<0.0 && neutron==false) {
      cout << endl;
      cout << "Neutron orbitals: " << endl;
      cout << "Index State      Energy       Nodes Deg." << endl;
      neutron=true;
    }
    cout.width(2);
    cout << i << " ";
    cout.width(12);
    cout.setf(ios::left);
    cout << rn.unocc_levels[i].state << " ";
    cout.unsetf(ios::left);
    cout << rn.unocc_levels[i].eigen << " "
	 << rn.unocc_levels[i].nodes << "     "
	 << rn.unocc_levels[i].twojp1 << endl;
  }
  cout << endl;

  cout << "Proton RMS radius       :  " << rn.rprms << " fm." << endl;
  cout << "Neutron RMS radius      :  " << rn.rnrms << " fm." << endl;
  cout << "Total energy per nucleon: " << rn.etot << " MeV." << endl;
  cout << endl;
  
  t.test_rel(-7.842551,rn.etot,1.0e-5,"Lead binding energy");

#ifdef O2SCL_HDF

  hdf_file hf;
  hdf_file hf2;
  hf.open_or_create("ex_rmf_nuc_prof.o2");
  hf2.open_or_create("ex_rmf_nuc_chden.o2");
  o2_shared_ptr<table_units<> >::type profiles=rn.get_profiles();
  o2_shared_ptr<table_units<> >::type charge_dens=rn.get_chden();
  hdf_output(hf,*profiles,"profiles");
  hdf_output(hf2,*charge_dens,"charge_densities");
  hf.close();
  hf2.close();

  //profiles.line_of_names(((string)"r rhop rhon sig ome rho ")+
  //"coul chden rhosp rhosn ");
  //chden_table.line_of_names("x chden1 chdenc");

#endif

  t.report();

  return 0;
}
// End of example
