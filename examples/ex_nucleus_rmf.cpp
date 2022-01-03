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

// sphinx-example-start
/* Example: ex_nucleus_rmf.cpp
   -------------------------------------------------------------------

   This example uses the NL3 interaction to compute the structure
   of Lead-208
*/
#include <o2scl/nucleus_rmf.h>
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

int lead_chden_exp(std::shared_ptr<table_units<> > profiles);

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  nucleus_rmf rn;
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

  std::shared_ptr<table_units<> > profiles=rn.get_profiles();
  lead_chden_exp(profiles);
  
#ifdef O2SCL_HDF

  hdf_file hf;
  hdf_file hf2;
  hf.open_or_create("ex_nucleus_rmf_prof.o2");
  hf2.open_or_create("ex_nucleus_rmf_chden.o2");
  std::shared_ptr<table_units<> > charge_dens=rn.get_chden();
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

int lead_chden_exp(std::shared_ptr<table_units<> > profiles) {
  
  // The experimental charge density for Lead 208 from de Vries, et
  // al. At. Data Nucl. Data Tables 36 (1987) 495 using the
  // Sum-of-Gaussians method
  double rr[12]={0.1,0.7,1.6,2.1,2.7,3.5,4.2,
                 5.1,6.0,6.6,7.6,8.7};
  double qq[12]={0.003845,0.009724,0.033093,0.000120,
                 0.083107,0.080869,0.139957,0.260892,
                 0.336013,0.0033637,0.018729,0.000020};
  double g=1.7/sqrt(1.5);
  double a[12];
  for(size_t i=0;i<12;i++) {
    a[i]=82.0*qq[i]/2.0/pow(o2scl_const::pi,1.5)/
      pow(g,3.0)/(1.0+2.0*rr[i]*rr[i]/g/g);
  }
  
  // Add experimental profile to table
  profiles->new_column("chden_exp");
  for(size_t i=0;i<profiles->get_nlines();i++) {
    double val=0.0;
    for(size_t j=0;j<12;j++) {
      val+=a[j]*(exp(-pow((profiles->get("r",i)-rr[j])/g,2.0))+
		 exp(-pow((profiles->get("r",i)+rr[j])/g,2.0)));
    }
    profiles->set("chden_exp",i,val);
  }
  return 0;
}
