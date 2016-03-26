/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

/* Example: ex_nucmass.cpp
   -------------------------------------------------------------------
   Demonstrate nuclear mass formulas by comparing them with 
   the Audi et al. (2012) atomic mass evaluation. This example
   computes the absolute deviation in the mass excess. 
*/
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/table_units.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucdist.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_wlw.h>
#include <o2scl/nucmass_ktuy.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  // ------------------------------------------------------
  // The most recent experimental masses from the 2012 AME as a baseline

  nucmass_ame_exp ame;
  o2scl_hdf::ame_load(ame,"12");

  // ------------------------------------------------------
  // Instantiate and load all of the nuclear mass objects. Some of
  // them require a separate HDF5 file

  nucmass_semi_empirical se;
  nucmass_mnmsk mnmsk;
  o2scl_hdf::mnmsk_load(mnmsk);
  nucmass_hfb hfb14;
  o2scl_hdf::hfb_load(hfb14,14);
  nucmass_hfb_sp hfb21;
  o2scl_hdf::hfb_sp_load(hfb21,21);
  nucmass_hfb_sp hfb27;
  o2scl_hdf::hfb_sp_load(hfb27,27);
  nucmass_ame_exp ame03;
  o2scl_hdf::ame_load(ame03,"03");
  nucmass_dz_table dz;
  nucmass_ktuy ktuy05;
  nucmass_dvi dvi;
  nucmass_wlw ws32("WS3.2");
  nucmass_wlw ws36("WS3.6");

  // ------------------------------------------------------
  // List of pointers to all masses for convenience

  static const size_t n_tables=11;
  nucmass *massp[n_tables]={&se,&mnmsk,&hfb14,&hfb21,&hfb27,
			    &ame03,&dz,&ktuy05,&dvi,&ws32,&ws36};

  // ------------------------------------------------------
  // Create a list of all of the experimental masses

  vector<nucleus> ame_dist;
  nucdist_set(ame_dist,ame,"N>7 && Z>7");

  // ------------------------------------------------------
  // Fit the semi-empirical and DvI (2009) mass formulas to 
  // the 2012 AME data

  static const size_t n_fits=2;
  nucmass_fit mf;
  // The default number of trials isn't enough for the DvI 
  // model, so we increase it
  mf.def_mmin.ntrial*=10;
  // Use the same list as above
  mf.dist=ame_dist;
  // The RMS deviation in the mass excess
  double res;
  // Fit both mass formulas
  nucmass_fit_base *fitp[n_fits]={&se,&dvi};
  for(size_t i=0;i<n_fits;i++) {
    mf.fit(*(fitp[i]),res);
  }

  // ------------------------------------------------------
  // Create a table to store the data

  table_units<> tu;
  tu.line_of_names(((string)"Z N ame se mnmsk hfb14 hfb21 ")+
		   "hfb27 ame03 dz96 ktuy05 dvi ws32 ws36");

  // ------------------------------------------------------
  // Fill the table

  for(size_t i=0;i<ame_dist.size();i++) {
    vector<double> line;
    line.push_back(ame_dist[i].Z);
    line.push_back(ame_dist[i].N);
    double ame_mass=ame.mass_excess(ame_dist[i].Z,ame_dist[i].N);
    line.push_back(ame_mass);
    for(size_t j=0;j<n_tables;j++) {
      if (massp[j]->is_included(ame_dist[i].Z,ame_dist[i].N)) {
	double val=ame_mass-
	  massp[j]->mass_excess(ame_dist[i].Z,ame_dist[i].N);
	line.push_back(val);
      } else {
	line.push_back(0.0);
      }
    }
    tu.line_of_data(line);
  }
  
  // ------------------------------------------------------
  // Output the table to a file

  hdf_file hf;
  hf.open_or_create("ex_nucmass_table.o2");
  hdf_output(hf,tu,"nuclear_masses");
  hf.close();
 
  t.report();
  return 0;
}
// End of example
