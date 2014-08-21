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

/* Example: ex_nucmass.cpp
   -------------------------------------------------------------------
   Demonstrate nuclear mass formulas by comparing them with 
   the Audi et al. (2003) atomic mass evaluation. This example
   computes the RMS deviation in the mass excess. 
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
#include <o2scl/nucmass_ktuy.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  // The most recent experimental masses as a baseline
  nucmass_ame_exp ame;
  o2scl_hdf::ame_load(ame,"12");

  // Instantiate and load all of the nuclear mass objects
  nucmass_semi_empirical sm;
  nucmass_mnmsk mn;
  o2scl_hdf::mnmsk_load(mn);
  nucmass_hfb hfb;
  o2scl_hdf::hfb_load(hfb,14);
  nucmass_ame_exp ame03;
  o2scl_hdf::ame_load(ame03,"03");
  nucmass_dz_table dz;
  nucmass_ktuy ktuy;

  // List of pointers to all masses
  size_t n_tables=6;
  nucmass *massp[7]={&sm,&mn,&hfb,&ame03,&dz,&ktuy};

  // Create a distribution with all of the experimental masses
  vector<nucleus> ame_dist;
  nucdist_set(ame_dist,ame);

  // Create a smaller distribution to fit to
  vector<nucleus> fit_dist;
  nucdist_set(ame_dist,ame,"N>7 & Z>7");

  // Fit to the experimental masses
  size_t n_fits=2;
  nucmass_fit mf;
  mf.set_dist(ame_dist);
  double res;
  nucmass_fit_base *fitp[2]={&sm,&mn};
  for(size_t i=0;i<n_fits;i++) {
    mf.fit(*fitp,res);
  }

  // Create a table to store the data
  table_units<> tu;
  tu.line_of_names("Z N ame sm mn hfb ame03 dz ktuy");

  // Create the table
  for(size_t i=0;i<ame_dist.size();i++) {
    vector<double> line;
    line.push_back(ame_dist[i].Z);
    line.push_back(ame_dist[i].N);
    double ame=ame.mass_excess(ame_dist[i].Z,ame_dist[i].N);
    line.push_back(ame);
    for(size_t j=0;j<n_tables;j++) {
      if (massp->is_included(ame_dist[i].Z,ame_dist[i].N)) {
	line.push_back(ame-massp->mass_excess(ame_dist[i].Z,ame_dist[i].N));
      } else {
	line.push_back(0.0);
      }
    }
    tu.line_of_data(line);
  }
  
  // Output the table to a file
  hdf_file hf;
  hf.open_or_create("ex_nucmass_table.o2");
  hdf_output(hf,tu,"nuclear_masses");
  hf.close();
 
  t.report();
  return 0;
}
// End of example
