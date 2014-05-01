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

  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

  nucmass_ame_exp ame;
  o2scl_hdf::ame_load(ame,"12");

  nucdist_full ame_fd(ame);
  
  nucmass_semi_empirical sm;
  nucmass_fit mf;
  double res;
  mf.set_exp_mass(ame);
  mf.fit(sm,res);
  mf.eval(sm,res);
  cout.width(15);
  cout << "Semi-empir:" << " " << res << endl;
  
  nucmass_mnmsk mn;
  o2scl_hdf::mnmsk_load(mn);
  mf.eval(mn,res);
  cout.width(15);
  cout << "FRDM:" << " " << res << endl;

  nucmass_hfb hfb;
  o2scl_hdf::hfb_load(hfb,14);
  mf.eval(hfb,res);
  cout.width(15);
  cout << "HFB:" << " " << res << endl;

  nucmass_ame_exp ame03;
  o2scl_hdf::ame_load(ame03,"03");
  mf.eval(ame03,res);
  cout.width(15);
  cout << "AME '03:" << " " << res << endl;

  nucmass_dz_table dz;
  mf.eval(dz,res);
  cout.width(15);
  cout << "Duf-zuk:" << " " << res << endl;

  nucmass_ktuy ktuy;
  mf.eval(ktuy,res);
  cout.width(15);
  cout << "Koura:" << " " << res << endl;
  
  table_units<> tu;
  tu.line_of_names("Z N ame sm mn hfb ame03 dz ktuy");

  typedef nucdist::iterator iter;

  iter ame_it=ame_fd.begin();
  for(iter it=ame_fd.begin();it!=ame_fd.end();it++) {
    double m_ame=ame.mass_excess(it->Z,it->N);
    double m_sm=0.0;
    if (it->Z>=8 && it->N>=8) {
      m_sm=sm.mass_excess(it->Z,it->N);
    }
    double m_mn=0.0;
    if (mn.is_included(it->Z,it->N)) {
      m_mn=mn.mass_excess(it->Z,it->N);
    }
    double m_hfb=0.0;
    if (hfb.is_included(it->Z,it->N)) {
      m_hfb=hfb.mass_excess(it->Z,it->N);
    }
    double m_ame03=0.0;
    if (ame03.is_included(it->Z,it->N)) {
      m_ame03=ame03.mass_excess(it->Z,it->N);
    }
    double m_dz=0.0;
    if (dz.is_included(it->Z,it->N)) {
      m_dz=dz.mass_excess(it->Z,it->N);
    }
    double m_ktuy=0.0;
    if (ktuy.is_included(it->Z,it->N)) {
      m_ktuy=ktuy.mass_excess(it->Z,it->N);
    }
    double line[9]={((double)it->Z),((double)it->N),m_ame,m_sm,
		    m_mn,m_hfb,m_ame03,m_dz,m_ktuy};
    tu.line_of_data(9,line);
  }

  hdf_file hf;
  hf.open_or_create("ex_nucmass_table.o2");
  hdf_output(hf,tu,"nuclear_masses");
  hf.close();
 
  t.report();
  return 0;
}
// End of example
