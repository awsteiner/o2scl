/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/cloud_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"20");
    nucdist_set(dist,ame);
    
    nucmass_ame2 ame2;
    ame2.load("20",false,1);
    nucdist_set(dist2,ame2);
    
    cout << dist.size() << " " << dist2.size() << endl;
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame2::entry en2=ame2.get_ZN(dist2[i].Z,dist2[i].N);
      
      if (en.NMZ!=en2.NMZ) exit(-1);
      if (en.N!=en2.N) exit(-1);
      if (en.Z!=en2.Z) exit(-1);
      if (en.A!=en2.A) exit(-1);
      if (((std::string)(&en.el[0]))!=((std::string)(&en2.el[0]))) exit(-1);
      if (((std::string)(&en.orig[0]))!=
          ((std::string)(&en2.orig[0]))) exit(-1);
      if (en.mass!=en2.mass) exit(-1);
      if (en.dmass!=en2.dmass) exit(-1);
      if (en.mass_acc!=en2.mass_acc) exit(-1);
      if (en.be!=en2.be) exit(-1);
      if (en.dbe!=en2.dbe) exit(-1);
      if (en.be_acc!=en2.be_acc) exit(-1);
      if (en.beoa!=en2.beoa) exit(-1);
      if (en.dbeoa!=en2.dbeoa) exit(-1);
      if (en.beoa_acc!=en2.beoa_acc) exit(-1);
      if (((std::string)(&en.bdmode[0]))!=
          ((std::string)(&en2.bdmode[0]))) exit(-1);
      if (en.bde!=en2.bde) exit(-1);
      if (en.dbde!=en2.dbde) exit(-1);
      if (en.bde_acc!=en2.bde_acc) exit(-1);
      if (en.A2!=en2.A2) exit(-1);
      if (en.amass!=en2.amass) exit(-1);
      if (en.damass!=en2.damass) exit(-1);
      if (en.amass_acc!=en2.amass_acc) exit(-1);
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"20round");
    nucdist_set(dist,ame);
    
    nucmass_ame2 ame2;
    ame2.load("20round",false,1);
    nucdist_set(dist2,ame2);
    
    cout << dist.size() << " " << dist2.size() << endl;
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame2::entry en2=ame2.get_ZN(dist2[i].Z,dist2[i].N);
      
      if (en.NMZ!=en2.NMZ) exit(-1);
      if (en.N!=en2.N) exit(-1);
      if (en.Z!=en2.Z) exit(-1);
      if (en.A!=en2.A) exit(-1);
      if (((std::string)(&en.el[0]))!=((std::string)(&en2.el[0]))) exit(-1);
      if (((std::string)(&en.orig[0]))!=
          ((std::string)(&en2.orig[0]))) exit(-1);
      if (en.mass!=en2.mass) exit(-1);
      if (en.dmass!=en2.dmass) exit(-1);
      if (en.mass_acc!=en2.mass_acc) exit(-1);
      if (en.be!=en2.be) exit(-1);
      if (en.dbe!=en2.dbe) exit(-1);
      if (en.be_acc!=en2.be_acc) exit(-1);
      if (en.beoa!=en2.beoa) exit(-1);
      if (en.dbeoa!=en2.dbeoa) exit(-1);
      if (en.beoa_acc!=en2.beoa_acc) exit(-1);
      if (((std::string)(&en.bdmode[0]))!=
          ((std::string)(&en2.bdmode[0]))) exit(-1);
      if (en.bde!=en2.bde) exit(-1);
      if (en.dbde!=en2.dbde) exit(-1);
      if (en.bde_acc!=en2.bde_acc) exit(-1);
      if (en.A2!=en2.A2) exit(-1);
      if (en.amass!=en2.amass) exit(-1);
      if (en.damass!=en2.damass) exit(-1);
      if (en.amass_acc!=en2.amass_acc) exit(-1);
    }

  }

  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"16");
    nucdist_set(dist,ame);
    
    nucmass_ame2 ame2;
    ame2.load("16",false,1);
    nucdist_set(dist2,ame2);
    
    cout << dist.size() << " " << dist2.size() << endl;
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame2::entry en2=ame2.get_ZN(dist2[i].Z,dist2[i].N);
      
      if (en.NMZ!=en2.NMZ) exit(-1);
      if (en.N!=en2.N) exit(-1);
      if (en.Z!=en2.Z) exit(-1);
      if (en.A!=en2.A) exit(-1);
      if (((std::string)(&en.el[0]))!=((std::string)(&en2.el[0]))) exit(-1);
      if (((std::string)(&en.orig[0]))!=
          ((std::string)(&en2.orig[0]))) exit(-1);
      if (en.mass!=en2.mass) exit(-1);
      if (en.dmass!=en2.dmass) exit(-1);
      if (en.mass_acc!=en2.mass_acc) exit(-1);
      if (en.be!=en2.be) exit(-1);
      if (en.dbe!=en2.dbe) exit(-1);
      if (en.be_acc!=en2.be_acc) exit(-1);
      if (en.beoa!=en2.beoa) exit(-1);
      if (en.dbeoa!=en2.dbeoa) exit(-1);
      if (en.beoa_acc!=en2.beoa_acc) exit(-1);
      if (((std::string)(&en.bdmode[0]))!=
          ((std::string)(&en2.bdmode[0]))) exit(-1);
      if (en.bde!=en2.bde) exit(-1);
      if (en.dbde!=en2.dbde) exit(-1);
      if (en.bde_acc!=en2.bde_acc) exit(-1);
      if (en.A2!=en2.A2) exit(-1);
      if (en.amass!=en2.amass) exit(-1);
      if (en.damass!=en2.damass) exit(-1);
      if (en.amass_acc!=en2.amass_acc) exit(-1);
    }

  }
  
  
  t.report();
  return 0;
}

