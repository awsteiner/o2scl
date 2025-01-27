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
  t.set_output_level(3);

#ifdef O2SCL_NEVER_DEFINED
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"20");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("20",false,1);
    nucdist_set(dist2,ame);

    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      //t.test_gen(((std::string)(&en.bdmode[0]))==
      //((std::string)(&en2.bdmode[0])),"16");
      //<< en.bde << " " << en2.bde << endl;
      //t.test_gen(en.bde==en2.bde,"17");
      //cout << en.N << " " 
      //<< (std::string)(&en.el[0]) << " "
      //<< en.dbde << " " << en2.dbde << endl;
      //t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      //t.test_gen(en.A2==en2.A2,"20");
      //std::cout << en.Z << " " << en.N << " "
      //<< (std::string)(&en.el[0]) << " "
      //<< en.damass << " " << en2.damass << std::endl;
      //t.test_gen(en.amass==en2.amass,"21");
      //t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;

    nucmass_ame ame;
    ame_load(ame,"20round");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("20round",false,2);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }

  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"16");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("16",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"16round");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("16round",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"12");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("12",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"03");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("03",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"03round");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("03round",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"95exp");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("95exp",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }
  
  if (true) {
    
    vector<nucleus> dist, dist2;
    
    nucmass_ame ame;
    ame_load(ame,"95rmd");
    nucdist_set(dist,ame);
    
    nucmass_ame ame;
    ame.load("95rmd",false,1);
    nucdist_set(dist2,ame);
    
    t.test_gen(dist.size()==dist2.size(),"size");
    
    for(size_t i=0;i<dist.size();i++) {
      nucmass_ame::entry en=ame.get_ZN(dist[i].Z,dist[i].N);
      nucmass_ame::entry en2=ame.get_ZN(dist2[i].Z,dist2[i].N);
      
      t.test_gen(en.NMZ==en2.NMZ,"1");
      t.test_gen(en.N==en2.N,"2");
      t.test_gen(en.Z==en2.Z,"3");
      t.test_gen(en.A==en2.A,"4");
      t.test_gen(((std::string)(&en.el[0]))==((std::string)(&en2.el[0])),"5");
      t.test_gen(((std::string)(&en.orig[0]))==
          ((std::string)(&en2.orig[0])),"6");
      t.test_gen(en.mass==en2.mass,"7");
      t.test_gen(en.dmass==en2.dmass,"8");
      t.test_gen(en.mass_acc==en2.mass_acc,"9");
      t.test_gen(en.be==en2.be,"10");
      t.test_gen(en.dbe==en2.dbe,"11");
      t.test_gen(en.be_acc==en2.be_acc,"12");
      t.test_gen(en.beoa==en2.beoa,"13");
      t.test_gen(en.dbeoa==en2.dbeoa,"14");
      t.test_gen(en.beoa_acc==en2.beoa_acc,"15");
      t.test_gen(((std::string)(&en.bdmode[0]))==
      ((std::string)(&en2.bdmode[0])),"16");
      t.test_gen(en.bde==en2.bde,"17");
      t.test_gen(en.dbde==en2.dbde,"18");
      t.test_gen(en.bde_acc==en2.bde_acc,"19");
      t.test_gen(en.A2==en2.A2,"20");
      t.test_gen(en.amass==en2.amass,"21");
      t.test_gen(en.damass==en2.damass,"22");
      t.test_gen(en.amass_acc==en2.amass_acc,"23");
    }

  }

#endif
  
  t.report();
  return 0;
}

