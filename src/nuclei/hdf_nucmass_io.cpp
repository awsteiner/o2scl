/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

void o2scl_hdf::ame_load(o2scl::ame_mass &ame, std::string file_name, 
			 std::string table_name) {

  size_t offset[23]={HOFFSET(o2scl::ame_mass::entry,NMZ),
		     HOFFSET(o2scl::ame_mass::entry,N),
		     HOFFSET(o2scl::ame_mass::entry,Z),
		     HOFFSET(o2scl::ame_mass::entry,A),
		     HOFFSET(o2scl::ame_mass::entry,el),
		     HOFFSET(o2scl::ame_mass::entry,orig),
		     HOFFSET(o2scl::ame_mass::entry,mass),
		     HOFFSET(o2scl::ame_mass::entry,dmass),
		     HOFFSET(o2scl::ame_mass::entry,mass_acc),
		     HOFFSET(o2scl::ame_mass::entry,be),
		     HOFFSET(o2scl::ame_mass::entry,dbe),
		     HOFFSET(o2scl::ame_mass::entry,be_acc),
		     HOFFSET(o2scl::ame_mass::entry,beoa),
		     HOFFSET(o2scl::ame_mass::entry,dbeoa),
		     HOFFSET(o2scl::ame_mass::entry,beoa_acc),
		     HOFFSET(o2scl::ame_mass::entry,bdmode),
		     HOFFSET(o2scl::ame_mass::entry,bde),
		     HOFFSET(o2scl::ame_mass::entry,dbde),
		     HOFFSET(o2scl::ame_mass::entry,bde_acc),
		     HOFFSET(o2scl::ame_mass::entry,A2),
		     HOFFSET(o2scl::ame_mass::entry,amass),
		     HOFFSET(o2scl::ame_mass::entry,damass),
		     HOFFSET(o2scl::ame_mass::entry,amass_acc)};
  o2scl::ame_mass::entry ae;
  size_t sizes[23]={sizeof(ae.NMZ),
		    sizeof(ae.N),
		    sizeof(ae.Z),
		    sizeof(ae.A),
		    sizeof(ae.el),
		    sizeof(ae.orig),
		    sizeof(ae.mass),
		    sizeof(ae.dmass),
		    sizeof(ae.mass_acc),
		    sizeof(ae.be),
		    sizeof(ae.dbe),
		    sizeof(ae.be_acc),
		    sizeof(ae.beoa),
		    sizeof(ae.dbeoa),
		    sizeof(ae.beoa_acc),
		    sizeof(ae.bdmode),
		    sizeof(ae.bde),
		    sizeof(ae.dbde),
		    sizeof(ae.bde_acc),
		    sizeof(ae.A2),
		    sizeof(ae.amass),
		    sizeof(ae.damass),
		    sizeof(ae.amass_acc)};

  hdf_file hf;
  hf.open(file_name);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in ame_load().",exc_efailed);
  }
  hf.gets_fixed("reference",reference);
    
  o2scl::ame_mass::entry *m=new o2scl::ame_mass::entry[nrecords];
  herr_t status=H5TBread_table(file,table_name.c_str(),
			       sizeof(o2scl::ame_mass::entry),
			       offset,sizes,m);
  ame.n=nrecords;
  ame.mass=m;
  ame.reference=reference;
  ame.last=nrecords/2;
    
  hf.close();

  return;
}

void o2scl_hdf::ame_load(o2scl::ame_mass &ame, std::string name) {
  
  std::string file_name, table_name;
  file_name=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  if (name=="95exp") {
    file_name+="/ame95exp.o2";
    table_name="ame95exp.o2";
  } else if (name=="95rmd") {
    file_name+="/ame95rmd.o2";
    table_name="ame95rmd.o2";
  } else if (name=="03round") {
    file_name+="/ame03round.o2";
    table_name="ame03round.o2";
  } else if (name=="03") {
    file_name+="/ame03.o2";
    table_name="ame03.o2";
  } else {
    file_name+="/ame12.o2";
    table_name="ame12.o2";
  }
  
  ame_load(ame,file_name,table_name);
  return;
}

void o2scl_hdf::mnmsk_load(o2scl::mnmsk_mass &mnmsk, string dir) {
    
  if (dir.size()==0) {
    dir=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }
    
  std::string fname, tname="/mnmsk.o2";
  fname=dir+tname;
  
  size_t offset[34]={HOFFSET(o2scl::mnmsk_mass_entry,N),
		     HOFFSET(o2scl::mnmsk_mass_entry,Z),
		     HOFFSET(o2scl::mnmsk_mass_entry,A),
		     HOFFSET(o2scl::mnmsk_mass_entry,eps2),
		     HOFFSET(o2scl::mnmsk_mass_entry,eps3),
		     HOFFSET(o2scl::mnmsk_mass_entry,eps4),
		     HOFFSET(o2scl::mnmsk_mass_entry,eps6),
		     HOFFSET(o2scl::mnmsk_mass_entry,eps6sym),
		     HOFFSET(o2scl::mnmsk_mass_entry,beta2),
		     HOFFSET(o2scl::mnmsk_mass_entry,beta3),
		     HOFFSET(o2scl::mnmsk_mass_entry,beta4),
		     HOFFSET(o2scl::mnmsk_mass_entry,beta6),
		     HOFFSET(o2scl::mnmsk_mass_entry,Emic),
		     HOFFSET(o2scl::mnmsk_mass_entry,Mth),
		     HOFFSET(o2scl::mnmsk_mass_entry,Mexp),
		     HOFFSET(o2scl::mnmsk_mass_entry,sigmaexp),
		     HOFFSET(o2scl::mnmsk_mass_entry,EmicFL),
		     HOFFSET(o2scl::mnmsk_mass_entry,MthFL),
		     HOFFSET(o2scl::mnmsk_mass_entry,spinp),
		     HOFFSET(o2scl::mnmsk_mass_entry,spinn),
		     HOFFSET(o2scl::mnmsk_mass_entry,gapp),
		     HOFFSET(o2scl::mnmsk_mass_entry,gapn),
		     HOFFSET(o2scl::mnmsk_mass_entry,be),
		     HOFFSET(o2scl::mnmsk_mass_entry,S1n),
		     HOFFSET(o2scl::mnmsk_mass_entry,S2n),
		     HOFFSET(o2scl::mnmsk_mass_entry,PA),
		     HOFFSET(o2scl::mnmsk_mass_entry,PAm1),
		     HOFFSET(o2scl::mnmsk_mass_entry,PAm2),
		     HOFFSET(o2scl::mnmsk_mass_entry,Qbeta),
		     HOFFSET(o2scl::mnmsk_mass_entry,Tbeta),
		     HOFFSET(o2scl::mnmsk_mass_entry,S1p),
		     HOFFSET(o2scl::mnmsk_mass_entry,S2p),
		     HOFFSET(o2scl::mnmsk_mass_entry,Qalpha),
		     HOFFSET(o2scl::mnmsk_mass_entry,Talpha)};

  o2scl::mnmsk_mass_entry me;
  size_t sizes[34]={sizeof(me.N),
		    sizeof(me.Z),
		    sizeof(me.A),
		    sizeof(me.eps2),
		    sizeof(me.eps3),
		    sizeof(me.eps4),
		    sizeof(me.eps6),
		    sizeof(me.eps6sym),
		    sizeof(me.beta2),
		    sizeof(me.beta3),
		    sizeof(me.beta4),
		    sizeof(me.beta6),
		    sizeof(me.Emic),
		    sizeof(me.Mth),
		    sizeof(me.Mexp),
		    sizeof(me.sigmaexp),
		    sizeof(me.EmicFL),
		    sizeof(me.MthFL),
		    sizeof(me.spinp),
		    sizeof(me.spinn),
		    sizeof(me.gapp),
		    sizeof(me.gapn),
		    sizeof(me.be),
		    sizeof(me.S1n),
		    sizeof(me.S2n),
		    sizeof(me.PA),
		    sizeof(me.PAm1),
		    sizeof(me.PAm2),
		    sizeof(me.Qbeta),
		    sizeof(me.Tbeta),
		    sizeof(me.S1p),
		    sizeof(me.S2p),
		    sizeof(me.Qalpha),
		    sizeof(me.Talpha)};

  hdf_file hf;
  hf.open(fname);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in mnmsk_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  o2scl::mnmsk_mass_entry *m=new o2scl::mnmsk_mass_entry[nrecords];
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::mnmsk_mass_entry),offset,sizes,m);
    
  mnmsk.set_data(nrecords,m,reference);
    
  hf.close();

  return;
}

void o2scl_hdf::hfb_load(o2scl::hfb_mass &hfb, size_t model, 
			string dir) {
    
  if (dir.size()==0) {
    dir=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }
    
  std::string fname, tname;
  if (model==2) {
    tname="/hfb2.o2";
  } else if (model==8) {
    tname="/hfb8.o2";
  } else {
    tname="/hfb14.o2";
  }
  fname=dir+tname;
  
  size_t offset[12]={HOFFSET(o2scl::hfb_mass_entry,N),
		     HOFFSET(o2scl::hfb_mass_entry,Z),
		     HOFFSET(o2scl::hfb_mass_entry,A),
		     HOFFSET(o2scl::hfb_mass_entry,bet2),
		     HOFFSET(o2scl::hfb_mass_entry,bet4),
		     HOFFSET(o2scl::hfb_mass_entry,Rch),
		     HOFFSET(o2scl::hfb_mass_entry,def_wig),
		     HOFFSET(o2scl::hfb_mass_entry,Sn),
		     HOFFSET(o2scl::hfb_mass_entry,Sp),
		     HOFFSET(o2scl::hfb_mass_entry,Qbet),
		     HOFFSET(o2scl::hfb_mass_entry,Mcal),
		     HOFFSET(o2scl::hfb_mass_entry,Err)};
    
  o2scl::hfb_mass_entry he;

  size_t sizes[12]={sizeof(he.N),
		    sizeof(he.Z),
		    sizeof(he.A),
		    sizeof(he.bet2),
		    sizeof(he.bet4),
		    sizeof(he.Rch),
		    sizeof(he.def_wig),
		    sizeof(he.Sn),
		    sizeof(he.Sp),
		    sizeof(he.Qbet),
		    sizeof(he.Mcal),
		    sizeof(he.Err)};
    
  hdf_file hf;
  hf.open(fname);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in hfb_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  o2scl::hfb_mass_entry *m=new o2scl::hfb_mass_entry[nrecords];
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::hfb_mass_entry),offset,sizes,m);
    
  hfb.set_data(nrecords,m,reference);
    
  hf.close();

  return;
}

void o2scl_hdf::hfb_sp_load(hfb_sp_mass &hfb, size_t model, string dir) {
  
  if (dir.size()==0) {
    dir=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }
    
  std::string fname, tname;
  if (model==17) {
    tname="/hfb17.o2";
  } else {
    tname="/hfb21.o2";
  }
  fname=dir+tname;
  
  size_t offset[16]={HOFFSET(o2scl::hfb_sp_mass_entry,N),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Z),
		     HOFFSET(o2scl::hfb_sp_mass_entry,A),
		     HOFFSET(o2scl::hfb_sp_mass_entry,bet2),
		     HOFFSET(o2scl::hfb_sp_mass_entry,bet4),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Rch),
		     HOFFSET(o2scl::hfb_sp_mass_entry,def_wig),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Sn),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Sp),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Qbet),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Mcal),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Err),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Jexp),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Jth),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Pexp),
		     HOFFSET(o2scl::hfb_sp_mass_entry,Pth)};
    
  o2scl::hfb_sp_mass_entry he;

  size_t sizes[16]={sizeof(he.N),
		    sizeof(he.Z),
		    sizeof(he.A),
		    sizeof(he.bet2),
		    sizeof(he.bet4),
		    sizeof(he.Rch),
		    sizeof(he.def_wig),
		    sizeof(he.Sn),
		    sizeof(he.Sp),
		    sizeof(he.Qbet),
		    sizeof(he.Mcal),
		    sizeof(he.Err),
		    sizeof(he.Jexp),
		    sizeof(he.Jth),
		    sizeof(he.Pexp),
		    sizeof(he.Pth)};
    
  hdf_file hf;
  hf.open(fname);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in ame_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  o2scl::hfb_sp_mass_entry *m=new o2scl::hfb_sp_mass_entry[nrecords];
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::hfb_sp_mass_entry),offset,sizes,m);
    
  hfb.set_data(nrecords,m,reference);
    
  hf.close();

  return;
}

