/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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

#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

void o2scl_hdf::mnmsk_load(o2scl::nucmass_mnmsk &mnmsk, std::string model,
			   string filename) {

  if (model=="mnmsk97") {
    if (filename.size()==0) {
      filename=o2scl::o2scl_settings.get_data_dir()+"/nucmass/mnmsk.o2";
    }
  } else {
    model="msis16";
    if (filename.size()==0) {
      filename=o2scl::o2scl_settings.get_data_dir()+"/nucmass/msis16.o2";
    }
  }

  if (model=="mnmsk97") {
    size_t offset[34]={HOFFSET(o2scl::nucmass_mnmsk::entry,N),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Z),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,A),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,eps2),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,eps3),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,eps4),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,eps6),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,eps6sym),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,beta2),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,beta3),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,beta4),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,beta6),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Emic),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Mth),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Mexp),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,sigmaexp),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,EmicFL),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,MthFL),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,spinp),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,spinn),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,gapp),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,gapn),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,be),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,S1n),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,S2n),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,PA),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,PAm1),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,PAm2),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Qbeta),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Tbeta),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,S1p),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,S2p),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Qalpha),
		       HOFFSET(o2scl::nucmass_mnmsk::entry,Talpha)};

    o2scl::nucmass_mnmsk::entry me;
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
    hf.open(filename);
    hid_t file=hf.get_current_id();

    int nrecords;
    std::string reference;
    hf.geti("nrecords",nrecords);
    if (nrecords<=0) {
      O2SCL_ERR("Number of records <= 0 in mnmsk_load().",exc_efailed);
    }
    hf.gets("reference",reference);

    std::vector<o2scl::nucmass_mnmsk::entry> m(nrecords);
    //=new o2scl::nucmass_mnmsk::entry[nrecords];
    herr_t status=H5TBread_table
      (file,"mnmsk.o2",sizeof(o2scl::nucmass_mnmsk::entry),offset,sizes,
       &(m[0]));
    
    mnmsk.set_data(nrecords,m,reference);
    
    hf.close();

  } else {

    hdf_file hf;
    hf.open(filename);
    std::string name;
    table_units<> tab;
    hdf_input_n(hf,tab,name);
    hf.close();
    size_t nr=tab.get_nlines();
    std::vector<o2scl::nucmass_mnmsk::entry> m(nr);
    for(size_t j=0;j<nr;j++) {
      m[j].N=tab.get("N",j);
      m[j].Z=tab.get("Z",j);
      m[j].A=tab.get("A",j);
      m[j].eps2=tab.get("eps2",j);
      m[j].eps3=tab.get("eps3",j);
      m[j].eps4=tab.get("eps4",j);
      m[j].eps6=tab.get("eps6",j);
      m[j].eps6sym=0.0;
      m[j].beta2=tab.get("beta2",j);
      m[j].beta3=tab.get("beta3",j);
      m[j].beta4=tab.get("beta4",j);
      m[j].beta6=tab.get("beta6",j);
      m[j].Emic=tab.get("Emic",j);
      m[j].Mth=tab.get("Mth",j);
      m[j].Mexp=tab.get("Mexp",j);
      m[j].sigmaexp=tab.get("sigmaexp",j);
      m[j].EmicFL=tab.get("EFLmic",j);
      m[j].MthFL=tab.get("MFLth",j);
      m[j].spinp[0]='N';
      m[j].spinp[1]='A';
      m[j].spinp[2]='\0';
      m[j].spinn[0]='N';
      m[j].spinn[1]='A';
      m[j].spinn[2]='\0';
      m[j].gapp=1.0e99;
      m[j].gapn=1.0e99;
      m[j].be=1.0e99;
      m[j].S1n=1.0e99;
      m[j].S2n=1.0e99;
      m[j].PA=1.0e99;
      m[j].PAm1=1.0e99;
      m[j].PAm2=1.0e99;
      m[j].Qbeta=1.0e99;
      m[j].Tbeta=1.0e99;
      m[j].S1p=1.0e99;
      m[j].S2p=1.0e99;
      m[j].Qalpha=1.0e99;
      m[j].Talpha=1.0e99;
    }
    string reference=((std::string)"Möller, Sierk, Ichikawa, and ")+
      "Sagawa, At. Data and Nucl. Data Tables 109 (2016), 1.";
    mnmsk.set_data(nr,m,reference);
  }

  return;
}

void o2scl_hdf::hfb_load(o2scl::nucmass_hfb &hfb, size_t model, 
			 string filename) {
    
  if (filename.size()==0) {
    filename=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }
    
  std::string tname;
  if (model==2) {
    tname="/hfb2.o2";
  } else if (model==8) {
    tname="/hfb8.o2";
  } else if (model==14) {
    tname="/hfb14.o2";
  } else {
    tname="/hfb14_v0.o2";
  }
  filename=filename+tname;
  
  size_t offset[12]={HOFFSET(o2scl::nucmass_hfb::entry,N),
		     HOFFSET(o2scl::nucmass_hfb::entry,Z),
		     HOFFSET(o2scl::nucmass_hfb::entry,A),
		     HOFFSET(o2scl::nucmass_hfb::entry,bet2),
		     HOFFSET(o2scl::nucmass_hfb::entry,bet4),
		     HOFFSET(o2scl::nucmass_hfb::entry,Rch),
		     HOFFSET(o2scl::nucmass_hfb::entry,def_wig),
		     HOFFSET(o2scl::nucmass_hfb::entry,Sn),
		     HOFFSET(o2scl::nucmass_hfb::entry,Sp),
		     HOFFSET(o2scl::nucmass_hfb::entry,Qbet),
		     HOFFSET(o2scl::nucmass_hfb::entry,Mcal),
		     HOFFSET(o2scl::nucmass_hfb::entry,Err)};
    
  o2scl::nucmass_hfb::entry he;

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
  hf.open(filename);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in hfb_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  std::vector<o2scl::nucmass_hfb::entry> m(nrecords);
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::nucmass_hfb::entry),offset,
     sizes,&(m[0]));
    
  hfb.set_data(nrecords,m,reference);
    
  hf.close();

  return;
}

void o2scl_hdf::hfb_sp_load(nucmass_hfb_sp &hfb, size_t model,
			    string filename) {
  
  if (filename.size()==0) {
    filename=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }
    
  std::string tname;
  if (model==17) {
    tname="/hfb17.o2";
  } else if (model==21) {
    tname="/hfb21.o2";
  } else if (model==22) {
    tname="/hfb22.o2";
  } else if (model==23) {
    tname="/hfb23.o2";
  } else if (model==24) {
    tname="/hfb24.o2";
  } else if (model==25) {
    tname="/hfb25.o2";
  } else if (model==26) {
    tname="/hfb26.o2";
  } else if (model==27) {
    tname="/hfb27.o2";
  } else if (model==28) {
    tname="/hfb28.o2";
  } else if (model==29) {
    tname="/hfb29.o2";
  } else if (model==30) {
    tname="/hfb30.o2";
  } else if (model==31) {
    tname="/hfb31.o2";
  } else {
    tname="/hfb32.o2";
  }
  filename=filename+tname;
  
  size_t offset[16]={HOFFSET(o2scl::nucmass_hfb_sp::entry,N),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Z),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,A),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,bet2),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,bet4),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Rch),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,def_wig),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Sn),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Sp),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Qbet),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Mcal),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Err),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Jexp),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Jth),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Pexp),
		     HOFFSET(o2scl::nucmass_hfb_sp::entry,Pth)};
    
  o2scl::nucmass_hfb_sp::entry he;

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
  hf.open(filename);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in ame_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  std::vector<o2scl::nucmass_hfb_sp::entry> m(nrecords);
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::nucmass_hfb_sp::entry),offset,
     sizes,&(m[0]));

  hfb.set_data(nrecords,m,reference);

  hf.close();

  return;
}

void o2scl_hdf::bskg_load(nucmass_hfb_sp &hfb, size_t model,
                          string filename) {

  if (filename.size()==0) {
    filename=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
  }

  std::string tname;
  if (model==1) {
    tname="/bskg1.o2";
  } else if (model==2) {
    tname="/bskg2.o2";
  } else if (model==3) {
    tname="/bskg3.o2";
  } else if (model==4) {
    tname="/bskg4.o2";
  } else {
    O2SCL_ERR("Invalid model in bskg_load().",exc_einval);
  }
  filename=filename+tname;

  // BSkG1, BSkG2, and BSkG4 fill in a different set of extra
  // fields than BSkG3 (see the documentation in hdf_nucmass_io.h),
  // so they are stored in a table with 8 additional columns beyond
  // the 30 used for BSkG3.
  if (model==1 || model==2 || model==4) {

    size_t offset2[38]={HOFFSET(o2scl::nucmass_hfb_sp::entry,N),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Z),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,A),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,bet2),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,bet4),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Rch),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Sn),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Sp),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Qbet),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Mcal),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Err),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Jexp),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Jth),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Pexp),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Pth),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,gamma),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,beta20),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,beta22),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,beta30),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,beta32),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,S2n),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,S2p),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,delta3n),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,delta3p),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,delta5n),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,delta5p),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,rc4),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,I1),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,I2),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,I3),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,Erot),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,avgap_n),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,avgap_p),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,rc_exp),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,rc_err),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,MOI),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,par_p),
                        HOFFSET(o2scl::nucmass_hfb_sp::entry,par_n)};

    o2scl::nucmass_hfb_sp::entry he2;

    size_t sizes2[38]={sizeof(he2.N),
                       sizeof(he2.Z),
                       sizeof(he2.A),
                       sizeof(he2.bet2),
                       sizeof(he2.bet4),
                       sizeof(he2.Rch),
                       sizeof(he2.Sn),
                       sizeof(he2.Sp),
                       sizeof(he2.Qbet),
                       sizeof(he2.Mcal),
                       sizeof(he2.Err),
                       sizeof(he2.Jexp),
                       sizeof(he2.Jth),
                       sizeof(he2.Pexp),
                       sizeof(he2.Pth),
                       sizeof(he2.gamma),
                       sizeof(he2.beta20),
                       sizeof(he2.beta22),
                       sizeof(he2.beta30),
                       sizeof(he2.beta32),
                       sizeof(he2.S2n),
                       sizeof(he2.S2p),
                       sizeof(he2.delta3n),
                       sizeof(he2.delta3p),
                       sizeof(he2.delta5n),
                       sizeof(he2.delta5p),
                       sizeof(he2.rc4),
                       sizeof(he2.I1),
                       sizeof(he2.I2),
                       sizeof(he2.I3),
                       sizeof(he2.Erot),
                       sizeof(he2.avgap_n),
                       sizeof(he2.avgap_p),
                       sizeof(he2.rc_exp),
                       sizeof(he2.rc_err),
                       sizeof(he2.MOI),
                       sizeof(he2.par_p),
                       sizeof(he2.par_n)};

    hdf_file hf2;
    hf2.open(filename);
    hid_t file2=hf2.get_current_id();

    int nrecords2;
    std::string reference2;
    hf2.geti("nrecords",nrecords2);
    if (nrecords2<=0) {
      O2SCL_ERR("Number of records <= 0 in bskg_load().",exc_efailed);
    }
    hf2.gets("reference",reference2);

    std::vector<o2scl::nucmass_hfb_sp::entry> m2(nrecords2);
    herr_t status2=H5TBread_table
      (file2,tname.c_str(),sizeof(o2scl::nucmass_hfb_sp::entry),offset2,
       sizes2,&(m2[0]));

    hfb.set_data(nrecords2,m2,reference2);

    hf2.close();

    return;
  }

  size_t offset[30]={HOFFSET(o2scl::nucmass_hfb_sp::entry,N),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Z),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,A),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,bet2),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,bet4),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Rch),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Sn),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Sp),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Qbet),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Mcal),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Err),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Jexp),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Jth),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Pexp),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,Pth),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,gamma),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,beta20),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,beta22),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,beta30),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,beta32),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,S2n),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,S2p),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,delta3n),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,delta3p),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,delta5n),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,delta5p),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,rc4),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,I1),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,I2),
                     HOFFSET(o2scl::nucmass_hfb_sp::entry,I3)};

  o2scl::nucmass_hfb_sp::entry he;

  size_t sizes[30]={sizeof(he.N),
                    sizeof(he.Z),
                    sizeof(he.A),
                    sizeof(he.bet2),
                    sizeof(he.bet4),
                    sizeof(he.Rch),
                    sizeof(he.Sn),
                    sizeof(he.Sp),
                    sizeof(he.Qbet),
                    sizeof(he.Mcal),
                    sizeof(he.Err),
                    sizeof(he.Jexp),
                    sizeof(he.Jth),
                    sizeof(he.Pexp),
                    sizeof(he.Pth),
                    sizeof(he.gamma),
                    sizeof(he.beta20),
                    sizeof(he.beta22),
                    sizeof(he.beta30),
                    sizeof(he.beta32),
                    sizeof(he.S2n),
                    sizeof(he.S2p),
                    sizeof(he.delta3n),
                    sizeof(he.delta3p),
                    sizeof(he.delta5n),
                    sizeof(he.delta5p),
                    sizeof(he.rc4),
                    sizeof(he.I1),
                    sizeof(he.I2),
                    sizeof(he.I3)};

  hdf_file hf;
  hf.open(filename);
  hid_t file=hf.get_current_id();

  int nrecords;
  std::string reference;
  hf.geti("nrecords",nrecords);
  if (nrecords<=0) {
    O2SCL_ERR("Number of records <= 0 in bskg_load().",exc_efailed);
  }
  hf.gets("reference",reference);

  std::vector<o2scl::nucmass_hfb_sp::entry> m(nrecords);
  herr_t status=H5TBread_table
    (file,tname.c_str(),sizeof(o2scl::nucmass_hfb_sp::entry),offset,
     sizes,&(m[0]));

  hfb.set_data(nrecords,m,reference);

  hf.close();

  return;
}

