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

#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl;
using namespace o2scl_hdf;

void o2scl_hdf::gogny_load(o2scl::eos_had_gogny &ge, std::string model, 
			   bool external) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();

  hdf_file hf;

  if (external) {
    fname=model;
    hf.open(fname);
    std::string table_name;
    hdf_input(hf,ge.t3d,table_name);
    hf.close();
  } else {
    fname=dir+"/gogny.o2";
    hf.open(fname);
    cout << "Reading." << endl;
    hdf_input(hf,ge.t3d,model);
    cout << "Done reading." << endl;
    hf.close();
  }
  
  return;
}

void o2scl_hdf::rmf_load(o2scl::eos_had_rmf &rmf, std::string model, 
			bool external) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    fname=dir+"/rmfdata/"+model+".o2";
  }
  
  hdf_file hf;
  hf.open(fname);

  double gs=0.0, gw=0.0, gr=0.0;
  int itmp;
  bool oakstyle, tokistyle, standardstyle;

  hf.geti_def("oakstyle",0,itmp);
  oakstyle=itmp;
  hf.geti_def("tokistyle",0,itmp);
  tokistyle=itmp;
  hf.geti_def("standardstyle",0,itmp);
  standardstyle=itmp;
  
  hf.getd("ms",rmf.ms);
  hf.getd("mw",rmf.mw);
  hf.getd("mr",rmf.mr);
  hf.getd("mnuc",rmf.mnuc);
  rmf.ms/=o2scl_const::hc_mev_fm; 
  rmf.mw/=o2scl_const::hc_mev_fm; 
  rmf.mr/=o2scl_const::hc_mev_fm; 
  rmf.mnuc/=o2scl_const::hc_mev_fm;
  
  if (standardstyle==true) {
    hf.getd("Cs2",rmf.cs);
    hf.getd("Cw2",rmf.cw);
    hf.getd("Cr2",rmf.cr);
    hf.getd("b",rmf.b);
    hf.getd("c",rmf.c);
    rmf.cs=sqrt(rmf.cs)/rmf.mnuc;
    rmf.cw=sqrt(rmf.cw)/rmf.mnuc;
    rmf.cr=sqrt(rmf.cr)/rmf.mnuc;
  } else if (oakstyle==true || tokistyle==true) {
    hf.getd("gs",gs);
    hf.getd("gw",gw);
    hf.getd("gr",gr);
    hf.getd("g2",rmf.b);
    hf.getd("g3",rmf.c);
    rmf.b/=-rmf.mnuc*pow(fabs(gs),3.0);
    rmf.c/=pow(gs,4.0);
    gr*=2.0;
    rmf.cs=gs/rmf.ms;
    rmf.cw=gw/rmf.mw;
    rmf.cr=gr/rmf.mr;
  } else {
    hf.getd("cs",rmf.cs);
    hf.getd("cw",rmf.cw);
    hf.getd("cr",rmf.cr);
    hf.getd("b",rmf.b);
    hf.getd("c",rmf.c);
  }
  
  if (tokistyle==true) {
    hf.getd_def("zeta",0.0,rmf.zeta);
    rmf.zeta=rmf.zeta/gw/gw/gw/gw*6.0;
  } else {
    hf.getd_def("zeta",0.0,rmf.zeta);
  }
    
  hf.getd_def("xi",0.0,rmf.xi);
  hf.getd_def("a1",0.0,rmf.a1);
  hf.getd_def("a2",0.0,rmf.a2);
  hf.getd_def("a3",0.0,rmf.a3);
  hf.getd_def("a4",0.0,rmf.a4);
  hf.getd_def("a5",0.0,rmf.a5);
  hf.getd_def("a6",0.0,rmf.a6);
  hf.getd_def("b1",0.0,rmf.b1);
  hf.getd_def("b2",0.0,rmf.b2);
  hf.getd_def("b3",0.0,rmf.b3);

  return;
}
  
void o2scl_hdf::skyrme_load(o2scl::eos_had_skyrme &sk, std::string model, 
			   bool external) {

  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    fname=dir+"/skdata/"+model+".o2";
  }

  hdf_file hf;
  hf.open(fname);

  hf.getd("t0hc",sk.t0);
  hf.getd("t2hc",sk.t2);
  hf.getd("t3hc",sk.t3);
  sk.t0/=o2scl_const::hc_mev_fm; 
  sk.t2/=o2scl_const::hc_mev_fm; 
  sk.t3/=o2scl_const::hc_mev_fm; 

  hf.getd("x0",sk.x0);
  hf.getd("x2",sk.x2);
  hf.getd("x3",sk.x3);
  
  int itmp;
  hf.geti_def("dpfix",0,itmp);
  bool dpfix=itmp;

  if (dpfix==true) {
    sk.t1=-sk.t2/3.0*(5.0+4.0*sk.x2);
    sk.x1=-(4.0+5.0*sk.x2)/(5.0+4.0*sk.x2);
    sk.alpha=1.0/3.0;
    sk.a=1.0;
    sk.b=0.0;
  } else {
    hf.getd_def("t1hc",0.0,sk.t1);
    sk.t1/=o2scl_const::hc_mev_fm;
    hf.getd_def("x1",0.0,sk.x1);
    hf.getd_def("a",1.0,sk.a);
    hf.getd_def("b",0.0,sk.b);
    hf.getd_def("alpha",0.0,sk.alpha);
  }

  hf.geti_def("pdmode",0,itmp);
  bool pdmode=itmp;
  if (pdmode==true) {
    double pfp, pfn;
    hf.getd("pairfp",pfp);
    hf.getd("pairfn",pfn);
    sk.W0=(pfp+pfn)/4.0/o2scl_const::hc_mev_fm;
  } else {
    hf.getd("W0hc",sk.W0);
    sk.W0/=o2scl_const::hc_mev_fm;
  }

  hf.gets_def("reference","",sk.reference);
  hf.getd_def("b4",0.0,sk.b4);
  hf.getd_def("b4p",0.0,sk.b4p);

  return;
}

#ifdef O2SCL_NEVER_DEFINED

/*
  This is a draft of a new version designed to name the 
  fields a bit more sensibly. It's not finished yet and
  we need to redo the Skyrme model data files
*/

void o2scl_hdf::skyrme_load(hdf_file &hf, o2scl::eos_had_skyrme &sk, 
			    std::string name) {

  hf.getd("t0_hc",sk.t0);
  hf.getd("t1_hc",sk.t1);
  hf.getd("t2_hc",sk.t2);
  hf.getd("t3_hc",sk.t3);
  sk.t0/=o2scl_const::hc_mev_fm; 
  sk.t1/=o2scl_const::hc_mev_fm;
  sk.t2/=o2scl_const::hc_mev_fm; 
  sk.t3/=o2scl_const::hc_mev_fm; 

  hf.getd("x0",sk.x0);
  hf.getd("x1",sk.x1);
  hf.getd("x2",sk.x2);
  hf.getd("x3",sk.x3);

  hf.getd("alpha",sk.alpha);

  hf.getd_def("a",1.0,sk.a);
  hf.getd_def("b",0.0,sk.b);
  
  hf.getd_def("delta_n_hc",sk.W0);
  sk.W0/=o2scl_const::hc_mev_fm;
  //hf.getd_def("delta_p_hc",sk.W0);
  //sk.W0/=o2scl_const::hc_mev_fm;

  hf.gets_def("reference","",sk.reference);
  hf.getd_def("b4_hc",0.0,sk.b4);
  sk.b4/=o2scl_const::hc_mev_fm;
  hf.getd_def("b4p_hc",0.0,sk.b4p);
  sk.b4p/=o2scl_const::hc_mev_fm;

  return;
}

#endif

void o2scl_hdf::skyrme_write(o2scl::eos_had_skyrme &sk, std::string model) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  fname=dir+"/skdata/"+model+".o2";
  
  hdf_file hf;
  hf.open(fname);
  skyrme_write(hf,sk,model);
  hf.close();

  return;
}
  
void o2scl_hdf::skyrme_write(hdf_file &hf, o2scl::eos_had_skyrme &sk, 
			    std::string name) {

  /*
  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
      
  // Add typename
  hf.sets_fixed("o2scl_type","eos_had_skyrme");
  */

  // Write data
  hf.setd("t0hc",sk.t0*o2scl_const::hc_mev_fm);
  hf.setd("t1hc",sk.t1*o2scl_const::hc_mev_fm);
  hf.setd("t2hc",sk.t2*o2scl_const::hc_mev_fm);
  hf.setd("t3hc",sk.t3*o2scl_const::hc_mev_fm);
  hf.setd("x0",sk.x0);
  hf.setd("x1",sk.x1);
  hf.setd("x2",sk.x2);
  hf.setd("x3",sk.x3);
  hf.setd("a",sk.a);
  hf.setd("b",sk.b);
  hf.setd("alpha",sk.alpha);
  hf.setd("delta_n_hc",sk.W0*o2scl_const::hc_mev_fm);
  hf.setd("delta_p_hc",sk.W0*o2scl_const::hc_mev_fm);
  hf.setd("b4",sk.b4);
  hf.setd("b4p",sk.b4p);
  hf.setd("W0hc",0.0);
  hf.sets("reference",sk.reference);

  /*
  // Close eos_had_skyrme group
  hf.close_group(group);
  
  // Return location to previous value
  hf.set_current_id(top);
  */

  return;
}
  
