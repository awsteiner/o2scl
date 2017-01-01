/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2017, Andrew W. Steiner

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
#include <o2scl/eos_had_potential.h>

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

void o2scl_hdf::eos_had_strings_list() {
  vector<string> list;
  list.push_back("apr");
  list.push_back("gogny");
  list.push_back("skyrme a");
  list.push_back("skyrme b");
  list.push_back("skyrme BSk10");
  list.push_back("skyrme BSk11");
  list.push_back("skyrme BSk12");
  list.push_back("skyrme BSk13");
  list.push_back("skyrme BSk14");
  list.push_back("skyrme BSk16");
  list.push_back("skyrme BSk1");
  list.push_back("skyrme BSk2");
  list.push_back("skyrme BSk2p");
  list.push_back("skyrme BSk3");
  list.push_back("skyrme BSk4");
  list.push_back("skyrme BSk5");
  list.push_back("skyrme BSk6");
  list.push_back("skyrme BSk7");
  list.push_back("skyrme BSk8");
  list.push_back("skyrme BSk9");
  list.push_back("skyrme E");
  list.push_back("skyrme Es");
  list.push_back("skyrme FitA");
  list.push_back("skyrme FitB");
  list.push_back("skyrme FitK");
  list.push_back("skyrme FitKs");
  list.push_back("skyrme FitL");
  list.push_back("skyrme Gs");
  list.push_back("skyrme KDE0v1");
  list.push_back("skyrme KDE0v");
  list.push_back("skyrme LNS");
  list.push_back("skyrme Ly5");
  list.push_back("skyrme MSk1");
  list.push_back("skyrme MSk2");
  list.push_back("skyrme MSk3");
  list.push_back("skyrme MSk4");
  list.push_back("skyrme MSk5");
  list.push_back("skyrme MSk5s");
  list.push_back("skyrme MSk6");
  list.push_back("skyrme MSk7");
  list.push_back("skyrme MSk8");
  list.push_back("skyrme MSk9");
  list.push_back("skyrme MSkA");
  list.push_back("skyrme mst0.81");
  list.push_back("skyrme mst0.90");
  list.push_back("skyrme mst1");
  list.push_back("skyrme NRAPR2");
  list.push_back("skyrme NRAPR");
  list.push_back("skyrme PeEVs");
  list.push_back("skyrme PeHF");
  list.push_back("skyrme PeSIs");
  list.push_back("skyrme QMC1");
  list.push_back("skyrme QMC2");
  list.push_back("skyrme QMC3");
  list.push_back("skyrme RATP");
  list.push_back("skyrme Rs");
  list.push_back("skyrme SGII");
  list.push_back("skyrme SGI");
  list.push_back("skyrme SIII");
  list.push_back("skyrme SIIIs");
  list.push_back("skyrme SII");
  list.push_back("skyrme SI");
  list.push_back("skyrme SIp");
  list.push_back("skyrme SIV");
  list.push_back("skyrme SK255");
  list.push_back("skyrme SK272");
  list.push_back("skyrme SkI1");
  list.push_back("skyrme SkI2");
  list.push_back("skyrme SkI3");
  list.push_back("skyrme SkI4");
  list.push_back("skyrme SkI5");
  list.push_back("skyrme SkI6");
  list.push_back("skyrme SkkT8");
  list.push_back("skyrme SkM1");
  list.push_back("skyrme SkMDIx0");
  list.push_back("skyrme SkMDIx1");
  list.push_back("skyrme SkMDIxm1");
  list.push_back("skyrme SkMDIxm2");
  list.push_back("skyrme SkM");
  list.push_back("skyrme SkMP");
  list.push_back("skyrme SkMs");
  list.push_back("skyrme SkNF1");
  list.push_back("skyrme SkNF2");
  list.push_back("skyrme SkO");
  list.push_back("skyrme SkOp");
  list.push_back("skyrme SkP");
  list.push_back("skyrme SKRA");
  list.push_back("skyrme SkSC10");
  list.push_back("skyrme SkSC11");
  list.push_back("skyrme SkSC14");
  list.push_back("skyrme SkSC15");
  list.push_back("skyrme SkSC1");
  list.push_back("skyrme SkSC2");
  list.push_back("skyrme SkSC3");
  list.push_back("skyrme SkSC4");
  list.push_back("skyrme SkSC4o");
  list.push_back("skyrme SkSC5");
  list.push_back("skyrme SkSC6");
  list.push_back("skyrme SkT1");
  list.push_back("skyrme SkT1s");
  list.push_back("skyrme SkT2");
  list.push_back("skyrme SkT3");
  list.push_back("skyrme SkT3s");
  list.push_back("skyrme SkT4");
  list.push_back("skyrme SkT5");
  list.push_back("skyrme SkT6");
  list.push_back("skyrme SkT7");
  list.push_back("skyrme SkT8");
  list.push_back("skyrme SkT9");
  list.push_back("skyrme SkTK");
  list.push_back("skyrme SkT");
  list.push_back("skyrme SkXce");
  list.push_back("skyrme SkXm");
  list.push_back("skyrme SkX");
  list.push_back("skyrme Skxs15");
  list.push_back("skyrme Skxs20");
  list.push_back("skyrme Skxs25");
  list.push_back("skyrme Skyrme1p");
  list.push_back("skyrme SKz0");
  list.push_back("skyrme SKz1");
  list.push_back("skyrme SKz2");
  list.push_back("skyrme SKz3");
  list.push_back("skyrme SKz4");
  list.push_back("skyrme SKzm1");
  list.push_back("skyrme SLy0");
  list.push_back("skyrme SLy10");
  list.push_back("skyrme SLy1");
  list.push_back("skyrme SLy230a");
  list.push_back("skyrme SLy230b");
  list.push_back("skyrme SLy2");
  list.push_back("skyrme SLy3");
  list.push_back("skyrme SLy4");
  list.push_back("skyrme SLy5");
  list.push_back("skyrme SLy6");
  list.push_back("skyrme SLy7");
  list.push_back("skyrme SLy8");
  list.push_back("skyrme SLy9");
  list.push_back("skyrme SV-bas");
  list.push_back("skyrme SVII");
  list.push_back("skyrme SVI");
  list.push_back("skyrme SV-K218");
  list.push_back("skyrme SV-K226");
  list.push_back("skyrme SV-K241");
  list.push_back("skyrme SV-kap00");
  list.push_back("skyrme SV-kap02");
  list.push_back("skyrme SV-kap06");
  list.push_back("skyrme SV-mas07");
  list.push_back("skyrme SV-mas08");
  list.push_back("skyrme SV-mas10");
  list.push_back("skyrme SV-min");
  list.push_back("skyrme SV");
  list.push_back("skyrme SV-sym28");
  list.push_back("skyrme SV-sym32");
  list.push_back("skyrme SV-sym34");
  list.push_back("skyrme SV-tls");
  list.push_back("skyrme T");
  list.push_back("skyrme UNEDF0");
  list.push_back("skyrme UNEDF1");
  list.push_back("skyrme UNEDF2");
  list.push_back("skyrme v070");
  list.push_back("skyrme v075");
  list.push_back("skyrme v080");
  list.push_back("skyrme v090");
  list.push_back("skyrme v100");
  list.push_back("skyrme v105");
  list.push_back("skyrme v110");
  list.push_back("skyrme Z");
  list.push_back("skyrme Zs");
  list.push_back("skyrme Zss");
  list.push_back("rmf BMPII");
  list.push_back("rmf BMPI");
  list.push_back("rmf es25n15");
  list.push_back("rmf es25new");
  list.push_back("rmf es25");
  list.push_back("rmf es25small");
  list.push_back("rmf es275n15");
  list.push_back("rmf es275new");
  list.push_back("rmf es275");
  list.push_back("rmf es30n15");
  list.push_back("rmf es30new");
  list.push_back("rmf es30");
  list.push_back("rmf es325n15");
  list.push_back("rmf es325");
  list.push_back("rmf es35n15");
  list.push_back("rmf es35");
  list.push_back("rmf FPWC");
  list.push_back("rmf FSUGold");
  list.push_back("rmf IUFSU");
  list.push_back("rmf L1");
  list.push_back("rmf L2");
  list.push_back("rmf L3");
  list.push_back("rmf L-BF");
  list.push_back("rmf L-HS");
  list.push_back("rmf L-W");
  list.push_back("rmf L-Z");
  list.push_back("rmf NL-065");
  list.push_back("rmf NL-06");
  list.push_back("rmf NL-075");
  list.push_back("rmf NL-07");
  list.push_back("rmf NL1");
  list.push_back("rmf NL2");
  list.push_back("rmf NL3");
  list.push_back("rmf NL4");
  list.push_back("rmf NL-B1");
  list.push_back("rmf NL-B2");
  list.push_back("rmf NL-SH");
  list.push_back("rmf NL-Z");
  list.push_back("rmf PL-40");
  list.push_back("rmf PL-Z");
  list.push_back("rmf RAPRhdp");
  list.push_back("rmf RAPR");
  list.push_back("rmf S271");
  list.push_back("rmf SFHo");
  list.push_back("rmf SFHx");
  list.push_back("rmf SR1");
  list.push_back("rmf SR2");
  list.push_back("rmf SR3");
  list.push_back("rmf TM1");
  list.push_back("rmf TM2");
  list.push_back("rmf Z271");
  list.push_back("pot MDI0");
  list.push_back("pot MDI1");
  list.push_back("pot BGBD_das");
  list.push_back("pot PAL11");
  list.push_back("pot PAL12");
  list.push_back("pot PAL13");
  list.push_back("pot PAL21");
  list.push_back("pot PAL22");
  list.push_back("pot PAL23");
  list.push_back("pot PAL31");
  list.push_back("pot PAL32");
  list.push_back("pot PAL33");
  list.push_back("pot BPALb11");
  list.push_back("pot BPALb12");
  list.push_back("pot BPALb13");
  list.push_back("pot BPALb21");
  list.push_back("pot BPALb22");
  list.push_back("pot BPALb23");
  list.push_back("pot BPALb31");
  list.push_back("pot SL12");
  list.push_back("pot GBD0");
  list.push_back("pot GBD1");
  list.push_back("pot CKLxm2");
  list.push_back("pot CKLxm1");
  list.push_back("pot CKLx0");
  list.push_back("pot CKLx1");

  vector<string> reformatted;
  o2scl::screenify(list.size(),list,reformatted);
  cout << "Current list (" << list.size() << ") of EOSs:" << endl;
  for(size_t i=0;i<reformatted.size();i++) {
    cout << reformatted[i] << endl;
  }
  return;
}

eos_had_base *o2scl_hdf::eos_had_strings(std::string type,
					 std::string name) {
  if (type=="skyrme") {
    eos_had_skyrme *sk=new eos_had_skyrme;
    skyrme_load(*sk,name);
    return sk;
  } else if (type=="apr") {
    eos_had_apr *apr=new eos_had_apr;
    return apr;
  } else if (type=="gogny") {
    eos_had_gogny *gogny=new eos_had_gogny;
    return gogny;
  } else if (type=="pot") {
    eos_had_potential *pot=new eos_had_potential;
    if (name=="MDI0") {
      pot->Au=-95.98/o2scl_const::hc_mev_fm;
      pot->Al=-120.57/o2scl_const::hc_mev_fm;
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->x=0.0;
      pot->form=pot->mdi_form;
    } else if (name=="MDI1") {
      pot->Au=-187.27/o2scl_const::hc_mev_fm;
      pot->Al=-29.28/o2scl_const::hc_mev_fm;
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->x=1.0;
      pot->form=pot->mdi_form;
    } else if (name=="BGBD_das") {
      pot->Au=-192.0/o2scl_const::hc_mev_fm;
      pot->Al=-96.0/o2scl_const::hc_mev_fm;
      pot->B=203.3/o2scl_const::hc_mev_fm;
      pot->Cu=-84.53/o2scl_const::hc_mev_fm;
      pot->Cl=-65.472/o2scl_const::hc_mev_fm;
      pot->sigma=7.0/6.0;
      pot->x=1.0/15.0;
      pot->form=pot->gbd_form;
    } else if (name=="PAL11") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=1;
      pot->form=pot->pal_form;
    } else if (name=="PAL12") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=2;
      pot->form=pot->pal_form;
    } else if (name=="PAL13") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=3;
      pot->form=pot->pal_form;
    } else if (name=="PAL21") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=1;
      pot->form=pot->pal_form;
    } else if (name=="PAL22") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=2;
      pot->form=pot->pal_form;
    } else if (name=="PAL23") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=3;
      pot->form=pot->pal_form;
    } else if (name=="PAL31") {
      pot->A=-46.65/o2scl_const::hc_mev_fm;
      pot->B=39.54/o2scl_const::hc_mev_fm;
      pot->Bp=0.3;
      pot->sigma=1.663;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=1;
      pot->form=pot->pal_form;
    } else if (name=="PAL32") {
      pot->A=-46.65/o2scl_const::hc_mev_fm;
      pot->B=39.54/o2scl_const::hc_mev_fm;
      pot->Bp=0.3;
      pot->sigma=1.663;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=2;
      pot->form=pot->pal_form;
    } else if (name=="PAL33") {
      pot->A=-46.65/o2scl_const::hc_mev_fm;
      pot->B=39.54/o2scl_const::hc_mev_fm;
      pot->Bp=0.3;
      pot->sigma=1.663;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.0;
      pot->x3=0.0;
      pot->z1=0.0;
      pot->z2=0.0;
      pot->sym_index=3;
      pot->form=pot->pal_form;
    } else if (name=="BPALb11") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=-1.361;
      pot->x3=-0.244;
      pot->z1=-13.91/o2scl_const::hc_mev_fm;
      pot->z2=16.69/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb12") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=-1.361;
      pot->x3=-0.244;
      pot->z1=-13.91/o2scl_const::hc_mev_fm;
      pot->z2=16.69/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb13") {
      pot->A=75.94/o2scl_const::hc_mev_fm;
      pot->B=-30.880/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.498;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=-1.903;
      pot->x3=-1.056;
      pot->z1=-1.83/o2scl_const::hc_mev_fm;
      pot->z2=5.09/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb21") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.086;
      pot->x3=0.561;
      pot->z1=-18.4/o2scl_const::hc_mev_fm;
      pot->z2=46.27/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb22") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.086;
      pot->x3=0.561;
      pot->z1=-18.4/o2scl_const::hc_mev_fm;
      pot->z2=46.27/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb23") {
      pot->A=440.94/o2scl_const::hc_mev_fm;
      pot->B=-213.41/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.927;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.086;
      pot->x3=0.561;
      pot->z1=-18.4/o2scl_const::hc_mev_fm;
      pot->z2=46.27/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="BPALb31") {
      pot->A=-46.65/o2scl_const::hc_mev_fm;
      pot->B=39.45/o2scl_const::hc_mev_fm;
      pot->Bp=0.3;
      pot->sigma=1.663;
      pot->C1=-83.84/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=0.376;
      pot->x3=0.246;
      pot->z1=-12.23/o2scl_const::hc_mev_fm;
      pot->z2=-2.98/o2scl_const::hc_mev_fm;
      pot->form=pot->bpal_form;
    } else if (name=="SL12") {
      pot->A=3.706/o2scl_const::hc_mev_fm;
      pot->B=-31.155/o2scl_const::hc_mev_fm;
      pot->Bp=0.0;
      pot->sigma=0.453;
      pot->C1=-41.28/o2scl_const::hc_mev_fm;
      pot->C2=23.0/o2scl_const::hc_mev_fm;
      pot->x0=-3.548;
      pot->x3=-0.5;
      pot->z1=-13.355/o2scl_const::hc_mev_fm;
      pot->z2=2.789/o2scl_const::hc_mev_fm;
      pot->form=pot->sl_form;
    } else if (name=="GBD0") {
      pot->Au=-109.85/o2scl_const::hc_mev_fm;
      pot->Al=-191.30/o2scl_const::hc_mev_fm;
      pot->B=205.66/o2scl_const::hc_mev_fm;
      pot->Cu=-118.80/o2scl_const::hc_mev_fm;
      pot->Cl=-26.26/o2scl_const::hc_mev_fm;
      pot->sigma=7.0/6.0;
      pot->x=0.0;
      pot->form=pot->gbd_form;
    } else if (name=="GBD1") {
      pot->Au=-299.69/o2scl_const::hc_mev_fm;
      pot->Al=-1.46/o2scl_const::hc_mev_fm;
      pot->B=205.66/o2scl_const::hc_mev_fm;
      pot->Cu=-118.80/o2scl_const::hc_mev_fm;
      pot->Cl=-26.26/o2scl_const::hc_mev_fm;
      pot->sigma=7.0/6.0;
      pot->x=1.0;
    } else if (name=="CKLxm2") {
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->form=pot->mdi_form;
      pot->x=-2.0;
      pot->Au=-95.98/o2scl_const::hc_mev_fm-2.0*pot->B*
	pot->x/(pot->sigma+1.0);
      pot->Al=-120.75/o2scl_const::hc_mev_fm+2.0*pot->B*
	pot->x/(pot->sigma+1.0);
    } else if (name=="CKLxm1") {
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->form=pot->mdi_form;
      pot->x=-2.0;
      pot->Au=-95.98/o2scl_const::hc_mev_fm-2.0*pot->B*
	pot->x/(pot->sigma+1.0);
      pot->Al=-120.75/o2scl_const::hc_mev_fm+2.0*pot->B*
	pot->x/(pot->sigma+1.0);
    } else if (name=="CKLx0") {
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->form=pot->mdi_form;
      pot->x=-1.0;
      pot->Au=-95.98/o2scl_const::hc_mev_fm-2.0*pot->B*
	pot->x/(pot->sigma+1.0);
      pot->Al=-120.75/o2scl_const::hc_mev_fm+2.0*pot->B*
	pot->x/(pot->sigma+1.0);
    } else if (name=="CKLx1") {
      pot->B=106.35/o2scl_const::hc_mev_fm;
      pot->sigma=4.0/3.0;
      pot->Cu=-103.40/o2scl_const::hc_mev_fm;
      pot->Cl=-11.70/o2scl_const::hc_mev_fm;
      pot->form=pot->mdi_form;
      pot->x=1.0;
      pot->Au=-95.98/o2scl_const::hc_mev_fm-2.0*pot->B*
	pot->x/(pot->sigma+1.0);
      pot->Al=-120.75/o2scl_const::hc_mev_fm+2.0*pot->B*
	pot->x/(pot->sigma+1.0);
    }
    pot->rho0=0.16;
    pot->Lambda=1.5*cbrt(1.5*o2scl_const::pi2*pot->rho0);
    pot->Lambda2=3.0*cbrt(1.5*o2scl_const::pi2*pot->rho0);
    return pot;
  } else if (type=="rmf") {
    eos_had_rmf *rmf=new eos_had_rmf;
    rmf_load(*rmf,name);
    return rmf;
  }
  O2SCL_ERR("Type not understood in eos_had_strings().",exc_einval);
  return 0;
}
