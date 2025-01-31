/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2021-2025, Andrew W. Steiner
  
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
#include <o2scl/vec_stats.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/nucleus_bin.h>
#include <o2scl/eos_had_skyrme.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

nucleus_bin::nucleus_bin() {

  verbose=1;
  precision=4;

  //std::cout << "Reading nuclear mass tables." << std::endl;
    
  ame20exp.load("20");
  ame20round.load("20round");
  
  o2scl_hdf::mnmsk_load(m16,"msis16");
    
  o2scl_hdf::hfb_sp_load(hfb21,21);
  o2scl_hdf::hfb_sp_load(hfb22,22);
  o2scl_hdf::hfb_sp_load(hfb23,23);
  o2scl_hdf::hfb_sp_load(hfb24,24);
  o2scl_hdf::hfb_sp_load(hfb25,25);
  o2scl_hdf::hfb_sp_load(hfb26,26);
  o2scl_hdf::hfb_sp_load(hfb27,27);

  wlw1.load("WS3.2");
  wlw2.load("WS3.3");
  wlw3.load("WS3.6");
  wlw4.load("WS3_RBF");
  wlw5.load("WS4_RBF");
    
  ddme2.load_be("ddme2","E",1.0);
  ddmed.load_be("ddmed","E",1.0);
  ddpc1.load_be("ddpc1","E",1.0);
  nl3s.load_be("nl3s","E",1.0);
  sly4.load_be("sly4_all","Binding_Energy__MeV_",1.0);
  skms.load_be("skms_all","Binding_Energy__MeV_",1.0);
  skp.load_be("skp_all","Binding_Energy__MeV_",1.0);
  sv_min.load_be("sv-min_all","Binding_Energy__MeV_",1.0);
  unedf0.load_be("unedf0_all","Binding_Energy__MeV_",1.0);
  unedf1.load_be("unedf1_all","Binding_Energy__MeV_",1.0);
  
  //std::cout << "Done reading nuclear mass tables." << std::endl;

  nmd={&ame20exp,&ame20round,&m16,
       &hfb21,&hfb22,&hfb23,
       &hfb24,&hfb25,&hfb26,&hfb27,
       &wlw1,&wlw2,&wlw3,&wlw4,&wlw5,
       &dz,&ddme2,&ddmed,&ddpc1,&nl3s,&sly4,&skms,
       &skp,&sv_min,&unedf0,&unedf1};
       
  n_tables=nmd.size();

  table_names={"AME exp 20","AME rnd 20","MSIS 16",
	       "HFB21","HFB22","HFB23",
	       "HFB24","HFB25","HFB26","HFB27",
	       "WLW 10","WLLW 10","LWDW 11","WL 11","WLWM 14",
	       "DZ 95","DDME2","DDMED","DDPC1","NL3S",
               "SLy4","SKM*","SkP","SV-min","UNEDF0","UNEDF1"};
  if (table_names.size()!=nmd.size()) {
    O2SCL_ERR2("Table list size problem in ",
              "nucleus_bin::nucleus_bin().",o2scl::exc_esanity);
  }

  nmfd={&se,&frdm,&dzf,&dzf33,&frdm_shell,&ldrop_shell};
  n_fits=nmfd.size();
    
  fit_names={"Semi-empirical","FRDM","DZ fit 10","DZ fit 33","FRDM shell",
             "Liq drop shell"};

  o2scl_hdf::skyrme_load(sk,"SLy4");
  ldrop_shell.set_eos_had_temp_base(sk);
  
  if (true) {
    ubvector p(5);
    p[0]=1.534012717750970e+01;
    p[1]=2.261032660770839e+01;
    p[2]=1.617898540187606e+01;
    p[3]=6.941932590175626e-01;
    p[4]=1.162846797031895e+01;
    se.fit_fun(5,p);
  }
  if (true) {
    ubvector p(10);
    p[0]=3.613692309944690e-01;
    p[1]=9.216607474662888e-01;
    p[2]=5.170200810075372e+01;
    p[3]=1.947090687326399e+01;
    p[4]=4.084492571046316e+01;
    p[5]=4.380015938017932e+01;
    p[6]=2.092293766732548e+01;
    p[7]=9.025789227326884e-01;
    p[8]=1.618679100791941e+00;
    p[9]=2.905166877962070e-01;
    frdm.fit_fun(10,p);
  }
  if (true) {
    ubvector p(10);
    p[0]=7.067453190240083e-01;
    p[1]=1.777839083682484e+01;
    p[2]=1.633713770563290e+01;
    p[3]=3.791636389950607e+01;
    p[4]=5.483429204783442e+01;
    p[5]=3.496320427597641e-01;
    p[6]=1.561301751231576e+00;
    p[7]=1.942514988061515e-02;
    p[8]=4.392327118808018e+01;
    p[9]=6.218195946263423e+00;
    dzf.fit_fun(10,p);
  }
  if (true) {
    ubvector p(33);
    p[0]=9.088827842830527e+00;
    p[1]=6.511564581987122e+00;
    p[2]=4.462956108688138e+00;
    p[3]=2.074222936534796e+01;
    p[4]=1.754662224554880e+00;
    p[5]=7.740340448804256e+00;
    p[6]=-4.369191221566433e+00;
    p[7]=-3.418368107217408e+01;
    p[8]=-3.578727861463152e-01;
    p[9]=7.269557856660328e-01;
    p[10]=-7.501870792510353e-01;
    p[11]=-3.771736481691423e+00;
    p[12]=-1.779481045580030e-01;
    p[13]=-9.057241132560844e-01;
    p[14]=3.952971227188688e-01;
    p[15]=1.800359380599373e+00;
    p[16]=2.357449926816584e-01;
    p[17]=1.051823513035010e+00;
    p[18]=8.973655324903573e+00;
    p[19]=5.612237450969798e+01;
    p[20]=1.824848865708848e+01;
    p[21]=7.426985283747025e+01;
    p[22]=-2.690428328513786e+01;
    p[23]=-1.280799866769777e+02;
    p[24]=-4.169614442419837e+00;
    p[25]=-2.917812343383905e+01;
    p[26]=-3.792308695125627e+01;
    p[27]=-5.385348566074595e+01;
    p[28]=1.582609566178451e+00;
    p[29]=5.568450543337717e+00;
    p[30]=7.055984771856043e-01;
    p[31]=6.196805198023444e+00;
    p[32]=1.997170958626300e+01;
    dzf33.fit_fun(33,p);
  }
  if (true) {
    ubvector p(14);
    p[0]=2.862652716272743e+00;
    p[1]=1.129807927418693e+00;
    p[2]=3.559392452233342e+01;
    p[3]=1.657931612140274e+01;
    p[4]=2.369337715147342e+01;
    p[5]=3.468919264623499e+01;
    p[6]=2.635234363179206e+01;
    p[7]=4.266240233807385e-01;
    p[8]=1.257558181549950e+03;
    p[9]=3.250151576606326e+00;
    p[10]=-1.460177890508211e+00;
    p[11]=1.961924905272880e-02;
    p[12]=2.010431359292397e-03;
    p[13]=8.873157422094327e-02;
    frdm_shell.fit_fun(14,p);
  }
  if (true) {
    ubvector p(11);
    p[0]=8.993854058157809e-01;
    p[1]=9.692983376594209e-01;
    p[2]=9.129367034093190e-01;
    p[3]=9.709856845570888e-01;
    p[4]=-1.004780591364136e-02;
    p[5]=1.392239609972227e-01;
    p[6]=1.096713554301717e+01;
    p[7]=-1.498528880855952e+00;
    p[8]=1.433936004772412e-02;
    p[9]=1.666488914701955e-03;
    p[10]=1.137565990431053e-01;
    ldrop_shell.fit_fun(11,p);
  }

  older_tables=false;
}

void nucleus_bin::update_older_tables() {

  if (older_tables==true && ame16.get_nentries()==0) {

    ame95rmd.load("95rmd");
    ame95exp.load("95exp");
    ame03round.load("03round");
    ame03.load("03");
    ame12.load("12");
    ame16.load("16");

    o2scl_hdf::mnmsk_load(m95,"mnmsk97");
  
    o2scl_hdf::hfb_load(hfb2,2);
    o2scl_hdf::hfb_load(hfb8,8);
    o2scl_hdf::hfb_load(hfb14,14);
    o2scl_hdf::hfb_load(hfb14_v0,15);
    o2scl_hdf::hfb_sp_load(hfb17,17);
    
    kt.load("04");
    kt2.load("05");
    
    sdnp1.load("sdnp03");
    sdnp2.load("sd_skp_04");
    sdnp3.load("sd_sly4_04");

    nmd.push_back(&ame95rmd);
    nmd.push_back(&ame95exp);
    nmd.push_back(&ame03round);
    nmd.push_back(&ame03);
    nmd.push_back(&ame12);
    nmd.push_back(&ame16);

    nmd.push_back(&m95);
    
    nmd.push_back(&kt);
    nmd.push_back(&kt2);
    
    nmd.push_back(&sdnp1);
    nmd.push_back(&sdnp2);
    nmd.push_back(&sdnp3);
    
    nmd.push_back(&hfb2);
    nmd.push_back(&hfb8);
    nmd.push_back(&hfb14);
    nmd.push_back(&hfb14_v0);
    nmd.push_back(&hfb17);
    
    table_names.push_back("AME rmd 95");
    table_names.push_back("AME 95 exp");
    table_names.push_back("AME rnd 03");
    table_names.push_back("AME 03");
    table_names.push_back("AME 12");
    table_names.push_back("AME 16");
    
    table_names.push_back("MNMSK 95");
    
    table_names.push_back("KTUY 04");
    table_names.push_back("KTUY 05");
    
    table_names.push_back("SDNP1 05");
    table_names.push_back("SDNP2 05");
    table_names.push_back("SDNP3 05");
    
    table_names.push_back("HFB2");
    table_names.push_back("HFB8");
    table_names.push_back("HFB14");
    table_names.push_back("HFB14_v0");
    table_names.push_back("HFB17");

  }

  if (older_tables==false && ame16.get_nentries()>0) {

    ame95rmd.clear();
    ame95exp.clear();
    ame03round.clear();
    ame03.clear();
    ame12.clear();
    ame16.clear();

    m95.clear();
    
    hfb2.clear();
    hfb8.clear();
    hfb14.clear();
    hfb14_v0.clear();
    hfb17.clear();
    
    kt.clear();
    kt2.clear();
    
    sdnp1.clear();
    sdnp2.clear();
    sdnp3.clear();

    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame95rmd),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame95exp),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame03round),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame03),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame12),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&ame16),nmd.end());
    
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&m95),nmd.end());
    
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&kt),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&kt2),nmd.end());
    
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&sdnp1),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&sdnp2),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&sdnp3),nmd.end());
    
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&hfb2),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&hfb8),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&hfb14),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&hfb14_v0),nmd.end());
    nmd.erase(std::remove(nmd.begin(),nmd.end(),&hfb17),nmd.end());

    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME rmd 95"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME 95 exp"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME rnd 03"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME 03"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME 12"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "AME 16"),table_names.end());
    
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "MNMSK 95"),table_names.end());
    
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "KTUY 04"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "KTUY 05"),table_names.end());
    
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "SDNP1 05"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "SDNP2 05"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "SDNP3 05"),table_names.end());
    
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "HFB2"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "HFB8"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "HFB14"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "HFB14_v0"),table_names.end());
    table_names.erase(std::remove(table_names.begin(),
                                  table_names.end(),
                                  "HFB17"),table_names.end());
  }

  n_tables=nmd.size();
  
  if (table_names.size()!=nmd.size()) {
    O2SCL_ERR2("Table list size problem in ",
              "nucleus_bin::nucleus_bin().",o2scl::exc_esanity);
  }
  
  return;
}

int nucleus_bin::get(std::vector<std::string> &sv, bool itive_com) {

  update_older_tables();
  
  if (sv.size()<3) {
    cerr << "Get function needs Z and N." << endl;
    return 1;
  }

  size_t Z=o2scl::stoszt(sv[1]);
  size_t N=o2scl::stoszt(sv[2]);

  cout << "Z=" << Z << " N=" << N << " A=" << Z+N
       << " symbol=" << nmi.Ztoel(Z)
       << " name=" << nmi.Ztoname(Z) << endl;
  cout << endl;

  size_t left_column=14;
    
  cout.width(left_column+1);
  cout << "Model ";
  cout.width(6+precision);
  cout << "mass ex. MeV ";
  cout.width(6+precision);
  cout << "BE/A MeV " << "  ";
  cout.width(left_column+1);
  cout << "Model  ";
  cout.width(6+precision);
  cout << "mass ex. MeV ";
  cout.width(6+precision);
  cout << "BE/A MeV " << endl;
  cout << endl;
    
  nucleus nuc;

  cout.precision(precision);
  cout.setf(ios::showpos);
  
  int n_out=0, last_endl=0;
  for(size_t i=0;i<n_tables;i+=2) {
    if (nmd[i]->is_included(Z,N)) {
      int ret=nmd[i]->get_nucleus(Z,N,nuc);
      cout.width(left_column);
      cout << table_names[i] << " "
	   << nuc.mex*o2scl_const::hc_mev_fm << " "
	   << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << "  ";
      n_out++;
    }
    // End line after every pair of outputs
    if (n_out-last_endl==2) {
      cout << endl;
      last_endl=n_out;
    }
    if (i+1<n_tables) {
      if (nmd[i+1]->is_included(Z,N)) {
        int ret=nmd[i+1]->get_nucleus(Z,N,nuc);
        cout.width(left_column);
        cout << table_names[i+1] << " "
             << nuc.mex*o2scl_const::hc_mev_fm << " "
             << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << "  ";
        n_out++;
      }
    }
    // End line after every pair of outputs
    if (n_out-last_endl==2) {
      cout << endl;
      last_endl=n_out;
    }
  }
  // Additional end line if there is a left-over output
  if (n_out>last_endl) cout << endl;
  // No extra end line if no tables had results
  if (n_out>0) {
    cout << endl;
  }

  n_out=0;
  last_endl=0;
  for(size_t i=0;i<n_fits;i+=2) {
    if (nmfd[i]->is_included(Z,N)) {
      int ret=nmfd[i]->get_nucleus(Z,N,nuc);
      cout.width(left_column);
      cout << fit_names[i] << " "
	   << nuc.mex*o2scl_const::hc_mev_fm << " "
	   << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << "  ";
      n_out++;
    }
    // End line after every pair of outputs
    if (n_out-last_endl==2) {
      cout << endl;
      last_endl=n_out;
    }
    if (i+1<n_fits) {
      if (nmfd[i+1]->is_included(Z,N)) {
        int ret=nmfd[i+1]->get_nucleus(Z,N,nuc);
        cout.width(left_column);
        cout << fit_names[i+1] << " "
             << nuc.mex*o2scl_const::hc_mev_fm << " "
             << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << "  ";
        n_out++;
      }
    }
    // End line after every pair of outputs
    if (n_out-last_endl==2) {
      cout << endl;
      last_endl=n_out;
    }
  }
  // Additional end line if there is a left-over output
  if (n_out>last_endl) cout << endl;
  // No extra end line if no tables had results
  if (n_out>0) {
    cout << endl;
  }
  
  cout.unsetf(ios::showpos);

  if (ame20exp.is_included(Z,N)) {
    nucmass_ame::entry en=ame20exp.get_ZN(Z,N);
    cout << "From AME 2020:" << endl;
    cout << "data origin: "
         << ((std::string)&en.orig[0]) << std::endl;
    cout << "mass excess unc.: " << en.dmass << " keV" << endl;
    cout << "beta decay energy: " << en.bde << " keV" << endl;
    cout << "beta decay energy unc.: " << en.dbde << " keV" << endl;
    cout << "beta-decay mode: " 
         << ((std::string)&en.bdmode[0]) << std::endl;
    cout << "atomic mass: " << en.amass << " keV" << endl;
    cout << "atomic mass unc.: " << en.damass << " keV" << endl;
    cout << endl;
  }

  if (hfb27.is_included(Z,N)) {
    nucmass_hfb_sp::entry en=hfb27.get_ZN(Z,N);
    if (en.Jexp>1.0e98 && en.Pexp>98) {
      cout << "HFB27 spin and parity, theory: "
           << en.Jth << " " << en.Pth << endl;
    } else {
      cout << "HFB27 spin and parity, expt.: "
           << en.Jexp << " " << en.Pexp << ", theory: "
           << en.Jth << " " << en.Pth << endl;
    }
    cout << endl;
  }
  
  if (ame20exp.is_included(Z,N)) {
    nucmass_ame::entry en=ame20exp.get_ZN(Z,N);
    cout.setf(ios::left);
    cout << "Nubase entries:" << endl;
    cout << "Z i ";
    cout.width(precision+8);
    cout << " mass ";
    cout.width(precision+8);
    cout << " dmass ";
    cout << "a ";
    cout.width(precision+8);
    cout << " exc_energy ";
    cout.width(precision+8);
    cout << " dexc_en. ";
    cout << "a ";
    cout << "ori";
    cout << "iu";
    cout << "ii";
    cout.width(precision+8);
    cout << " hlife ";
    cout << "a ";
    cout << "hun ";
    cout.width(9);
    cout << "dhlife ";
    cout.width(16);
    cout << "spin parity ";
    cout << "ENyr ";
    cout << "disc ";
    cout << "decay intensity" << endl;
    
    for(size_t j=0;j<en.props.size();j++) {
      cout << en.props[j].Znote << " ";
      if (en.props[j].isomer=='\0') {
        cout << "_ ";
      } else {
        cout << en.props[j].isomer << " ";
      }
      cout.setf(ios::showpos);
      cout << en.props[j].mass << " ";
      cout << en.props[j].dmass << " ";
      cout.unsetf(ios::showpos);
      cout << en.props[j].mass_acc << " ";
      cout.setf(ios::showpos);
      cout << en.props[j].exc_energy << " ";
      cout << en.props[j].dexc_energy << " ";
      cout.unsetf(ios::showpos);
      cout << en.props[j].exc_energy_acc << " ";
      cout.width(2);
      cout << ((std::string)(&(en.props[j].origin[0]))) << " ";
      if (en.props[j].isomer_unc=='\0') {
        cout << "_ ";
      } else {
        cout << en.props[j].isomer_unc << " ";
      }
      if (en.props[j].isomer_inv=='\0') {
        cout << "_ ";
      } else {
        cout << en.props[j].isomer_inv << " ";
      }
      cout.setf(ios::showpos);
      cout << en.props[j].hlife << " ";
      cout.unsetf(ios::showpos);
      cout << en.props[j].hlife_acc << " ";
      cout.width(3);
      cout << ((std::string)(&(en.props[j].hl_unit[0]))) << " ";
      cout.width(8);
      cout << ((std::string)(&(en.props[j].dhlife[0]))) << " ";
      cout.width(15);
      cout << ((std::string)(&(en.props[j].spinp[0]))) << " ";
      cout.width(4);
      cout << en.props[j].ENSDF_year << " ";
      cout.width(4);
      cout << en.props[j].discovery << " ";
      if (en.props[j].decay_intensity[0]=='\0') {
        cout << "_" << endl;
      } else {
        cout << ((std::string)(&(en.props[j].decay_intensity[0]))) << endl;
      }
    }
    cout.unsetf(ios::left);
    
  }
  
  cout.precision(6);
    
  return 0;
}

int nucleus_bin::tables(std::vector<std::string> &sv, bool itive_com) {

  update_older_tables();
  
  size_t left_column=18;
    
  cout << "Number of entries in nuclear mass tables: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    cout.width(2);
    cout << i << " ";
    cout.width(left_column);
    cout << table_names[i] << ": " << nmd[i]->get_nentries() << endl;
  }
    
  return 0;
}

int nucleus_bin::refs(std::vector<std::string> &sv, bool itive_com) {

  update_older_tables();
  
  size_t left_column=18;

  cout << "References: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    cout.width(left_column);
    cout << table_names[i] << ": " << nmd[i]->reference << endl;
  }
    
  return 0;
}

int nucleus_bin::cdist(std::vector<std::string> &sv, bool itive_com) {

  update_older_tables();
  
  // Set a large distribution
  nucdist_set(moller_dist,m95);

  // Find all nucleus in all the tables
  size_t min_N=400, min_Z=400, max_N=0, max_Z=0;
  for(size_t i=0;i<moller_dist.size();i++) {
    bool included=true;
    // Skip the experimental tables and some theory tables
    // which have only very limited data
    for(size_t j=7;j<n_tables;j++) {
      if (j!=25 && j!=26 && j!=27) {
	if (nmd[j]->is_included(moller_dist[i].Z,
				moller_dist[i].N)==false) {
	  included=false;
	}
      }
    }
    if (included==true) {
      common_dist.push_back(moller_dist[i]);
      if (moller_dist[i].Z<((int)min_Z)) min_Z=moller_dist[i].Z;
      if (moller_dist[i].N<((int)min_N)) min_N=moller_dist[i].N;
      if (moller_dist[i].Z>((int)max_Z)) max_Z=moller_dist[i].Z;
      if (moller_dist[i].N>((int)max_N)) max_N=moller_dist[i].N;
    }
  }
  cout << "Size of common distribution: " << common_dist.size() << endl;

  // Create a table object which has a list of the common
  // distribution, along with the standard deviation in
  // the mass excess and the standard deviation divided by
  // mass number
  
  cout << "Storing in table." << endl;
  table<> tab;
  tab.line_of_names("N Z sd sdoa");
  
  vector<double> mass_temp;
  for(size_t i=0;i<common_dist.size();i++) {
    mass_temp.clear();
    for(size_t j=7;j<n_tables;j++) {
      if (j!=25 && j!=26 && j!=27) {
	mass_temp.push_back(nmd[j]->mass_excess(common_dist[i].Z,
						common_dist[i].N));
      }
    }
    double stddev=o2scl::vector_stddev(mass_temp.size(),mass_temp);
    if (i%1000==0) {
      cout.width(4);
      cout << i << " ";
      cout.width(3);
      cout << common_dist[i].N << " ";
      cout.width(2);
      cout << common_dist[i].Z << " " << stddev << endl;
    }
    if (fabs(stddev)<1.0e8 && (common_dist[i].Z!=21 ||
			       common_dist[i].N!=15)) {
      double line[4]={((double)common_dist[i].N),
		      ((double)common_dist[i].Z),
		      stddev,stddev/common_dist[i].A};
      tab.line_of_data(4,line);
    }
  }

  // Convert the table to a table3d object for plotting
  table3d t3d;
  t3d.read_table(tab,"N","Z",0.0,0);

  // Output to file
  hdf_file hf;
  hf.open_or_create("cdist.o2");
  hdf_output(hf,(const table3d &)t3d,"t3d");
  hdf_output(hf,tab,"tab");
  hf.close();
  
  return 0;
}

int nucleus_bin::fit_method(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()>=2) {
    size_t method=o2scl::stoszt(sv[1]);
    if (method<7) {
      fitter.fit_method=method;
    }
  }
  return 0;
}
  
int nucleus_bin::compare(std::vector<std::string> &sv, bool itive_com) {

  update_older_tables();
  
  size_t left_column=18;

  double res;

  cout << "Tables: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    nucdist_pair_set(fitter.dist,ame20exp,*nmd[i]);
    cout.width(left_column);
    fitter.eval(*nmd[i],res);
    cout << table_names[i] << ": ";
    cout.width(4);
    cout << fitter.dist.size() << " " << res << endl;
  }
  
  nucdist_set(fitter.dist,ame20exp);
  
  cout << "Formulas: " << endl;
  for(size_t i=0;i<n_fits;i++) {
    cout.width(left_column);
    fitter.eval(*nmfd[i],res);
    cout << fit_names[i] << ": " << res << endl;
  }
  cout << endl;

  return 0;
}

int nucleus_bin::fit(std::vector<std::string> &sv, bool itive_com) {

  fitter.def_mmin.ntrial*=100;
  
  double res;
  kwargs kw;
  if (sv.size()>=3) {
    kw.set(sv[2]);
  }
  fitter.def_mmin.tol_abs=kw.get_double("tol_abs",1.0e-4);
  
  nucdist_set(fitter.dist,ame20exp);
  ldrop_shell.large_vals_unphys=true;

  fitter.def_mmin.verbose=kw.get_int("verbose",0);

  size_t ix_start=0;
  size_t ix_end=n_fits;
  if (sv.size()>=2) {
    ix_start=o2scl::stoszt(sv[1]);
    ix_end=ix_start+1;
  }
  for(size_t i=ix_start;i<ix_end;i++) {
    fitter.eval(*nmfd[i],res);
    cout << fit_names[i] << ": Before fit: " << res << endl;
    fitter.fit(*nmfd[i],res);
    fitter.eval(*nmfd[i],res);
    cout << fit_names[i] << ": After fit: " << res << endl;
    ubvector p(nmfd[i]->nfit);
    nmfd[i]->guess_fun(nmfd[i]->nfit,p);
    for(size_t k=0;k<nmfd[i]->nfit;k++) {
      cout << "    p[" << k << "]=" << dtos(p[k],0) << ";" << endl;
    }
  }
    
  return 0;
}

void nucleus_bin::setup_cli(o2scl::cli &cl) {

  static const int nopt=7;
  o2scl::comm_option_s options[nopt]={
    {0,"ZN","Information for a nucleus given Z and N.",
     2,2,"<Z> <N>",((std::string)"The 'ZN' command outputs ")+
     "the binding energy for a specified nucleus for all "+
     "tables and models.",new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::get),o2scl::cli::comm_option_both},
    {0,"tables","List available tables and number of nuclei",
     0,0,"",((std::string)"The 'tables' command lists all ")+
     "available tables (both experimental and theoretical).",
     new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::tables),o2scl::cli::comm_option_both},
    {0,"compare","Compare theoretical mass models to experiment",
     0,0,"",((std::string)"The 'compare' command compares all of the ")+
     "mass tables to experiment (currently AME 2020).",
     new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::compare),o2scl::cli::comm_option_both},
    {0,"fit","Fit theoretical mass model(s) to experiment",
     0,2,"[index] [kwargs]",
     ((std::string)"The 'fit' command adjusts the fit parameters ")+
     "for one or all of the models for which these parameters can be "+
     "varied and optimizes their fit with experiment. This takes "+
     "awhile, likely at least an hour.",
     new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::fit),o2scl::cli::comm_option_both},
    {0,"cdist","Create a distribution of nuclei common to several tables",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::cdist),o2scl::cli::comm_option_both},
    {0,"refs","List the references for all of the tables and models",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::refs),o2scl::cli::comm_option_both},
    {0,"fit-method","Desc",
     1,1,"","",new o2scl::comm_option_mfptr<nucleus_bin>
     (this,&nucleus_bin::fit_method),o2scl::cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);
    
  p_verbose.i=&verbose;
  p_verbose.help="Verbosity parameter (default 1)";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  p_precision.i=&precision;
  p_precision.help="Precision parameter (default 4)";
  cl.par_list.insert(make_pair("precision",&p_precision));

  p_older_tables.b=&older_tables;
  p_older_tables.help="Older tables (default false)";
  cl.par_list.insert(make_pair("older_tables",&p_older_tables));

  return;
}
