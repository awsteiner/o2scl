/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2021-2024, Andrew W. Steiner
  
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
#include "nucleus.h"

#include <o2scl/vec_stats.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

nucleus_class::nucleus_class() {

  verbose=1;
    
  std::cout << "Reading nuclear mass tables." << std::endl;
    
  o2scl_hdf::ame_load(ame95rmd,"95rmd");
  o2scl_hdf::ame_load(ame95exp,"95exp");
  o2scl_hdf::ame_load(ame03round,"03round");
  o2scl_hdf::ame_load(ame03,"03");
  o2scl_hdf::ame_load(ame12,"12");
  o2scl_hdf::ame_load(ame16,"16");
  o2scl_hdf::ame_load(ame20exp,"20");
  o2scl_hdf::ame_load(ame20round,"20round");
    
  o2scl_hdf::mnmsk_load(m95);
    
  o2scl_hdf::hfb_load(hfb2,2);
  o2scl_hdf::hfb_load(hfb8,8);
  o2scl_hdf::hfb_load(hfb14,14);
  o2scl_hdf::hfb_load(hfb14_v0,15);
  o2scl_hdf::hfb_sp_load(hfb17,17);
  o2scl_hdf::hfb_sp_load(hfb21,21);
  o2scl_hdf::hfb_sp_load(hfb22,22);
  o2scl_hdf::hfb_sp_load(hfb23,23);
  o2scl_hdf::hfb_sp_load(hfb24,24);
  o2scl_hdf::hfb_sp_load(hfb25,25);
  o2scl_hdf::hfb_sp_load(hfb26,26);
  o2scl_hdf::hfb_sp_load(hfb27,27);

  kt.load("04");
  kt2.load("05");

  wlw1.load("WS3.2");
  wlw2.load("WS3.3");
  wlw3.load("WS3.6");
    
  sdnp1.load("sdnp03");
  sdnp2.load("sd_skp_04");
  sdnp3.load("sd_sly4_04");
  std::cout << "Done reading nuclear mass tables." << std::endl;

  nmd={&ame95rmd,&ame95exp,&ame03round,&ame03,
       &ame12,&ame16,&ame20exp,&ame20round,&m95,&kt,
       &kt2,&hfb2,&hfb8,&hfb14,&hfb14_v0,
       &hfb17,&hfb21,&hfb22,&hfb23,
       &hfb24,&hfb25,&hfb26,&hfb27,
       &wlw1,&wlw2,&wlw3,&sdnp1,
       &sdnp2,&sdnp3,&dz,&ddme2,&ddmed,&ddpc1,
       &nl3s,&sly4,&skms,&skp,&sv_min,&unedf0,&unedf1};
       
  n_tables=nmd.size();

  table_names={"AME (1995) rmd.","AME (1995) exp.",
	       "AME (2003) round.","AME (2003)",
	       "AME (2012)","AME (2016)","AME (2020) exp.",
               "AME (2020) round.",
	       "MNMSK (1995)","KTUY (2004)","KTUY (2005)",
	       "HFB2","HFB8","HFB14","HFB14_v0",
	       "HFB17","HFB21","HFB22","HFB23",
	       "HFB24","HFB25","HFB26","HFB27",
	       "WLW1 (2005)","WLW2 (2005)","WLW3 (2005)",
	       "SDNP1 (2005)","SDNP2 (2005)","SDNP3 (2005)",
	       "DZ (1995)","DDME2","DDMED","DDPC1","NL3S",
               "SLy4","SKM*","SkP","SV-min","UNEDF0","UNEDF1"};

  nmfd={&se,&frdm,&dzf,&dzf33};
  n_fits=nmfd.size();
    
  fit_names={"Semi-empirical","FRDM","DZ fit 10","DZ fit 33"};
}

int nucleus_class::get(std::vector<std::string> &sv, bool itive_com) {

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

  size_t left_column=18;
    
  cout.width(left_column+2);
  cout << "Model  ";
  cout.width(10);
  cout << "mass ex. MeV ";
  cout.width(10);
  cout << "BE/A MeV " << endl;
  cout << endl;
    
  nucleus nuc;

  for(size_t i=0;i<n_tables;i++) {
    if (nmd[i]->is_included(Z,N)) {
      int ret=nmd[i]->get_nucleus(Z,N,nuc);
      cout.width(left_column);
      cout << table_names[i] << ": "
	   << nuc.mex*o2scl_const::hc_mev_fm << " "
	   << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << endl;
    }
  }
  cout << endl;
    
  for(size_t i=0;i<4;i++) {
    if (nmfd[i]->is_included(Z,N)) {
      int ret=nmfd[i]->get_nucleus(Z,N,nuc);
      cout.width(left_column);
      cout << fit_names[i] << ": "
	   << nuc.mex*o2scl_const::hc_mev_fm << " "
	   << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << endl;
    }
  }
    
  return 0;
}

int nucleus_class::tables(std::vector<std::string> &sv, bool itive_com) {

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

int nucleus_class::refs(std::vector<std::string> &sv, bool itive_com) {

  size_t left_column=18;

  cout << "References: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    cout.width(left_column);
    cout << table_names[i] << ": " << nmd[i]->reference << endl;
  }
    
  return 0;
}

int nucleus_class::cdist(std::vector<std::string> &sv, bool itive_com) {

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

int nucleus_class::fits(std::vector<std::string> &sv, bool itive_com) {

  size_t left_column=18;

  nucdist_set(fitter.dist,ame16);
  fitter.def_mmin.ntrial*=10;

  cout << "Tables: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    if (i!=12 && i!=13 && i!=14 && i!=15 && i!=16 && i!=17 &&
	i!=18 && i!=25 && i!=26 && i!=27) {
      cout.width(left_column);
      double res;
      fitter.eval(*nmd[i],res);
      cout << table_names[i] << ": " << res << endl;
    }
  }
  cout << endl;
  
  cout << "Before fit: " << endl;
  for(size_t i=0;i<n_fits;i++) {
    cout.width(left_column);
    double res;
    fitter.eval(*nmfd[i],res);
    cout << fit_names[i] << ": " << res << endl;
  }
  cout << endl;

  cout << "After fit: " << endl;
  for(size_t i=0;i<n_fits;i++) {
    cout.width(left_column);
    double res;
    if (i==3) fitter.def_mmin.verbose=1;
    fitter.fit(*nmfd[i],res);
    cout << fit_names[i] << ": " << res << endl;
  }
    
  return 0;
}

void nucleus_class::setup_cli(o2scl::cli &cl) {
  
  static const int nopt=5;
  o2scl::comm_option_s options[nopt]={
    {0,"ZN","Get by Z and N.",
     2,2,"","",new o2scl::comm_option_mfptr<nucleus_class>
     (this,&nucleus_class::get),o2scl::cli::comm_option_both},
    {0,"tables","",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_class>
     (this,&nucleus_class::tables),o2scl::cli::comm_option_both},
    {0,"fits","",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_class>
     (this,&nucleus_class::fits),o2scl::cli::comm_option_both},
    {0,"cdist","Create a common distribution of nuclei",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_class>
     (this,&nucleus_class::cdist),o2scl::cli::comm_option_both},
    {0,"refs","",
     0,0,"","",new o2scl::comm_option_mfptr<nucleus_class>
     (this,&nucleus_class::refs),o2scl::cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);
    
  p_verbose.i=&verbose;
  p_verbose.help="Verbose parameter (default 1)";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  return;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  nucleus_class nuc;
  cli cl;
  
  nuc.setup_cli(cl);

  cl.run_auto(argc,argv);
  
  return 0;
}
