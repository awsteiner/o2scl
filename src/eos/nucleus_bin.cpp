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

  //std::cout << "Reading nuclear mass tables." << std::endl;
    
  o2scl_hdf::ame_load(ame95rmd,"95rmd");
  o2scl_hdf::ame_load(ame95exp,"95exp");
  o2scl_hdf::ame_load(ame03round,"03round");
  o2scl_hdf::ame_load(ame03,"03");
  o2scl_hdf::ame_load(ame12,"12");
  o2scl_hdf::ame_load(ame16,"16");
  o2scl_hdf::ame_load(ame20exp,"20");
  o2scl_hdf::ame_load(ame20round,"20round");
    
  o2scl_hdf::mnmsk_load(m95,"mnmsk97");
  o2scl_hdf::mnmsk_load(m16,"msis16");
    
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
  wlw4.load("WS3_RBF");
  wlw5.load("WS4_RBF");
    
  sdnp1.load("sdnp03");
  sdnp2.load("sd_skp_04");
  sdnp3.load("sd_sly4_04");
  
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

  nmd={&ame95rmd,&ame95exp,&ame03round,&ame03,
       &ame12,&ame16,&ame20exp,&ame20round,&m95,&m16,&kt,
       &kt2,&hfb2,&hfb8,&hfb14,&hfb14_v0,
       &hfb17,&hfb21,&hfb22,&hfb23,
       &hfb24,&hfb25,&hfb26,&hfb27,
       &wlw1,&wlw2,&wlw3,&wlw4,&wlw5,&sdnp1,
       &sdnp2,&sdnp3,&dz,&ddme2,&ddmed,&ddpc1,
       &nl3s,&sly4,&skms,&skp,&sv_min,&unedf0,&unedf1};
       
  n_tables=nmd.size();

  table_names={"AME rmd 95","AME 95 exp",
	       "AME rnd 03","AME 03",
	       "AME 12","AME 16","AME exp 20","AME rnd 20",
	       "MNMSK 95","MSIS 16","KTUY 04","KTUY 05",
	       "HFB2","HFB8","HFB14","HFB14_v0",
	       "HFB17","HFB21","HFB22","HFB23",
	       "HFB24","HFB25","HFB26","HFB27",
	       "WLW 10","WLLW 10","LWDW 11","WL 11","WLWM 14",
	       "SDNP1 05","SDNP2 05","SDNP3 05",
	       "DZ 95","DDME2","DDMED","DDPC1","NL3S",
               "SLy4","SKM*","SkP","SV-min","UNEDF0","UNEDF1"};

  nmfd={&se,&frdm,&dzf,&dzf33,&frdm_shell,&ldrop_shell};
  n_fits=nmfd.size();
    
  fit_names={"Semi-empirical","FRDM","DZ fit 10","DZ fit 33","FRDM shell",
             "Ldrop shell"};

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
    p[0]=6.262656723970295e-01;
    p[1]=1.062118484793815e+00;
    p[2]=4.278509731327071e+01;
    p[3]=1.712865295142668e+01;
    p[4]=2.750315169290062e+01;
    p[5]=3.585974615126219e+01;
    p[6]=2.488088599031062e+01;
    p[7]=6.952310612512534e-01;
    p[8]=1.133367354754059e+00;
    p[9]=2.989576509602983e-01;
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
    p[0]=2.657399203505667e+00;
    p[1]=1.136316849511381e+00;
    p[2]=3.548541967676024e+01;
    p[3]=1.649009712939987e+01;
    p[4]=2.331520244072635e+01;
    p[5]=3.437294342444157e+01;
    p[6]=2.660294896136239e+01;
    p[7]=4.409661805078772e-01;
    p[8]=2.501490274580595e+01;
    p[9]=1.897879047364143e+00;
    p[10]=-1.460475719478483e+00;
    p[11]=1.928679183044334e-02;
    p[12]=1.901033655300475e-03;
    p[13]=8.970261351311239e-02;    
    frdm_shell.fit_fun(14,p);
  }
  if (true) {
    ubvector p(11);
    p[0]=8.994776301007761e-01;
    p[1]=9.679865078598426e-01;
    p[2]=8.751369188587536e-01;
    p[3]=9.710432146736609e-01;
    p[4]=-9.041294789462331e-03;
    p[5]=1.390547985659261e-01;
    p[6]=1.246579548574642e+01;
    p[7]=-1.493972115439528e+00;
    p[8]=1.419539065031770e-02;
    p[9]=1.659654542326672e-03;
    p[10]=1.136613448382515e-01;
    ldrop_shell.fit_fun(11,p);
  }
}

int nucleus_bin::get(std::vector<std::string> &sv, bool itive_com) {

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
  cout.width(10);
  cout << "mass ex. MeV ";
  cout.width(10);
  cout << "BE/A MeV " << "  ";
  cout.width(left_column+1);
  cout << "Model  ";
  cout.width(10);
  cout << "mass ex. MeV ";
  cout.width(10);
  cout << "BE/A MeV " << endl;
  cout << endl;
    
  nucleus nuc;

  cout.precision(4);
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
    
  for(size_t i=0;i<n_fits;i+=2) {
    if (nmfd[i]->is_included(Z,N)) {
      int ret=nmfd[i]->get_nucleus(Z,N,nuc);
      cout.width(left_column);
      cout << fit_names[i] << " "
	   << nuc.mex*o2scl_const::hc_mev_fm << " "
	   << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << " ";
    }
    if (i+1<n_fits) {
      if (nmfd[i+1]->is_included(Z,N)) {
        int ret=nmfd[i+1]->get_nucleus(Z,N,nuc);
        cout.width(left_column);
        cout << fit_names[i+1] << " "
             << nuc.mex*o2scl_const::hc_mev_fm << " "
             << nuc.be*o2scl_const::hc_mev_fm/(Z+N) << endl;
      }
    } else {
      cout << endl;
    }
  }
  cout.precision(6);
    
  return 0;
}

int nucleus_bin::tables(std::vector<std::string> &sv, bool itive_com) {

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

  size_t left_column=18;

  cout << "References: " << endl;
  for(size_t i=0;i<n_tables;i++) {
    cout.width(left_column);
    cout << table_names[i] << ": " << nmd[i]->reference << endl;
  }
    
  return 0;
}

int nucleus_bin::cdist(std::vector<std::string> &sv, bool itive_com) {

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

int nucleus_bin::compare(std::vector<std::string> &sv, bool itive_com) {

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

  std::cout << "1." << endl;
  fitter.def_mmin.ntrial*=100;
  
  double res;
  kwargs kw;
  if (sv.size()>=3) {
    kw.set(sv[2]);
  }
  std::cout << "2." << endl;
  
  nucdist_set(fitter.dist,ame20exp);
  ldrop_shell.large_vals_unphys=true;

  std::cout << "3." << endl;
  fitter.def_mmin.verbose=kw.get_int("verbose",0);

  size_t ix_start=0;
  size_t ix_end=n_fits;
  std::cout << "4." << endl;
  if (sv.size()>=2) {
    ix_start=o2scl::stoszt(sv[1]);
    ix_end=ix_start+1;
  }
  for(size_t i=ix_start;i<ix_end;i++) {
    std::cout << "5." << endl;
    cout << fit_names[i] << ": "<< endl;
    fitter.fit(*nmfd[i],res);
    cout << fit_names[i] << ": " << res << endl;
    ubvector p(nmfd[i]->nfit);
    nmfd[i]->guess_fun(nmfd[i]->nfit,p);
    for(size_t k=0;k<nmfd[i]->nfit;k++) {
      cout << "    p[" << k << "]=" << dtos(p[k],0) << ";" << endl;
    }
    std::cout << "6." << endl;
  }
    
  return 0;
}

void nucleus_bin::setup_cli(o2scl::cli &cl) {

  std::cout << "A." << std::endl;
  
  static const int nopt=6;
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
     (this,&nucleus_bin::refs),o2scl::cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);
    
  p_verbose.i=&verbose;
  p_verbose.help="Verbosity parameter (default 1)";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  std::cout << "B." << std::endl;
  
  return;
}
