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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/eos_sn.h>
#include <o2scl/cli.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/cloud_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

/** \brief Desc
 */
class ex_eos_sn {

protected:

  /// Directory for the data files
  string directory;

  /// Level of output
  int verbose;

  /// Name of EOS table
  string name;

  /// Generic EOS pointer
  o2scl::eos_sn_base *genp;

  /// \name The supernova EOS objects
  //@{
  o2scl::eos_sn_ls ls;  
  o2scl::eos_sn_sht sht;  
  o2scl::eos_sn_stos stos;  
  o2scl::eos_sn_hfsl hfsl;  
  o2scl::eos_sn_oo oo;  
  //@}

  /** \brief Desc
   */
  int slices(std::vector<std::string> &sv, bool itive_com) {

    table3d tab;
    vector<size_t> index(3);

    // For Ye=0.5
    index[1]=genp->A.lookup_grid(1,0.5);
    genp->A.copy_slice_align(0,2,index,tab,"A_Ye0.5");
    genp->mun.copy_slice_align(0,2,index,tab,"mun_Ye0.5");
    genp->mup.copy_slice_align(0,2,index,tab,"mup_Ye0.5");
    genp->Eint.copy_slice_align(0,2,index,tab,"Eint_Ye0.5");
    genp->E.copy_slice_align(0,2,index,tab,"E_Ye0.5");
    genp->Sint.copy_slice_align(0,2,index,tab,"Sint_Ye0.5");
    genp->S.copy_slice_align(0,2,index,tab,"S_Ye0.5");
    genp->Xn.copy_slice_align(0,2,index,tab,"Xn_Ye0.5");
    genp->Xp.copy_slice_align(0,2,index,tab,"Xp_Ye0.5");
    genp->Xalpha.copy_slice_align(0,2,index,tab,"Xalpha_Ye0.5");
    genp->Xnuclei.copy_slice_align(0,2,index,tab,"Xnuclei_Ye0.5");

    // For Ye=0.1
    index[1]=genp->A.lookup_grid(1,0.1);
    genp->A.copy_slice_align(0,2,index,tab,"A_Ye0.1");
    genp->mun.copy_slice_align(0,2,index,tab,"mun_Ye0.1");
    genp->mup.copy_slice_align(0,2,index,tab,"mup_Ye0.1");
    genp->Eint.copy_slice_align(0,2,index,tab,"Eint_Ye0.1");
    genp->E.copy_slice_align(0,2,index,tab,"E_Ye0.1");
    genp->Sint.copy_slice_align(0,2,index,tab,"Sint_Ye0.1");
    genp->S.copy_slice_align(0,2,index,tab,"S_Ye0.1");
    genp->Xn.copy_slice_align(0,2,index,tab,"Xn_Ye0.1");
    genp->Xp.copy_slice_align(0,2,index,tab,"Xp_Ye0.1");
    genp->Xalpha.copy_slice_align(0,2,index,tab,"Xalpha_Ye0.1");
    genp->Xnuclei.copy_slice_align(0,2,index,tab,"Xnuclei_Ye0.1");

    // Output to file
    string outfile=name+"_slices.o2";
    hdf_file hf;
    hf.open_or_create(outfile);
    hdf_output(hf,tab,"Atab");
    hf.close();

    return 0;
  }

  /// \name Functions to load the EOSs
  //@{
  /** \brief Desc
   */
  int ls_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;

    if (sv[1]=="ls") {
      fname+="ls.dat";
      name="ls_ls";
    } else if (sv[1]=="skm") {
      fname+="skm.dat";
      name="ls_skm";
    } else if (sv[1]=="ska") {
      fname+="ska.dat";
      name="ls_ska";
    } else if (sv[1]=="sk1") {
      fname+="sk1.dat";
      name="ls_sk1";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    ls.verbose=verbose;
    ls.load(fname);
    genp=&ls;
    //test_mgr t;
    //ls.check_eg(t);
    //ls.check_free_energy();
    
    return 0;

  }

  /** \brief Desc
   */
  int stos_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;

    if (sv[1]=="stos") {
      fname+="eos1.tab";
      name="stos";
    } else if (sv[1]=="stos") {
      fname+="eos2.tab";
      name="stos";
    } else if (sv[1]=="stos") {
      fname+="eos3.tab";
      name="stos";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    stos.verbose=verbose;
    stos.load(fname,0);
    genp=&stos;
    
    return 0;

  }

  /** \brief Desc
   */
  int sht_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;
    string fname2=directory;

    sht.verbose=verbose;
    if (sv[1]=="fsu17") {
      fname+="FSU1.7eos1.01.dat";
      fname2+="FSU1.7eosb1.01.dat";
      name="sht_fsu17";
      sht.load(fname,eos_sn_sht::mode_17);
      //cout << sht.A.interp(0.02,0.5,1.0) << " ";
      //cout << sht.E.interp(0.02,0.5,1.0) << endl;
      sht.load(fname2,eos_sn_sht::mode_17b);
      //cout << sht.A.interp(0.02,0.5,1.0) << " ";
      //cout << sht.E.interp(0.02,0.5,1.0) << " ";
      //cout << sht.Eint.interp(0.02,0.5,1.0) << endl;
      exit(-1);
    } else if (sv[1]=="fsu21") {
      fname+="FSU2.1eos1.01.dat";
      fname2+="FSU2.1eosb1.01.dat";
      name="sht_fsu21";
      sht.load(fname,eos_sn_sht::mode_21);
      sht.load(fname2,eos_sn_sht::mode_21b);
    } else if (sv[1]=="nl3") {
      fname+="NL3eos1.03.dat";
      fname2+="NL3eosb1.03.dat";
      name="sht_nl3";
      sht.load(fname,eos_sn_sht::mode_NL3);
      sht.load(fname2,eos_sn_sht::mode_NL3b);
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    genp=&sht;
    
    return 0;

  }

  /** \brief Desc
   */
  int hfsl_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;

    if (sv[1]=="sfho") {
      fname+="sfho_frdm_shen98_v1.03.tab";
      name="hfsl_sfho";
    } else if (sv[1]=="sfhx") {
      fname+="sfhx_frdm_shen98_v1.03.tab";
      name="hfsl_sfhx";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    hfsl.verbose=verbose;
    hfsl.load(fname);
    genp=&hfsl;
    
    return 0;

  }

  /** \brief Desc
   */
  int eg(std::vector<std::string> &sv, bool itive_com) {
    if (true) {
      genp->relf.upper_limit_fac=40.0;
      genp->relf.dit->tol_abs=1.0e-11;
      genp->relf.dit->tol_rel=1.0e-11;
      genp->relf.nit->tol_abs=1.0e-11;
      genp->relf.nit->tol_rel=1.0e-11;
    }
    genp->compute_eg();
    return 0;
  }
  
  /** \brief Desc
   */
  int interp(std::vector<std::string> &sv, bool itive_com) {
    
    if (sv.size()<4) {
      cerr << "Not enough arguments specified in 'interp'." << endl;
      return 1;
    }
    if (genp==0 || genp->is_loaded()==false) {
      cerr << "No EOS table loaded in 'interp'." << endl;
      return 1;
    }

    double nB=o2scl::stod(sv[1]);
    double Ye=o2scl::stod(sv[2]);
    double T=o2scl::stod(sv[3]);
    
    cout << "nB= " << nB << " 1/fm^3" << endl;
    cout << "Ye= " << Ye << endl;
    cout << "T= " << T << " MeV" << endl;
    if (genp->data_with_leptons()) {
      cout << "F= " << genp->F.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "E= " << genp->E.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "P= " << genp->P.interp_linear(nB,Ye,T) << " MeV/fm^3" << endl;
      cout << "S= " << genp->S.interp_linear(nB,Ye,T) << endl;
    }
    if (genp->data_baryons_only()) {
      cout << "Fint= " << genp->Fint.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "Eint= " << genp->Eint.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "Pint= " << genp->Pint.interp_linear(nB,Ye,T) << " MeV/fm^3" 
	   << endl;
      cout << "Sint= " << genp->Sint.interp_linear(nB,Ye,T) << endl;
    }
    cout << "mun= " << genp->mun.interp_linear(nB,Ye,T) << " MeV" << endl;
    cout << "mup= " << genp->mup.interp_linear(nB,Ye,T) << " MeV" << endl;
    cout << "Z= " << genp->Z.interp_linear(nB,Ye,T) << endl;
    cout << "A= " << genp->A.interp_linear(nB,Ye,T) << endl;
    cout << "Xn= " << genp->Xn.interp_linear(nB,Ye,T) << endl;
    cout << "Xp= " << genp->Xp.interp_linear(nB,Ye,T) << endl;
    cout << "Xalpha= " << genp->Xalpha.interp_linear(nB,Ye,T) << endl;
    cout << "Xnuclei= " << genp->Xnuclei.interp_linear(nB,Ye,T) << endl;
    cout << endl;
    
    for(size_t i=0;i<(genp->n_oth);i++) {
      cout << genp->oth_names[i] << "= " 
	   << genp->other[i].interp_linear(nB,Ye,T) << endl;
    }
    cout << endl;

    return 0;
  }

  /** \brief Desc
   */
  int oo_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;
    size_t mode;

    if (sv[1]=="fsu17") {
      cloud_file cf;
      cf.verbose=2;
      fname+="GShenFSU_1.7EOS_rho280_temp180_ye52_version_1.1_20120817.h5";
      std::string sha=((std::string)"57aec0f5011caf0333fc93ea818c786b0e2")+
	"180e975425e1f4d90a3458b46f131";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="fsu17";
    } else if (sv[1]=="fsu21") {
      cloud_file cf;
      cf.verbose=2;
      fname+="GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5";
      std::string sha=((std::string)"cdf857d69884f2661c857c6bcce501af75d")+
	"48e51bb14a5fab89872df5ed834f6";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="fsu21";
    } else if (sv[1]=="sht_nl3") {
      cloud_file cf;
      cf.verbose=2;
      fname+="GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5";
      std::string sha=((std::string)"982cd2249e08895a3580dc4969ab79add72")+
	"d1d7ce4464e946c18f9950edb7bfe";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="sht_nl3";
    } else if (sv[1]=="stos") {
      cloud_file cf;
      cf.verbose=2;
      fname+="HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5";
      std::string sha=((std::string)"3b7c598bf56ec12d734e13a97daf1eeb1f5")+
	"8f59849c5f65c4f9f72dd292b177c";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="stos";
      mode=eos_sn_oo::stos_mode;
    } else if (sv[1]=="stos_hyp") {
      fname+="HShen_HyperonEOS_rho220_temp180_ye65_version_1.1_20131007.h5";
      name="stos_hyp";
      mode=eos_sn_oo::stos_mode;
    } else if (sv[1]=="dd2") {
      fname+="Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="dd2";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="fsg") {
      fname+="Hempel_FSGEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="fsg";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="hfsl_nl3") {
      fname+="Hempel_NL3EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="hfsl_nl3";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="sfho") {
      fname+="Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5";
      name="sfho";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="sfhx") {
      fname+="Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="sfhx";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="tm1") {
      fname+="Hempel_TM1EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="tm1";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="tma") {
      fname+="Hempel_TMAEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
      name="tma";
      mode=eos_sn_oo::hfsl_mode;
    } else if (sv[1]=="ls180") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"3172f0f7b542fa1bd2e7a46f1b2e62c848f")+
	"0d9a979e546902ad3f3b6285e27ca";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="ls180";
      mode=eos_sn_oo::ls_mode;
    } else if (sv[1]=="ls220") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"d8c4d4f1315942a663e96fc6452f66d90fc")+
	"87f283e0ed552c8141d1ddba34c19";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS220_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS220_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="ls220";
      mode=eos_sn_oo::ls_mode;
    } else if (sv[1]=="ls375") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"13d31df4944b968f1de799e5fc37881eae9")+
	"8cd3f2d3bc14648698661cef35bdd";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS375_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       sha,((string)"https://isospin.roam.utk.edu/")+
		       "public/eos_tables/scollapse/LS375_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",fname,".o2scl_data");
      name="ls375";
      mode=eos_sn_oo::ls_mode;
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    oo.verbose=verbose;
    oo.load(fname,mode);
    genp=&oo;
    //test_mgr t;
    //ls.check_eg(t);
    //ls.check_free_energy();
    
    return 0;

  }
  //@}
  
public:
  
  ex_eos_sn() {
    directory=".";
    verbose=1;
    genp=0;
  }

  /** \brief Desc<
   */
  void run(int argc, char *argv[]) {

    // ---------------------------------------
    // Specify command-line option object
    
    cli cl;
    cl.prompt="ex_eos_sn> ";
    cl.gnu_intro=false;

    // ---------------------------------------
    // Set options
    
    static const int nopt=8;
    comm_option_s options[nopt]={
      {0,"ls","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::ls_fun),
       cli::comm_option_both},
      {0,"oo","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::oo_fun),
       cli::comm_option_both},
      {0,"sht","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::sht_fun),
       cli::comm_option_both},
      {0,"stos","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::stos_fun),
       cli::comm_option_both},
      {0,"hfsl","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::hfsl_fun),
       cli::comm_option_both},
      {0,"slices","short desc",
       0,0,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::slices),
       cli::comm_option_both},
      {0,"interp","short desc",
       3,3,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::interp),
       cli::comm_option_both},
      {0,"eg","short desc",
       0,0,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::eg),
       cli::comm_option_both}
    };
    cl.set_comm_option_vec(nopt,options);
    
    // ---------------------------------------
    // Set parameters

    cli::parameter_int p_verbose;
    p_verbose.i=&verbose;
    p_verbose.help=((string)"Desc.");
    cl.par_list.insert(make_pair("verbose",&p_verbose));

    cli::parameter_string p_directory;
    p_directory.str=&directory;
    p_directory.help=((string)"Desc.");
    cl.par_list.insert(make_pair("directory",&p_directory));

    // ---------------------------------------
    // Process command-line arguments and run
    
    cl.run_auto(argc,argv);

    return;
  }


};

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  ex_eos_sn gset;
  gset.run(argc,argv);
  
  return 0;
}
