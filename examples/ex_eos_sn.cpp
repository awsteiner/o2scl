/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

/** \brief A class for manipulating EOS tables 

    \note Highly experimental.
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
  o2scl::eos_sn_base nat;
  //@}

  /** \brief Desc
   */
  int slices(std::vector<std::string> &sv, bool itive_com) {

    table3d tab;
    vector<size_t> index(3);

    // For Ye=0.5
    index[1]=genp->A.lookup_grid(1,0.5);
    genp->A.copy_table3d_align(0,2,index,tab,"A_Ye0.5");
    genp->mun.copy_table3d_align(0,2,index,tab,"mun_Ye0.5");
    genp->mup.copy_table3d_align(0,2,index,tab,"mup_Ye0.5");
    genp->Eint.copy_table3d_align(0,2,index,tab,"Eint_Ye0.5");
    genp->E.copy_table3d_align(0,2,index,tab,"E_Ye0.5");
    genp->Sint.copy_table3d_align(0,2,index,tab,"Sint_Ye0.5");
    genp->S.copy_table3d_align(0,2,index,tab,"S_Ye0.5");
    genp->Xn.copy_table3d_align(0,2,index,tab,"Xn_Ye0.5");
    genp->Xp.copy_table3d_align(0,2,index,tab,"Xp_Ye0.5");
    genp->Xalpha.copy_table3d_align(0,2,index,tab,"Xalpha_Ye0.5");
    genp->Xnuclei.copy_table3d_align(0,2,index,tab,"Xnuclei_Ye0.5");

    // For Ye=0.1
    index[1]=genp->A.lookup_grid(1,0.1);
    genp->A.copy_table3d_align(0,2,index,tab,"A_Ye0.1");
    genp->mun.copy_table3d_align(0,2,index,tab,"mun_Ye0.1");
    genp->mup.copy_table3d_align(0,2,index,tab,"mup_Ye0.1");
    genp->Eint.copy_table3d_align(0,2,index,tab,"Eint_Ye0.1");
    genp->E.copy_table3d_align(0,2,index,tab,"E_Ye0.1");
    genp->Sint.copy_table3d_align(0,2,index,tab,"Sint_Ye0.1");
    genp->S.copy_table3d_align(0,2,index,tab,"S_Ye0.1");
    genp->Xn.copy_table3d_align(0,2,index,tab,"Xn_Ye0.1");
    genp->Xp.copy_table3d_align(0,2,index,tab,"Xp_Ye0.1");
    genp->Xalpha.copy_table3d_align(0,2,index,tab,"Xalpha_Ye0.1");
    genp->Xnuclei.copy_table3d_align(0,2,index,tab,"Xnuclei_Ye0.1");

    // Output to file
    string outfile=name+"_slices.o2";
    hdf_file hf;
    hf.open_or_create(outfile);
    hdf_output(hf,(const table3d &)tab,"Atab");
    hf.close();

    return 0;
  }

  /// \name Functions to load the EOSs
  //@{
  /** \brief Load a Lattimer-Swesty EOS
   */
  int ls_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=directory;

    if (sv[1]=="ls") {
      fname+="/ls.dat";
      name="ls_ls";
    } else if (sv[1]=="skm") {
      fname+="/skm.dat";
      name="ls_skm";
    } else if (sv[1]=="ska") {
      fname+="/ska.dat";
      name="ls_ska";
    } else if (sv[1]=="sk1") {
      fname+="/sk1.dat";
      name="ls_sk1";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    ls.verbose=verbose;
    ls.load(fname,0);
    genp=&ls;
    return 0;

  }

  /** \brief Desc
   */
  int native_fun(std::vector<std::string> &sv, bool itive_com) {

    string fname=sv[1];

    nat.verbose=verbose;
    cout << "Going to load." << endl;
    nat.load(fname,0);
    cout << "Done in load." << endl;
    genp=&nat;
    return 0;

  }

  /** \brief Desc
   */
  int create_ZoA(std::vector<std::string> &sv, bool itive_com) {

    if (sv.size()<2) {
      cerr << "Not enough arguments for create_ZoA." << endl;
      return 1;
    }

    vector<vector<double> > grid={genp->nB_grid,genp->Ye_grid,
				  genp->T_grid};
    tensor_grid3<vector<double>,vector<size_t> > ZoA;
    vector<size_t> sz={genp->n_nB,genp->n_Ye,genp->n_T};
    ZoA.resize(3,sz);
    ZoA.set_grid(grid);
    
    for(size_t iT=0;iT<genp->n_T;iT++) {
      cout << "iT: " << iT << endl;
      for(size_t iYe=0;iYe<genp->n_Ye;iYe++) {
	for(size_t inB=0;inB<genp->n_nB;inB++) {
	  if (genp->Xnuclei.get(inB,iYe,iT)>0.1) {
	    ZoA.get(inB,iYe,iT)=genp->Z.get(inB,iYe,iT)/
	      genp->A.get(inB,iYe,iT);
	  } else {
	    ZoA.get(inB,iYe,iT)=0.0;
	  }
	}
      }
    }

    hdf_file hf;
    hf.open_or_create(sv[1]);
    hdf_output(hf,ZoA,"ZoA");
    hf.close();

    return 0;

  }

  /** \brief Output to a file in native format
   */
  int output(std::vector<std::string> &sv, bool itive_com) {

    if (sv.size()<2) {
      cerr << "Output needs a filename." << endl;
      return 1;
    }
    
    genp->output(directory+"/"+sv[1]);
    
    return 0;
  }

  /** \brief Check the EOS
   */
  int check(std::vector<std::string> &sv, bool itive_com) {

    double v1=genp->check_eg();
    double v2;
    genp->check_free_energy(v2);
    cout << v1 << " " << v2 << endl;
    
    return 0;
  }

  /** \brief Load an H. Shen et al. EOS
   */
  int stos_fun(std::vector<std::string> &sv, bool itive_com) {

    stos.verbose=verbose;
    stos.load_auto(sv[1],directory);
    genp=&stos;

    /*
    string fname=directory;
    size_t mode=0;
    
    if (sv[1]=="stos") {
      fname+="eos1.tab";
      name="stos";
    } else if (sv[1]=="stos2") {
      fname+="eos2.tab";
      name="stos";
    } else if (sv[1]=="stos3") {
      fname+="eos3.tab";
      name="stos";
    } else if (sv[1]=="fyss") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"47d357600d875a2a24fbfb7b8064602")+
	"5434398a42113ffdf1f9121e32d9bdabb";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash
	("FYSS_ver_1_27.tab",
	 ((string)"https://isospin.roam.utk.edu/")+
	 "public_data/eos_tables/stos/"+
	 "FYSS_ver_1_27.tab",sha,directory);
      name="fyss";
      mode=eos_sn_stos::fyss_mode;
      fname=directory+"FYSS_ver_1_27.tab";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    stos.verbose=verbose;
    stos.load(fname,mode);
    genp=&stos;
    */
    
    return 0;

  }

  /** \brief Load a G. Shen et al. EOS
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
      sht.load(fname2,eos_sn_sht::mode_17b);
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

  /** \brief Load a Hempel et al. EOS
   */
  int hfsl_fun(std::vector<std::string> &sv, bool itive_com) {

    hfsl.verbose=verbose;
    hfsl.load_auto(sv[1],directory);
    genp=&hfsl;
    
    return 0;

  }

  /** \brief Compute the leptonic part of the EOS
   */
  int eg(std::vector<std::string> &sv, bool itive_com) {

    if (true) {
      // Make the electron EOS a bit more accurate
      genp->elep.frel.upper_limit_fac=40.0;
      genp->elep.frel.fri.dit.tol_abs=1.0e-11;
      genp->elep.frel.fri.dit.tol_rel=1.0e-11;
      genp->elep.frel.fri.nit.tol_abs=1.0e-11;
      genp->elep.frel.fri.nit.tol_rel=1.0e-11;
      genp->elep.frel.density_root->tol_rel=1.0e-10;
    }
    genp->compute_eg();
    return 0;
  }
  
  /** \brief Compute the EOS at one point
   */
  int point(std::vector<std::string> &sv, bool itive_com) {
    
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
    cout << "Baryons only EOS: " << genp->data_baryons_only() << endl;
    cout << "EOS with leptons: " << genp->data_with_leptons() << endl;
    cout << "Muons: " << genp->include_muons << endl;
    cout << endl;

    thermo th;
    if (genp->data_with_leptons()==false ||
	genp->data_baryons_only()==false) {
      
      genp->elep.frel.upper_limit_fac=40.0;
      genp->elep.frel.fri.dit.tol_abs=1.0e-11;
      genp->elep.frel.fri.dit.tol_rel=1.0e-11;
      genp->elep.frel.fri.nit.tol_abs=1.0e-11;
      genp->elep.frel.fri.nit.tol_rel=1.0e-11;
      genp->elep.frel.density_root->tol_rel=1.0e-10;

      double mue;
      genp->compute_eg_point(nB,Ye,T,th,mue);
      cout << "Automatically computing electron-photon EOS at "
	   << "specified point\n." << endl;
    }

    if (genp->data_with_leptons()==false) {
      cout << "F= " << genp->Fint.interp_linear(nB,Ye,T)+
	(th.ed*hc_mev_fm-T*th.en)/nB << " MeV" << endl;
      cout << "E= " << genp->Eint.interp_linear(nB,Ye,T)+
	th.ed/nB*hc_mev_fm << " MeV" << endl;
      cout << "P= " << genp->Pint.interp_linear(nB,Ye,T)+th.pr*hc_mev_fm
	   << " MeV/fm^3" << endl;
      cout << "S= " << genp->Sint.interp_linear(nB,Ye,T)+th.en/nB << endl;
    } else {
      cout << "F= " << genp->F.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "E= " << genp->E.interp_linear(nB,Ye,T) << " MeV" << endl;
      cout << "P= " << genp->P.interp_linear(nB,Ye,T) << " MeV/fm^3" << endl;
      cout << "S= " << genp->S.interp_linear(nB,Ye,T) << endl;
    }
    if (genp->data_baryons_only()==false) {
      cout << "Fint= " << genp->F.interp_linear(nB,Ye,T)-
	(th.ed*hc_mev_fm-T*th.en)/nB << " MeV" << endl;
      cout << "Eint= " << genp->E.interp_linear(nB,Ye,T)-
	th.ed/nB*hc_mev_fm << " MeV" << endl;
      cout << "Pint= " << genp->P.interp_linear(nB,Ye,T)-th.pr*hc_mev_fm
	   << " MeV/fm^3" << endl;
      cout << "Sint= " << genp->S.interp_linear(nB,Ye,T)-th.en/nB << endl;
    } else {
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

  /** \brief Load an EOS in the O'Connor-Ott format
   */
  int oo_fun(std::vector<std::string> &sv, bool itive_com) {

    oo.verbose=verbose;
    oo.load_auto(sv[1],directory);
    genp=&oo;
    
    return 0;

  }
  //@}
  
public:
  
  ex_eos_sn() {
    directory="~/.o2scl_data";
    verbose=1;
    genp=0;
  }

  /** \brief Main executable wrapper
   */
  void run(int argc, char *argv[]) {

    // ---------------------------------------
    // Specify command-line option object
    
    cli cl;
    cl.prompt="ex_eos_sn> ";
    cl.gnu_intro=false;

    // ---------------------------------------
    // Set options
    
    static const int nopt=12;
    comm_option_s options[nopt]={
      {0,"ls","Load an EOS in the Lattimer-Swesty format.",
       1,1,"<model>",((string)"Models \"ls\", \"skm\", \"ska\", ")+
       "or \"sk1\".",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::ls_fun),
       cli::comm_option_both},
      {0,"native","Load an EOS in the native format.",
       1,1,"<filename>","",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::native_fun),
       cli::comm_option_both},
      {0,"oo","Load an EOS in the O'Connor-Ott format.",
       1,1,"<model>",((string)"Models \"fsu17\", \"fsu21\", \"sht_nl3\", ")+
       "\"stos\", \"stos_hyp\", \"dd2\", \"fsg\", \"hfsl_nl3\", "+
       "\"sfho\", \"sfhx\", \"tm1\", \"tma\", \"ls180\", \"ls220\", "+
       "or \"ls375\".",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::oo_fun),
       cli::comm_option_both},
      {0,"sht","Load an EOS in the Shen-Horowitz-Teige format.",
       1,1,"<model>",((string)"Models \"fsu17\", \"fsu21\", ")+
       "or \"nl3\".",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::sht_fun),
       cli::comm_option_both},
      {0,"stos","Load an EOS in the Shen et al. format.",
       1,1,"<model>","Models are \"stos1\", \"stos2\" or \"stos3\".",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::stos_fun),
       cli::comm_option_both},
      {0,"hfsl","Load an EOS in the Hempel et al. format.",
       1,1,"<model>","Models are \"sfho\" and \"sfhx\"",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::hfsl_fun),
       cli::comm_option_both},
      {0,"slices","Construct slices.",
       0,0,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::slices),
       cli::comm_option_both},
      {0,"point","Interpolate the EOS at a specified (n_B,Y_e,T) point.",
       3,3,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::point),
       cli::comm_option_both},
      {0,"eg","Compute the electron-photon part of the EOS.",
       0,0,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::eg),
       cli::comm_option_both},
      {0,"output","Output to a file.",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::output),
       cli::comm_option_both},
      {0,"ZoA","",
       -1,-1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::create_ZoA),
       cli::comm_option_both},
      {0,"check","",
       0,0,"",((string)"long ")+"desc.",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::check),
       cli::comm_option_both}
    };
    cl.set_comm_option_vec(nopt,options);
    
    // ---------------------------------------
    // Set parameters

    cli::parameter_int p_verbose;
    p_verbose.i=&verbose;
    p_verbose.help=((string)"Verbosity parameter");
    cl.par_list.insert(make_pair("verbose",&p_verbose));

    cli::parameter_string p_directory;
    p_directory.str=&directory;
    p_directory.help=((string)"Directory for EOS table storage");
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
