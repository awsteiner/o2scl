/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

class eos_sn_compose : public eos_sn_base {
  
public:
  
  eos_sn_compose() {
  }
  
  /// Load table from filename \c fname with mode \c mode
  virtual void load() {
    
    //wordexp_single_file(fname);

    vector<double> grid;
    
    ifstream fin("eos.nb");
    int n_nB;
    // the first entry is ignored
    fin >> n_nB >> n_nB;
    nB_grid.resize(n_nB);
    for(int j=0;j<n_nB;j++) {
      fin >> nB_grid[j];
      grid.push_back(nB_grid[j]);
    }
    fin.close();
    
    ifstream fin2("eos.yq");
    int n_Ye;
    // the first entry is ignored
    fin2 >> n_Ye >> n_Ye;
    Ye_grid.resize(n_Ye);
    for(int j=0;j<n_Ye;j++) {
      fin2 >> Ye_grid[j];
      grid.push_back(Ye_grid[j]);
    }
    fin2.close();
    
    ifstream fin3("eos.t");
    int n_T;
    // the first entry is ignored
    fin3 >> n_T >> n_T;
    T_grid.resize(n_T);
    for(int j=0;j<n_T;j++) {
      fin3 >> T_grid[j];
      grid.push_back(T_grid[j]);
    }
    fin3.close();

    alloc();
    
    for(size_t i=0;i<n_base+n_oth;i++) {
      arr[i]->set_grid_packed(grid);
    }

    ifstream fin4("eos.thermo");
    fin4 >> m_neut;
    fin4 >> m_prot;
    
    double dtemp, dtemp2;
    for(int m=0;m<n_T;m++) {
      for(int k=0;k<n_Ye;k++) {
	for(int j=0;j<n_nB;j++) {
	  
	  // Skip the grid points since we know them already
	  fin4 >> dtemp;
	  fin4 >> dtemp;
	  fin4 >> dtemp;
	  
	  fin4 >> dtemp;
	  P.set(j,k,m,dtemp*nB_grid[j]);
	  fin4 >> dtemp;
	  S.set(j,k,m,dtemp);
	  
	  fin4 >> dtemp;
	  mun.set(j,k,m,(dtemp+1.0)*m_neut);
	  fin4 >> dtemp2;
	  mup.set(j,k,m,dtemp2*m_neut+(dtemp+1.0)*m_neut);
	  // Skip the lepton chemical potential
	  fin4 >> dtemp;
	  
	  fin4 >> dtemp;
	  F.set(j,k,m,(dtemp+1.0)*nB_grid[j]*m_neut);
	  fin4 >> dtemp;
	  E.set(j,k,m,(dtemp+1.0)*nB_grid[j]*m_neut);

	  // Skip the last column
	  fin4 >> dtemp;
	}
      }
    }
    fin4.close();

    ifstream fin5("eos.compo");
    
    for(int m=0;m<n_T;m++) {
      for(int k=0;k<n_Ye;k++) {
	for(int j=0;j<n_nB;j++) {
	  
	  // Skip the grid points since we know them already
	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  fin5 >> dtemp;

	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  fin5 >> dtemp;
	  
	  // This isn't right yet
	  fin5 >> dtemp;
	  A.set(j,k,m,dtemp);
	  fin5 >> dtemp;
	  Z.set(j,k,m,dtemp);
	  
	  // Skip the last column
	  fin5 >> dtemp;
	}
      }
    }
    fin5.close();

    // Loaded must be set to true before calling set_interp()
    n_oth=0;
    loaded=true;
    with_leptons_loaded=true;
    baryons_only_loaded=false;
    
    if (n_oth!=oth_names.size()) {
      O2SCL_ERR2("Number of names does not match number of data sets ",
		 "in eos_sn_oo::load().",exc_efailed);
    }
    
    // It is important that 'loaded' is set to true before the call to
    // set_interp_type().
    set_interp_type(itp_linear);
    
    if (verbose>0) {
      std::cout << "Done in eos_sn_compose::load()." << std::endl;
    }
    
    return;
  }
  
};

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
    ls.load(fname,0);
    genp=&ls;
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

  /** \brief Compute the leptonic part of the EOS
   */
  int eg(std::vector<std::string> &sv, bool itive_com) {

    if (true) {
      // Make the electron EOS a bit more accurate
      genp->relf.upper_limit_fac=40.0;
      genp->relf.dit->tol_abs=1.0e-11;
      genp->relf.dit->tol_rel=1.0e-11;
      genp->relf.nit->tol_abs=1.0e-11;
      genp->relf.nit->tol_rel=1.0e-11;
      genp->relf.density_root->tol_rel=1.0e-10;
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

    thermo th;
    if (genp->data_with_leptons()==false ||
	genp->data_baryons_only()==false) {
      
      genp->relf.upper_limit_fac=40.0;
      genp->relf.dit->tol_abs=1.0e-11;
      genp->relf.dit->tol_rel=1.0e-11;
      genp->relf.nit->tol_abs=1.0e-11;
      genp->relf.nit->tol_rel=1.0e-11;
      genp->relf.density_root->tol_rel=1.0e-10;

      double mue;
      genp->compute_eg_point(nB,Ye,T,th,mue);
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

    string fname;

    // Initialize to zero to prevent uninitialized variable errors
    size_t mode=0;

    if (sv[1]=="fsu17") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"57aec0f5011caf0333fc93ea818c786b0e2")+
	"180e975425e1f4d90a3458b46f131";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash
	("GShenFSU_1.7EOS_rho280_temp180_ye52_version_1.1_20120817.h5",
	 ((string)"https://isospin.roam.utk.edu/")+
	 "public_data/eos_tables/scollapse/GShenFSU_1.7EOS_rho280_"+
	 "temp180_ye52_version_1.1_20120817.h5",sha,directory);
      name="fsu17";
      mode=eos_sn_oo::sht_mode;
      fname=directory+
	"/GShenFSU_1.7EOS_rho280_temp180_ye52_version_1.1_20120817.h5";
    } else if (sv[1]=="fsu21") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"cdf857d69884f2661c857c6bcce501af75d")+
	"48e51bb14a5fab89872df5ed834f6";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash(((string)"GShenFSU_2.1EOS_rho280_temp180_")+
		       "ye52_version_1.1_20120824.h5",
		       ((string)"https://isospin.roam.utk.edu/public")+
		       "_data/eos_tables/scollapse/GShenFSU_2.1EOS_"+
		       "rho280_temp180_ye52_version_1.1_20120824.h5",
		       sha,directory);
      name="fsu21";
      mode=eos_sn_oo::sht_mode;
      fname=directory+
	"/GShenFSU_2.1EOS_rho280_temp180_ye52_version_1.1_20120824.h5";
    } else if (sv[1]=="sht_nl3") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"982cd2249e08895a3580dc4969ab79add72")+
	"d1d7ce4464e946c18f9950edb7bfe";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash(((string)"GShen_NL3EOS_rho280_temp180_ye52_")+
		       "version_1.1_20120817.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/GShen_NL3EOS_"+
		       "rho280_temp180_ye52_version_1.1_20120817.h5",
		       sha,directory);
      name="sht_nl3";
      mode=eos_sn_oo::sht_mode;
      fname=directory+
	"/GShen_NL3EOS_rho280_temp180_ye52_version_1.1_20120817.h5";
    } else if (sv[1]=="stos") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"3b7c598bf56ec12d734e13a97daf1eeb1f5")+
	"8f59849c5f65c4f9f72dd292b177c";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/HShenEOS_rho220_"+
		       "temp180_ye65_version_1.1_20120817.h5",
		       sha,directory);
      name="stos";
      mode=eos_sn_oo::stos_mode;
      fname=directory+"/HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5";
    } else if (sv[1]=="stos_hyp") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="stos_hyp";
      mode=eos_sn_oo::stos_mode;
      fname=directory+
	"/HShen_HyperonEOS_rho220_temp180_ye65_version_1.1_20131007.h5";
    } else if (sv[1]=="dd2") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="dd2";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="fsg") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="fsg";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_FSGEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="hfsl_nl3") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="hfsl_nl3";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_NL3EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="sfho") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"82a4acd670189917800567f6b75bb2a")+
	"3605f6ae7f9068215a1eec0acf924cb3d";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("_",((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/Hempel_SFHo"+
		       "EOS_rho222_temp180_ye60_version_1.1_20120817.h5",
		       sha,directory);
      name="sfho";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="sfhx") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"8651770ee78fb3dede5af1fe0cec33d6")+
	"bfc86ef2bd8505ab99db4d31f236fc44";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash(((string)"Hempel_SFHxEOS_rho234_temp180_ye60_")+
		       "version_1.1_20120817.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/Hempel_SFHx"+
		       "EOS_rho234_temp180_ye60_version_1.1_20120817.h5",
		       sha,directory);
      name="sfhx";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="tm1") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="tm1";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_TM1EOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="tma") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"");
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/",
		       sha,directory);
      name="tma";
      mode=eos_sn_oo::hfsl_mode;
      fname=directory+
	"/Hempel_TMAEOS_rho234_temp180_ye60_version_1.1_20120817.h5";
    } else if (sv[1]=="ls180") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"3172f0f7b542fa1bd2e7a46f1b2e62c848f")+
	"0d9a979e546902ad3f3b6285e27ca";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS180_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/LS180_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",sha,directory);
      name="ls180";
      mode=eos_sn_oo::ls_mode;
      fname=directory+"/LS180_234r_136t_50y_analmu_20091212_SVNr26.h5";
    } else if (sv[1]=="ls220") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"d8c4d4f1315942a663e96fc6452f66d90fc")+
	"87f283e0ed552c8141d1ddba34c19";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS220_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/LS220_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",sha,directory);
      name="ls220";
      mode=eos_sn_oo::ls_mode;
      fname=directory+"/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5";
    } else if (sv[1]=="ls375") {
      cloud_file cf;
      cf.verbose=2;
      std::string sha=((std::string)"13d31df4944b968f1de799e5fc37881eae9")+
	"8cd3f2d3bc14648698661cef35bdd";
      cf.hash_type=cloud_file::sha256;
      cf.get_file_hash("LS375_234r_136t_50y_analmu_20091212_SVNr26.h5",
		       ((string)"https://isospin.roam.utk.edu/")+
		       "public_data/eos_tables/scollapse/LS375_234r_136t_50y_"+
		       "analmu_20091212_SVNr26.h5",sha,directory);
      name="ls375";
      mode=eos_sn_oo::ls_mode;
      fname=directory+"/LS375_234r_136t_50y_analmu_20091212_SVNr26.h5";
    } else if (sv[1]=="acmp_apr_sna") {
      mode=eos_sn_oo::sht_mode;
      fname="APR_0000_rho393_temp133_ye66_gitM180edd5_20190225.h5";
    } else if (sv[1]=="acmp_apr_nse") {
      mode=eos_sn_oo::sht_mode;
      fname="APR_0000_rho393_temp133_ye66_gitM180edd5_20190225.h5";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    oo.verbose=verbose;
    oo.load(fname,mode);
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
    
    static const int nopt=10;
    comm_option_s options[nopt]={
      {0,"ls","Load an EOS in the Lattimer-Swesty format.",
       1,1,"<model>",((string)"Models \"ls\", \"skm\", \"ska\", ")+
       "or \"sk1\".",
       new comm_option_mfptr<ex_eos_sn>(this,&ex_eos_sn::ls_fun),
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
