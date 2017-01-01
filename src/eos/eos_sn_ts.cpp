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
#include <o2scl/eos_sn_gen.h>
#include <o2scl/cli.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

class test_class {

public:

  /// Desc
  string home_dir;

  /// Desc
  int verbose;

  /// Desc
  eos_sn_gen *genp;

  /// Desc
  string name;

  test_class() {
    home_dir="";
    verbose=1;
    genp=0;
  }

  /// Desc
  int ls_fun(std::vector<std::string> &sv, bool itive_com) {

    test_mgr t;
    
    string fname=home_dir;

    if (sv[1]=="ls") {
      fname+="/pkgs/ls_eos2/ls.dat";
      name="ls_ls";
    } else if (sv[1]=="skm") {
      fname+="/pkgs/ls_eos2/skm.dat";
      name="ls_skm";
    } else if (sv[1]=="ska") {
      fname+="/pkgs/ls_eos2/ska.dat";
      name="ls_ska";
    } else if (sv[1]=="sk1") {
      fname+="/pkgs/ls_eos2/sk1.dat";
      name="ls_sk1";
    } else {
      O2SCL_ERR("Need EOS type.",exc_efailed);
    }
    
    ls_eos ls;  
    ls.verbose=verbose;
    ls.load(fname);
    genp=&ls;
    //ls.check_eg(t);
    //ls.check_free_energy();

    if (false) {
      
      double cf=(o2scl_cgs::mass_neutron+o2scl_cgs::mass_proton)/2.0*1.0e39;
      double nb, ye, T;
      
      cout << "Look at the equilibrium nuclei for nearly zero T and Ye=0.5:" 
	   << endl;
      cout << "rho       Ye        T         "
	   << "Z         A         N_skin: " << endl;
      for(double rho=1.0e9;rho<=1.0e15;rho*=2.0) {
	cout << rho << " " << 0.5 << " " << 0.32 << " "
	     << ls.Z.interp(rho/cf,0.5,0.32) << " "
	     << ls.A.interp(rho/cf,0.5,0.32) << " "
	     << ls.Nskin.interp(rho/cf,0.5,0.32) << endl;
      }
      cout << endl;
      
      for(nb=1.0e-4;nb<=1.01;nb*=sqrt(10.0)) {
	cout << nb << " " << ls.E.interp(nb,0.5,0.03) << " "
	     << ls.P.interp(nb,0.5,0.03) << endl;
      }
    }
    
    return 0;

  }

  int old_fun(std::vector<std::string> &sv, bool itive_com) {

    string mode;
    string fname;
    double cf;
    test_mgr t;

    if (mode=="oo") {
      
      oo_eos oo;  
      oo.verbose=2;
      
      /*
	if (((string)argc[3])=="ls") {
	oo.load(fname,0);
	} else {
	oo.load(fname,1);
	}
      */

      for(double nb=1.0e-4;nb<=1.01;nb*=sqrt(10.0)) {
	cout << nb << " " << oo.E.interp(nb,0.5,0.03) << " "
	     << oo.P.interp(nb,0.5,0.03) << endl;
      }

    } else if (mode=="stos") {

      stos_eos stos;
      stos.load(fname);

      double nb, ye, T;
    
      cout << "Look at the equilibrium nuclei for nearly zero T and Ye=0.5:" 
	   << endl;
      cout << "rho       Ye        T         "
	   << "Z         A         X_nuclei: " << endl;
      for(double rho=1.0e9;rho<=1.0e15;rho*=2.0) {
	cout << rho << " " << 0.5 << " " << 0.32 << " "
	     << stos.Z.interp(rho/cf,0.5,0.32) << " "
	     << stos.A.interp(rho/cf,0.5,0.32) << " "
	     << stos.Xnuclei.interp(rho/cf,0.5,0.32) << endl;
      }
      cout << endl;

      stos.compute_eg();

      for(nb=1.0e-4;nb<=1.01;nb*=sqrt(10.0)) {
	cout << nb << " " << stos.E.interp(nb,0.5,0.03) << " "
	     << stos.P.interp(nb,0.5,0.03) << endl;
      }

    } else if (mode=="sht") {

      string fname2="";//argc[3];
      sht_eos sht, sht2;
      int imode=0;//o2scl::stoi(argc[4]);
      sht.load(fname,imode);
      sht2.load(fname2,imode+2);

      //eos_sn_gen_ts sht $(HOME)/pkgs/sht_eos/FSU2.1eos1.01.dat 
      //$(HOME)/pkgs/sht_eos/FSU2.1eosb1.01.dat 1		 
      //>> gen_sn_test.scr

      double nb, ye, T;
    
      if (false) {
	cout << "Look at the equilibrium nuclei for nearly zero T and Ye=0.5:" 
	     << endl;
	cout << "rho       Ye        T         "
	     << "Z         A         X_nuclei: " << endl;
	for(double rho=1.0e9;rho<=1.0e15;rho*=2.0) {
	  cout << rho << " " << 0.5 << " " << 0.32 << " "
	       << sht.Z.interp(rho/cf,0.5,0.32) << " "
	       << sht.A.interp(rho/cf,0.5,0.32) << " "
	       << sht.Xnuclei.interp(rho/cf,0.5,0.32) << endl;
	}
	cout << endl;
	
	cout << "Look at the equilibrium nuclei for nearly zero T and Ye=0.5:" 
	     << endl;
	cout << "rho       Ye        T         "
	     << "Z         A         X_nuclei: " << endl;
	for(double rho=1.0e9;rho<=1.0e15;rho*=2.0) {
	  cout << rho << " " << 0.5 << " " << 0.32 << " "
	       << sht2.Z.interp(rho/cf,0.5,0.32) << " "
	       << sht2.A.interp(rho/cf,0.5,0.32) << " "
	       << sht2.Xnuclei.interp(rho/cf,0.5,0.32) << endl;
	}
	cout << endl;
      }

      {
	cout << "Checking electrons and photons:" << endl;
	
	int i=10, j=10, k=10;
	
	t.set_output_level(0);
	
	for(size_t ell=0;ell<40;ell++) {
	  if (ell%3==0) {
	    i+=49;
	    i=i%sht.n_nB;
	  }
	  if (ell%3==1) {
	    j+=49;
	    j=j%sht.n_Ye;
	  }
	  if (ell%3==2) {
	    k+=49;
	    k=k%sht.n_T;
	  }
	  double nb1, ye1, T1;
	  nb1=sht.E.get_grid(0,i);
	  ye1=sht.E.get_grid(1,j);
	  T1=sht.E.get_grid(2,k);

	  if (ye1>0.0) {

	    sht.electron.n=nb1*ye1;
	    cout << nb1 << " " << ye1 << " " << T1 << " " << std::flush;
	    sht.relf.pair_density(sht.electron,T1/hc_mev_fm);
	    sht.photon.massless_calc(T1/hc_mev_fm);
	    
	    double E_eg=(sht.electron.ed+sht.photon.ed)/nb1*hc_mev_fm;
	    double val=sht.E.interp(nb1,ye1,T1)-sht2.Eint.interp(nb1,ye1,T1);
	    double P_eg=(sht.electron.pr+sht.photon.pr)*hc_mev_fm;
	    double val2=sht.P.interp(nb1,ye1,T1)-sht2.Pint.interp(nb1,ye1,T1);
	    double S_eg=(sht.electron.en+sht.photon.en)/nb1;
	    double val3=sht.S.interp(nb1,ye1,T1)-sht2.Sint.interp(nb1,ye1,T1);

	    // Currently, the entropy doesn't match so well at low
	    // temperatures and high densities and the pressure
	    // doesn't match so well at high temperature and low
	    // densities. 
	    cout << fabs(val-E_eg)/fabs(val) << " ";
	    cout << fabs(val2-P_eg)/fabs(val2) << " ";
	    cout << fabs(val3-S_eg)/fabs(val3) << endl;
	    if (fabs(val2-P_eg)/fabs(val2)>1.0e-2) {
	      cout << "\tP " << sht.P.interp(nb1,ye1,T1) << " " 
		   << sht2.Pint.interp(nb1,ye1,T1) << " " << P_eg << endl;
	    }
	    if (fabs(val3-S_eg)/fabs(val3)>1.0e-2) {
	      cout << "\tS " << sht.S.interp(nb1,ye1,T1) << " " 
		   << sht2.Sint.interp(nb1,ye1,T1) << " " << S_eg << endl;
	    }
	    t.test_rel(fabs(val-E_eg)/fabs(val),0.0,1.0e-2,"e and g");
	    t.test_rel(fabs(val2-P_eg)/fabs(val2),0.0,4.0e-2,"e and g");
	    t.test_rel(fabs(val3-S_eg)/fabs(val3),0.0,1.0e-1,"e and g");
	  }
	       
	}
	
	if (t.report()==false) {
	  O2SCL_ERR("Function sht_eos::load() failed.",exc_efailed);
	}
	
      }
      
    } else if (mode=="hfsl") {

      hfsl_eos hfsl;
      hfsl.load(fname);

      double nb, ye, T;
    
      cout << "Look at the equilibrium nuclei for nearly zero T and Ye=0.5:" 
	   << endl;
      cout << "rho       Ye        T         "
	   << "Z         A         X_nuclei: " << endl;
      for(double rho=1.0e9;rho<=1.0e15;rho*=2.0) {
	cout << rho << " " << 0.5 << " " << 0.32 << " "
	     << hfsl.Z.interp(rho/cf,0.5,0.32) << " "
	     << hfsl.A.interp(rho/cf,0.5,0.32) << " "
	     << hfsl.Xnuclei.interp(rho/cf,0.5,0.32) << endl;
      }
      cout << endl;

    }

    t.report();

    return 0;
  }

  void run(int argc, char *argv[]) {

    home_dir=getenv("HOME");
    
    // ---------------------------------------
    // Specify command-line option object
    
    cli cl;
    cl.prompt="eos_sn_gen_ts> ";
    cl.gnu_intro=false;
    
    // ---------------------------------------
    // Set options
    
    static const int nopt=1;
    comm_option_s options[nopt]={
      {0,"ls","short desc",
       1,1,"",((string)"long ")+"desc.",
       new comm_option_mfptr<test_class>(this,&test_class::ls_fun),
       cli::comm_option_both}
    };
    cl.set_comm_option_vec(nopt,options);
    
    // ---------------------------------------
    // Set parameters

    cli::parameter_int p_verbose;
    p_verbose.i=&verbose;
    p_verbose.help=((string)"Desc.");
    cl.par_list.insert(make_pair("verbose",&p_verbose));

    // ---------------------------------------
    // Process command-line arguments and run
    
    cl.run_auto(argc,argv);

    return;
  }


};

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_class tc;
  tc.run(argc,argv);

  return 0;
}
