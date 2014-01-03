/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/gen_sn_eos.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

gen_sn_eos::gen_sn_eos() : cu(o2scl_settings.get_convert_units()) {
  n_nB=0;
  n_Ye=0;
  n_T=0;
  n_oth=0;

  arr[0]=&F;
  arr[1]=&Fint;
  arr[2]=&E;
  arr[3]=&Eint;
  arr[4]=&P;
  arr[5]=&Pint;
  arr[6]=&S;
  arr[7]=&Sint;
  arr[8]=&mun;
  arr[9]=&mup;
  arr[10]=&Z;
  arr[11]=&A;
  arr[12]=&Xn;
  arr[13]=&Xp;
  arr[14]=&Xalpha;
  arr[15]=&Xnuclei;
  for(size_t i=0;i<20;i++) arr[16+i]=&other[i];

  photon.init(0.0,2.0);
  
  electron.init(cu.convert("kg","1/fm",o2scl_mks::mass_electron),2.0);
  muon.init(cu.convert("kg","1/fm",o2scl_mks::mass_muon),2.0);
  include_muons=false;

  verbose=1;

  loaded=false;
  with_leptons_loaded=false;
  baryons_only_loaded=false;
}

gen_sn_eos::~gen_sn_eos() {
  if (loaded) free();
}

void gen_sn_eos::alloc() {
  size_t dim[3]={n_nB,n_Ye,n_T};
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->resize(3,dim);
  }
  return;
}

void gen_sn_eos::free() {
  if (loaded) {
    std::vector<size_t> tmp;
    for(size_t i=0;i<n_base+n_oth;i++) {
      if (arr[i]->total_size()>0) {
	arr[i]->resize(0,tmp);
      }
    }
  }
  loaded=false;
  return;
}

void gen_sn_eos::set_interp_type(size_t interp_type) {
  if (!loaded) {
    O2SCL_ERR("File not loaded in gen_sn_eos::set_interp().",
		  exc_einval);
  }
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_interp_type(interp_type);
  }
  return;
}

int gen_sn_eos::compute_eg() {

  if (verbose>0) {
    cout << "Adding automatically computed electrons and photons." << endl;
  }
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded (loaded=false) in ",
	       "gen_sn_eos::compute_eg().",exc_einval);
  }

  for(size_t i=0;i<n_nB;i++) {
    if (verbose>0 && i%5==0) {
      cout << (i+1) << "/" << n_nB << endl;
    }
    for(size_t j=0;j<n_Ye;j++) {
      for(size_t k=0;k<n_T;k++) {

	double nb1, ye1, T1;
	nb1=E.get_grid(0,i);
	ye1=E.get_grid(1,j);
	T1=E.get_grid(2,k);

	double deg=nb1/pow(T1/hc_mev_fm,3.0);

	photon.massless_calc(T1/hc_mev_fm);
	electron.n=nb1*ye1;
	if (deg>1.0e-5) {
	  relf.density_root->err_nonconv=false;
	  relf.pair_density(electron,T1/hc_mev_fm);
	  relf.density_root->err_nonconv=true;
	} else {
	  electron.ed=0.0;
	  electron.pr=0.0;
	  electron.en=0.0;
	}
	if (include_muons) {
	  muon.mu=electron.mu;
	  relf.pair_mu(muon,T1/hc_mev_fm);
	}
	
	double E_eg=(electron.ed+photon.ed)/nb1*hc_mev_fm;
	double P_eg=(electron.pr+photon.pr)*hc_mev_fm;
	double S_eg=(electron.en+photon.en)/nb1;
	double F_eg=E_eg-T1*S_eg;

	if (include_muons) {
	  E_eg+=muon.ed/nb1*hc_mev_fm;
	  P_eg+=muon.pr*hc_mev_fm;
	  S_eg+=muon.en/nb1;
	  F_eg+=muon.ed/nb1*hc_mev_fm-T1*muon.en/nb1;
	}

	if (baryons_only_loaded==true) {
	  E.set(i,j,k,Eint.get(i,j,k)+E_eg);
	  P.set(i,j,k,Pint.get(i,j,k)+P_eg);
	  S.set(i,j,k,Sint.get(i,j,k)+S_eg);
	  F.set(i,j,k,Fint.get(i,j,k)+F_eg);
	} else {
	  Eint.set(i,j,k,E.get(i,j,k)-E_eg);
	  Pint.set(i,j,k,P.get(i,j,k)-P_eg);
	  Sint.set(i,j,k,S.get(i,j,k)-S_eg);
	  Fint.set(i,j,k,F.get(i,j,k)-F_eg);
	}
      }
    }
  }

  if (baryons_only_loaded==true) {
    with_leptons_loaded=true;
  } else {
    baryons_only_loaded=true;
  }
  
  return 0;
}

void gen_sn_eos::check_free_energy() {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded (loaded=false) in ",
	       "gen_sn_eos::check_free_energy().",exc_einval);
  }

  std::cout << "Checking free energy. " << std::endl;
  
  int i=10, j=10, k=10;
  
  double sum=0.0;
  int count=0;
  
  for(size_t ell=0;ell<400;ell++) {
    if (ell%10==0) cout << ell << "/400" << endl;
    if (ell%3==0) {
      i+=49;
      i=i%n_nB;
    }
    if (ell%3==1) {
      j+=49;
      j=j%n_Ye;
    }
    if (ell%3==2) {
      k+=49;
      k=k%n_T;
    }
    double nb1, ye1, T1;
    nb1=E.get_grid(0,i);
    ye1=E.get_grid(1,j);
    T1=E.get_grid(2,k);

    if (baryons_only_loaded) {
      sum+=fabs(Eint.get(i,j,k)-
		Fint.get(i,j,k)-T1*Sint.get(i,j,k))/fabs(Eint.get(i,j,k));
    } else {
      sum+=fabs(E.get(i,j,k)-
		F.get(i,j,k)-T1*S.get(i,j,k))/fabs(E.get(i,j,k));
    }
    count++;
  }
  
  cout << "loaded: " << loaded << endl;
  cout << "baryons_only_loaded: " << baryons_only_loaded << endl;
  cout << "with_leptons_loaded: " << with_leptons_loaded << endl;
  cout << "sum/count: " << sum/count << endl;
  cout << endl;
  
  return;
}

void gen_sn_eos::beta_eq_sfixed(size_t i, double entr,
				double &nb, double &E_beta, 
				double &P_beta, double &Ye_beta,
				double &Z_beta, double &A_beta,
				double &T_beta) {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded in ",
	       "gen_sn_eos::beta_eq_s4().",exc_einval);
  }
  if (i>=n_nB) {
    O2SCL_ERR2("Too high for baryon grid in ",
	       "gen_sn_eos::beta_eq_s4().",exc_einval);
  }
  if (with_leptons_loaded==false) {
    compute_eg();
  }
  // Get baryon density from grid
  nb=E.get_grid(0,i);
  double Emin=0.0;
  // The electron fraction grid point corresponding to the minimum
  size_t j_found=0;
  // The temperature grid point for s=4 at the Ye grid point j_found
  size_t k_min_j=0;
  for(size_t j=0;j<n_Ye;j++) {
    // For each electron fraction, we need to find the 
    // grid point surrounding s=4
    bool found=false;
    size_t k_found=0;
    for(size_t k=0;k<n_T-1;k++) {
      if (S.get(i,j,k)<entr && S.get(i,j,k+1)>=entr) {
	k_found=k; 
	k=n_T;
	found=true;
      }
    }
    if (found==false) {
      if (S.get(i,j,0)<entr) {
	k_found=n_T-2;
      }
    }
    if (j==0) {
      Emin=E.get(i,j,k_found)+(E.get(i,j,k_found+1)-E.get(i,j,k_found))*
	(entr-S.get(i,j,k_found))/(S.get(i,j,k_found+1)-
				   S.get(i,j,k_found));
    } else {
      double Ethis=E.get(i,j,k_found)+
	(E.get(i,j,k_found+1)-E.get(i,j,k_found))*
	(entr-S.get(i,j,k_found))/(S.get(i,j,k_found+1)-
				   S.get(i,j,k_found));
      if (Ethis<Emin) {
	j_found=j;
	Emin=Ethis;
	k_min_j=k_found;
      }
    }
  }

  // Interpolate final results
  E_beta=Emin;
  P_beta=P.get(i,j_found,k_min_j)+
    (P.get(i,j_found,k_min_j+1)-P.get(i,j_found,k_min_j))*
    (entr-S.get(i,j_found,k_min_j))/
    (S.get(i,j_found,k_min_j+1)-S.get(i,j_found,k_min_j));
  Z_beta=Z.get(i,j_found,k_min_j)+
    (Z.get(i,j_found,k_min_j+1)-Z.get(i,j_found,k_min_j))*
    (entr-S.get(i,j_found,k_min_j))/
    (S.get(i,j_found,k_min_j+1)-S.get(i,j_found,k_min_j));
  A_beta=A.get(i,j_found,k_min_j)+
    (A.get(i,j_found,k_min_j+1)-A.get(i,j_found,k_min_j))*
    (entr-S.get(i,j_found,k_min_j))/
    (S.get(i,j_found,k_min_j+1)-S.get(i,j_found,k_min_j));
  Ye_beta=E.get_grid(1,j_found);
  T_beta=E.get_grid(2,k_min_j)+
    (E.get_grid(2,k_min_j+1)-E.get_grid(2,k_min_j))*
    (entr-S.get(i,j_found,k_min_j))/
    (S.get(i,j_found,k_min_j+1)-S.get(i,j_found,k_min_j));
	
  return;
}

void ls_eos::load(std::string fname) {
  
  double dtemp;

  if (loaded) free();

  std::ifstream fin;
  fin.open(fname.c_str());
      
  if (verbose>0) {
    std::cout << "Loading EOS from file '" << fname << "'." << std::endl;
  }

  // Read grid size and allocate memory

  fin >> n_nB >> n_T >> n_Ye;
  n_oth=11;

  double *nb=new double[n_nB];
  double *ye=new double[n_Ye];
  double *temp=new double[n_T];
  std::vector<double> grid, tgrid;

  alloc();
  
  // Read and set grid
  for(size_t i=0;i<n_nB;i++) {
    fin >> nb[i];
    grid.push_back(nb[i]);
  }
  for(size_t i=0;i<n_T;i++) {
    fin >> temp[i];
    tgrid.push_back(temp[i]);
  }
  for(size_t i=0;i<n_Ye;i++) {
    fin >> ye[i];
    grid.push_back(ye[i]);
  }
  // Reorder Ye and T
  for(size_t i=0;i<n_T;i++) {
    grid.push_back(tgrid[i]);
  }
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_grid_packed(grid);
  }

  // Read data into tensor objects
  
  for(size_t l=0;l<26;l++) {

    if (verbose>0) {
      std::cout << "Reading data section (" << l+1 << "/26)" << std::endl;
    }

    for(size_t k=0;k<n_Ye;k++) {
      for(size_t j=0;j<n_T;j++) {
	for(size_t i=0;i<n_nB;i++) {
	  fin >> dtemp;
	  if (l==0) P.set(i,k,j,dtemp);
	  else if (l==1) F.set(i,k,j,dtemp);
	  else if (l==2) S.set(i,k,j,dtemp);
	  else if (l==3) E.set(i,k,j,dtemp);
	  else if (l==4) mun.set(i,k,j,dtemp);
	  else if (l==5) mup.set(i,k,j,dtemp);
	  else if (l==6) Xn.set(i,k,j,dtemp);
	  else if (l==7) Xp.set(i,k,j,dtemp);
	  else if (l==8) {
	    Xalpha.set(i,k,j,dtemp);
	    Xnuclei.set(i,k,j,1.0-Xn.get(i,k,j)-
			Xp.get(i,k,j)-Xalpha.get(i,k,j));
	  }
	  else if (l==9) Pint.set(i,k,j,dtemp);
	  else if (l==10) Fint.set(i,k,j,dtemp);
	  else if (l==11) Sint.set(i,k,j,dtemp);
	  else if (l==12) {
	    double dtemp2=dtemp-Eint.get_grid(1,k)*
	      (cu.convert("kg","1/fm",o2scl_mks::mass_neutron)-
	       cu.convert("kg","1/fm",o2scl_mks::mass_proton))*hc_mev_fm;
	    Eint.set(i,k,j,dtemp2);
	  }
	  else if (l==13) other[0].set(i,k,j,dtemp);
	  else if (l==14) other[1].set(i,k,j,dtemp);
	  else if (l==15) other[2].set(i,k,j,dtemp);
	  else if (l==16) other[3].set(i,k,j,dtemp);
	  else if (l==17) other[4].set(i,k,j,dtemp);
	  else if (l==18) other[5].set(i,k,j,dtemp);
	  else if (l==19) other[6].set(i,k,j,dtemp);
	  else if (l==20) A.set(i,k,j,dtemp);
	  else if (l==21) {
	    Z.set(i,k,j,dtemp*A.get(i,k,j));
	  } else if (l==22) other[7].set(i,k,j,dtemp);
	  else if (l==23) other[8].set(i,k,j,dtemp);
	  else if (l==24) other[9].set(i,k,j,dtemp);
	  else if (l==25) other[10].set(i,k,j,dtemp);
	  if (i>=n_nB || k>=n_Ye || j>=n_T) {
	    loaded=false;
	    O2SCL_ERR2("Index problem in ",
		       "ls_eos::load().",exc_einval);
	  }
	}
      }
    }
  }
      
  fin.close();

  delete[] nb;
  delete[] ye;
  delete[] temp;

  if (verbose>0) {
    std::cout << "Done." << std::endl;
  }

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=true;
  baryons_only_loaded=true;
  
  set_interp_type(itp_linear);

  return;
}

int ls_eos::check_eg(test_mgr &tm) {
  
  if (!baryons_only_loaded || !with_leptons_loaded) {
    O2SCL_ERR("Not enough data loaded in check_eg().",exc_efailed);
  }

  // Double check lepton and photon contribution
  if (verbose>0) {
    std::cout << "Checking leptons and photons. " << std::endl;
  }
  
  int i=10, j=10, k=10;
  
  double sum=0.0;
  int count=0;

  for(size_t ell=0;ell<400;ell++) {
    if (ell%10==0) cout << ell << "/400" << endl;
    if (ell%3==0) {
      i+=49;
      i=i%n_nB;
    }
    if (ell%3==1) {
      j+=49;
      j=j%n_Ye;
    }
    if (ell%3==2) {
      k+=49;
      k=k%n_T;
    }
    double nb1, ye1, T1;
    nb1=E.get_grid(0,i);
    ye1=E.get_grid(1,j);
    T1=E.get_grid(2,k);
    electron.n=nb1*ye1;

    cout.precision(3);
    cout.width(3);
    cout << ell << " ";
    cout << nb1 << " " << ye1 << " " << T1 << " " << flush;

    relf.pair_density(electron,T1/hc_mev_fm);
    photon.massless_calc(T1/hc_mev_fm);
    
    double E_eg=(electron.ed+photon.ed)/nb1*hc_mev_fm-1.3*ye1;
    double P_eg=(electron.pr+photon.pr)*hc_mev_fm;
    double S_eg=(electron.en+photon.en)/nb1;
    double F_eg=E_eg-T1*S_eg;

    cout.setf(ios::showpos);
    sum+=fabs(E_eg-(E.interp(nb1,ye1,T1)-Eint.interp(nb1,ye1,T1)))/
      fabs(E.interp(nb1,ye1,T1)-Eint.interp(nb1,ye1,T1));
    cout << (E_eg-(E.interp(nb1,ye1,T1)-Eint.interp(nb1,ye1,T1)))/
      fabs(E.interp(nb1,ye1,T1)-Eint.interp(nb1,ye1,T1)) << " ";
    sum+=fabs(P_eg-(P.interp(nb1,ye1,T1)-Pint.interp(nb1,ye1,T1)))/
      fabs(P.interp(nb1,ye1,T1)-Pint.interp(nb1,ye1,T1));
    cout << (P_eg-(P.interp(nb1,ye1,T1)-Pint.interp(nb1,ye1,T1)))/
      fabs(P.interp(nb1,ye1,T1)-Pint.interp(nb1,ye1,T1)) << " ";
    sum+=fabs(S_eg-(S.interp(nb1,ye1,T1)-Sint.interp(nb1,ye1,T1)))/
      fabs(S.interp(nb1,ye1,T1)-Sint.interp(nb1,ye1,T1));
    cout << (S_eg-(S.interp(nb1,ye1,T1)-Sint.interp(nb1,ye1,T1)))/
      fabs(S.interp(nb1,ye1,T1)-Sint.interp(nb1,ye1,T1)) << " ";
    sum+=fabs(F_eg-(F.interp(nb1,ye1,T1)-Fint.interp(nb1,ye1,T1)))/
      fabs(F.interp(nb1,ye1,T1)-Fint.interp(nb1,ye1,T1));
    cout << (F_eg-(F.interp(nb1,ye1,T1)-Fint.interp(nb1,ye1,T1)))/
      fabs(F.interp(nb1,ye1,T1)-Fint.interp(nb1,ye1,T1)) << endl;
    cout.unsetf(ios::showpos);

    if (ell==6 || ell==32 || ell==13) {
      cout.precision(8);
      cout << nb1 << " " << ye1 << " " << T1 << " " 
	   << E.interp(nb1,ye1,T1) << " " << Eint.interp(nb1,ye1,T1) << endl;
      cout << E.interp(nb1,ye1,T1)-Eint.interp(nb1,ye1,T1) << " "
	   << E_eg << " " << 1.3*ye1 << " "
	   << electron.n*electron.m/nb1*hc_mev_fm << " " << electron.mu << endl;
      cout.precision(3);
    }
    
    count+=4;
  }
  
  cout << sum/count << endl;
  cout << endl;
  tm.test_gen((sum/count)<1.0e-1,"electron and photons");
  
  return 0;
}

/*
  This function reads the HDF5 EOS tables generated by E. O'Connor
  and C. Ott in \ref OConnor10. The tables are available from 
  
  http://stellarcollapse.org/equationofstate
  
  and are available under a creative commons
  attribution-noncommercial-share alike license. This \o2 code to
  read those tables is licensed (along with all \o2 code) under
  the GPLv3 license (with permission from Evan O'Connor).
*/
void oo_eos::load(std::string fname, size_t mode) {
  
  double dtemp;
  
  if (loaded) free();

  std::ifstream fin;
  fin.open(fname.c_str());
      
  if (verbose>0) {
    std::cout << "Loading EOS from file '" << fname << "'." << std::endl;
  }

  // Read grid size and allocate memory

  hdf_file hf;
  hf.open(fname);

  int inb, it, iye;
  hf.geti("pointsrho",inb);
  hf.geti("pointstemp",it);
  hf.geti("pointsye",iye);
  hf.getd("energy_shift",energy_shift);

  n_nB=inb;
  n_T=it;
  n_Ye=iye;
  if (verbose>0) {
    cout << n_nB << " " << n_T << " " << n_Ye << " " 
	 << energy_shift << endl;
  }

  n_oth=8;
  if (hfsl_mode) n_oth+=3;

  alloc();

  std::vector<double> ye_grid, t_grid, grid;

  hf.getd_vec_copy("logrho",rho);
  hf.getd_vec_copy("logtemp",t_grid);
  hf.getd_vec_copy("ye",ye_grid);
  
  for(size_t i=0;i<rho.size();i++) {
    // Undo the log
    rho[i]=pow(10.0,rho[i]);
    // Convert from g/cm^3 to baryon density through the 
    // atomic mass unit
    double nb=rho[i]/o2scl_cgs::unified_atomic_mass/1.0e39;
    grid.push_back(nb);
  }
  for(size_t i=0;i<ye_grid.size();i++) grid.push_back(ye_grid[i]);
  for(size_t i=0;i<t_grid.size();i++) grid.push_back(pow(10.0,t_grid[i]));

  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_grid_packed(grid);
  }
  
  size_t nloop=19;
  vector<string> names;
  vector<size_t> indices;

  names.push_back("Abar");
  names.push_back("Xa");
  names.push_back("Xh");
  names.push_back("Xn");
  names.push_back("Xp");
  names.push_back("Zbar");
  names.push_back("cs2");
  names.push_back("dedt");
  names.push_back("dpderho");
  names.push_back("dpdrhoe");
  names.push_back("entropy");
  names.push_back("gamma");
  names.push_back("logenergy");
  names.push_back("logpress");
  names.push_back("mu_e");
  names.push_back("mu_n");
  names.push_back("mu_p");
  names.push_back("muhat");
  names.push_back("munu");
  indices.push_back(11);
  indices.push_back(14);
  indices.push_back(15);
  indices.push_back(12);
  indices.push_back(13);
  indices.push_back(10);
  indices.push_back(16);
  indices.push_back(17);
  indices.push_back(18);
  indices.push_back(19);
  indices.push_back(6);
  indices.push_back(20);
  indices.push_back(2);
  indices.push_back(4);
  indices.push_back(21);
  indices.push_back(8);
  indices.push_back(9);
  indices.push_back(22);
  indices.push_back(23);
  
  if (mode==hfsl_mode) {
    nloop+=3;
    names.push_back("X3he");
    names.push_back("X4li");
    names.push_back("Xt");
    indices.push_back(24);
    indices.push_back(25);
    indices.push_back(26);
  }
		  
  for(size_t i=0;i<nloop;i++) {
    if (verbose>0) {
      cout << "Loading data section " << i 
	   << " for quantity " << names[i] << endl;
    }
    tensor3<> dat;
    hf.getd_ten(names[i],dat);
    for(size_t j=0;j<n_nB;j++) {
      for(size_t k=0;k<n_Ye;k++) {
	for(size_t m=0;m<n_T;m++) {
	  if (i==12) {
	    // For log energy, first undo the log and add the shift
	    double energy=-energy_shift+pow(10.0,dat.get(k,m,j));
	    // Multiply by atomic mass unit to get erg
	    energy*=o2scl_cgs::unified_atomic_mass;
	    // Then convert to MeV
	    energy*=1.0e-6/o2scl_cgs::electron_volt;
	    // Set the new value
	    arr[indices[i]]->set(j,k,m,energy);
	  } else if (i==13) {
	    // For log pressure, first undo the log
	    double press=pow(10.0,dat.get(k,m,j));
	    // Then convert to MeV/fm^3
	    press*=1.0e-45/o2scl_cgs::electron_volt;
	    arr[indices[i]]->set(j,k,m,press);
	  } else if (i==15) {
	    // Neutron chemical potential
	    if (mode==ls_mode) {
	      arr[indices[i]]->set(j,k,m,dat.get(k,m,j));
	    } else {
	      arr[indices[i]]->set
		(j,k,m,dat.get(k,m,j)+938.0-cu.convert
		 ("kg","MeV",o2scl_mks::mass_neutron));
	    }
	  } else if (i==16) {
	    // Proton chemical potential
	    if (mode==ls_mode) {
	      arr[indices[i]]->
		set(j,k,m,dat.get(k,m,j)+
		    cu.convert("kg","MeV",o2scl_mks::mass_neutron)-
		    cu.convert("kg","MeV",o2scl_mks::mass_proton));
	    } else {
	      arr[indices[i]]->set
		(j,k,m,dat.get(k,m,j)+938.0-
		 cu.convert("kg","MeV",o2scl_mks::mass_proton));
	    }
	  } else {
	    arr[indices[i]]->set(j,k,m,dat.get(k,m,j));
	  }
	}
      }
    }
  }

  hf.close();

  if (verbose>0) {
    std::cout << "Done." << std::endl;
  }

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=true;
  baryons_only_loaded=false;
  
  set_interp_type(itp_linear);

  return;
}

void stos_eos::load(std::string fname, size_t mode) {

  size_t ndat=17;
  if (mode==quark_mode) ndat++;

  if (loaded) free();

  std::ifstream fin;
  std::string tstr;

  double dtemp;

  fin.open(fname.c_str());

  if (verbose>0) {
    std::cout << "Loading EOS from file '" << fname << "'." << std::endl;
  }
  
  n_nB=104;
  n_T=31;
  n_Ye=71;
  n_oth=5;
  if (mode==quark_mode) n_oth++;

  alloc();

  double *ye=new double[n_Ye];
  std::vector<double> grid;

  for(size_t k=0;k<n_nB;k++) {
    // We set the density grid later, fill with zeros for now
    grid.push_back(0.0);
  }

  for(size_t j=0;j<n_Ye;j++) {
    ye[j]=pow(10.0,((double)j)*0.025-2.0);
    grid.push_back(ye[j]);
  }
  // Temperature grid from Matthias Hempel
  double temp[31]={0.1,0.12,0.15,0.2,0.25,0.32,0.4,0.5,0.63,0.8,
		   1.0,1.2,1.5,2.0,2.5,3.2,4.0,5.0,6.3,8.0,10.0,
		   12.0,15.0,20.0,25.0,32.0,40.0,50.0,63.0,80.0,100.0};
  for(size_t i=0;i<n_T;i++) {
    grid.push_back(temp[i]);
  }

  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_grid_packed(grid);
  }

  // Set columns not available to zero
  E.set_all(0.0);
  P.set_all(0.0);
  S.set_all(0.0);
  F.set_all(0.0);

  // Create a table to store each section
  table<> t(n_nB);
  t.set_interp_type(itp_linear);
  t.set_nlines(n_nB);

  // Dummy column names which aren't used
  std::string cols;
  for(size_t ell=0;ell<ndat;ell++) {
    cols+="c"+itos(ell)+" ";
  }
  t.line_of_names(cols);

  // Read the data
  for(size_t i=0;i<n_T;i++) {

    if (i!=0) {
      getline(fin,tstr);
    }
    getline(fin,tstr);
    fin >> tstr;
    fin >> temp[i];

    if (verbose>0) {
      std::cout << "Reading section for temperature=" << temp[i] << "." 
		<< std::endl;
    }

    getline(fin,tstr);
    for(size_t j=0;j<n_Ye;j++) {
      
      for(size_t k=0;k<n_nB;k++) {
	for(size_t ell=0;ell<ndat;ell++) {
	  fin >> dtemp;
	  t.set(ell,k,dtemp);
	}
      }

      if (i==0 && j==0) {
	
	// Set the density grid from the first section
	for(size_t k=0;k<n_nB;k++) {
	  grid[k]=t[1][k];
	}
	for(size_t k=0;k<n_base+n_oth;k++) {
	  arr[k]->set_grid_packed(grid);
	}
	    
      } 

      for(size_t ell=0;ell<ndat;ell++) {
	for(size_t k=0;k<n_nB;k++) {

	  if (ell==0) {
	    // The first n_nb elements of grid already contains the
	    // baryon density, so use grid[k] for interpolation.
	    // other[0] is log(rho)
	    other[0].set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==1) {
	    // Baryon density
	    other[1].set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==2) {
	    // Log of proton fraction
	    other[2].set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==3) {
	    // Proton fraction
	    other[3].set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==4) {
	    // The free energy in the table is stored with respect
	    // to the proton mass, so we rescale it
	    double Fint_tmp=t.interp(1,grid[k],ell);
	    double Ye_tmp=Fint.get_grid(1,j);
	    double Fint_new=Fint_tmp+o2scl_const::hc_mev_fm*
	      (cu.convert("kg","1/fm",o2scl_mks::mass_proton)-Ye_tmp*
	       cu.convert("kg","1/fm",o2scl_mks::mass_proton)-
	       (1.0-Ye_tmp)*cu.convert("kg","MeV",o2scl_mks::mass_neutron));
	    Fint.set(k,j,i,Fint_new);
	  } else if (ell==5) {
	    // The internal energy in the table is stored with respect
	    // to the atomic mass, so we rescale it
	    double Eint_tmp=t.interp(1,grid[k],ell);
	    double Ye_tmp=Eint.get_grid(1,j);
	    double Eint_new=Eint_tmp+
	      (cu.convert("kg","MeV",o2scl_mks::unified_atomic_mass)-
	       Ye_tmp*cu.convert("kg","MeV",o2scl_mks::mass_proton)-
	       (1.0-Ye_tmp)*cu.convert("kg","MeV",o2scl_mks::mass_neutron));
	    Eint.set(k,j,i,Eint_new);
	  } else if (ell==6) {
	    Sint.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==7) {
	    A.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==8) {
	    Z.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==9) {
	    // Nucleon effective mass
	    other[4].set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==10) {
	    Xn.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==11) {
	    Xp.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==12) {
	    Xalpha.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==13) {
	    Xnuclei.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==14) {
	    Pint.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==15) {
	    mun.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==16) {
	    mup.set(k,j,i,t.interp(1,grid[k],ell));
	  } else if (ell==16) {
	    other[5].set(k,j,i,t.interp(1,grid[k],ell));
	  } 
	  if (k>=n_nB || j>=n_Ye || i>=n_T) {
	    loaded=false;
	    O2SCL_ERR2("Index problem in ",
		       "stos_eos::load().",exc_einval);
	  }
	}
      }

      // Read the blank line between sections
      getline(fin,tstr);
    }
  }

  fin.close();

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=false;
  baryons_only_loaded=true;

  set_interp_type(itp_linear);

  {
    // Double check the grid 
    if (verbose>0) {
      std::cout << "Checking grid. " << std::endl;
    }

    int i=10, j=10, k=10;

    test_mgr tm;
    tm.set_output_level(0);

    for(size_t ell=0;ell<400;ell++) {
      if (ell%3==0) {
	i+=49;
	i=i%n_nB;
      }
      if (ell%3==1) {
	j+=49;
	j=j%n_Ye;
      }
      if (ell%3==2) {
	k+=49;
	k=k%n_T;
      }
      double nb1, ye1, T1;
      nb1=E.get_grid(0,i);
      ye1=E.get_grid(1,j);
      tm.test_rel(nb1,nB.get(i,j,k),1.0e-10,"nb1");
      tm.test_rel(ye1,Yp.get(i,j,k),1.0e-6,"ye1");
    }

    if (tm.report()==false) {
      O2SCL_ERR("Function stos_eos::load() failed.",exc_efailed);
    }

  }

  if (verbose>0) {
    std::cout << "Done." << std::endl;
  }

  return;
}

void sht_eos::load(std::string fname, size_t mode) {

  if (verbose>0) {
    cout << "In sht_eos::load(), loading EOS from file\n\t'" 
	 << fname << "'." << endl;
  }

  // Commenting this out so we can load both 17 and 17b
  //if (loaded) free();
  
  std::ifstream fin;
  fin.open(fname.c_str());
      
  if (mode==mode_17 || mode==mode_17b) {
    n_nB=336;
    n_T=109;
    n_Ye=53;
    n_oth=5;
  } else {
    // This works for both FSU2.1 and NL3
    n_nB=328;
    n_T=109;
    n_Ye=53;
    n_oth=5;
  }

  tensor_grid3 *Eptr, *Fptr, *Sptr, *Pptr;
  if (mode==mode_17b || mode==mode_21b || mode==mode_NL3b) {
    Eptr=&Eint;
    Pptr=&Pint;
    Sptr=&Sint;
    Fptr=&Fint;
  } else {
    Eptr=&E;
    Pptr=&P;
    Sptr=&S;
    Fptr=&F;
  }

  double dtemp;
  double *nb=new double[n_nB];
  double *ye=new double[n_Ye];
  double *temp=new double[n_T];
  std::vector<double> grid;

  alloc();
  
  for(size_t i=0;i<n_nB;i++) {
    nb[i]=pow(10.0,-8.0+((double)i)*0.025);
    grid.push_back(nb[i]);
  }
  ye[0]=0.0;
  grid.push_back(ye[0]);
  for(size_t i=1;i<n_Ye;i++) {
    ye[i]=0.05+((double)(i-1)*0.01);
    grid.push_back(ye[i]);
  }
  temp[0]=0.0;
  grid.push_back(temp[0]);
  for(size_t i=1;i<n_T;i++) {
    temp[i]=pow(10.0,-0.8+((double)(i-1))*0.025);
    grid.push_back(temp[i]);
  }
  for(size_t k=0;k<n_base+n_oth;k++) {
    arr[k]->set_grid_packed(grid);
  }
  
  string stmp;

  for(size_t j=0;j<n_T;j++) {
    if (verbose>0) {
      cout << "Loading section for temperature=" << temp[j] << "." << endl;
    }
    for(size_t k=0;k<n_Ye;k++) {
      for(size_t i=0;i<n_nB;i++) {
	for(size_t ell=0;ell<16;ell++) {
	  fin >> stmp;
	  dtemp=o2scl::stod(stmp);
	  if (ell==0) {
	    // Temperature
	    other[0].set(i,k,j,dtemp);
	  } else if (ell==1) {
	    // Proton fraction
	    other[1].set(i,k,j,dtemp);
	    // Baryon density
	  } else if (ell==2) {
	    other[2].set(i,k,j,dtemp);
	  } else if (ell==3) {
	    // Rescale it relative to Ye*mp+(1-Ye)*mn
	    double Ye_tmp=Fptr->get_grid(1,k);
	    double dtemp2=dtemp+939.0-o2scl_const::hc_mev_fm*
	      (Ye_tmp*cu.convert("kg","1/fm",o2scl_mks::mass_proton)+
	       (1.0-Ye_tmp)*cu.convert("kg","1/fm",o2scl_mks::mass_neutron));
	    Fptr->set(i,k,j,dtemp2);
	  } else if (ell==4) {
	    Pptr->set(i,k,j,dtemp);
	  } else if (ell==5) {
	    Sptr->set(i,k,j,dtemp);
	    Eptr->set(i,k,j,Fptr->get(i,k,j)+temp[j]*Sptr->get(i,k,j));
	  } else if (ell==6) {
	    mun.set(i,k,j,dtemp+939-
		    cu.convert("kg","1/fm",o2scl_mks::mass_neutron)*
		    hc_mev_fm);
	  } else if (ell==7) {
	    mup.set(i,k,j,dtemp+939-
		    cu.convert("kg","1/fm",o2scl_mks::mass_proton)*
		    hc_mev_fm);
	  } else if (ell==8) {
	    // Electron chemical potential (this already includes
	    // the electron mass)
	    other[3].set(i,k,j,dtemp);
	  } else if (ell==9) {
	    A.set(i,k,j,dtemp);
	  } else if (ell==10) {
	    Z.set(i,k,j,dtemp);
	  } else if (ell==11) {
	    Xn.set(i,k,j,dtemp);
	  } else if (ell==12) {
	    Xp.set(i,k,j,dtemp);
	  } else if (ell==13) {
	    Xalpha.set(i,k,j,dtemp);
	  } else if (ell==14) {
	    Xnuclei.set(i,k,j,dtemp);
	  } else {
	    // Nucleon effective mass
	    other[4].set(i,k,j,dtemp);
	  }
	  if (i>=n_nB || j>=n_T || k>=n_Ye) {
	    loaded=false;
	    O2SCL_ERR2("Index problem in ",
		       "sht_eos::load().",exc_einval);
	  }
	}
      }
    }
    
  }
      
  fin.close();

  if (check_grid) {
    // Double check the grid 
    if (verbose>0) {
      cout << "Checking grid in sht_eos::load()." << endl;
    }

    int i=10, j=10, k=10;

    test_mgr tm;
    tm.set_output_level(0);

    for(size_t ell=0;ell<400;ell++) {
      if (ell%3==0) {
	i+=49;
	i=i%n_nB;
      }
      if (ell%3==1) {
	j+=49;
	j=j%n_Ye;
      }
      if (ell%3==2) {
	k+=49;
	k=k%n_T;
      }
      double nb1, ye1, T1;
      nb1=E.get_grid(0,i);
      ye1=E.get_grid(1,j);
      T1=E.get_grid(2,k);
      tm.test_rel(nb1,other[2].get(i,j,k),4.0e-6,"nb1");
      tm.test_rel(ye1,other[1].get(i,j,k),1.0e-6,"ye1");
      tm.test_rel(T1,other[0].get(i,j,k),1.0e-6,"T1");
    }

    if (tm.report()==false) {
      O2SCL_ERR("Grid check in sht_eos::load() failed.",exc_efailed);
    }

  }

  loaded=true;
  if (mode==mode_17b || mode==mode_21b || mode==mode_NL3b) {
    with_leptons_loaded=false;
    baryons_only_loaded=true;
  } else {
    with_leptons_loaded=true;
    baryons_only_loaded=false;
  }

  if (verbose>0) {
    std::cout << "Done in sht_eos::load()." << std::endl;
  }

  return;
}

void hfsl_eos::load(std::string fname) {

  if (loaded) free();

  std::ifstream fin;
  std::string tstr;

  if (verbose>0) {
    cout << "Loading EOS from file '" << fname << "'." << endl;
  }
  
  fin.open(fname.c_str());

  n_nB=326;
  n_T=81;
  n_Ye=60;
  n_oth=7;

  double dtemp;
  double *nb=new double[n_nB];
  double *ye=new double[n_Ye];
  double *temp=new double[n_T];
  std::vector<double> grid;

  alloc();

  for(size_t i=0;i<n_nB;i++) {
    nb[i]=pow(10.0,((double)i)*0.04-12);
    grid.push_back(nb[i]);
  }
  for(size_t i=0;i<n_Ye;i++) {
    ye[i]=0.01*((double)(i+1));
    grid.push_back(ye[i]);
  }
  for(size_t i=0;i<n_T;i++) {
    temp[i]=pow(10.0,((double)i)*0.04-1.0);
    grid.push_back(temp[i]);
  }
  for(size_t ell=0;ell<n_base+n_oth;ell++) {
    arr[ell]->set_grid_packed(grid);
  }

  for(size_t i=0;i<n_T;i++) {
    if (i!=0) {
      getline(fin,tstr);
    }
    getline(fin,tstr);
    fin >> tstr;
    fin >> temp[i];

    if (verbose>0) {
      cout << "Reading section for temperature=" 
	   << temp[i] << "." << endl;
    }

    getline(fin,tstr);
    for(size_t j=0;j<n_Ye;j++) {
      for(size_t k=0;k<n_nB;k++) {
        for(size_t ell=0;ell<19;ell++) {
	  fin >> dtemp;
	  if (ell==0) {
	    other[0].set(k,j,i,dtemp);
	  } else if (ell==1) {
	    other[1].set(k,j,i,dtemp);
	  } else if (ell==2) {
	    other[2].set(k,j,i,dtemp);
	  } else if (ell==3) {
	    other[3].set(k,j,i,dtemp);
	  } else if (ell==4) {
	    // The free energy in the table is stored with respect
	    // to the proton mass, so we rescale it
	    double Ye_tmp=Fint.get_grid(1,j);
	    double dtemp2=dtemp+o2scl_const::hc_mev_fm*
	      (cu.convert("kg","1/fm",o2scl_mks::mass_proton)-
	       Ye_tmp*cu.convert("kg","1/fm",o2scl_mks::mass_proton)-
	       (1.0-Ye_tmp)*cu.convert("kg","1/fm",o2scl_mks::mass_neutron));
	    Fint.set(k,j,i,dtemp2);
	  } else if (ell==5) {
	    // The internal energy in the table is stored with respect
	    // to the atomic mass, so we rescale it
	    double Ye_tmp=Eint.get_grid(1,j);
	    double dtemp2=dtemp+o2scl_const::hc_mev_fm*
	      (cu.convert("kg","1/fm",o2scl_mks::unified_atomic_mass)-
	       Ye_tmp*cu.convert("kg","1/fm",o2scl_mks::mass_proton)-
	       (1.0-Ye_tmp)*cu.convert("kg","1/fm",o2scl_mks::mass_neutron));
	    Eint.set(k,j,i,dtemp2);
	  } else if (ell==6) {
	    Sint.set(k,j,i,dtemp);
	  } else if (ell==7) {
	    A.set(k,j,i,dtemp);
	  } else if (ell==8) {
	    Z.set(k,j,i,dtemp);
	  } else if (ell==9) {
	    other[4].set(k,j,i,dtemp);
	  } else if (ell==10) {
	    Xn.set(k,j,i,dtemp);
	  } else if (ell==11) {
	    Xp.set(k,j,i,dtemp);
	  } else if (ell==12) {
	    Xalpha.set(k,j,i,dtemp);
	  } else if (ell==13) {
	    Xnuclei.set(k,j,i,dtemp);
	  } else if (ell==14) {
	    Pint.set(k,j,i,dtemp);
	  } else if (ell==15) {
	    mun.set(k,j,i,dtemp);
	  } else if (ell==16) {
	    mup.set(k,j,i,dtemp);
	  } else if (ell==17) {
	    other[5].set(k,j,i,dtemp);
	  } else if (ell==18) {
	    other[6].set(k,j,i,dtemp);
	  } 
          if (k>=n_nB || j>=n_Ye || i>=n_T) {
            loaded=false;
            O2SCL_ERR("Index problem in hfsl_eos::load().",exc_einval);
          }
        }
      }
      getline(fin,tstr);
    }
  }

  fin.close();

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=false;
  baryons_only_loaded=true;

  set_interp_type(itp_linear);

  if (check_grid) {

    // Double check the grid 
    if (verbose>0) {
      cout << "Checking grid. " << endl;
    }

    int i=10, j=10, k=10;

    test_mgr tm;
    tm.set_output_level(2);

    for(size_t ell=0;ell<400;ell++) {
      if (ell%3==0) {
	i+=49;
	i=i%n_nB;
      }
      if (ell%3==1) {
	j+=49;
	j=j%n_Ye;
      }
      if (ell%3==2) {
	k+=49;
	k=k%n_T;
      }
      double nb1, ye1, T1;
      nb1=other[0].get_grid(0,i);
      ye1=other[0].get_grid(1,j);
      T1=other[0].get_grid(2,k);
      double nb2=other[1].get(i,j,k);
      double ye2=other[3].get(i,j,k);
      if (nb2>0.0) {
	tm.test_rel(nb1,nb2,1.0e-6,"nb1");
	tm.test_rel(ye1,ye2,1.0e-6,"ye1");
      }
      if (verbose>1) {
	cout << nb1 << " " << nb2 << " " << fabs(nb1-nb2)/nb2 << " " 
	     << ye1 << " " << ye2 << " " << fabs(ye1-ye2)/ye2 << endl;
      }
    }

    if (tm.report()==false) {
      O2SCL_ERR("Function hfsl_eos::load() failed.",exc_efailed);
    }

  }

  if (verbose>0) {
    std::cout << "Done." << std::endl;
  }

  return;
}
