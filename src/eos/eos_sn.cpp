/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <o2scl/eos_sn.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>
#include <o2scl/cloud_file.h>
#include <o2scl/classical.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

eos_sn_base::eos_sn_base() : cu(o2scl_settings.get_convert_units()) {
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
  for(size_t i=0;i<30;i++) arr[16+i]=&other[i];

  photon.init(0.0,2.0);
  
  electron.init(cu.convert("kg","1/fm",o2scl_mks::mass_electron),2.0);
  muon.init(cu.convert("kg","1/fm",o2scl_mks::mass_muon),2.0);
  include_muons=false;

  verbose=1;

  loaded=false;
  with_leptons_loaded=false;
  baryons_only_loaded=false;

  m_neut=o2scl_mks::mass_neutron*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0)*
    o2scl_const::hc_mev_fm;
  m_prot=o2scl_mks::mass_proton*
    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0)*
    o2scl_const::hc_mev_fm;
}

eos_sn_base::~eos_sn_base() {
  if (loaded) free();
}

void eos_sn_base::output(std::string fname) {

  if (loaded==false) {
    O2SCL_ERR("Not loaded in eos_sn_base::output().",
	      exc_efailed);
  }
  
  if (verbose>0) {
    cout << "eos_sn_base::output(): Output to file named '"
	 << fname << "'." << endl;
  }
  
  hdf_file hf;

  wordexp_single_file(fname);

  hf.open_or_create(fname);

  // Grid
  hf.set_szt("n_nB",n_nB);
  hf.set_szt("n_Ye",n_Ye);
  hf.set_szt("n_T",n_T);
  hf.setd_vec("nB_grid",nB_grid);
  hf.setd_vec("Ye_grid",Ye_grid);
  hf.setd_vec("T_grid",T_grid);

  // Flags
  hf.seti("baryons_only",data_baryons_only());
  hf.seti("with_leptons",data_with_leptons());
  if (include_muons) {
    hf.seti("include_muons",1);
  } else {
    hf.seti("include_muons",0);
  }

  if (with_leptons_loaded) {
    hdf_output(hf,F,"F");
    hdf_output(hf,E,"E");
    hdf_output(hf,P,"P");
    hdf_output(hf,S,"S");
  }
  if (baryons_only_loaded) {
    hdf_output(hf,Fint,"Fint");
    hdf_output(hf,Eint,"Eint");
    hdf_output(hf,Pint,"Pint");
    hdf_output(hf,Sint,"Sint");
  }

  // Chemical potentials
  hdf_output(hf,mun,"mun");
  hdf_output(hf,mup,"mup");

  // Composition
  hdf_output(hf,Z,"Z");
  hdf_output(hf,A,"A");
  hdf_output(hf,Xn,"Xn");
  hdf_output(hf,Xp,"Xp");
  hdf_output(hf,Xnuclei,"Xnuclei");
  hdf_output(hf,Xalpha,"Xalpha");
  
  // Nucleon masses
  hf.setd("m_neut",m_neut);
  hf.setd("m_prot",m_prot);

  // Other data 
  hf.set_szt("n_oth",n_oth);
  if (n_oth>0) {
    if (oth_names.size()!=n_oth ||
	oth_units.size()!=n_oth) {
      cout << oth_names.size() << " "
	   << oth_units.size() << " " << n_oth << endl;
      O2SCL_ERR("Lists oth_names or oth_units don't match n_oth.",
		o2scl::exc_einval);
    }
    hf.sets_vec("oth_names",oth_names);
    hf.sets_vec("oth_units",oth_units);
    for(size_t i=0;i<n_oth;i++) {
      hdf_output(hf,other[i],oth_names[i]);
    }
  }

  hf.close();

  if (verbose>0) {
    cout << "eos_sn_base::output(): Done with output." << endl;
  }

  return;
}

void eos_sn_base::load(std::string fname, size_t mode) {

  wordexp_single_file(fname);
  
  if (verbose>0) {
    cout << "In eos_sn_base::load(), loading EOS from file\n\t'"
	 << fname << "'." << endl;
  }

  if (loaded) free();
  
  hdf_file hf;

  hf.open(fname);

  // Grid
  hf.get_szt("n_nB",n_nB);
  hf.get_szt("n_Ye",n_Ye);
  hf.get_szt("n_T",n_T);
  hf.getd_vec("nB_grid",nB_grid);
  hf.getd_vec("Ye_grid",Ye_grid);
  hf.getd_vec("T_grid",T_grid);

  // Flags
  int itmp;
  hf.geti("baryons_only",itmp);
  if (itmp==1) baryons_only_loaded=true;
  else baryons_only_loaded=false;
  hf.geti("with_leptons",itmp);
  if (itmp==1) with_leptons_loaded=true;
  else with_leptons_loaded=false;
  hf.geti("include_muons",itmp);
  if (itmp==1) include_muons=true;
  else include_muons=false;

  if (with_leptons_loaded) {
    hdf_input(hf,F,"F");
    hdf_input(hf,E,"E");
    hdf_input(hf,P,"P");
    hdf_input(hf,S,"S");
  }
  if (baryons_only_loaded) {
    hdf_input(hf,Fint,"Fint");
    hdf_input(hf,Eint,"Eint");
    hdf_input(hf,Pint,"Pint");
    hdf_input(hf,Sint,"Sint");
  }

  // Chemical potentials
  hdf_input(hf,mun,"mun");
  hdf_input(hf,mup,"mup");

  // Composition
  hdf_input(hf,Z,"Z");
  hdf_input(hf,A,"A");
  hdf_input(hf,Xn,"Xn");
  hdf_input(hf,Xp,"Xp");
  hdf_input(hf,Xnuclei,"Xnuclei");
  hdf_input(hf,Xalpha,"Xalpha");
  
  // Nucleon masses
  hf.getd("m_neut",m_neut);
  hf.getd("m_prot",m_prot);

  // Other data
  hf.get_szt("n_oth",n_oth);
  if (n_oth>0) {
    hf.gets_vec("oth_names",oth_names);
    hf.gets_vec("oth_units",oth_units);
    for(size_t i=0;i<n_oth;i++) {
      hdf_input(hf,other[i],oth_names[i]);
    }
  }
  
  hf.close();

  loaded=true;

  // It is important that 'loaded' is set to true before the call to
  // set_interp_type().
  set_interp_type(itp_linear);

  if (verbose>0) {
    cout << "Done in eos_sn_base::load()." << endl;
  }

  return;
}

void eos_sn_base::alloc() {
  size_t dim[3]={n_nB,n_Ye,n_T};
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->resize(3,dim);
  }
  return;
}

void eos_sn_base::free() {
  if (loaded) {
    std::vector<size_t> tmp;
    for(size_t i=0;i<n_base+n_oth;i++) {
      if (arr[i]->total_size()>0) {
	arr[i]->resize(0,tmp);
      }
    }
  }
  loaded=false;
  oth_names.clear();
  oth_units.clear();
  return;
}

void eos_sn_base::set_interp_type(size_t interp_type) {
  if (!loaded) {
    O2SCL_ERR("File not loaded in eos_sn_base::set_interp().",
	      exc_einval);
  }
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_interp_type(interp_type);
  }
  return;
}

void eos_sn_base::compute_eg_point(double nB, double Ye, double T,
				   thermo &th, double &mue) {
  
  photon.massless_calc(T/hc_mev_fm);
  electron.n=nB*Ye;
  
  // Provide the initial guess for the electron
  // chemical potential
  double guess2=mue;
  electron.mu=mue;

  relf.err_nonconv=false;
  
  // For lower densities, pair_density() has problems computing the
  // chemical potential from the density, so we loosen the tolerances
  // a little bit.
  relf.alt_solver.tol_rel=1.0e-6;
  
  int retx=relf.pair_density(electron,T/hc_mev_fm);
  if (retx!=0) {

    std::cout << "Failed. Trying again." << std::endl;
    relf.pair_mu(electron,T/hc_mev_fm);
    cout << electron.n << endl;
    relf.calc_mu_ndeg(electron,T/hc_mev_fm,1.0e-18,false);
    cout << electron.n << endl;
    relf.calc_mu_ndeg(electron,T/hc_mev_fm,1.0e-18,true);
    cout << electron.n << endl;
    exit(-1);
    
    double best_mue=guess2;
    double diff=1.0;
    for(double fact=0.2;fact>1.0e-6;fact/=2.0) {
      for(mue=guess2/(1.0+fact);mue<guess2*(1.0+fact);mue*=1.0+fact/20.0) {
        electron.mu=mue;
        relf.pair_mu(electron,T/hc_mev_fm);
        if (fabs(electron.n-nB*Ye)<diff) {
          diff=fabs(electron.n-nB*Ye);
          best_mue=mue;
        }
        cout << mue << " " << electron.n << " " << nB*Ye << endl;
      }
      cout << "fgb: " << fact << " " << guess2 << " " << best_mue << endl;
      guess2=best_mue;
    }
    
    cout << endl;

    if (true) {
      relf.verbose=2;
      //electron.mu=guess2;
      electron.mu=best_mue;
      electron.n=nB*Ye;
      int retx2=relf.pair_density(electron,T/hc_mev_fm);
      cout << "retx2: " << retx2 << endl;
      //char ch;
      //cin >> ch;
      if (retx2!=0) {
        O2SCL_ERR2("Function fermion_rel::pair_density() failed in ",
                   "eos_sn_base::compute_eg_point().",o2scl::exc_efailed);
      }
    } else {
      electron.mu=best_mue;
      relf.pair_mu(electron,T/hc_mev_fm);
      electron.n=nB*Ye;
    }
  }
  
  if (include_muons) {
    muon.mu=electron.mu;
    relf.pair_mu(muon,T/hc_mev_fm);
  }

  th.ed=electron.ed+photon.ed;
  th.pr=electron.pr+photon.pr;
  th.en=electron.en+photon.en;
  
  if (include_muons) {
    th.ed+=muon.ed;
    th.en+=muon.en;
    th.pr+=muon.pr;
  }

  mue=electron.mu;

  return;
}

void eos_sn_base::compute_eg() {

  if (verbose>0) {
    cout << "Adding automatically computed electrons and photons." << endl;
  }
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded (loaded=false) in ",
	       "eos_sn_base::compute_eg().",exc_einval);
  }

  for(int i=n_nB-1;i>=0;i--) {
    if (verbose>0 && i%5==0) {
      cout << (i+1) << "/" << n_nB << endl;
    }
    for(size_t j=0;j<n_Ye;j++) {
      for(size_t k=0;k<n_T;k++) {

	double nb1, ye1, T1;
	nb1=E.get_grid(0,i);
	ye1=E.get_grid(1,j);
	T1=E.get_grid(2,k);

	thermo th;
	double mue2;
	compute_eg_point(nb1,ye1,T1,th,mue2);
	
	double E_eg=th.ed/nb1*hc_mev_fm;
	double P_eg=th.pr*hc_mev_fm;
	double S_eg=th.en/nb1;
	double F_eg=E_eg-T1*S_eg;

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
  
  return;
}

void eos_sn_base::check_free_energy(double &avg) {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded (loaded=false) in ",
	       "eos_sn_base::check_free_energy().",exc_einval);
  }

  if (verbose>0) {
    std::cout << "Checking free energy. " << std::endl;
  }
  
  int i=10, j=10, k=10;
  
  double sum=0.0;
  int count=0;
  
  for(size_t ell=0;ell<400;ell++) {
    if (verbose>0) {
      if (ell%10==0) cout << ell << "/400" << endl;
    }
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
    }
    if (with_leptons_loaded) {
      sum+=fabs(E.get(i,j,k)-
		F.get(i,j,k)-T1*S.get(i,j,k))/fabs(E.get(i,j,k));
    }
    count++;
  }
  
  avg=sum/count;

  if (verbose>0) {
    cout << "loaded: " << loaded << endl;
    cout << "baryons_only_loaded: " << baryons_only_loaded << endl;
    cout << "with_leptons_loaded: " << with_leptons_loaded << endl;
    cout << "sum/count: " << sum/count << endl;
    cout << endl;
  }
  
  return;
}

void eos_sn_base::check_composition(double &max1, double &max2) {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded (loaded=false) in ",
	       "eos_sn_base::check_composition().",exc_einval);
  }

  if (verbose>0) {
    std::cout << "Checking composition. " << std::endl;
  }
  
  int i=10, j=10, k=10;
  
  max1=0.0;
  max2=0.0;
  
  for(size_t ell=0;ell<400;ell++) {
    if (verbose>0) {
      if (ell%10==0) cout << ell << "/400" << endl;
    }
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
    
    double tot=Xp.get(i,j,k)+Xnuclei.get(i,j,k)+
      Xn.get(i,j,k)+Xalpha.get(i,j,k);
    if (fabs(tot-1.0)>max1) max1=fabs(tot-1.0);
    
    double np=nb1*Xp.get(i,j,k);
    double nh=nb1*Xnuclei.get(i,j,k)/A.get(i,j,k);
    double na=nb1*Xalpha.get(i,j,k)/4.0;
    double Ych=(np+Z.get(i,j,k)*nh+2.0*na)/nb1;
    if (fabs(Ych-ye1)>max2) max2=fabs(Ych-ye1);
  }

  if (verbose>0) {
    cout << "Maximum deviation for sum of mass fractions: " 
	 << max1 << endl;
    cout << "Maximum deviation for charge fraction relative to Y_e: " 
	 << max2 << endl;
  }
  
  return;
}

void eos_sn_base::beta_eq_sfixed(double nB, double entr,
				 double &Ye, double &T) {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded in ",
	       "eos_sn_base::beta_eq_sfixed().",exc_einval);
  }
  if (with_leptons_loaded==false) {
    compute_eg();
  }

  // Linear interpolation for the temperature
  interp<> it(itp_linear);

  // Create a vector for the temperature and entropy
  ubvector Tvec(n_T), Svec(n_T);
  for(size_t k=0;k<n_T;k++) {
    Tvec[k]=S.get_grid(2,k);
  }

  // Create vectors to interpolate for the minimum free energy
  vector<double> Yevec, Fvec;
  for(size_t j=0;j<n_Ye;j++) {
    double Ye=F.get_grid(1,j);
    Yevec.push_back(Ye);

    // Create vectors to find the temperature corresponding to the
    // specified entropy
    for(size_t k=0;k<n_T;k++) {
      Svec[k]=S.interp_linear(nB,Ye,Tvec[k]);
    }
    T=it.eval(entr,n_T,Svec,Tvec);

    // Add the free energy at this temperature for this electron 
    // fraction
    Fvec.push_back(F.interp_linear(nB,Ye,T));
  }
  
  // Use those vectors to get Ye
  size_t ix=vector_min_index<vector<double>,double>(n_Ye,Fvec);
  if (ix==0 || ix==n_Ye-1) {
    Ye=Yevec[ix];
    return;
  }
  Ye=quadratic_extremum_x(Yevec[ix-1],Yevec[ix],Yevec[ix+1],
			  Fvec[ix-1],Fvec[ix],Fvec[ix+1]);
  
  // Now compute the temperature for this particular Ye
  for(size_t k=0;k<n_T;k++) {
    Svec[k]=S.interp_linear(nB,Ye,Tvec[k]);
  }
  T=it.eval(entr,n_T,Svec,Tvec);

  return;
}

void eos_sn_base::beta_eq_Tfixed(double nB, double T, double &Ye) {
  
  if (loaded==false) {
    O2SCL_ERR2("No data loaded in ",
	       "eos_sn_base::beta_eq_Tfixed().",exc_einval);
  }
  if (with_leptons_loaded==false) {
    compute_eg();
  }

  // Create vectors to interpolate for the minimum free energy
  vector<double> Yevec, Fvec;
  for(size_t j=0;j<n_Ye;j++) {
    double Ye=F.get_grid(1,j);
    Yevec.push_back(Ye);
    Fvec.push_back(F.interp_linear(nB,Ye,T));
  }
  
  // Use those vectors to get Ye
  size_t ix=vector_min_index<vector<double>,double>(n_Ye,Fvec);
  if (ix==0 || ix==n_Ye-1) {
    Ye=Yevec[ix];
    return;
  }
  Ye=quadratic_extremum_x(Yevec[ix-1],Yevec[ix],Yevec[ix+1],
			  Fvec[ix-1],Fvec[ix],Fvec[ix+1]);
  
  return;
}

double eos_sn_base::check_eg() {
  
  if (!baryons_only_loaded) {
    O2SCL_ERR("Baryons only EOS not loaded in check_eg().",exc_efailed);
  }
  if (!with_leptons_loaded) {
    O2SCL_ERR("With leptons EOS not loaded in check_eg().",exc_efailed);
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

    if (ye1>1.0e-20) {
      electron.n=nb1*ye1;

      cout.precision(3);
      cout.width(3);
      cout << ell << " ";
      cout << nb1 << " " << ye1 << " " << T1 << " " << flush;

      relf.pair_density(electron,T1/hc_mev_fm);
      photon.massless_calc(T1/hc_mev_fm);
    
      double E_eg=(electron.ed+photon.ed)/nb1*hc_mev_fm;
      double P_eg=(electron.pr+photon.pr)*hc_mev_fm;
      double S_eg=(electron.en+photon.en)/nb1;
      double F_eg=E_eg-T1*S_eg;

      cout.setf(ios::showpos);

      double diff, rat;

      // -----------------------------------------------------
      // Compare the energy

      double val_E=E.get(i,j,k);
      double val_Eint=Eint.get(i,j,k);
      diff=E_eg-(val_E-val_Eint);

      // Rescale difference by the maximum of X and Xint
      rat=diff;
      if (fabs(val_E)>fabs(val_Eint)) rat/=val_E;
      else rat/=val_Eint;
    
      sum+=fabs(rat);
      cout << rat << " ";

      // -----------------------------------------------------
      // Compare the pressure

      double val_P=P.get(i,j,k);
      double val_Pint=Pint.get(i,j,k);
      diff=P_eg-(val_P-val_Pint);

      // Rescale difference by the maximum of X and Xint
      rat=diff;
      if (fabs(val_P)>fabs(val_Pint)) rat/=val_P;
      else rat/=val_Pint;
    
      sum+=fabs(rat);
      cout << rat << " ";

      // -----------------------------------------------------
      // Compare the entropy

      if (T1>0.0) {
	double val_S=S.get(i,j,k);
	double val_Sint=Sint.get(i,j,k);
	diff=S_eg-(val_S-val_Sint);
	
	// Rescale difference by the maximum of X and Xint
	rat=diff;
	if (fabs(val_S)>fabs(val_Sint)) rat/=val_S;
	else rat/=val_Sint;
	
	sum+=fabs(rat);
      } else {
	rat=0.0;
      }
      cout << rat << " ";

      // -----------------------------------------------------
      // Compare the free energy

      double val_F=F.get(i,j,k);
      double val_Fint=Fint.get(i,j,k);
      diff=F_eg-(val_F-val_Fint);

      // Rescale difference by the maximum of X and Xint
      rat=diff;
      if (fabs(val_F)>fabs(val_Fint)) rat/=val_F;
      else rat/=val_Fint;
    
      sum+=fabs(rat);
      cout << rat << endl;
    
      // -----------------------------------------------------

      cout.unsetf(ios::showpos);
    
      count+=4;

    }
  }

  if (verbose>0) {
    cout << "Check electrons photons result: " << sum/count << endl;
  }
  
  return sum/count;
}

void eos_sn_ls::load(std::string fname, size_t mode) {
  
  wordexp_single_file(fname);

  if (verbose>0) {
    cout << "In eos_sn_ls::load(), loading EOS from file\n\t'" 
	 << fname << "'." << endl;
  }

  double dtemp;

  if (loaded) free();

  std::ifstream fin;
  fin.open(fname.c_str());
      
  // Read grid size and allocate memory

  fin >> n_nB >> n_T >> n_Ye;
  n_oth=11;
  size_t ndat=26;
  alloc();

  // Read and set grid
  std::vector<double> grid, tgrid;
  for(size_t i=0;i<n_nB;i++) {
    fin >> dtemp;
    grid.push_back(dtemp);
  }
  for(size_t i=0;i<n_T;i++) {
    fin >> dtemp;
    tgrid.push_back(dtemp);
  }
  for(size_t i=0;i<n_Ye;i++) {
    fin >> dtemp;
    grid.push_back(dtemp);
  }
  // Reorder Ye and T
  for(size_t i=0;i<n_T;i++) {
    grid.push_back(tgrid[i]);
  }
  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_grid_packed(grid);
  }

  // Read data into tensor objects
  
  for(size_t l=0;l<ndat;l++) {

    if (verbose>0) {
      std::cout << "Reading data section (" << l+1 << "/" << ndat 
		<< ")" << std::endl;
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
	  // filling factor
	  else if (l==13) other[0].set(i,k,j,dtemp);
	  // baryon density inside nuclei
	  else if (l==14) other[1].set(i,k,j,dtemp);
	  // dPdn
	  else if (l==15) other[2].set(i,k,j,dtemp);
	  // dPdT
	  else if (l==16) other[3].set(i,k,j,dtemp);
	  // dPdY
	  else if (l==17) other[4].set(i,k,j,dtemp);
	  // dsdT
	  else if (l==18) other[5].set(i,k,j,dtemp);
	  // dsdY
	  else if (l==19) other[6].set(i,k,j,dtemp);
	  else if (l==20) A.set(i,k,j,dtemp);
	  else if (l==21) {
	    Z.set(i,k,j,dtemp*A.get(i,k,j));
	  } 
	  // Number of neutrons in skin
	  else if (l==22) other[7].set(i,k,j,dtemp);
	  // baryon density outside nuclei
	  else if (l==23) other[8].set(i,k,j,dtemp);
	  // x_out
	  else if (l==24) other[9].set(i,k,j,dtemp);
	  // mu
	  else if (l==25) other[10].set(i,k,j,dtemp);
	  if (i>=n_nB || k>=n_Ye || j>=n_T) {
	    loaded=false;
	    O2SCL_ERR2("Index problem in ",
		       "eos_sn_ls::load().",exc_einval);
	  }
	}
      }
    }
  }
  fin.close();

  oth_names.push_back("fill");
  oth_names.push_back("nb_in");
  oth_names.push_back("dPdn");
  oth_names.push_back("dPdT");
  oth_names.push_back("dPdY");
  oth_names.push_back("dsdT");
  oth_names.push_back("dsdY");
  oth_names.push_back("Nskin");
  oth_names.push_back("nb_out");
  oth_names.push_back("x_out");
  oth_names.push_back("mu");

  oth_units.push_back("");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("MeV");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("MeV/fm^3");
  oth_units.push_back("1/(MeV*fm^3)");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("");
  oth_units.push_back("MeV");

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=true;
  baryons_only_loaded=true;
  
  if (n_oth!=oth_names.size()) {
    O2SCL_ERR2("Number of names does not match number of data sets ",
	       "in eos_sn_ls::load().",exc_efailed);
  }

  // It is important that 'loaded' is set to true before the call to
  // set_interp_type().
  set_interp_type(itp_linear);

  if (verbose>0) {
    std::cout << "Done in eos_sn_ls::load()." << std::endl;
  }

  return;
}

double eos_sn_ls::check_eg() {
  
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
    sum+=fabs(E_eg-(E.interp_linear(nb1,ye1,T1)-
		    Eint.interp_linear(nb1,ye1,T1)))/
      fabs(E.interp_linear(nb1,ye1,T1)-Eint.interp_linear(nb1,ye1,T1));
    cout << (E_eg-(E.interp_linear(nb1,ye1,T1)-
		   Eint.interp_linear(nb1,ye1,T1)))/
      fabs(E.interp_linear(nb1,ye1,T1)-Eint.interp_linear(nb1,ye1,T1)) << " ";
    sum+=fabs(P_eg-(P.interp_linear(nb1,ye1,T1)-
		    Pint.interp_linear(nb1,ye1,T1)))/
      fabs(P.interp_linear(nb1,ye1,T1)-Pint.interp_linear(nb1,ye1,T1));
    cout << (P_eg-(P.interp_linear(nb1,ye1,T1)-
		   Pint.interp_linear(nb1,ye1,T1)))/
      fabs(P.interp_linear(nb1,ye1,T1)-Pint.interp_linear(nb1,ye1,T1)) << " ";
    sum+=fabs(S_eg-(S.interp_linear(nb1,ye1,T1)-
		    Sint.interp_linear(nb1,ye1,T1)))/
      fabs(S.interp_linear(nb1,ye1,T1)-Sint.interp_linear(nb1,ye1,T1));
    cout << (S_eg-(S.interp_linear(nb1,ye1,T1)-
		   Sint.interp_linear(nb1,ye1,T1)))/
      fabs(S.interp_linear(nb1,ye1,T1)-Sint.interp_linear(nb1,ye1,T1)) << " ";
    sum+=fabs(F_eg-(F.interp_linear(nb1,ye1,T1)-
		    Fint.interp_linear(nb1,ye1,T1)))/
      fabs(F.interp_linear(nb1,ye1,T1)-Fint.interp_linear(nb1,ye1,T1));
    cout << (F_eg-(F.interp_linear(nb1,ye1,T1)-
		   Fint.interp_linear(nb1,ye1,T1)))/
      fabs(F.interp_linear(nb1,ye1,T1)-Fint.interp_linear(nb1,ye1,T1)) << endl;
    cout.unsetf(ios::showpos);
    
    count+=4;
  }
  
  if (verbose>0) {
    cout << "Check electrons photons result: " << sum/count << endl;
  }
  
  return sum/count;
}

void eos_sn_oo::load_auto(std::string model, std::string directory) {

  string fname, name;

  // Initialize to zero to prevent uninitialized variable errors
  size_t mode=0;

  if (model=="fsu17") {
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
  } else if (model=="fsu21") {
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
  } else if (model=="sht_nl3") {
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
  } else if (model=="stos") {
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
  } else if (model=="stos_hyp") {
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
  } else if (model=="dd2") {
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
  } else if (model=="fsg") {
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
  } else if (model=="hfsl_nl3") {
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
  } else if (model=="sfho") {
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
  } else if (model=="sfhx") {
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
  } else if (model=="tm1") {
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
  } else if (model=="tma") {
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
  } else if (model=="ls180") {
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
  } else if (model=="ls220") {
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
  } else if (model=="ls375") {
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
  } else if (model=="sro_NRAPR") {
    cloud_file cf;
    cf.verbose=2;
    std::string sha=((std::string)"f2f211bb93ecc5587c029a741f01bee0")+
      "3e2bc399e188d1d14354c35eed70ad14";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash("NRAPR_3335_rho391_temp163_ye66.h5.bz2",
		     ((string)"https://isospin.roam.utk.edu/")+
		     "public_data/eos_tables/sro/"+
		     "NRAPR_3335_rho391_temp163_ye66.h5.bz2",sha,directory);
    name="sro_NRAPR";
    mode=eos_sn_oo::sht_mode;
    fname=directory+"/NRAPR_3335_rho391_temp163_ye66.h5";
    std::string cmd=((string)"bunzip2 -k ")+directory+"/"+
      "NRAPR_3335_rho391_temp163_ye66.h5.bz2";
    cout << "Decompressing with: " << cmd << endl;
    int iret=system(cmd.c_str());
    if (iret!=0) {
      cout << "Decompression failed." << endl;
    }
  } else if (model=="acmp_apr_sna") {
    mode=eos_sn_oo::sht_mode;
    fname="APR_0000_rho393_temp133_ye66_gitM180edd5_20190225.h5";
  } else if (model=="acmp_apr_nse") {
    mode=eos_sn_oo::sht_mode;
    fname="APR_0000_rho393_temp133_ye66_gitM180edd5_20190225.h5";
  } else {
    O2SCL_ERR("Need EOS type.",exc_efailed);
  }
    
  load(fname,mode);
    
  return;

}

void eos_sn_oo::load(std::string fname, size_t mode) {
  
  wordexp_single_file(fname);

  if (verbose>0) {
    cout << "In eos_sn_oo::load(), loading EOS from file\n\t'" 
	 << fname << "'." << endl;
  }

  double dtemp;
  
  if (loaded) free();

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
  if (verbose>1) {
    cout << "n_nB, n_T, n_Ye, energy_shift: " 
	 << n_nB << " " << n_T << " " << n_Ye << " " 
	 << energy_shift << endl;
  }

  n_oth=8;
  if (mode==hfsl_mode) n_oth+=4;

  if (mode==sht_mode) {
    m_neut=939.0;
    m_prot=939.0;
  } else if (mode==stos_mode) {
    m_neut=938.0;
    m_prot=938.0;
  } else if (mode==hfsl_mode) {
    m_neut=939.565346;
    m_prot=938.272013;
  } else if (mode==ls_mode) {
    m_neut=939.0;
    m_prot=939.0;
  } else {
    O2SCL_ERR2("Invalid value of 'mode' in ",
	       "eos_sn_oo::load().",exc_einval);
  }

  alloc();

  std::vector<double> ye_grid, t_grid, grid;

  hf.getd_vec("logrho",rho);
  hf.getd_vec("logtemp",t_grid);
  hf.getd_vec("ye",ye_grid);

  nB_grid.clear();
  for(size_t i=0;i<rho.size();i++) {
    // Undo the log
    rho[i]=pow(10.0,rho[i]);
    // Convert from g/cm^3 to baryon density through the 
    // atomic mass unit
    double nb=rho[i]/o2scl_cgs::unified_atomic_mass/1.0e39;
    grid.push_back(nb);
    nB_grid.push_back(nb);
  }
  Ye_grid.clear();
  for(size_t i=0;i<ye_grid.size();i++) {
    grid.push_back(ye_grid[i]);
    Ye_grid.push_back(ye_grid[i]);
  }
  T_grid.clear();
  for(size_t i=0;i<t_grid.size();i++) {
    grid.push_back(pow(10.0,t_grid[i]));
    T_grid.push_back(pow(10.0,t_grid[i]));
  }

  for(size_t i=0;i<n_base+n_oth;i++) {
    arr[i]->set_grid_packed(grid);
  }
  
  size_t ndat=19;
  // Names of sections in the HDF5 file
  vector<string> names;
  vector<size_t> indices;

  // 0-4
  names.push_back("Abar");
  indices.push_back(11);
  names.push_back("Xa");
  indices.push_back(14);
  names.push_back("Xh");
  indices.push_back(15);
  names.push_back("Xn");
  indices.push_back(12);
  names.push_back("Xp");
  indices.push_back(13);

  // 5-9
  names.push_back("Zbar");
  indices.push_back(10);
  names.push_back("cs2");
  indices.push_back(16);
  names.push_back("dedt");
  indices.push_back(17);
  names.push_back("dpderho");
  indices.push_back(18);
  names.push_back("dpdrhoe");
  indices.push_back(19);

  // 10-14
  names.push_back("entropy");
  indices.push_back(6);
  names.push_back("gamma");
  indices.push_back(20);
  names.push_back("logenergy");
  indices.push_back(2);
  names.push_back("logpress");
  indices.push_back(4);
  names.push_back("mu_e");
  indices.push_back(21);

  // 15-18
  names.push_back("mu_n");
  indices.push_back(8);
  names.push_back("mu_p");
  indices.push_back(9);
  names.push_back("muhat");
  indices.push_back(22);
  names.push_back("munu");
  indices.push_back(23);

  oth_names.push_back("cs2");
  oth_names.push_back("dedt");
  oth_names.push_back("dpderho");
  oth_names.push_back("dpdrhoe");
  oth_names.push_back("gamma");
  oth_names.push_back("mu_e");
  oth_names.push_back("muhat");
  oth_names.push_back("munu");

  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("MeV");
  oth_units.push_back("MeV");
  oth_units.push_back("MeV");

  if (mode==hfsl_mode) {
    ndat+=4;
    names.push_back("X3he");
    indices.push_back(24);
    names.push_back("X4li");
    indices.push_back(25);
    names.push_back("Xt");
    indices.push_back(26);
    names.push_back("Xd");
    indices.push_back(27);

    oth_names.push_back("X3he");
    oth_names.push_back("X4li");
    oth_names.push_back("Xt");
    oth_names.push_back("Xd");

    oth_units.push_back("");
    oth_units.push_back("");
    oth_units.push_back("");
    oth_units.push_back("");
  }
		  
  for(size_t i=0;i<ndat;i++) {
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
	    arr[indices[i]]->set(j,k,m,dat.get(k,m,j));

	    //if (mode==ls_mode) {
	    //} else {
	    //arr[indices[i]]->set
	    //(j,k,m,dat.get(k,m,j)+938.0-cu.convert
	    //("kg","MeV",o2scl_mks::mass_neutron));
	    //}
	  } else if (i==16) {
	    // Proton chemical potential
	    arr[indices[i]]->set(j,k,m,dat.get(k,m,j));

	    //if (mode==ls_mode) {
	    //arr[indices[i]]->
	    //set(j,k,m,dat.get(k,m,j)+
	    //cu.convert("kg","MeV",o2scl_mks::mass_neutron)-
	    //cu.convert("kg","MeV",o2scl_mks::mass_proton));
	    //} else {
	    //arr[indices[i]]->set
	    //(j,k,m,dat.get(k,m,j)+938.0-
	    //cu.convert("kg","MeV",o2scl_mks::mass_proton));
	    //}
	  } else {
	    arr[indices[i]]->set(j,k,m,dat.get(k,m,j));
	  }
	  // Set the free energy per baryon from the energy per 
	  // baryon and entropy per baryon
	  if (i==ndat-1) {
	    double T=F.get_grid(2,m);
	    F.set(j,k,m,E.get(j,k,m)-T*S.get(j,k,m));
	  }
	}
      }
    }
  }

  hf.close();

  // Loaded must be set to true before calling set_interp()
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
    std::cout << "Done in eos_sn_oo::load()." << std::endl;
  }

  return;
}

void eos_sn_stos::load_auto(std::string model, std::string directory) {
  size_t mode=orig_mode;
  std::string fname;
  if (model=="stos") {
    cloud_file cf;
    cf.verbose=2;
    //std::string sha=((std::string)"a7c34b355df12adc7aa8ca369e1")+
    //"6d40cbd7e04a0faf6d80d58e788d1610e45bf";
    std::string sha=((std::string)"fe0d7dc0940c9a2d5b7cc200f5d6")+
      "2c49a8e37f42cb51101e7fb7adf466ab2c9a";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("eos1.tab",
       ((string)"https://isospin.roam.utk.edu/")+
       "public_data/eos_tables/stos/"+
       "eos1.tab",sha,directory);
    mode=orig_mode;
    fname=directory+"/eos1.tab";
  } else if (model=="fyss") {
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
    mode=fyss_mode;
    fname=directory+"/FYSS_ver_1_27.tab";
  } else {
    O2SCL_ERR("Need EOS type.",exc_efailed);
  }
  
  load(fname,mode);
  return;
}

void eos_sn_stos::load(std::string fname, size_t mode) {

  wordexp_single_file(fname);

  if (verbose>0) {
    cout << "In eos_sn_stos::load(), loading EOS from file\n\t'" 
	 << fname << "'." << endl;
  }

  size_t ndat=17;
  if (mode==quark_mode) ndat++;
  if (mode==fyss_mode) ndat=25;

  if (loaded) free();

  std::ifstream fin;
  std::string tstr;

  double dtemp;

  fin.open(fname.c_str());

  n_nB=104;
  n_T=31;
  n_Ye=71;
  n_oth=5;
  oth_names.clear();
  oth_names.push_back("oth1");
  oth_names.push_back("oth2");
  oth_names.push_back("oth3");
  oth_names.push_back("oth4");
  oth_names.push_back("oth5");
  oth_units.clear();
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("");
  if (mode==quark_mode) {
    n_oth++;
    oth_names.push_back("oth6");
    oth_units.push_back("");
  }
  if (mode==fyss_mode) {
    n_nB=110;
    n_Ye=65;
    n_T=91;
  }

  alloc();

  std::vector<double> grid;

  nB_grid.clear();
  if (mode==fyss_mode) {
    for(size_t k=0;k<n_nB;k++) {
      double gcm3=pow(10.0,5.1+((double)k)*0.1);
      // AWS 7/28/2020: Computed with "acol -set precision 10
      // -get-conv g/cm^3 MeV/fm^3"
      double MeVfm3=gcm3*5.6095886038e-13;
      grid.push_back(MeVfm3/939.0);
      nB_grid.push_back(MeVfm3/939.0);
      //cout << k << " " << n_nB << " " << gcm3 << endl;
    }
  } else {
    for(size_t k=0;k<n_nB;k++) {
      // We set the density grid later, fill with zeros for now
      grid.push_back(0.0);
    }
  }

  Ye_grid.clear();
  if (mode==fyss_mode) {
    for(size_t j=0;j<n_Ye;j++) {
      grid.push_back(0.01*((double)j)+0.01);
      //cout << j << " " << n_Ye << " " << 0.01*((double)j)+0.01 << endl;
      Ye_grid.push_back(0.01*((double)j)+0.01);
    }
  } else {
    for(size_t j=0;j<n_Ye;j++) {
      double ye_temp=pow(10.0,((double)j)*0.025-2.0);
      grid.push_back(ye_temp);
      Ye_grid.push_back(ye_temp);
    }
  }

  T_grid.clear();
  if (mode==fyss_mode) {
    for(size_t i=0;i<n_T;i++) {
      double log10T=-1.0+((double)i)*0.04;
      grid.push_back(pow(10.0,log10T));
      //cout << i << " " << n_T << " " << log10T << " "
      //<< pow(10.0,log10T) << endl;
      T_grid.push_back(pow(10.0,log10T));
    }
  } else {
    // Temperature grid from Matthias Hempel
    double hempel_temp[31]={0.1,0.12,0.15,0.2,0.25,0.32,0.4,0.5,0.63,0.8,
			    1.0,1.2,1.5,2.0,2.5,3.2,4.0,5.0,6.3,8.0,10.0,
			    12.0,15.0,20.0,25.0,32.0,40.0,50.0,63.0,80.0,
			    100.0};
    for(size_t i=0;i<n_T;i++) {
      grid.push_back(hempel_temp[i]);
    }
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

    double tptr;
    if (mode==fyss_mode) {
      fin >> tstr;
      //cout << "1. " << tstr << endl;
      fin >> tstr;
      //cout << "2. " << tstr << endl;
      fin >> tstr;
      //cout << "3. " << tstr << endl;
      fin >> tstr;
      //cout << "4. " << tstr << endl;
      fin >> tptr;
    } else {
      if (i!=0) {
	getline(fin,tstr);
      }
      getline(fin,tstr);
      fin >> tstr;
      fin >> tptr;
    }

    if (verbose>0) {
      std::cout << "Reading section for temperature=" << tptr << "." 
		<< std::endl;
    }

    getline(fin,tstr);
    
    for(size_t j=0;j<n_Ye;j++) {
      //cout << "Herej: " << j << endl;
      
      for(size_t k=0;k<n_nB;k++) {
	//cout << "Herek: " << k << endl;
	if (mode==fyss_mode) {
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  Fint.set(k,j,i,dtemp);
	  fin >> dtemp;
	  Eint.set(k,j,i,dtemp);
	  fin >> dtemp;
	  Sint.set(k,j,i,dtemp);
	  fin >> dtemp;
	  A.set(k,j,i,dtemp);
	  fin >> dtemp;
	  Z.set(k,j,i,dtemp);
	  /*
	    cout << nB_grid[k] << " ";
	    cout << Ye_grid[j] << " ";
	    cout << T_grid[i] << " ";
	    cout << A.get(k,j,i) << " ";
	    cout << Z.get(k,j,i) << endl;
	    char ch;
	    cin >> ch;
	  */
	  fin >> dtemp;
	  fin >> dtemp;
	  Xn.set(k,j,i,dtemp);
	  fin >> dtemp;
	  Xp.set(k,j,i,dtemp);
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  Pint.set(k,j,i,dtemp);
	  fin >> dtemp;
	  mun.set(k,j,i,dtemp);
	  fin >> dtemp;
	  mup.set(k,j,i,dtemp);
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  fin >> dtemp;
	  Xalpha.set(k,j,i,dtemp);
	  fin >> dtemp;
	  Xnuclei.set(k,j,i,dtemp);
	  fin >> dtemp;
	  fin >> dtemp;
	} else {
	  for(size_t ell=0;ell<ndat;ell++) {
	    fin >> dtemp;
	    t.set(ell,k,dtemp);
	  }
	}
      }

      if (i==0 && j==0 && mode!=fyss_mode) {
	
	// Set the density grid from the first section
	for(size_t k=0;k<n_nB;k++) {
	  grid[k]=t[1][k];
	  nB_grid[k]=grid[k];
	}
	for(size_t k=0;k<n_base+n_oth;k++) {
	  arr[k]->set_grid_packed(grid);
	}
	    
      } 
      
      if (mode!=fyss_mode) {
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
			 "eos_sn_stos::load().",exc_einval);
	    }
	  }
	}
      }

      if (mode!=fyss_mode) {
	// Read the blank line between sections
	getline(fin,tstr);
      }
    }
  }

  fin.close();

  // Loaded must be set to true before calling set_interp()
  loaded=true;
  with_leptons_loaded=false;
  baryons_only_loaded=true;

  // It is important that 'loaded' is set to true before the call to
  // set_interp_type().
  set_interp_type(itp_linear);

  // Double check the grid 
  if (mode!=fyss_mode && check_grid) {
    if (verbose>0) {
      std::cout << "Checking grid in eos_sn_stos::load(). " << std::endl;
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
      O2SCL_ERR("Check grid failed in eos_sn_stos::load().",exc_efailed);
    }

  }

  if (verbose>0) {
    std::cout << "Done in eos_sn_stos::load()." << std::endl;
  }

  return;
}

void eos_sn_sht::load(std::string fname, size_t mode) {

  wordexp_single_file(fname);

  if (verbose>0) {
    cout << "In eos_sn_sht::load(), loading EOS from file\n\t'" 
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
    // These are valid for both FSU2.1 and NL3
    n_nB=328;
    n_T=109;
    n_Ye=53;
    n_oth=5;
  }
  size_t ndat=16;

  tensor_grid3<> *Eptr, *Fptr, *Sptr, *Pptr;
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
  std::vector<double> grid;

  alloc();
  
  for(size_t i=0;i<n_nB;i++) {
    double nb_temp=pow(10.0,-8.0+((double)i)*0.025);
    grid.push_back(nb_temp);
  }
  double ye_temp=0.0;
  grid.push_back(ye_temp);
  for(size_t i=1;i<n_Ye;i++) {
    ye_temp=0.05+((double)(i-1)*0.01);
    grid.push_back(ye_temp);
  }
  double temp_temp=0.0;
  grid.push_back(temp_temp);
  for(size_t i=1;i<n_T;i++) {
    temp_temp=pow(10.0,-0.8+((double)(i-1))*0.025);
    grid.push_back(temp_temp);
  }
  for(size_t k=0;k<n_base+n_oth;k++) {
    arr[k]->set_grid_packed(grid);
  }
  
  string stmp;

  for(size_t j=0;j<n_T;j++) {
    if (verbose>0) {
      cout << "Loading section for temperature=" << A.get_grid(2,j) 
	   << "." << endl;
    }
    for(size_t k=0;k<n_Ye;k++) {
      for(size_t i=0;i<n_nB;i++) {
	for(size_t ell=0;ell<ndat;ell++) {
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
	    // Free energy per baryon
	    Fptr->set(i,k,j,dtemp);
	  } else if (ell==4) {
	    Pptr->set(i,k,j,dtemp);
	  } else if (ell==5) {
	    // Entropy per baryon
	    Sptr->set(i,k,j,dtemp);
	    // Also compute energy per baryon
	    double T=Eptr->get_grid(2,j);
	    Eptr->set(i,k,j,Fptr->get(i,k,j)+T*Sptr->get(i,k,j));
	  } else if (ell==6) {
	    // Neutron chemical potential
	    mun.set(i,k,j,dtemp);
	  } else if (ell==7) {
	    // Proton chemical potential
	    mup.set(i,k,j,dtemp);
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
		       "eos_sn_sht::load().",exc_einval);
	  }
	}
      }
    }
    
  }
      
  fin.close();

  oth_names.clear();
  oth_names.push_back("T");
  oth_names.push_back("Yp");
  oth_names.push_back("nB");
  oth_names.push_back("mue");
  oth_names.push_back("M_star");

  oth_units.clear();
  oth_units.push_back("MeV");
  oth_units.push_back("");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("MeV");
  oth_units.push_back("MeV");
  
  if (n_oth!=oth_names.size()) {
    O2SCL_ERR2("Number of names does not match number of data sets ",
	       "in eos_sn_sht::load().",exc_efailed);
  }

  if (check_grid) {
    // Double check the grid 
    if (verbose>0) {
      cout << "Checking grid in eos_sn_sht::load()." << endl;
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
      O2SCL_ERR("Grid check in eos_sn_sht::load() failed.",exc_efailed);
    }

  }

  // Here, we make sure not to reset the values of the flags if an EOS
  // has already been loaded. This is important so we can load both
  // the EOS table with lepton contributions and the table without.
  if (mode==mode_17b || mode==mode_21b || mode==mode_NL3b) {
    if (loaded==false) with_leptons_loaded=false;
    baryons_only_loaded=true;
  } else {
    with_leptons_loaded=true;
    if (loaded==false) baryons_only_loaded=false;
  }
  loaded=true;

  if (verbose>0) {
    std::cout << "Done in eos_sn_sht::load()." << std::endl;
  }

  return;
}

void eos_sn_hfsl::load(std::string fname, size_t mode) {

  wordexp_single_file(fname);

  if (verbose>0) {
    cout << "In eos_sn_hfsl::load(), loading EOS from file\n\t'" 
	 << fname << "'." << endl;
  }

  if (loaded) free();

  std::ifstream fin;
  std::string tstr;

  fin.open(fname.c_str());

  n_nB=326;
  n_T=81;
  n_Ye=60;
  n_oth=7;
  size_t ndat=19;

  double dtemp;
  std::vector<double> grid;

  alloc();

  for(size_t i=0;i<n_nB;i++) {
    double nb_temp=pow(10.0,((double)i)*0.04-12);
    grid.push_back(nb_temp);
    nB_grid.push_back(nb_temp);
  }
  for(size_t i=0;i<n_Ye;i++) {
    double ye_temp=0.01*((double)(i+1));
    grid.push_back(ye_temp);
    Ye_grid.push_back(ye_temp);
  }
  double temp_temp;
  for(size_t i=0;i<n_T;i++) {
    temp_temp=pow(10.0,((double)i)*0.04-1.0);
    grid.push_back(temp_temp);
    T_grid.push_back(temp_temp);
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
    fin >> temp_temp;

    if (verbose>0) {
      cout << "Reading section for temperature=" 
	   << temp_temp << "." << endl;
    }

    getline(fin,tstr);
    for(size_t j=0;j<n_Ye;j++) {
      for(size_t k=0;k<n_nB;k++) {
        for(size_t ell=0;ell<ndat;ell++) {
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
	    // The free energy in the Shen '98 format is stored with respect
	    // to the proton mass, so we rescale it
	    double Ye_tmp=Fint.get_grid(1,j);
	    double dtemp2=dtemp+938.0-Ye_tmp*m_prot-(1.0-Ye_tmp)*m_neut;
	    Fint.set(k,j,i,dtemp2);
	  } else if (ell==5) {
	    // The internal energy in the table is stored with respect
	    // to the atomic mass, so we rescale it
	    double Ye_tmp=Eint.get_grid(1,j);
	    double dtemp2=dtemp+m_amu-Ye_tmp*m_prot-(1.0-Ye_tmp)*m_neut;
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
            O2SCL_ERR("Index problem in eos_sn_hfsl::load().",exc_einval);
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

  // It is important that 'loaded' is set to true before the call to
  // set_interp_type().
  set_interp_type(itp_linear);

  oth_names.push_back("log_rho");
  oth_names.push_back("nB");
  oth_names.push_back("log_Y");
  oth_names.push_back("Yp");
  oth_names.push_back("M_star");
  oth_names.push_back("A_light");
  oth_names.push_back("Z_light");
  
  oth_units.push_back("");
  oth_units.push_back("1/fm^3");
  oth_units.push_back("");
  oth_units.push_back("");
  oth_units.push_back("MeV");
  oth_units.push_back("");
  oth_units.push_back("");
  
  if (n_oth!=oth_names.size()) {
    O2SCL_ERR2("Number of names does not match number of data sets ",
	       "in eos_sn_hfsl::load().",exc_efailed);
  }

  if (check_grid) {

    // Double check the grid 
    if (verbose>0) {
      cout << "Checking grid. " << endl;
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
      O2SCL_ERR("Function eos_sn_hfsl::load() failed.",exc_efailed);
    }

  }

  if (verbose>0) {
    std::cout << "Done in eos_sn_hfsl::load()." << std::endl;
  }

  return;
}
