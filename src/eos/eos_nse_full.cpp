/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014, Andrew W. Steiner
  
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
#include <o2scl/eos_nse_full.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

dense_matter::dense_matter() {
  n.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  p.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  n.inc_rest_mass=false;
  p.inc_rest_mass=false;
  n.non_interacting=false;
  p.non_interacting=false;
  
  e.init(o2scl_settings.get_convert_units().convert
         ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  e.inc_rest_mass=true;
  e.non_interacting=true;

  mu.init(o2scl_settings.get_convert_units().convert
	  ("kg","1/fm",o2scl_mks::mass_muon),2.0);
  mu.inc_rest_mass=true;
  mu.non_interacting=true;

  photon.init(0.0,2.0);
}

double dense_matter::average_a() {
  double ntot=0.0;
  for (size_t i=0;i<dist.size();i++) {
    ntot+=dist[i].n;
  }
  if (ntot==0.0) return 0.0;
  return cbrt(3.0/4.0/o2scl_const::pi/ntot);
}

double dense_matter::average_A() {
    
  double ntot=0.0;
  double Antot=0.0;
    
  for (size_t i=0;i<dist.size();i++) {
    ntot+=dist[i].n;
    Antot+=dist[i].A*dist[i].n;
  }
  if (ntot==0.0) return 0.0;
    
  return Antot/ntot;
}
  
double dense_matter::average_Z() {
    
  double ntot=0.0;
  double Zntot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    ntot+=dist[i].n;
    Zntot+=dist[i].Z*dist[i].n;
  }
  if (ntot==0.0) return 0.0;
    
  return Zntot/ntot;
}

double dense_matter::average_N() {
    
  double ntot=0.0;
  double Nntot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    ntot+=dist[i].n;
    Nntot+=dist[i].N*dist[i].n;
  }
  if (ntot==0.0) return 0.0;
    
  return Nntot/ntot;
}

double dense_matter::nB_nuclei() {
    
  double nBtot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    nBtot+=dist[i].n*(dist[i].Z+dist[i].N);
  }
  
  return nBtot;
}

double dense_matter::impurity() {
    
  // First compute average proton number
  double ave_Z=average_Z();

  double sum=0.0, ntot=0.0;
    
  for (size_t i=0;i<dist.size();i++) {
    sum+=dist[i].n*pow(dist[i].Z-ave_Z,2.0);
    ntot+=dist[i].n;
  }
  
  if (ntot==0.0) return 0.0;

  return sum/ntot;
}

nucmass_densmat::nucmass_densmat() {
  massp=&ame;
  o2scl_hdf::ame_load(ame,"12");
}

void nucmass_densmat::set_mass(nucmass &nm) {
  massp=&nm;
  return;
}

void nucmass_densmat::test_derivatives(double eps, double &t1, double &t2, 
				       double &t3, double &t4) {
      
  double Z=26.0;
  double N=30.0;

  // None of these can be zero because we divide by them in 
  // the tests below
  double npout=0.005;
  double nnout=0.02;
  double nneg=0.01;
  double T=0.01;

  double E2, E1, dEdnp, dEdnn, dEdnneg, dEdT;
  double temp1, temp2, temp3, temp4;
  binding_energy_densmat_derivs(Z,N,npout,nnout,nneg,T,E1,
				dEdnp,dEdnn,dEdnneg,dEdT);
      
  binding_energy_densmat_derivs(Z,N,npout*(1.0+eps),nnout,nneg,T,E2,
				temp1,temp2,temp3,temp4);
  if (fabs(dEdnp)<1.0e-20) {
    t1=fabs(dEdnp-(E2-E1)/(npout*eps));
  } else {
    t1=fabs(dEdnp-(E2-E1)/(npout*eps))/fabs(dEdnp);
  }
      
  binding_energy_densmat_derivs(Z,N,npout,nnout*(1.0+eps),nneg,T,E2,
				temp1,temp2,temp3,temp4);
  if (fabs(dEdnp)<1.0e-20) {
    t2=fabs(dEdnn-(E2-E1)/(nnout*eps));
  } else {
    t2=fabs(dEdnn-(E2-E1)/(nnout*eps))/fabs(dEdnn);
  }
      
  binding_energy_densmat_derivs(Z,N,npout,nnout,nneg*(1.0+eps),T,E2,
				temp1,temp2,temp3,temp4);
  if (fabs(dEdnp)<1.0e-20) {
    t3=fabs(dEdnneg-(E2-E1)/(nneg*eps));
  } else {
    t3=fabs(dEdnneg-(E2-E1)/(nneg*eps))/fabs(dEdnneg);
  }
      
  binding_energy_densmat_derivs(Z,N,npout,nnout,nneg,T*(1.0+eps),E2,
				temp1,temp2,temp3,temp4);
  if (fabs(dEdnp)<1.0e-20) {
    t4=fabs(dEdT-(E2-E1)/(T*eps));
  } else {
    t4=fabs(dEdT-(E2-E1)/(T*eps))/fabs(dEdT);
  }
      
  return;
}

void nucmass_densmat::binding_energy_densmat_derivs
(double Z, double N, double npout, double nnout, 
 double nneg, double T, double &E, double &dEdnp, double &dEdnn,
 double &dEdnneg, double &dEdT) {

  // Half saturation density
  double n0o2=0.08;

  if (nneg<npout) {
    O2SCL_ERR2("Not enough negative charges in nucmass_densmat::",
	       "binding_energy_densmat_derivs().",exc_einval);
  }
  if (npout>n0o2) {
    O2SCL_ERR2("Too many protons in nucmass_densmat::",
	       "binding_energy_densmat_derivs().",exc_einval);
  }

  // Radii
  double R_p_3=3.0*Z/4.0/o2scl_const::pi/(n0o2-npout);
  double R_n_3=3.0*N/4.0/o2scl_const::pi/(n0o2-nnout);
  double R_p=cbrt(R_p_3), R_n=cbrt(R_n_3);
  double R_WS_3=R_p_3*(n0o2-npout)/(nneg-npout);
  double R_WS=cbrt(R_WS_3);

  if (R_p>R_WS) {
    O2SCL_ERR2("Proton radius larger than cell in nucmass_densmat::",
	       "binding_energy_densmat_derivs().",exc_einval);
  }
  if (R_n>R_WS) {
    O2SCL_ERR2("Neutron radius larger than cell in nucmass_densmat::",
	       "binding_energy_densmat_derivs().",exc_einval);
  }

  // Volume fractions
  double chi_p=R_p_3/R_WS_3;
  double chi_n=R_n_3/R_WS_3;
      
  // Add the finite-size part of the Coulomb energy
  double fdu=0.2*chi_p-0.6*cbrt(chi_p);
  double coul=(Z+N)*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure*R_p*R_p*
    pow(fabs(n0o2-npout),2.0)/0.16*fdu;

  // Total binding energy
  E=massp->mass_excess_d(Z,N)+coul+(Z+N)*massp->m_amu-Z*massp->m_elec-
    N*massp->m_neut-Z*massp->m_prot;
      
  // Derivatives
  double dfof=(0.2-0.2*pow(chi_p,-2.0/3.0))/fdu;
  double dchi_dnp=-(n0o2-nneg)/pow(n0o2-npout,2.0);
  double dchi_dnneg=1.0/(n0o2-npout);

  dEdnp=-4.0/3.0*coul/(n0o2-npout)+coul*dfof*dchi_dnp;
  dEdnneg=coul*dfof*dchi_dnneg;
  dEdT=0.0;
  dEdnn=0.0;

  return;
}

eos_nse_full::eos_nse_full() {
  invalid_config=-10;
  inc_lept_phot=true;

  // Ensure fermion_rel doesn't throw exceptions for convergence
  // errors for electrons
  relf.density_root->err_nonconv=false;
  relf.err_nonconv=false;

  // Load Skyrme EOS
  o2scl_hdf::skyrme_load(sk,"SLy4");
  
  ehtp=&sk;
  
  massp=&nuc_dens;
  
  nucdist_set(def_dist,nuc_dens.ame);
  ad=&def_dist;

  def_mmin.ntrial*=100;
  def_mmin.tol_rel/=1000.0;
  def_mmin.tol_abs/=1000.0;
}

double eos_nse_full::free_energy_nr(const ubvector &n_nuc, dense_matter &dm) {

  double nB_nuc=0.0, np_nuc=0.0;

  for(size_t i=0;i<n_nuc.size();i++) {
    dm.dist[i].n=n_nuc[i];
    nB_nuc+=n_nuc[i]*(dm.dist[i].Z+dm.dist[i].N);
    np_nuc+=n_nuc[i]*dm.dist[i].Z;
  }

  dm.p.n=dm.nB*dm.Ye-np_nuc;
  dm.n.n=dm.nB-nB_nuc-dm.p.n;

  // Return a large value if necessary
  //dm.e.n=dm.nB*dm.Ye;
  //if (dm.e.n<dm.p.n || dm.n.n<0.0 || dm.p.n<0.0) return 1.0e6;

  int ret=calc_density_noneq_nr(dm);
  if (ret==invalid_config) return 1.0e4;

  /*
    cout << "fen: ";
    cout.precision(4);
    vector_out(cout,n_nuc);
    cout.precision(6);
    cout << " " << dm.th.ed-dm.T*dm.th.en << endl;
  */

  return dm.th.ed-dm.T*dm.th.en;
}

int eos_nse_full::calc_density_fixcomp_nr(dense_matter &dm, int verbose) {

  ubvector n_nuc(dm.dist.size()), n_nuc2(dm.dist.size());
  for(size_t i=0;i<n_nuc.size();i++) {
    n_nuc[i]=dm.dist[i].n;
    n_nuc2[i]=n_nuc[i]*1.01;
  }

  multi_funct11 mf=std::bind
    (std::mem_fn<double(const ubvector &,dense_matter &)>
     (&eos_nse_full::free_energy_nr),
     this,std::placeholders::_2,std::ref(dm));

  double fr_min=0.0;

  //fr_min=free_energy_nr(n_nuc,dm);
  //cout << fr_min << endl;
  //vector_out(cout,n_nuc,true);

  def_mmin.mmin_twovec(n_nuc.size(),n_nuc,n_nuc2,fr_min,mf);

  // Perform a final function evaluation (important to set the final
  // nuclear densities in the 'dm' parameter)
  fr_min=free_energy_nr(n_nuc,dm);

  return 0;
}

int eos_nse_full::calc_density_noneq_nr(dense_matter &dm, int verbose) {
  
  // -----------------------------------------------------------
  // Sanity checks

  if (!o2scl::is_finite(dm.n.n) || 
      !o2scl::is_finite(dm.p.n)) {
    O2SCL_ERR2("Neutron or proton density not finite in ",
	       "eos_nse_full::calc_density_noneq_nr().",exc_esanity);
  }
  if (dm.n.m<0.0 || dm.p.m<0.0) {
    O2SCL_ERR2("Mass negative in ",
	       "eos_nse_full::calc_density_noneq_nr().",exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].inc_rest_mass==true) {
      O2SCL_ERR("Wrong inc_rest_mass for nuclei in noneq_nr().",
		exc_esanity);
    }
  }
  if (dm.n.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for neutrons in noneq_nr().",
	      exc_esanity);
  }
  if (dm.p.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for protons in noneq_nr().",
	      exc_esanity);
  }
  if (dm.e.inc_rest_mass==false) {
    O2SCL_ERR("Wrong inc_rest_mass for electrons in noneq_nr().",
	      exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].non_interacting==false) {
      O2SCL_ERR("Wrong non_interacting for nuclei in noneq_nr().",
		exc_esanity);
    }
  }
  if (dm.n.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for neutrons in noneq_nr().",
	      exc_esanity);
  }
  if (dm.p.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for protons in noneq_nr().",
	      exc_esanity);
  }
  if (dm.e.non_interacting==false) {
    O2SCL_ERR("Wrong non_interacting for electrons in noneq_nr().",
	      exc_esanity);
  }

  if (dm.n.n<0.0 || dm.p.n<0.0) {
    if (verbose>0) {
      cout << "Neutron (" << dm.n.n << ") or proton ("
	   << dm.p.n << ") density negative." << endl;
    }
    return invalid_config;
  }
  
  // -----------------------------------------------------------
  // Compute properties of homogeneous matter

  int ret=ehtp->calc_temp_e(dm.n,dm.p,dm.T,dm.drip_th);
  if (ret!=0) return ret;

  // Add contribution from dripped nucleons
  dm.th=dm.drip_th;

  // Compute Ye, nB and electron density
  double nB, Ye;
  nB=dm.n.n+dm.p.n;
  Ye=dm.p.n;
  for(size_t i=0;i<dm.dist.size();i++) {
    nB+=dm.dist[i].n*((double)dm.dist[i].A);
    Ye+=dm.dist[i].n*((double)dm.dist[i].Z);
  }
  Ye/=nB;
  dm.e.n=Ye*nB;

  // -----------------------------------------------------------
  // Verbose output

  if (verbose>0) {
    cout.setf(ios::left);
    cout.setf(ios::showpos);
    cout << "--------------------------------------"
	 << "--------------------------------------" << endl;
    cout << "noneq_nr():" << endl;
    cout << endl;
    cout << "nB, Ye, T (MeV): " << nB << " " << Ye << " "
	 << dm.T*hc_mev_fm << endl;
    cout << endl;
    cout << "                                n            mu          "
	 << "  ed            en" << endl;
    cout.width(20);
    cout << "Neutrons: " << dm.n.n << " " << dm.n.mu << endl;
    cout.width(20);
    cout << "Protons: " << dm.p.n << " " << dm.p.mu << endl;
    cout.width(48);
    cout << "Dripped nucleons: " 
	 << dm.drip_th.ed << " " << dm.drip_th.en << endl;
  }

  // -----------------------------------------------------------
  // Leptons and photons

  if (inc_lept_phot) {

    // Add electrons
    dm.e.n=Ye*nB;
    if (dm.e.n<dm.p.n) {
      if (verbose>0) {
	cout << "Electron density too small to match proton density."
	     << endl;
	cout << "np: " << dm.p.n << " ne: " << dm.e.n << endl;
      }
      return invalid_config;
    }
    ret=relf.calc_density(dm.e,dm.T);
    if (ret!=0) return ret;
    dm.th.ed+=dm.e.ed;
    dm.th.en+=dm.e.en;

    if (verbose>0) {
      cout.width(20);
      cout << "Electrons: " 
	   << dm.e.n << " " << dm.e.mu << " " << dm.e.ed << " "
	   << dm.e.en << endl;
    }

    if (false) {
      // Add muons
      dm.mu.mu=dm.e.mu;
      relf.calc_mu(dm.mu,dm.T);
      dm.th.ed+=dm.mu.ed;
      dm.th.en+=dm.mu.en;
      
      if (verbose>0) {
	cout.width(20);
	cout << "Muons: " 
	     << dm.mu.n << " " << dm.mu.mu << " " << dm.mu.ed << " "
	     << dm.mu.en << endl;
      }
    }

    // Add photons
    dm.photon.massless_calc(dm.T);
    dm.th.ed+=dm.photon.ed;
    dm.th.en+=dm.photon.en;

    if (verbose>0) {
      cout.width(20);
      cout << "Photons: " 
	   << dm.photon.n << " " << 0.0 << " " << dm.photon.ed << " "
	   << dm.photon.en << endl;
    }
  }

  // Vectors for derivatives of nuclear binding energy
  ubvector vec_dEdne(dm.dist.size()), vec_dEdnp(dm.dist.size());

  // Compute the properties of the nuclei
  for(size_t i=0;i<dm.dist.size();i++) {
    
    if (dm.dist[i].n>0.0) {

      // Check that the proton density isn't too large in
      // the presence of nuclei
      if (dm.p.n>0.08) {
	if (verbose>0) {
	  cout << "External proton density (" << dm.p.n 
	       << ") too large." << endl;
	}
	return invalid_config;
      }
      
      // Check that the electron density isn't too large
      // in the presence of nuclei
      double fac=(0.08-dm.p.n)/(dm.e.n-dm.p.n);
      if (1.0>fac) {
	if (verbose>0) {
	  cout << "Proton radius negative or larger than cell size." << endl;
	  cout << "fac: " << fac << endl;
	}
	return invalid_config;
      }
      
      // Create a reference for this nucleus
      nucleus &nuc=dm.dist[i];
      
      // Compute nuclear binding energy and total mass
      double dEdnp, dEdnn, dEdne, dEdT;
      massp->binding_energy_densmat_derivs
	(nuc.Z,nuc.N,dm.p.n,dm.n.n,dm.e.n,dm.T,nuc.be,dEdnp,dEdnn,dEdne,dEdT);
      nuc.be/=hc_mev_fm;
      nuc.m=nuc.Z*dm.p.m+nuc.N*dm.n.m+nuc.be;
      vec_dEdnp[i]=dEdnp;
      vec_dEdne[i]=dEdne;

      // Translational energy
      cla.calc_density(nuc,dm.T);

      // Update thermo object with information from nucleus
      dm.th.ed+=nuc.be*nuc.n+nuc.ed;
      dm.th.en+=nuc.en;

      if (verbose>0) {
	string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	  itos(((int)(nuc.N+1.0e-8)))+"): ";
	cout.width(20);
	cout << s << nuc.n << " " << nuc.mu << " " 
	     << nuc.ed+nuc.be*nuc.n << " " << nuc.en << endl;
      }
      
    } else if (dm.dist[i].n<0.0) {
      if (verbose>0) {
	cout << "Density of nucleus: " << i << " negative." << endl;
      }
      return invalid_config;
    }

  }

  // -----------------------------------------------------------
  // Compute etas

  // Ensure the eta vector has the correct size
  if (dm.eta_nuc.size()!=dm.dist.size()) {
    dm.eta_nuc.resize(dm.dist.size());
  }

  dm.eta_n=dm.n.mu;
  dm.eta_p=dm.p.mu+dm.e.mu;

  for(size_t i=0;i<dm.dist.size();i++) {

    double dmudm_i=-1.5*dm.T/dm.dist[i].m;
    double dfdm_i=dm.dist[i].n*dmudm_i;
    dm.eta_nuc[i]=dm.dist[i].be+dm.dist[i].mu+dm.dist[i].Z*dm.e.mu;
    dm.eta_p+=(dm.dist[i].n+dfdm_i)*(vec_dEdnp[i]+vec_dEdne[i])/hc_mev_fm;
    
    for(size_t j=0;j<dm.dist.size();j++) {
      
      if (dm.dist[i].n>0.0) {
	
	double dmudm=-1.5*dm.T/dm.dist[j].m;
	double dfdm=dm.dist[j].n*dmudm;
	
	dm.eta_nuc[i]+=dm.dist[i].Z*(dm.dist[j].n+dfdm)*
	  vec_dEdne[j]/hc_mev_fm;
      }
    }
  }
      

  // -----------------------------------------------------------
  // Computation of pressure

  dm.th.pr=-dm.th.ed+dm.n.n*dm.eta_n+dm.p.n*dm.eta_p+
    dm.T*dm.th.en;
  for(size_t i=0;i<dm.dist.size();i++) {
    dm.th.pr+=dm.dist[i].n*dm.eta_nuc[i];
  }

  // -----------------------------------------------------------
  // Verbose output

  if (verbose>0) {
    cout.width(48);
    cout << "Total: " << dm.th.ed << " " << dm.th.en << endl;
    cout << endl;
    cout << "Free energy: " << dm.th.ed-dm.T*dm.th.en << endl;
    cout << endl;
    cout << "Contributions to pressure:" << endl;
    cout << "                                n           eta          "
	 << "  pr" << endl;
    cout.width(48);
    cout << "- Energy: " << -dm.th.ed << endl;
    cout.width(48);
    cout << "T * Entropy: " << dm.T*dm.th.en << endl;
    cout.width(20);
    cout << "Neutrons: " << dm.n.n << " " << dm.eta_n << " "
	 << dm.n.n*dm.eta_n << endl;
    cout.width(20);
    cout << "Protons: " << dm.p.n << " " << dm.eta_p << " "
	 << dm.p.n*dm.eta_p << endl;
    for(size_t i=0;i<dm.dist.size();i++) {
      nucleus &nuc=dm.dist[i];
      string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	itos(((int)(nuc.N+1.0e-8)))+"): ";
      cout.width(20);
      cout << s << nuc.n << " " << dm.eta_nuc[i] << " "
	   << nuc.n*dm.eta_nuc[i] << endl;
    }
    cout.width(48);
    cout << "Total pressure: " << dm.th.pr << endl;
    cout << "--------------------------------------"
	 << "--------------------------------------" << endl;
    cout.unsetf(ios::left);
    cout.unsetf(ios::showpos);
  }
  
  return success;
}

int eos_nse_full::density_match_nr(dense_matter &dm) {
  
  double nn_fix=(1.0-dm.Ye)*dm.nB;
  double np_fix=dm.Ye*dm.nB;

  // Fix negative densities
  if (dm.n.n<0.0) dm.n.n=0.0;
  if (dm.p.n<0.0) dm.p.n=0.0;
  for(vector<nucleus>::iterator it=dm.dist.begin();it!=dm.dist.end();it++) {
    if (it->n<0.0) {
      dm.dist.erase(it);
      it=dm.dist.begin();
    }
  }

  // Initial evaluation
  int ret=calc_density_noneq_nr(dm);
  
  // If it's invalid, try to fix
  if (ret==invalid_config) {

    double shift=dm.nB*1.0e-9;

    while(ret==invalid_config && shift<dm.nB) {

      // Use current shift to create a new configuration
      double nB_tot=0.0;
      for(size_t i=0;i<dm.dist.size();i++) {
	if (dm.dist[i].n>shift) dm.dist[i].n-=shift;
	nB_tot+=(dm.dist[i].Z+dm.dist[i].N)*dm.dist[i].n;
      }
      double nB_corr=dm.nB-nB_tot;
      dm.p.n+=dm.Ye*nB_corr;
      dm.p.n+=(1.0-dm.Ye)*nB_corr;

      // Evaluate new 
      ret=calc_density_noneq_nr(dm);
      shift*=10.0;
    }
  }

  // If we couldn't fix, throw
  if (ret==invalid_config) {
    O2SCL_ERR2("Could not find valid configuration in ",
	       "eos_nse_full::density_match_nr().",exc_efailed);
  } else if (ret!=0) {
    // Or if there is some other failure, return
    return ret;
  }

  // Update for density match
  if (fabs((1.0-dm.Ye)*dm.nB-nn_fix)/nn_fix>1.0e-6 ||
      fabs(dm.Ye*dm.nB-np_fix)/np_fix>1.0e-6) {
    
    // First make sure we don't over estimate the particle number
    double rat_n=dm.nB*(1.0-dm.Ye)/nn_fix;
    double rat_p=dm.nB*dm.Ye/np_fix;
    double ratio;
    if (rat_n>rat_p) ratio=rat_n;
    else ratio=rat_p;
    cout << "Adjusting by ratio: " << ratio << endl;
    dm.n.n/=ratio;
    dm.p.n/=ratio;
    for(size_t i=0;i<dm.dist.size();i++) {
      dm.dist[i].n/=ratio;
    }
    
    // Now add free neutrons and protons to match
    calc_density_noneq_nr(dm);
    cout << "Adjusting neutrons by: " << nn_fix-(1.0-dm.Ye)*dm.nB << endl;
    dm.n.n+=nn_fix-(1.0-dm.Ye)*dm.nB;
    cout << "Adjusting protons by: " << np_fix-dm.Ye*dm.nB << endl;
    dm.p.n+=np_fix-dm.Ye*dm.nB;
  }
  
  ret=calc_density_noneq_nr(dm);

  if (ret==invalid_config) {
    O2SCL_ERR2("Did not produce valid configuration in ",
	       "eos_nse_full::density_match_nr().",exc_efailed);
  }

  if (false) {
    cout << "Density match: " << endl;
    cout << (1.0-dm.Ye)*dm.nB << " " << nn_fix << endl;
    cout << dm.Ye*dm.nB << " " << np_fix << endl;
    cout << endl;
  }

  if (fabs((1.0-dm.Ye)*dm.nB-nn_fix)/nn_fix>1.0e-6 ||
      fabs(dm.Ye*dm.nB-np_fix)/np_fix>1.0e-6) {
    O2SCL_ERR2("Density match failed in ",
	       "eos_nse_full::density_match_nr().",exc_esanity);
  }

  return ret;
}

void eos_nse_full::output(dense_matter &dm, int verbose) {
  cout << "nB=" << dm.nB << " fm^{-3}, Ye=" << dm.Ye << ", T="
       << dm.T*hc_mev_fm << " MeV" << endl;
  double nBnuc=0.0;
  double Xa=0.0, Xnuclei=0.0;
  for(size_t i=0;i<dm.dist.size();i++) {
    nBnuc+=dm.dist[i].n*(dm.dist[i].N+dm.dist[i].Z);
    if (dm.dist[i].Z==2 && dm.dist[i].N==2) Xa=dm.dist[i].n*4.0;
    else Xnuclei+=dm.dist[i].n*(dm.dist[i].Z+dm.dist[i].N);
    if (verbose>=2) {
      cout << "Z,N,n: " << dm.dist[i].Z << " " << dm.dist[i].N << " "
	   << dm.dist[i].n << endl;
    }
  }
  Xa/=dm.nB;
  Xnuclei/=dm.nB;
  cout << "nn,np,nBnuc: " << dm.n.n << " " << dm.p.n << " " 
       << nBnuc << " fm^{-3}" << endl;
  cout << "N,<Z>,<N>,<Q>: " << dm.dist.size() << " " 
       << dm.average_Z() << " "
       << dm.average_N() << " " << dm.impurity() << endl;
  cout << "fr: " << dm.th.ed-dm.T*dm.th.en << " fm^{-4}" << endl;
  if (verbose>=1) {
    if (inc_lept_phot) {
      cout << "F: " << (dm.th.ed-dm.T*dm.th.en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "E: " << (dm.th.ed)/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "S: " << (dm.th.en)/dm.nB << endl;
      double ed=dm.th.ed-dm.e.ed-dm.photon.ed;
      double en=dm.th.en-dm.e.en-dm.photon.en;
      cout << "Fint: " << (ed-dm.T*en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "Eint: " << (ed)/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "Sint: " << (en)/dm.nB << " MeV" << endl;
    } else {
      cout << "Fint: " << (dm.th.ed-dm.T*dm.th.en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "Eint: " << (dm.th.ed)/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "Sint: " << (dm.th.en)/dm.nB << endl;
    }
    cout << "<A>: " << dm.average_A() << endl;
    cout << "<Z>: " << dm.average_Z() << endl;
    cout << "Xalpha: " << Xa << endl;
    cout << "Xn: " << dm.n.n/dm.nB << endl;
    cout << "Xp: " << dm.p.n/dm.nB << endl;
    cout << "Xnuclei: " << Xnuclei << endl;
  }
  return;
}

int eos_nse_full::calc_density_nr(dense_matter &dm, int verbose) {

  double factor=1.0e-10;

  // Temporary storage for a distribution index
  size_t index;

  // Temporary storage for the baryon density in nuclei
  double nB_nuc;
  
  char ch;
  int ret;

  // Output
  calc_density_noneq_nr(dm,1);
  cout << "Initial guess." << endl;
  cin >> ch;
  
  // Initial minimization
  ret=calc_density_fixcomp_nr(dm);
  cout << "ret: " << ret << endl;

  // Output
  calc_density_noneq_nr(dm,1);
  cout << "Post initial minimization." << endl;
  cin >> ch;
  
  // Main loop
  bool done=false;
  while (done==false) {

    // Record guess distribution
    dense_matter dm2=dm;
    
    // Add new nuclei
    nB_nuc=dm.nB_nuclei();
    vector<nucleus> new_nuclei;
    if (false) {
      for(size_t i=0;i<dm.dist.size();i++) {
	if (dm.dist[i].n>nB_nuc*factor) {
	  int Z=dm.dist[i].Z;
	  int N=dm.dist[i].N;
	  for(int iz=-1;iz<=1;iz++) {
	    for(int in=-1;in<=1;in++) {
	      if (!dm.nuc_in_dist(Z+iz,N+in,index)) {
		nucleus nuc;
		nuc.Z=Z+iz;
		nuc.N=N+in;
		nuc.A=Z+N;
		nuc.n=nB_nuc*factor;
		nuc.g=2.0;
		new_nuclei.push_back(nuc);
		dm.dist.push_back(nuc);
	      }
	    }
	  }
	}
      }
    } else {
      nucleus nuc;
      nuc.Z=25;
      nuc.N=25;
      nuc.A=50;
      nuc.n=nB_nuc/100.0;
      nuc.g=2.0;
      new_nuclei.push_back(nuc);
      dm.dist.push_back(nuc);
    }

    // Readjust densities
    ret=density_match_nr(dm);
    cout << "ret: " << ret << endl;

    // Output
    ret=calc_density_noneq_nr(dm,1);
    cout << "retx: " << ret << endl;
    cout << "Post density match." << endl;
    cin >> ch;

    for(size_t i=0;i<5;i++) {

      // Perform new minimization
      cout << "Going to fixcomp: " << endl;
      ret=calc_density_fixcomp_nr(dm);
      cout << "ret: " << ret << endl;
      
      // Output
      calc_density_noneq_nr(dm,1);
      cout << "Post iterative minimization: " << i << endl;
      cin >> ch;
    }

    // Prune distribution
    nB_nuc=dm.nB_nuclei();
    for(vector<nucleus>::iterator it=dm.dist.begin();
	it!=dm.dist.end();it++) {
      if (it->n<nB_nuc*factor) {
	dm.dist.erase(it);
	it=dm.dist.begin();
      }
    }
    
    // Output
    calc_density_noneq_nr(dm,1);
    cout << "Post prune." << endl;
    cin >> ch;

    // Test to see if distribution has changed
    done=true;
    for(size_t i=0;i<new_nuclei.size();i++) {
      if (dm.nuc_in_dist(new_nuclei[i].Z,new_nuclei[i].N,index)) done=false;
    }

    // If the distribution has changed, perform another iteration
    cout << "Done: " << done << endl;
  }

  return 0;
}

