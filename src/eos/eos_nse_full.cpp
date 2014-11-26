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
#include <o2scl/root_brent_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

eos_nse_full::eos_nse_full() {
  inc_lept_phot=true;

  // Ensure fermion_rel doesn't throw exceptions for convergence
  // errors for electrons
  err_nonconv=true;

  // It's important not to automatically load masses from
  // HDF5 by default because this causes issues instantiating
  // this class with many processors
  ehtp=0;
  
  massp=&nuc_dens;
  ad=0;

  // Make the default minimizer more accurate
  def_mmin.ntrial=10000;
  def_mmin.tol_rel=1.0e-7;
  def_mmin.tol_abs=1.0e-7;

  inc_prot_coul=true;
  include_muons=false;

  verbose=1;
}

double eos_nse_full::free_energy(const ubvector &n_nuc, dense_matter &dm) {

  double nB_nuc=0.0, np_nuc=0.0;

  if (n_nuc.size()!=dm.dist.size()) {
    O2SCL_ERR2("Sizes incommensurate in eos_nse_full::",
	       "free_energy().",exc_einval);
  }

  for(size_t i=0;i<n_nuc.size();i++) {
    dm.dist[i].n=n_nuc[i];
    nB_nuc+=n_nuc[i]*(dm.dist[i].Z+dm.dist[i].N);
    np_nuc+=n_nuc[i]*dm.dist[i].Z;
  }

  dm.p.n=dm.nB*dm.Ye-np_nuc;
  dm.n.n=dm.nB-nB_nuc-dm.p.n;

  // Return a large value if necessary

  int ret=calc_density_noneq(dm);
  if (ret==invalid_config) return 1.0e4;

  return dm.th.ed-dm.T*dm.th.en;
}

int eos_nse_full::calc_density_by_min(dense_matter &dm) {

  ubvector n_nuc(dm.dist.size()), n_nuc2(dm.dist.size());
  for(size_t i=0;i<n_nuc.size();i++) {
    n_nuc[i]=dm.dist[i].n;
    n_nuc2[i]=n_nuc[i]*1.01;
  }

  multi_funct11 mf=std::bind
    (std::mem_fn<double(const ubvector &,dense_matter &)>
     (&eos_nse_full::free_energy),
     this,std::placeholders::_2,std::ref(dm));

  double fr_min=0.0;

  int ret=def_mmin.mmin_twovec(n_nuc.size(),n_nuc,n_nuc2,fr_min,mf);
  if (ret!=success) {
    O2SCL_CONV2_RET("Minimizer failed in eos_nse_full::",
		    "calc_density_saha().",exc_ebadfunc,err_nonconv);
  }

  // Perform a final function evaluation (important to set the final
  // nuclear densities in the 'dm' parameter)
  fr_min=free_energy(n_nuc,dm);

  return 0;
}

int eos_nse_full::calc_density_saha(dense_matter &dm) {

  // Check user-specified conditions
  if (dm.nB<0.0 || dm.Ye<0.0 || dm.T<0.0) {
    O2SCL_ERR2("Cannot use negative densities or temperatures in ",
	       "eos_nse_full::calc_density_saha().",exc_einval);
  }
  
  if (dm.nB*dm.Ye>=0.08 || dm.T*hc_mev_fm>30.0) {
    dm.n.n=dm.nB*(1.0-dm.Ye);
    dm.p.n=dm.nB*dm.Ye;
    dm.e.n=dm.p.n;
    for(size_t i=0;i<dm.dist.size();i++) {
      dm.dist[i].n=0.0;
    }
    return calc_density_noneq(dm);
  }

  // Vector for initial guess and equations
  ubvector x(2), y(2);
  x[0]=dm.n.n;
  x[1]=dm.p.n;

  // Adjust if initial guesses are negative
  if (x[0]<0.0) x[0]=dm.nB*(1.0-dm.Ye);
  if (x[1]<0.0) x[1]=dm.nB*dm.Ye;

  // Adjust proton density if it's larger than electron density.
  // It's important to make the adjustment large enough so that
  // small steps from the solver don't make the configuration 
  // invalid. 
  if (dm.nB*dm.Ye<x[1]) {
    x[1]=dm.nB*dm.Ye*(1.0-1.0e-2);
  }
  
  // Perform initial function evaluation 
  int ret;
  ret=solve_fixnp(2,x,y,dm);
  if (ret!=success) {
    if (verbose>0) {
      cout << "Initial point failed in "
	   << "eos_nse_full::calc_density_saha()." << endl;
    }
    O2SCL_CONV2_RET("Initial point failed in eos_nse_full::",
		    "calc_density_saha().",exc_ebadfunc,err_nonconv);
  }

  // Iterate to fix density for initial guess if necessary
  size_t it=0;
  while (it<100 && dm.baryon_density()>10.0) {
    x[0]/=1.5;
    x[1]/=1.5;
    ret=solve_fixnp(2,x,y,dm);
    it++;
  }

  // Call solver
  mm_funct11 mf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,dense_matter &)>
     (&eos_nse_full::solve_fixnp),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(dm));
  ret=def_mroot.msolve(2,x,mf);
  if (ret!=success) {
    if (verbose>0) {
      cout << "Solver failed in eos_nse_full::calc_density_saha()." << endl;
    }
    O2SCL_CONV2_RET("Solver failed in eos_nse_full::",
		    "calc_density_saha().",exc_ebadfunc,err_nonconv);
  }

  // Final function evaluation
  ret=solve_fixnp(2,x,y,dm);
  if (ret!=success) {
    cout << "Final function evaluation failed in "
	 << "eos_nse_full::calc_density_saha()." << endl;
    O2SCL_CONV2_RET("Final function evaluation failed in eos_nse_full::",
		    "calc_density_saha().",exc_ebadfunc,err_nonconv);
  }
  return 0;
}

int eos_nse_full::solve_fixnp(size_t n, const ubvector &x, ubvector &y,
			      dense_matter &dm) {
  dm.n.n=x[0];
  dm.p.n=x[1];
  if (x[0]<0.0 || x[1]<0.0 || !o2scl::is_finite(x[0]) || 
      !o2scl::is_finite(x[1])) {
    if (verbose>0) {
      cout << "Not positive and finite point in "
	   << "eos_nse_full::solve_fixnp() " << x[0] << " " << x[1] << endl;
    }
    return exc_ebadfunc;
  }
  int ret=calc_density_fixnp(dm);
  if (ret!=0) return ret;
  y[0]=2.0*(dm.baryon_density()-dm.nB)/(dm.baryon_density()+dm.nB);
  y[1]=2.0*(dm.electron_fraction()-dm.Ye)/(dm.electron_fraction()+dm.Ye);
  if (fabs(y[0])>=2.0 || fabs(y[1])>=2.0) {
    return exc_ebadfunc;
  }
  if (!o2scl::is_finite(y[0]) || !o2scl::is_finite(y[1])) {
    if (verbose>0) {
      cout << "Not finite point returning from "
	   << "eos_nse_full::solve_fixnp()." << endl;
    }
    return exc_ebadfunc;
  }
  return 0;
}

int eos_nse_full::calc_density_fixnp(dense_matter &dm) {

  if (ehtp==0) {
    O2SCL_ERR2("Homogeneous matter EOS not specified in ",
	       "eos_nse_full::calc_density_fixnp().",exc_efailed);
  }

  // -----------------------------------------------------------
  // Sanity checks

  if (!o2scl::is_finite(dm.n.n) || 
      !o2scl::is_finite(dm.p.n)) {
    O2SCL_ERR2("Neutron or proton density not finite in ",
	       "eos_nse_full::calc_density_fixnp().",exc_esanity);
  }
  if (dm.n.m<0.0 || dm.p.m<0.0) {
    O2SCL_ERR2("Mass negative in ",
	       "eos_nse_full::calc_density_fixnp().",exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].inc_rest_mass==true) {
      O2SCL_ERR("Wrong inc_rest_mass for nuclei in fixnp().",
		exc_esanity);
    }
  }
  if (dm.n.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for neutrons in fixnp().",
	      exc_esanity);
  }
  if (dm.p.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for protons in fixnp().",
	      exc_esanity);
  }
  if (dm.e.inc_rest_mass==false) {
    O2SCL_ERR("Wrong inc_rest_mass for electrons in fixnp().",
	      exc_esanity);
  }
  if (dm.mu.inc_rest_mass==false) {
    O2SCL_ERR("Wrong inc_rest_mass for muons in fixnp().",
	      exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].non_interacting==false) {
      O2SCL_ERR("Wrong non_interacting for nuclei in fixnp().",
		exc_esanity);
    }
  }
  if (dm.n.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for neutrons in fixnp().",
	      exc_esanity);
  }
  if (dm.p.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for protons in fixnp().",
	      exc_esanity);
  }
  if (dm.e.non_interacting==false) {
    O2SCL_ERR("Wrong non_interacting for electrons in fixnp().",
	      exc_esanity);
  }
  if (dm.mu.non_interacting==false) {
    O2SCL_ERR("Wrong non_interacting for muons in fixnp().",
	      exc_esanity);
  }
  if (inc_prot_coul==true) {
    O2SCL_ERR2("Including protons in Coulomb not implemented in ",
	       "eos_nse_full::calc_density_fixnp().",exc_eunimpl);
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
  if (ret!=success) {
    if (verbose>0) {
      cout << "Homogeneous nucleonic EOS failed in "
	   << "eos_nse_full::calc_density_fixnp()." << endl;
    }
    O2SCL_CONV2_RET("Homogeneous nucleon EOS failed in eos_nse_full::",
		    "calc_density_fixnp().",exc_einval,err_nonconv);
  }

  // Add contribution from dripped nucleons
  dm.th=dm.drip_th;

  // -----------------------------------------------------------
  // Verbose output

  if (verbose>1) {
    cout.setf(ios::left);
    cout.setf(ios::showpos);
    cout << "--------------------------------------"
	 << "--------------------------------------" << endl;
    cout << "fixnp():" << endl;
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
    dm.e.n=dm.Ye*dm.nB;
    ret=relf.pair_density(dm.e,dm.T);
    if (ret!=success) {
      if (verbose>0) {
	cout << "Electron EOS failed in "
	     << "eos_nse_full::calc_density_fixnp()." << endl;
      }
      O2SCL_CONV2_RET("Electron EOS failed in eos_nse_full::",
		      "calc_density_fixnp().",exc_einval,err_nonconv);
    }
    dm.th.ed+=dm.e.ed;
    dm.th.en+=dm.e.en;

    if (verbose>1) {
      cout.width(20);
      cout << "Electrons: " 
	   << dm.e.n << " " << dm.e.mu << " " << dm.e.ed << " "
	   << dm.e.en << endl;
    }

    if (include_muons) {
      // Add muons
      dm.mu.mu=dm.e.mu;
      relf.pair_mu(dm.mu,dm.T);
      dm.th.ed+=dm.mu.ed;
      dm.th.en+=dm.mu.en;
      
      if (verbose>1) {
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

    if (verbose>1) {
      cout.width(20);
      cout << "Photons: " 
	   << dm.photon.n << " " << 0.0 << " " << dm.photon.ed << " "
	   << dm.photon.en << endl;
    }
  }

  // Negative charge density
  double n_neg=dm.Ye*dm.nB;
  if (include_muons) n_neg+=dm.mu.n;

  // -----------------------------------------------------------
  // Some invalid configuration checks

  if (n_neg<dm.p.n) {
    if (verbose>0) {
      cout << "Negative charge density too small to match proton density."
	   << "\n\tineos_nse_full::calc_density_fixnp()." << endl;
      cout << "np: " << dm.p.n << " ne: " << dm.e.n 
	   << " nmu: " << dm.mu.n << endl;
    }
    return invalid_config;
  }

  if (dm.dist.size()>0) {

    // Check that the proton density isn't too large in
    // the presence of nuclei
    if (dm.p.n>0.08) {
      if (verbose>0) {
	cout << "External proton density (" << dm.p.n 
	     << ") too large." << endl;
      }
      return invalid_config;
    }
    
    // Check that the negative charge density isn't too large
    // in the presence of nuclei
    double fac=(0.08-dm.p.n)/(n_neg-dm.p.n);
    if (1.0>=fac) {
      if (verbose>0) {
	cout << "Proton radius negative or larger than cell size in "
	     << "\n\teos_nse_full::calc_density_fixnp()." << endl;
	cout << "fac: " << fac << " " << n_neg << " " << dm.p.n << endl;
      }
      return invalid_config;
    }
  }

  // -----------------------------------------------------------
  // Properties of nuclear distribution

  // Vectors for derivatives of nuclear binding energy
  ubvector vec_dEdnneg(dm.dist.size());

  // Stepsize for verbose output
  size_t i_out=0, out_step=dm.dist.size()/10;
  if (out_step==0) out_step=1;

  // Main loop
  for(size_t i=0;i<dm.dist.size();i++) {

    // Create a reference for this nucleus
    nucleus &nuc=dm.dist[i];

    double condition;
    if (inc_prot_coul) {
      condition=nuc.N*(n_neg-dm.p.n)/nuc.Z/(0.08-dm.n.n);
    } else {
      condition=nuc.N*(n_neg)/nuc.Z/(0.08);
    }

    // If this nucleus is unphysical because R_n > R_{WS}, 
    // set it's density to zero and continue
    if (condition>=1.0) {
      
      nuc.n=0.0;
      nuc.ed=0.0;
      nuc.en=0.0;
      vec_dEdnneg[i]=0.0;

    } else {
    
      // Compute nuclear binding energy and total mass
      double dEdnp, dEdnn, dEdnneg, dEdT;
      // Include protons in the Coulomb energy
      if (inc_prot_coul) {
	massp->binding_energy_densmat_derivs
	  (nuc.Z,nuc.N,dm.p.n,dm.n.n,n_neg,dm.T,nuc.be,dEdnp,dEdnn,
	   dEdnneg,dEdT);
      } else {
	massp->binding_energy_densmat_derivs
	  (nuc.Z,nuc.N,0.0,0.0,n_neg,dm.T,nuc.be,dEdnp,dEdnn,dEdnneg,dEdT);
      }
      nuc.be/=hc_mev_fm;
      nuc.m=nuc.Z*dm.p.m+nuc.N*dm.n.m+nuc.be;
      vec_dEdnneg[i]=dEdnneg;
      
      // Use NSE to compute the chemical potential
      nuc.mu=nuc.Z*dm.p.mu+nuc.N*dm.n.mu-nuc.be;
    
      // Translational energy
      cla.calc_mu(nuc,dm.T);

      // Update thermo object with information from nucleus
      dm.th.ed+=nuc.be*nuc.n+nuc.ed;
      dm.th.en+=nuc.en;

      if ((verbose>1 && i==i_out) || verbose>2) {
	string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	  itos(((int)(nuc.N+1.0e-8)))+"): ";
	cout.width(20);
	cout << s << nuc.n << " " << nuc.mu << " " 
	     << nuc.ed+nuc.be*nuc.n << " " << nuc.en << endl;
	i_out+=out_step;
      }

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
    if (dm.dist[i].n>0.0) {
      dm.eta_nuc[i]=dm.dist[i].be+dm.dist[i].mu+dm.dist[i].Z*dm.e.mu;
    } else {
      dm.eta_nuc[i]=0.0;
    }
    // In eta_p, we don't include dEdnp terms which are zero
    dm.eta_p+=(dm.dist[i].n+dfdm_i)*(vec_dEdnneg[i])/hc_mev_fm;
    
    for(size_t j=0;j<dm.dist.size();j++) {
      
      if (dm.dist[i].n>0.0) {
	
	double dmudm=-1.5*dm.T/dm.dist[j].m;
	double dfdm=dm.dist[j].n*dmudm;
	
	dm.eta_nuc[i]+=dm.dist[i].Z*(dm.dist[j].n+dfdm)*
	  vec_dEdnneg[j]/hc_mev_fm;
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

  if (verbose>1) {
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
    i_out=0;
    for(size_t i=0;i<dm.dist.size();i++) {
      nucleus &nuc=dm.dist[i];
      if (i==i_out || verbose>2) {
	string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	  itos(((int)(nuc.N+1.0e-8)))+"): ";
	cout.width(20);
	cout << s << nuc.n << " " << dm.eta_nuc[i] << " "
	     << nuc.n*dm.eta_nuc[i] << endl;
	i_out+=out_step;
      }
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

int eos_nse_full::calc_density_noneq(dense_matter &dm) {
  
  if (ehtp==0) {
    O2SCL_ERR2("Homogeneous matter EOS not specified in ",
	       "eos_nse_full::calc_density_noneq().",exc_efailed);
  }

  // -----------------------------------------------------------
  // Sanity checks

  if (!o2scl::is_finite(dm.n.n) || 
      !o2scl::is_finite(dm.p.n)) {
    O2SCL_ERR2("Neutron or proton density not finite in ",
	       "eos_nse_full::calc_density_noneq().",exc_esanity);
  }
  if (dm.n.m<0.0 || dm.p.m<0.0) {
    O2SCL_ERR2("Mass negative in ",
	       "eos_nse_full::calc_density_noneq().",exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].inc_rest_mass==true) {
      O2SCL_ERR("Wrong inc_rest_mass for nuclei in noneq().",
		exc_esanity);
    }
  }
  if (dm.n.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for neutrons in noneq().",
	      exc_esanity);
  }
  if (dm.p.inc_rest_mass==true) {
    O2SCL_ERR("Wrong inc_rest_mass for protons in noneq().",
	      exc_esanity);
  }
  if (dm.e.inc_rest_mass==false) {
    O2SCL_ERR("Wrong inc_rest_mass for electrons in noneq().",
	      exc_esanity);
  }
  if (dm.mu.inc_rest_mass==false) {
    O2SCL_ERR("Wrong inc_rest_mass for muons in noneq().",
	      exc_esanity);
  }
  for(size_t i=0;i<dm.dist.size();i++) {
    if (dm.dist[i].non_interacting==false) {
      O2SCL_ERR("Wrong non_interacting for nuclei in noneq().",
		exc_esanity);
    }
  }
  if (dm.n.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for neutrons in noneq().",
	      exc_esanity);
  }
  if (dm.p.non_interacting==true) {
    O2SCL_ERR("Wrong non_interacting for protons in noneq().",
	      exc_esanity);
  }
  if (dm.e.non_interacting==false) {
    O2SCL_ERR("Wrong non_interacting for electrons in noneq().",
	      exc_esanity);
  }
  if (dm.mu.non_interacting==false) {
    O2SCL_ERR("Wrong non_interacting for muons in noneq().",
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
  if (ret!=success) {
    O2SCL_CONV2_RET("Homogeneous nucleon EOS failed in eos_nse_full::",
		    "calc_density_noneq().",exc_einval,err_nonconv);
  }

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

  if (verbose>1) {
    cout.setf(ios::left);
    cout.setf(ios::showpos);
    cout << "--------------------------------------"
	 << "--------------------------------------" << endl;
    cout << "noneq():" << endl;
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
    ret=relf.pair_density(dm.e,dm.T);
    if (ret!=success) {
      if (verbose>0) {
	cout << "Electron EOS failed in "
	     << "eos_nse_full::calc_density_noneq()." << endl;
      }
      O2SCL_CONV2_RET("Electron EOS failed in eos_nse_full::",
		      "calc_density_noneq().",exc_einval,err_nonconv);
    }
    dm.th.ed+=dm.e.ed;
    dm.th.en+=dm.e.en;

    if (verbose>1) {
      cout.width(20);
      cout << "Electrons: " 
	   << dm.e.n << " " << dm.e.mu << " " << dm.e.ed << " "
	   << dm.e.en << endl;
    }
    
    if (include_muons) {
      // Add muons
      dm.mu.mu=dm.e.mu;
      relf.pair_mu(dm.mu,dm.T);
      dm.th.ed+=dm.mu.ed;
      dm.th.en+=dm.mu.en;
      
      if (verbose>1) {
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

    if (verbose>1) {
      cout.width(20);
      cout << "Photons: " 
	   << dm.photon.n << " " << 0.0 << " " << dm.photon.ed << " "
	   << dm.photon.en << endl;
    }
  }

  double n_neg=dm.e.n;
  if (include_muons) n_neg+=dm.mu.n;

  // -----------------------------------------------------------
  // Finite precision errors can cause n_neg<dm.p.n, so we make a
  // small adjustment
  
  if (n_neg<dm.p.n) {
    dm.e.n*=(1.0+1.0e-12);
    n_neg=dm.e.n;
    if (include_muons) n_neg+=dm.mu.n;
  }

  // -----------------------------------------------------------
  // Some invalid configuration checks

  if (n_neg<dm.p.n) {
    if (verbose>0) {
      cout << "Negative charge density too small to match proton density"
	   << "\n\tin calc_density_noneq()." << endl;
      cout << "np: " << dm.p.n << " ne: " << dm.e.n 
           << " nmu: " << dm.mu.n << " n_neg: " << n_neg << " " 
	   << n_neg-dm.p.n << endl;
    }
    return invalid_config;
  }  

  // -----------------------------------------------------------
  // Properties of nuclear distribution

  // Vectors for derivatives of nuclear binding energy
  ubvector vec_dEdnneg(dm.dist.size()), vec_dEdnp(dm.dist.size());

  // Step size for verbose output
  size_t i_out=0, out_step=dm.dist.size()/10;
  if (out_step==0) out_step=1;

  // Main loop
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
      double fac=(0.08-dm.p.n)/(n_neg-dm.p.n);
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
      double dEdnp, dEdnn, dEdnneg, dEdT;
      if (inc_prot_coul) {
	// Include protons in the Coulomb energy
	massp->binding_energy_densmat_derivs
	  (nuc.Z,nuc.N,dm.p.n,dm.n.n,n_neg,dm.T,nuc.be,
	   dEdnp,dEdnn,dEdnneg,dEdT);
	nuc.be/=hc_mev_fm;
	nuc.m=nuc.Z*dm.p.m+nuc.N*dm.n.m+nuc.be;
	vec_dEdnp[i]=dEdnp;
	vec_dEdnneg[i]=dEdnneg;
      } else {
	// Don't include protons in the Coulomb energy
	massp->binding_energy_densmat_derivs
	  (nuc.Z,nuc.N,0.0,0.0,n_neg,dm.T,nuc.be,
	   dEdnp,dEdnn,dEdnneg,dEdT);
	nuc.be/=hc_mev_fm;
	nuc.m=nuc.Z*dm.p.m+nuc.N*dm.n.m+nuc.be;
	vec_dEdnp[i]=0.0;
	vec_dEdnneg[i]=dEdnneg;
      }

      // Translational energy
      cla.calc_density(nuc,dm.T);

      // Update thermo object with information from nucleus
      dm.th.ed+=nuc.be*nuc.n+nuc.ed;
      dm.th.en+=nuc.en;

      if (verbose>2 || (verbose>1 && i==i_out)) {
	string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	  itos(((int)(nuc.N+1.0e-8)))+"): ";
	cout.width(20);
	cout << s << nuc.n << " " << nuc.mu << " " 
	     << nuc.ed+nuc.be*nuc.n << " " << nuc.en << endl;
	i_out+=out_step;
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

    if (dm.dist[i].n>0.0) {
      
      double dmudm_i=-1.5*dm.T/dm.dist[i].m;
      double dfdm_i=dm.dist[i].n*dmudm_i;
      dm.eta_nuc[i]=dm.dist[i].be+dm.dist[i].mu+dm.dist[i].Z*dm.e.mu;
      dm.eta_p+=(dm.dist[i].n+dfdm_i)*(vec_dEdnp[i]+vec_dEdnneg[i])/hc_mev_fm;
      
      for(size_t j=0;j<dm.dist.size();j++) {
	
	double dmudm=-1.5*dm.T/dm.dist[j].m;
	double dfdm=dm.dist[j].n*dmudm;
	
	dm.eta_nuc[i]+=dm.dist[i].Z*(dm.dist[j].n+dfdm)*
	  vec_dEdnneg[j]/hc_mev_fm;
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

  if (verbose>1) {
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
    i_out=0;
    for(size_t i=0;i<dm.dist.size();i++) {
      nucleus &nuc=dm.dist[i];
      string s="Nucleus ("+itos(((int)(nuc.Z+1.0e-8)))+","+
	itos(((int)(nuc.N+1.0e-8)))+"): ";
      if (i==i_out || verbose>2) {
	cout.width(20);
	cout << s << nuc.n << " " << dm.eta_nuc[i] << " "
	     << nuc.n*dm.eta_nuc[i] << endl;
	i_out+=out_step;
      }
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

int eos_nse_full::density_match(dense_matter &dm) {
  
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
  int ret=calc_density_noneq(dm);
  
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
      ret=calc_density_noneq(dm);
      shift*=10.0;
    }
  }

  // If we couldn't fix, throw
  if (ret==invalid_config) {
    O2SCL_ERR2("Could not find valid configuration in ",
	       "eos_nse_full::density_match().",exc_efailed);
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
    calc_density_noneq(dm);
    cout << "Adjusting neutrons by: " << nn_fix-(1.0-dm.Ye)*dm.nB << endl;
    dm.n.n+=nn_fix-(1.0-dm.Ye)*dm.nB;
    cout << "Adjusting protons by: " << np_fix-dm.Ye*dm.nB << endl;
    dm.p.n+=np_fix-dm.Ye*dm.nB;
  }
  
  ret=calc_density_noneq(dm);

  if (ret==invalid_config) {
    O2SCL_ERR2("Did not produce valid configuration in ",
	       "eos_nse_full::density_match().",exc_efailed);
  }

  if (verbose>1) {
    cout << "Density match: " << endl;
    cout << (1.0-dm.Ye)*dm.nB << " " << nn_fix << endl;
    cout << dm.Ye*dm.nB << " " << np_fix << endl;
    cout << endl;
  }

  if (fabs((1.0-dm.Ye)*dm.nB-nn_fix)/nn_fix>1.0e-6 ||
      fabs(dm.Ye*dm.nB-np_fix)/np_fix>1.0e-6) {
    O2SCL_ERR2("Density match failed in ",
	       "eos_nse_full::density_match().",exc_esanity);
  }

  return ret;
}

void eos_nse_full::output(dense_matter &dm, int output_level) {
  cout << "nB=" << dm.nB << " fm^{-3}, Ye=" << dm.Ye << ", T="
       << dm.T*hc_mev_fm << " MeV" << endl;

  // Output nuclear densities and compute Xalpha and Xnuclei
  double nBnuc=0.0;
  double Xa=0.0, Xnuclei=0.0;
  size_t i_out=0, out_step=dm.dist.size()/10;
  if (out_step==0) out_step=1;
  for(size_t i=0;i<dm.dist.size();i++) {
    nBnuc+=dm.dist[i].n*(dm.dist[i].N+dm.dist[i].Z);
    if (dm.dist[i].Z==2 && dm.dist[i].N==2) Xa=dm.dist[i].n*4.0;
    else Xnuclei+=dm.dist[i].n*(dm.dist[i].Z+dm.dist[i].N);
    if (output_level>1 || (output_level==1 && i==i_out)) {
      cout << "Z,N,n: " << dm.dist[i].Z << " " << dm.dist[i].N << " "
	   << dm.dist[i].n << endl;
      i_out+=out_step;
    }
  }
  Xa/=dm.nB;
  Xnuclei/=dm.nB;

  // Main output
  cout << "nn,np,nBnuc: " << dm.n.n << " " << dm.p.n << " " 
       << nBnuc << " fm^{-3}" << endl;
  cout << "N,<Z>,<N>,<Q>: " << dm.dist.size() << " " 
       << dm.average_Z() << " "
       << dm.average_N() << " " << dm.impurity() << endl;
  cout << "fr: " << dm.th.ed-dm.T*dm.th.en << " fm^{-4}" << endl;
  if (output_level>=1) {
    if (inc_lept_phot) {
      cout << "F: " << (dm.th.ed-dm.T*dm.th.en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "E: " << (dm.th.ed)/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "P: " << dm.th.pr*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "S: " << (dm.th.en)/dm.nB << endl;
      double ed=dm.th.ed-dm.e.ed-dm.photon.ed;
      double pr=dm.th.pr-dm.e.pr-dm.photon.pr;
      double en=dm.th.en-dm.e.en-dm.photon.en;
      if (include_muons) {
	ed-=dm.mu.ed;
	pr-=dm.mu.pr;
	en-=dm.mu.en;
      }
      cout << "Fint: " << (ed-dm.T*en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "Eint: " << ed/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "Pint: " << pr*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "Sint: " << en/dm.nB << endl;
    } else {
      cout << "Fint: " << (dm.th.ed-dm.T*dm.th.en)/dm.nB*hc_mev_fm 
	   << " MeV" << endl;
      cout << "Eint: " << (dm.th.ed)/dm.nB*hc_mev_fm << " MeV" << endl;
      cout << "Pint: " << dm.th.pr*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "Sint: " << (dm.th.en)/dm.nB << endl;
    }
    cout << "<A>: " << dm.average_A() << endl;
    cout << "<Z>: " << dm.average_Z() << endl;
    cout << "Xalpha: " << Xa << endl;
    cout << "Xn: " << dm.n.n/dm.nB << endl;
    cout << "Xp: " << dm.p.n/dm.nB << endl;
    cout << "Xnuclei: " << Xnuclei << endl;
    cout << "mu_n: " << dm.n.mu*hc_mev_fm << " MeV" << endl;
    cout << "mu_p: " << dm.p.mu*hc_mev_fm << " MeV" << endl;
    cout << "eta_n: " << dm.eta_n*hc_mev_fm << " MeV" << endl;
    cout << "eta_p: " << dm.eta_p*hc_mev_fm << " MeV" << endl;
    if (include_muons) {
      cout << "Xmu: " << dm.mu.n/dm.nB << endl;
    }
  }
  return;
}

#ifdef O2SCL_NEVER_DEFINED

int eos_nse_full::calc_density(dense_matter &dm) {

  double factor=1.0e-10;

  // Temporary storage for a distribution index
  size_t index;

  // Temporary storage for the baryon density in nuclei
  double nB_nuc;
  
  char ch;
  int ret;

  // Output
  calc_density_noneq(dm,1);
  cout << "Initial guess." << endl;
  cin >> ch;
  
  // Initial minimization
  ret=calc_density_by_min(dm);
  cout << "ret: " << ret << endl;

  // Output
  calc_density_noneq(dm,1);
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
    ret=density_match(dm);
    cout << "ret: " << ret << endl;

    // Output
    ret=calc_density_noneq(dm,1);
    cout << "retx: " << ret << endl;
    cout << "Post density match." << endl;
    cin >> ch;

    for(size_t i=0;i<5;i++) {

      // Perform new minimization
      cout << "Going to by_min: " << endl;
      ret=calc_density_by_min(dm);
      cout << "ret: " << ret << endl;
      
      // Output
      calc_density_noneq(dm,1);
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
    calc_density_noneq(dm,1);
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

#endif
