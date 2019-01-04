/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <o2scl/nucmass_densmat.h>

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

#ifdef O2SCL_NEVER_DEFINED

void dense_matter::dist_properties(double &ntot, double &A, double &Z,
				   double &N, double &a) {
				   

  ntot=0.0;
  A=0.0;
  Z=0.0;
  N=0.0;
    
  for (size_t i=0;i<dist.size();i++) {
    double den_i=dist[i].n;
    ntot+=den_i;
    A+=dist[i].A*den_i;
    Z+=dist[i].Z*den_i;
    N+=dist[i].N*den_i;
  }

  A/=ntot;
  Z/=ntot;
  N/=ntot;
  a=cbrt(3.0/4.0/o2scl_const::pi/ntot);

  return;
}

double dense_matter::impurity(double ntot, double ave_Z) {
    
  double sum=0.0;
  
  for (size_t i=0;i<dist.size();i++) {
    sum+=dist[i].n*pow(dist[i].Z-ave_Z,2.0);
  }
  
  if (ntot==0.0) return 0.0;

  return sum/ntot;
}

void dense_matter::densities(double &nn, double &np) {
  double nn=n.n;
  double np=p.n;
  for (size_t i=0;i<dist.size();i++) {
    nn+=dist[i].n*dist[i].N;
    np+=dist[i].n*dist[i].Z;
  }
  return;
}

#endif

void dense_matter::output(std::ostream &out, int verbose) {

  out.setf(ios::scientific);
  out.setf(ios::left);
  out.setf(ios::showpos);
  out << "--------------------------------------"
      << "--------------------------------------" << endl;
  out << "       nB= " << nB << " Ye= " << Ye << " T= "
      << T*hc_mev_fm << endl;
  out << "Check: nB= " << baryon_density() << " Ye= " << electron_fraction()
      << endl;
  cout << "Nuclei: A= " << average_A() << "  Z= "
       << average_Z() << " Q= " << impurity() << endl;
  out << endl;
  out << "                                n            mu          "
      << "  ed            en" << endl;
  out.width(20);
  out << "Neutrons: " << n.n << " " << n.mu << endl;
  out.width(20);
  out << "Protons: " << p.n << " " << p.mu << endl;
  out.width(48);
  out << "Dripped nucleons: " 
      << drip_th.ed << " " << drip_th.en << endl;
  out.width(20);
  out << "Electrons: " 
      << e.n << " " << e.mu << " " << e.ed << " "
      << e.en << endl;
  out.width(20);
  out << "Muons: " 
      << mu.n << " " << mu.mu << " " << mu.ed << " "
      << mu.en << endl;
  out.width(20);
  out << "Photons: " 
      << photon.n << " " << 0.0 << " " << photon.ed << " "
      << photon.en << endl;
  
  size_t i_out;
  
  if (verbose>0) {
    
    i_out=0;
  
    for(size_t i=0;i<dist.size();i++) {
      o2scl::nucleus &nuc=dist[i];
    
      if ((verbose==1 && nuc.n>0.0 && i_out<10) ||
	  (verbose==2 && nuc.n>0.0) ||
	  verbose>=3) {
	string s="Nucleus ("+o2scl::itos(((int)(nuc.Z+1.0e-8)))+","+
	  o2scl::itos(((int)(nuc.N+1.0e-8)))+"): ";
	out.width(20);
	out << s << nuc.n << " " << nuc.mu << " " 
	    << nuc.ed+nuc.be*nuc.n << " " << nuc.en << endl;
	i_out++;
      }
    
    }

  } else {
    
    double sum_e=0.0, sum_s=0.0;
    for(size_t i=0;i<dist.size();i++) {
      o2scl::nucleus &nuc=dist[i];
      if (nuc.n>0.0) {
	sum_e+=nuc.ed+nuc.be*nuc.n;
	sum_s+=nuc.en;
      }
    }
    out << "Nuclei:                                         "
	<< sum_e << " " << sum_s << endl;
  }
  
  out.width(48);
  out << "Total: " << th.ed << " " << th.en << endl;
  out << "Free energy: " << th.ed-T*th.en << endl;
  out << endl;
  out.width(48);
  out << "Contributions to pressure:" << " pr" << endl;
  out.width(48);
  out << "- Energy: " << -th.ed << endl;
  out.width(48);
  out << "T * Entropy:         n             eta"
      << T*th.en << endl;
  out.width(20);
  out << "Neutrons: " << n.n << " " << eta_n << " "
      << n.n*eta_n << endl;
  out.width(20);
  out << "Protons: " << p.n << " " << eta_p << " "
      << p.n*eta_p << endl;

  if (verbose>0) {
    
    i_out=0;
    for(size_t i=0;i<dist.size();i++) {
      o2scl::nucleus &nuc=dist[i];
      
      if ((verbose==1 && nuc.n>0.0 && i_out<10) ||
	  (verbose==2 && nuc.n>0.0) ||
	  verbose>=3) {
	string s="Nucleus ("+o2scl::itos(((int)(nuc.Z+1.0e-8)))+","+
	  o2scl::itos(((int)(nuc.N+1.0e-8)))+"): ";
	out.width(20);
	out << s << nuc.n << " " << eta_nuc[i] << " "
	    << nuc.n*eta_nuc[i] << endl;
	i_out++;
      }
    }

  } else {
    double sum_p=0.0;
    for(size_t i=0;i<dist.size();i++) {
      o2scl::nucleus &nuc=dist[i];
      if (nuc.n>0.0) {
	sum_p+=nuc.n*eta_nuc[i];
      }
    }
    out << "Nuclei:                                         "
	<< sum_p << endl;
  }
  
  out.width(48);
  out << "Total pressure: " << th.pr << endl;
  out << "--------------------------------------"
      << "--------------------------------------" << endl;
  out.unsetf(ios::left);
  out.unsetf(ios::showpos);

  return;
}

double dense_matter::average_a() {
  double ntot=0.0;
  for (size_t i=0;i<dist.size();i++) {
    ntot+=dist[i].n;
  }
  if (ntot==0.0) return 0.0;
  return cbrt(3.0/4.0/o2scl_const::pi/ntot);
}

double dense_matter::nn_tot() {
  double nn_tot_val=n.n;
  for (size_t i=0;i<dist.size();i++) {
    nn_tot_val+=dist[i].n*dist[i].N;
  }
  return nn_tot_val;
}

double dense_matter::np_tot() {
  double np_tot_val=p.n;
  for (size_t i=0;i<dist.size();i++) {
    np_tot_val+=dist[i].n*dist[i].Z;
  }
  return np_tot_val;
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

double dense_matter::baryon_density_nuclei() {
    
  double nBtot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    nBtot+=dist[i].n*(dist[i].Z+dist[i].N);
  }
  
  return nBtot;
}

double dense_matter::baryon_density() {
    
  double nBtot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    nBtot+=dist[i].n*(dist[i].Z+dist[i].N);
  }
  
  return nBtot+n.n+p.n;
}

double dense_matter::electron_fraction() {
    
  double nptot=0.0;

  for (size_t i=0;i<dist.size();i++) {
    nptot+=dist[i].n*dist[i].Z;
  }
  
  return (nptot+p.n)/baryon_density();
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

bool dense_matter::nuc_in_dist(int Z, int N, size_t &index) {
  for(size_t i=0;i<dist.size();i++) {
    if (dist[i].Z==Z && dist[i].N==N) {
      index=i;
      return true;
    }
  }
  return false;
}

void dense_matter::prune_distribution(double factor) {
  
  double nB_nuc=baryon_density_nuclei();
  for(vector<nucleus>::iterator it=dist.begin();it!=dist.end();it++) {
    if (it->n<nB_nuc*factor) {
      dist.erase(it);
      it=dist.begin();
    }
  }
  
  return;
}

void dense_matter::copy_densities_from(dense_matter &dm2) {
  for(size_t i=0;i<dist.size();i++) {
    dist[i].n=0.0;
  }
  for(size_t i=0;i<dm2.dist.size();i++) {
    for(size_t j=0;j<dist.size();j++) {
      if (dist[j].Z==dm2.dist[i].Z &&
	  dist[j].N==dm2.dist[i].N) {
	dist[j].n=dm2.dist[i].n;
	j=dist.size();
      }
    }
  }
  return;
}

nucmass_densmat::nucmass_densmat() {
  // It's important not to automatically load masses from
  // HDF5 by default because this causes issues instantiating
  // this class with many processors
  massp=0;
}

void nucmass_densmat::set_mass(nucmass &nm) {
  massp=&nm;
  return;
}

void nucmass_densmat::test_derivatives(double eps, double &t1, double &t2, 
				       double &t3, double &t4) {

  if (massp==0) {
    O2SCL_ERR("Masses not specified in nucmass_densmat::test_derivatives().",
	      exc_efailed);
  }
      
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

void nucmass_densmat::binding_energy_densmat
(double Z, double N, double npout, double nnout, 
 double nneg, double T, double &E) {
  double dEdnn;
  double dEdnp;
  double dEdnneg;
  double dEdT;
  return binding_energy_densmat_derivs
    (Z,N,nnout,npout,nneg,T,E,dEdnp,dEdnn,dEdnneg,dEdT);
}

void nucmass_densmat::binding_energy_densmat_derivs
(double Z, double N, double npout, double nnout, 
 double nneg, double T, double &E, double &dEdnp, double &dEdnn,
 double &dEdnneg, double &dEdT) {

  if (massp==0) {
    O2SCL_ERR("Masses not specified in nucmass_densmat::test_derivatives().",
	      exc_efailed);
  }

  if (!massp->is_included(Z+1.0e-10,N+1.0e-10)) {
    O2SCL_ERR((((string)"Mass with Z=")+o2scl::dtos(Z)+" and N="+
	       o2scl::dtos(N)+" not included in nucmass_densmat"+
	       "::binding_energy_densmat_derivs().").c_str(),exc_einval);
  }

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

  // Allow a small error from finite precision
  if (R_p>R_WS) {
    R_p=R_WS*(1.0-1.0e-8);
    R_p_3=R_p*R_p*R_p;
  }
  if (R_p>R_WS) {
    cout << "Z,N,np,nn,ne: " << Z << " " << N << " " 
	 << npout << " " << nnout << " " << nneg << endl;
    cout << "Rn,Rp,RWS: " << R_n << " "<< R_p << " " << R_WS << endl;
    O2SCL_ERR2("Proton radius larger than cell in nucmass_densmat::",
	       "binding_energy_densmat_derivs().",exc_einval);
  }

  // Allow a small error from finite precision
  if (R_n>R_WS) {
    R_n=R_WS*(1.0-1.0e-8);
    R_n_3=R_n*R_n*R_n;
  }
  if (R_n>R_WS) {
    cout << "Z,N,np,nn,ne: " << Z << " " << N << " " 
	 << npout << " " << nnout << " " << nneg << endl;
    cout << "Rn,Rp,RWS: " << R_n << " "<< R_p << " " << R_WS << endl;
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
