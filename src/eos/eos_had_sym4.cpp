/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_had_sym4.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_sym4_base::eos_had_sym4_base() {
  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_electron),2.0);
}
  
int eos_had_sym4_base::calc_e_alpha(fermion &ne, fermion &pr, thermo &lth,
				double &alphak, double &alphap, double &alphat,
				double &diff_kin, double &diff_pot,
				double &ed_kin_nuc, double &ed_pot_nuc) {
  double nn=ne.n, np=pr.n;
  double eden, pres, mun, mup;
  int r1, r2, r3;
    
  double ed_kin, ed_pot, mu_n_kin, mu_p_kin, mu_n_pot, mu_p_pot;
    
  /// Compute the properties of neutron matter
    
  ne.n=nn+np;
  pr.n=0.0;
  r1=calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		mu_n_pot,mu_p_pot);
  double ed_pot_neut=ed_pot;
  double ed_kin_neut=ed_kin;
  double ed_neut=ed_pot+ed_kin;

  /// Compute the properties of nuclear matter
    
  ne.n=(nn+np)/2.0;
  pr.n=(nn+np)/2.0;
  r2=calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		mu_n_pot,mu_p_pot);
  ed_pot_nuc=ed_pot;
  ed_kin_nuc=ed_kin;
  double ed_nuc=ed_pot+ed_kin;

  /// Compute the properties of neutron-rich matter

  ne.n=(nn+np)*3.0/4.0;
  pr.n=(nn+np)/4.0;
  r3=calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		mu_n_pot,mu_p_pot);
  lth.ed=ed_kin+ed_pot;
  ne.mu=mu_n_kin+mu_n_pot;
  pr.mu=mu_p_kin+mu_p_pot;
  lth.pr=-lth.ed+ne.mu*ne.n+pr.mu*pr.n;
    
  alphak=(ed_kin_neut-ed_kin)/(ed_kin-ed_kin_nuc);
  alphap=(ed_pot_neut-ed_pot)/(ed_pot-ed_pot_nuc);
  alphat=(ed_neut-ed_pot-ed_kin)/(ed_pot+ed_kin-ed_nuc);

  diff_kin=ed_kin_neut-ed_kin_nuc;
  diff_pot=ed_pot_neut-ed_pot_nuc;

  if (r1!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }
  if (r2!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }
  if (r3!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }

  return 0;
}
    
double eos_had_sym4_base::calc_muhat(fermion &ne, fermion &pr) {
  double ed_kin, ed_pot, munk, mupk, munp, mupp;

  calc_e_sep(ne,pr,ed_kin,ed_pot,munk,mupk,munp,mupp);

  e.n=pr.n;
  fzt2.calc_density_zerot(e);
      
  return munk+munp-mupk-mupp-e.mu;
}

int eos_had_sym4_rmf::calc_e_sep
(fermion &ne, fermion &pr, double &ed_kin, 
 double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
 double &mu_n_pot, double &mu_p_pot) {

  thermo lth;
    
  int ret=eos_had_rmf::calc_e(ne,pr,lth);

  ed_kin=ne.ed+pr.ed;
  ed_pot=lth.ed-ed_kin;
  mu_n_kin=ne.nu;
  mu_p_kin=pr.nu;
  mu_n_pot=ne.mu-ne.nu;
  mu_p_pot=pr.mu-pr.nu;

  return ret;
}

int eos_had_sym4_apr::calc_e_sep
(fermion &ne, fermion &pr, double &ed_kin, 
 double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
 double &mu_n_pot, double &mu_p_pot) {

  thermo lth;

  int ret=eos_had_apr::calc_e(ne,pr,lth);
    
  double barn=ne.n+pr.n;
  double xp=pr.n/barn;
  double dxdnn=-xp/barn, dxdnp=(1.0-xp)/barn;
  double t1=barn*exp(-par[4]*barn);
  double t2=par[3]+par[5]*(1-xp);
  double t3=par[3]+xp*par[5];

  ed_kin=ne.ed+pr.ed;
  ed_pot=lth.ed-ed_kin;

  double dmsndnn=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-
				   t1*par[5]*dxdnn);
  double dmsndnp=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-
				   t1*par[5]*dxdnp);
  double dmspdnn=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+
				   t1*par[5]*dxdnn);
  double dmspdnp=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+
				   t1*par[5]*dxdnp);
  mu_n_kin=ne.nu-(ne.ed-ne.n*ne.m)/ne.ms*dmsndnn-
    (pr.ed-pr.n*pr.m)/pr.ms*dmspdnn;
  mu_p_kin=pr.nu-(pr.ed-pr.n*pr.m)/pr.ms*dmspdnp-
    (ne.ed-ne.n*ne.m)/ne.ms*dmsndnp;
    
  mu_n_pot=ne.mu-mu_n_kin;
  mu_p_pot=pr.mu-mu_p_kin;
    
  return ret;
}
  
int eos_had_sym4_skyrme::calc_e_sep
(fermion &ne, fermion &pr, double &ed_kin, 
 double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
 double &mu_n_pot, double &mu_p_pot) {

  thermo lth;
    
  int ret=eos_had_skyrme::calc_e(ne,pr,lth);

  ed_kin=ne.ed+pr.ed;
  ed_pot=lth.ed-ed_kin;
    
  double term=0.25*(t1*(1.0+x1/2.0)+t2*(1.0+x2/2.0));
  double term2=0.25*(t2*(0.5+x2)-t1*(0.5+x1));
  double gn=2.0*ne.ms*(ne.ed-ne.n*ne.m);
  double gp=2.0*pr.ms*(pr.ed-pr.n*pr.m);
    
  mu_n_kin=(gn+gp)*term+ne.nu+gn*term2;
  mu_p_kin=(gn+gp)*term+pr.nu+gp*term2;

  mu_n_pot=ne.mu-mu_n_kin;
  mu_p_pot=pr.mu-mu_p_kin;

  return ret;
}

double eos_had_sym4_mdi::energy_kin(double var) {
  double n, hamk, ham, ham1, ham2, ham3=0.0, xp;

  if (mode==nmode) neutron->n=var;
  else if (mode==pmode) proton->n=var;

  if (neutron->n<0.0) neutron->n=0.01;
  if (proton->n<0.0) proton->n=0.01;
  
  n=neutron->n+proton->n;
  xp=proton->n/n;
  double delta=1.0-2.0*xp;
  
  neutron->ms=neutron->m;
  proton->ms=proton->m;

  fzt2.calc_density_zerot(*neutron);
  fzt2.calc_density_zerot(*proton);

  hamk=neutron->ed+proton->ed;

  return hamk;
}

double eos_had_sym4_mdi::energy_pot(double var) {
  double n, hamk, ham, ham1, ham2, ham3=0.0, xp;
    
  if (mode==nmode) neutron->n=var;
  else if (mode==pmode) proton->n=var;

  if (neutron->n<0.0) neutron->n=0.01;
  if (proton->n<0.0) proton->n=0.01;
  
  n=neutron->n+proton->n;
  xp=proton->n/n;
  double delta=1.0-2.0*xp;
  
  neutron->ms=neutron->m;
  proton->ms=proton->m;

  fzt2.calc_density_zerot(*neutron);
  fzt2.calc_density_zerot(*proton);

  if (form==mdi_form || form==gbd_form) {

    ham1=Au/rho0*neutron->n*proton->n+
      Al/2.0/rho0*(neutron->n*neutron->n+proton->n*proton->n);
    ham2=B/(sigma+1.0)*pow(n,sigma+1.0)/pow(rho0,sigma)*
      (1.0-x*delta*delta);
    
    if (form==mdi_form) {
      ham3=(Cl*(mom_integral(neutron->kf,neutron->kf)+
		mom_integral(proton->kf,proton->kf))+
	    Cu*(mom_integral(proton->kf,neutron->kf)+
		mom_integral(neutron->kf,proton->kf)))/rho0;
    } else {
      
      double gn, gp;
      gn=Lambda*Lambda/pi2*(neutron->kf-Lambda*atan(neutron->kf/Lambda));
      gp=Lambda*Lambda/pi2*(proton->kf-Lambda*atan(proton->kf/Lambda));
      ham3=(Cl*(neutron->n*gn+proton->n*gp)+
	    Cu*(neutron->n*gp+proton->n*gn))/rho0;
      
    }
    
    ham=ham1+ham2+ham3;
      
  } else {
      
    ham1=2.0/3.0*A/rho0*((1.0+0.5*x0)*n*n-(0.5+x0)*
			 (neutron->n*neutron->n+proton->n*proton->n));
    double term=((1.0+0.5*x3)*n*n-
		 (0.5+x3)*(neutron->n*neutron->n+proton->n*proton->n))*
      pow(n,sigma-1.0);
    ham2=4.0/3.0*B/pow(rho0,sigma)*term/
      (1.0+4.0/3.0*Bp/pow(rho0,sigma-1.0)*term/n/n);
    double u=n/n0;

    if (form==bgbd_form) {

      double gn, gp;
      gn=Lambda*Lambda/pi2*(neutron->kf-Lambda*atan(neutron->kf/Lambda));
      gp=Lambda*Lambda/pi2*(proton->kf-Lambda*atan(proton->kf/Lambda));
      ham3=0.8/rho0*(C1+2.0*z1)*n*(gn+gp)+0.4/rho0*(C1-8.0*z1)*
	(neutron->n*gn+proton->n*gp);

    } else if (form==bpal_form) {

      double gn1, gp1, gn2, gp2;
      gn1=Lambda*Lambda/pi2*(neutron->kf-Lambda*atan(neutron->kf/Lambda));
      gp1=Lambda*Lambda/pi2*(proton->kf-Lambda*atan(proton->kf/Lambda));
      gn2=Lambda2*Lambda2/pi2*(neutron->kf-
			       Lambda2*atan(neutron->kf/Lambda2));
      gp2=Lambda2*Lambda2/pi2*(proton->kf-
			       Lambda2*atan(proton->kf/Lambda2));
      ham3=0.8/rho0*(C1+2.0*z1)*n*(gn1+gp1)+0.4/rho0*(C1-8.0*z1)*
	(neutron->n*gn1+proton->n*gp1);
      ham3+=0.8/rho0*(C2+2.0*z2)*n*(gn2+gp2)+0.4/rho0*(C2-8.0*z2)*
	(neutron->n*gn2+proton->n*gp2);

    } else if (form==sl_form) {

      double gn1, gp1, gn2, gp2;
      gn1=neutron->n-pow(neutron->kf,5.0)/5.0/pi2/Lambda/Lambda;
      gp1=proton->n-pow(proton->kf,5.0)/5.0/pi2/Lambda/Lambda;
      gn2=neutron->n-pow(neutron->kf,5.0)/5.0/pi2/Lambda2/Lambda2;
      gp2=proton->n-pow(proton->kf,5.0)/5.0/pi2/Lambda2/Lambda2;
      ham3=0.8/rho0*(C1+2.0*z1)*n*(gn1+gp1)+0.4/rho0*(C1-8.0*z1)*
	(neutron->n*gn1+proton->n*gp1);
      ham3+=0.8/rho0*(C2+2.0*z2)*n*(gn2+gp2)+0.4/rho0*(C2-8.0*z2)*
	(neutron->n*gn2+proton->n*gp2);

    }

    ham=ham1+ham2+ham3;

  }

  return ham;
}

int eos_had_sym4_mdi::calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			 double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			 double &mu_n_pot, double &mu_p_pot) {
  set_n_and_p(ne,pr);
    
  double tmp;
  funct_mfptr<eos_had_sym4_mdi> dfk(this,&eos_had_sym4_mdi::energy_kin);
  funct_mfptr<eos_had_sym4_mdi> dfp(this,&eos_had_sym4_mdi::energy_pot);
    
  mode=nmode;
  tmp=ne.n;
  mu_n_kin=mu_deriv_ptr->deriv(ne.n,dfk);
  mu_n_pot=mu_deriv_ptr->deriv(ne.n,dfp);
  ne.n=tmp;
    
  mode=pmode;
  tmp=pr.n;
  mu_p_kin=mu_deriv_ptr->deriv(pr.n,dfk);
  mu_p_pot=mu_deriv_ptr->deriv(pr.n,dfp);
  pr.n=tmp;
    
  mode=normal;
  ed_kin=energy_kin(pr.n);
  ed_pot=energy_pot(pr.n);
    
  /// Run the normal version to get the effective masses right
  thermo lth;
  int ret=eos_had_potential::calc_e(ne,pr,lth);

  return ret;
}

int eos_had_sym4_mdi::test_separation(fermion &ne, fermion &pr, test_mgr &t) {
  double ed_kin, ed_pot, mu_n_kin, mu_p_kin, 
    mu_n_pot, mu_p_pot;
    
  set_n_and_p(ne,pr);
    
  double tmp;
  funct_mfptr<eos_had_sym4_mdi> dfk(this,&eos_had_sym4_mdi::energy_kin);
  funct_mfptr<eos_had_sym4_mdi> dfp(this,&eos_had_sym4_mdi::energy_pot);
  int vpx=0;
    
  mode=nmode;
  tmp=ne.n;
  mu_n_kin=mu_deriv_ptr->deriv(ne.n,dfk);
  mu_n_pot=mu_deriv_ptr->deriv(ne.n,dfp);
  ne.n=tmp;
    
  mode=pmode;
  tmp=pr.n;
  mu_p_kin=mu_deriv_ptr->deriv(pr.n,dfk);
  mu_p_pot=mu_deriv_ptr->deriv(pr.n,dfp);
  pr.n=tmp;
    
  mode=normal;
  ed_kin=energy_kin(pr.n);
  ed_pot=energy_pot(pr.n);
    
  /// Run the normal version to get the effective masses right
  thermo lth;
  eos_had_potential::calc_e(ne,pr,lth);

  cout << ne.mu << " " << mu_n_kin+mu_n_pot << endl;
  cout << pr.mu << " " << mu_p_kin+mu_p_pot << endl;
  cout << lth.ed << " " << ed_kin+ed_pot << endl;
  t.test_rel(ne.mu,mu_n_kin+mu_n_pot,1.0e-5,"ts1");
  t.test_rel(pr.mu,mu_p_kin+mu_p_pot,1.0e-5,"ts2");
  t.test_rel(lth.ed,ed_kin+ed_pot,1.0e-5,"ts3");

  return 0;
}

int eos_had_sym4::set_base_eos(eos_had_sym4_base &seb) {
  sp=&seb;
  return 0;
}
  
int eos_had_sym4::test_eos(fermion &ne, fermion &pr, thermo &lth) {
  double nn=ne.n, np=pr.n;
  double eden, pres, mun, mup;
    
  double ed_kin, ed_pot, mu_n_kin, mu_p_kin, 
    mu_n_pot, mu_p_pot;
    
  sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		 mu_n_pot,mu_p_pot);

  eden=ed_kin;
  mun=mu_n_kin;
  mup=mu_p_kin;

  double edpot0=ed_pot;
  double munpot0=mu_n_pot;
  double muppot0=mu_p_pot;
  ne.n+=1.0e-4;
  sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		 mu_n_pot,mu_p_pot);
  cout << (ed_kin-eden)/1.0e-4 << " " << mun << endl;
  cout << (ed_pot-edpot0)/1.0e-4 << " " << munpot0 << endl;
  ne.n-=1.0e-4;
  pr.n+=1.0e-4;
  sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
		 mu_n_pot,mu_p_pot);
  cout << (ed_kin-eden)/1.0e-4 << " " << mup << endl;
  cout << (ed_pot-edpot0)/1.0e-4 << " " << muppot0 << endl;

  return 0;
}

int eos_had_sym4::calc_e(fermion &ne, fermion &pr, thermo &lth) {
  double nn=ne.n, np=pr.n;
  double eden, mun, mup;

  double ed_kin, ed_pot, mu_n_kin, mu_p_kin, mu_n_pot, mu_p_pot;
  double kf_n, kf_p;

  /// Compute kinetic energy part at the specified proton fraction

  int r1=sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
			mu_n_pot,mu_p_pot);

  eden=ed_kin;
  mun=mu_n_kin;
  mup=mu_p_kin;
  kf_n=ne.kf;
  kf_p=pr.kf;

  /// Compute the properties of neutron matter

  ne.n=nn+np;
  pr.n=0.0;
  int r2=sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
			mu_n_pot,mu_p_pot);
  double ed_pot_neut=ed_pot;
  double n_pot_mu_neut=mu_n_pot;
  // We don't use the proton chemical potential here for neutron
  // matter. Not sure if it makes sense.

  /// Compute the properties of nuclear matter

  ne.n=(nn+np)/2.0;
  pr.n=(nn+np)/2.0;
  int r3=sp->calc_e_sep(ne,pr,ed_kin,ed_pot,mu_n_kin,mu_p_kin,
			mu_n_pot,mu_p_pot);
  double ed_pot_nuc=ed_pot;
  double n_pot_mu_nuc=mu_n_pot;
  double p_pot_mu_nuc=mu_p_pot;

  // Combine everything fixing the isospin dependence of the 
  // symmetry energy
    
  double xp=nn/(nn+np);
  double xp2=xp*xp;
  double term=(-12.0*xp+76.0*xp2-128.0*xp2*xp+64.0*xp2*xp2+ 
	       28.0*xp*alpha/3.0-92.0*xp2*alpha/3.0+
	       128.0*xp2*xp*alpha/3.0-
	       64.0*xp2*xp2*alpha/3.0)/(1.0+alpha);
  double dtdx=(-12.0+152.0*xp-384.0*xp2+256.0*xp2*xp+ 
	       28.0*alpha/3.0-184.0*xp*alpha/3.0+
	       128.0*xp2*alpha-
	       256.0*xp2*xp*alpha/3.0)/(1.0+alpha);
    
  lth.ed=eden+ed_pot_nuc*term+ed_pot_neut*(1.0-term);
    
  ne.n=nn;
  pr.n=np;
    
  ne.mu=mun+n_pot_mu_nuc*term+n_pot_mu_neut*(1.0-term)+dtdx*
    (ed_pot_neut-ed_pot_nuc)*(xp-1.0)/(nn+np);
  pr.mu=mup+p_pot_mu_nuc*term+n_pot_mu_neut*(1.0-term)+dtdx*
    (ed_pot_neut-ed_pot_nuc)*xp/(nn+np);
    
  // Just use the thermodynamic identity for the pressure

  lth.pr=-lth.ed+ne.mu*ne.n+pr.mu*pr.n;

  // Reset the Fermi momenta

  ne.kf=kf_n;
  pr.kf=kf_p;
    
  if (r1!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }
  if (r2!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }
  if (r3!=0) {
    O2SCL_ERR_RET("Neutron matter failed in calc_e_alpha().",exc_efailed);
  }
  return 0;
}

