/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_had_rmf_hyp.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_rmf_hyp::eos_had_rmf_hyp() {
  xs=1.0;
  xw=1.0;
  xr=1.0;
  inc_cascade=true;

  def_lambda.init(mass_lambda_MeV_f<double>()/hc_mev_fm,2.0);
  def_sigma_p.init(mass_sigma_plus_MeV_f<double>()/hc_mev_fm,2.0);
  def_sigma_z.init(mass_sigma_zero_MeV_f<double>()/hc_mev_fm,2.0);
  def_sigma_m.init(mass_sigma_minus_MeV_f<double>()/hc_mev_fm,2.0);
  def_cascade_z.init(mass_cascade_zero_MeV_f<double>()/hc_mev_fm,2.0);
  def_cascade_m.init(mass_cascade_minus_MeV_f<double>()/hc_mev_fm,2.0);
  
  lambda=&def_lambda;
  sigma_p=&def_sigma_p;
  sigma_z=&def_sigma_z;
  sigma_m=&def_sigma_m;
  cascade_m=&def_cascade_m;
  cascade_z=&def_cascade_z;

  neutron->inc_rest_mass=true;
  neutron->non_interacting=false;
  proton->inc_rest_mass=true;
  proton->non_interacting=false;
  lambda->inc_rest_mass=true;
  lambda->non_interacting=false;
  sigma_p->inc_rest_mass=true;
  sigma_p->non_interacting=false;
  sigma_z->inc_rest_mass=true;
  sigma_z->non_interacting=false;
  sigma_m->inc_rest_mass=true;
  sigma_m->non_interacting=false;
  cascade_z->inc_rest_mass=true;
  cascade_z->non_interacting=false;
  cascade_m->inc_rest_mass=true;
  cascade_m->non_interacting=false;
}

int eos_had_rmf_hyp::calc_eq_hyp_p
(fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm, double sig, double ome, 
 double lrho, double &f1, double &f2, double &f3, thermo &lth) {
  
  double duds,us,fun,gs2,gr2,gw2,gs,gr,gw,sig2,ome2;
  double rho2,sig4,ome4,dfdw,fnn=0.0,fnp=0.0;

#if !O2SCL_NO_RANGE_CHECK
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10 ||
      fabs(lam.g-2.0)>1.0e-10 || fabs(sigp.g-2.0)>1.0e-10 ||
      fabs(sigz.g-2.0)>1.0e-10 || fabs(sigm.g-2.0)>1.0e-10 ||
      (inc_cascade==true &&
       (fabs(casz.g-2.0)>1.0e-10 || fabs(casm.g-2.0)>1.0e-10))) {
    O2SCL_ERR2("At least one of the baryon spin degeneracies ",
	       "is wrong in eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true ||
      lam.non_interacting==true || sigp.non_interacting==true ||
      sigz.non_interacting==true || sigm.non_interacting==true ||
      (inc_cascade==true && 
       (casz.non_interacting==true || casm.non_interacting==true))) {
    O2SCL_ERR2("At least one baryon is non-interacting in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
#endif
  
  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;

  double gss=xs*gs;
  double gss2=gss*gss;
  double gws=xw*gw;
  double gws2=gws*gws;
  double grs=xr*gr;
  double grs2=grs*grs;

  ne.ms=ne.m-gs*sig;
  pr.ms=pr.m-gs*sig;
  lam.ms=lam.m-gss*sig;
  sigp.ms=sigp.m-gss*sig;
  sigz.ms=sigz.m-gss*sig;
  sigm.ms=sigm.m-gss*sig;
  if (inc_cascade) {
    casm.ms=casm.m-gss*sig;
    casz.ms=casz.m-gss*sig;
  }

  if (ne.ms<0.0 || pr.ms<0.0 || lam.ms<0.0 || 
      sigp.ms<0.0 || sigz.ms<0.0 || sigm.ms<0.0 ||
      (inc_cascade && (casm.ms<0.0 || casz.ms<0.0))) {
    O2SCL_ERR2("Baryon mass negative in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",o2scl::exc_efailed);
  }
  
  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;
  lam.nu=lam.mu-gws*ome;
  sigp.nu=sigp.mu-gws*ome-grs*rho;
  sigz.nu=sigz.mu-gws*ome;
  sigm.nu=sigm.mu-gws*ome+grs*rho;
  if (inc_cascade) {
    casz.nu=casz.mu-gws*ome-0.5*grs*rho;
    casm.nu=casm.mu-gws*ome+0.5*grs*rho;
  }

  if (ne.nu<ne.ms || ne.ms<0.0) ne.kf=0.0;
  else ne.kf=sqrt(ne.nu*ne.nu-ne.ms*ne.ms);
  if (pr.nu<pr.ms || pr.ms<0.0) pr.kf=0.0;
  else pr.kf=sqrt(pr.nu*pr.nu-pr.ms*pr.ms);
  if (lam.nu<lam.ms || lam.ms<0.0) lam.kf=0.0;
  else lam.kf=sqrt(lam.nu*lam.nu-lam.ms*lam.ms);
  if (sigp.nu<sigp.ms || sigp.ms<0.0) sigp.kf=0.0;
  else sigp.kf=sqrt(sigp.nu*sigp.nu-sigp.ms*sigp.ms);
  if (sigz.nu<sigz.ms || sigz.ms<0.0) sigz.kf=0.0;
  else sigz.kf=sqrt(sigz.nu*sigz.nu-sigz.ms*sigz.ms);
  if (sigm.nu<sigm.ms || sigm.ms<0.0) sigm.kf=0.0;
  else sigm.kf=sqrt(sigm.nu*sigm.nu-sigm.ms*sigm.ms);
  if (inc_cascade) {
    if (casz.nu<casz.ms || casz.ms<0.0) casz.kf=0.0;
    else casz.kf=sqrt(casz.nu*casz.nu-casz.ms*casz.ms);
    if (casm.nu<casm.ms || casm.ms<0.0) casm.kf=0.0;
    else casm.kf=sqrt(casm.nu*casm.nu-casm.ms*casm.ms);
  }
  
  // We don't record error values, since these functions usually
  // always succeed
  fet->calc_mu_zerot(ne);
  fet->calc_mu_zerot(pr);
  fet->calc_mu_zerot(lam);
  fet->calc_mu_zerot(sigp);
  fet->calc_mu_zerot(sigz);
  fet->calc_mu_zerot(sigm);
  if (inc_cascade) {
    fet->calc_mu_zerot(casz);
    fet->calc_mu_zerot(casm);
  }

  sig2=sig*sig;
  sig4=sig2*sig2;
  ome2=ome*ome;
  ome4=ome2*ome2;
  rho2=lrho*lrho;

  duds=b*ne.m*gs2*gs*sig2+c*gs2*gs2*sig2*sig;
  
  fun=a1*sig+a2*sig2+a3*sig2*sig+a4*sig4+
    a5*sig4*sig+a6*sig4*sig2+b1*ome2+b2*ome4+b3*ome4*ome2;
  dfdw=2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome;
  
  us=b/3.0*ne.m*gs2*gs*sig2*sig+c/4.0*gs2*gs2*sig4;
  
  double ns_n=0.0;
  double ns_p=0.0;
  double ns_lam=0.0;
  double ns_sigp=0.0;
  double ns_sigz=0.0;
  double ns_sigm=0.0;
  double ns_casz=0.0;
  double ns_casm=0.0;
  if (ne.inc_rest_mass) {
    ns_n=1.0/ne.ms*(ne.ed-3.0*ne.pr);
    ns_p=1.0/pr.ms*(pr.ed-3.0*pr.pr);
    ns_lam=1.0/lam.ms*(lam.ed-3.0*lam.pr);
    ns_sigp=1.0/sigp.ms*(sigp.ed-3.0*sigp.pr);
    ns_sigz=1.0/sigz.ms*(sigz.ed-3.0*sigz.pr);
    ns_sigm=1.0/sigm.ms*(sigm.ed-3.0*sigm.pr);
    ns_casz=1.0/casz.ms*(casz.ed-3.0*casz.pr);
    ns_casm=1.0/casm.ms*(casm.ed-3.0*casm.pr);
  } else {
    ns_n=1.0/ne.ms*(ne.ed+ne.n*ne.m-3.0*ne.pr);
    ns_p=1.0/pr.ms*(pr.ed+pr.n*pr.m-3.0*pr.pr);
    ns_lam=1.0/lam.ms*(lam.ed+lam.n*lam.m-3.0*lam.pr);
    ns_sigp=1.0/sigp.ms*(sigp.ed+sigp.n*sigp.m-3.0*sigp.pr);
    ns_sigz=1.0/sigz.ms*(sigz.ed+sigz.n*sigz.m-3.0*sigz.pr);
    ns_sigm=1.0/sigm.ms*(sigm.ed+sigm.n*sigm.m-3.0*sigm.pr);
    ns_casz=1.0/casz.ms*(casz.ed+casz.n*casz.m-3.0*casz.pr);
    ns_casm=1.0/casm.ms*(casm.ed+casm.n*casm.m-3.0*casm.pr);
  }

  if (inc_cascade) {
    
    f1=ms*ms*sig-gs*(ns_n+ns_p)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig)-
      gss*(ns_lam+ns_sigp+ns_sigz+ns_sigm+ns_casz+ns_casm);
    f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
      (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome)-
      gws*(lam.n+sigp.n+sigz.n+sigm.n+casz.n+casm.n);
    f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
      2.0*gr2*lrho*fun-grs*(sigp.n-sigm.n+0.5*(casz.n-casm.n));
    
    lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+gr2*rho2*fun+
      ne.pr+pr.pr+lam.pr+sigp.pr+sigz.pr+sigm.pr+casz.pr+casm.pr;
    
    // Isospins
    // Neutron -1/2
    // Proton +1/2
    // Lambda 0
    // Sigma- -1
    // Sigma0 0
    // Sigma+ +1
    // Xi- -1/2
    // Xi0 +1/2

    lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed+casz.ed+casm.ed+
      zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
      gr2*rho2*(fun+ome*dfdw);
          
  } else {

    f1=ms*ms*sig-gs*(ns_n+ns_p)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig)-
      gss*(ns_lam+ns_sigp+ns_sigz+ns_sigm);
    f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
      (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome)-
      gws*(lam.n+sigp.n+sigz.n+sigm.n);
    f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
      2.0*gr2*lrho*fun-grs*0.5*(sigp.n-sigm.n);

    lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+gr2*rho2*fun+
      ne.pr+pr.pr+lam.pr+sigp.pr+sigz.pr+sigm.pr;
    
    // Isospins
    // Neutron -1/2
    // Proton +1/2
    // Lambda 0
    // Sigma- -1
    // Sigma0 0
    // Sigma+ +1
    // Xi- -1/2
    // Xi0 +1/2

    lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed+
      zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
      gr2*rho2*(fun+ome*dfdw);
    
  }

  return 0;
}

int eos_had_rmf_hyp::calc_eq_temp_hyp_p
(fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm, double T, double sig,
 double ome, double lrho, double &f1, double &f2, double &f3,
 thermo &lth) {
 
  
  double duds,us,fun,gs2,gr2,gw2,gs,gr,gw,sig2,ome2;
  double rho2,sig4,ome4,dfdw,fnn=0.0,fnp=0.0;

#if !O2SCL_NO_RANGE_CHECK
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10 ||
      fabs(lam.g-2.0)>1.0e-10 || fabs(sigp.g-2.0)>1.0e-10 ||
      fabs(sigz.g-2.0)>1.0e-10 || fabs(sigm.g-2.0)>1.0e-10 ||
      (inc_cascade==true &&
       (fabs(casz.g-2.0)>1.0e-10 || fabs(casm.g-2.0)>1.0e-10))) {
    O2SCL_ERR2("At least one of the baryon spin degeneracies ",
	       "is wrong in eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true ||
      lam.non_interacting==true || sigp.non_interacting==true ||
      sigz.non_interacting==true || sigm.non_interacting==true ||
      (inc_cascade==true && 
       (casz.non_interacting==true || casm.non_interacting==true))) {
    O2SCL_ERR2("At least one baryon is non-interacting in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",exc_einval);
  }
#endif
  
  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;

  double gss=xs*gs;
  double gss2=gss*gss;
  double gws=xw*gw;
  double gws2=gws*gws;
  double grs=xr*gr;
  double grs2=grs*grs;

  ne.ms=ne.m-gs*sig;
  pr.ms=pr.m-gs*sig;
  lam.ms=lam.m-gss*sig;
  sigp.ms=sigp.m-gss*sig;
  sigz.ms=sigz.m-gss*sig;
  sigm.ms=sigm.m-gss*sig;
  if (inc_cascade) {
    casm.ms=casm.m-gss*sig;
    casz.ms=casz.m-gss*sig;
  }

  if (ne.ms<0.0 || pr.ms<0.0 || lam.ms<0.0 || 
      sigp.ms<0.0 || sigz.ms<0.0 || sigm.ms<0.0 ||
      (inc_cascade && (casm.ms<0.0 || casz.ms<0.0))) {
    O2SCL_ERR2("Baryon mass negative in ",
	       "eos_had_rmf_hyp::calc_eq_hyp_p().",o2scl::exc_efailed);
  }
  
  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;
  lam.nu=lam.mu-gws*ome;
  sigp.nu=sigp.mu-gws*ome-grs*rho;
  sigz.nu=sigz.mu-gws*ome;
  sigm.nu=sigm.mu-gws*ome+grs*rho;
  if (inc_cascade) {
    casz.nu=casz.mu-gws*ome-0.5*grs*rho;
    casm.nu=casm.mu-gws*ome+0.5*grs*rho;
  }

  if (ne.nu<ne.ms || ne.ms<0.0) ne.kf=0.0;
  else ne.kf=sqrt(ne.nu*ne.nu-ne.ms*ne.ms);
  if (pr.nu<pr.ms || pr.ms<0.0) pr.kf=0.0;
  else pr.kf=sqrt(pr.nu*pr.nu-pr.ms*pr.ms);
  if (lam.nu<lam.ms || lam.ms<0.0) lam.kf=0.0;
  else lam.kf=sqrt(lam.nu*lam.nu-lam.ms*lam.ms);
  if (sigp.nu<sigp.ms || sigp.ms<0.0) sigp.kf=0.0;
  else sigp.kf=sqrt(sigp.nu*sigp.nu-sigp.ms*sigp.ms);
  if (sigz.nu<sigz.ms || sigz.ms<0.0) sigz.kf=0.0;
  else sigz.kf=sqrt(sigz.nu*sigz.nu-sigz.ms*sigz.ms);
  if (sigm.nu<sigm.ms || sigm.ms<0.0) sigm.kf=0.0;
  else sigm.kf=sqrt(sigm.nu*sigm.nu-sigm.ms*sigm.ms);
  if (inc_cascade) {
    if (casz.nu<casz.ms || casz.ms<0.0) casz.kf=0.0;
    else casz.kf=sqrt(casz.nu*casz.nu-casz.ms*casz.ms);
    if (casm.nu<casm.ms || casm.ms<0.0) casm.kf=0.0;
    else casm.kf=sqrt(casm.nu*casm.nu-casm.ms*casm.ms);
  }
  
  // We don't record error values, since these functions usually
  // always succeed
  fet->pair_mu(ne,T);
  fet->pair_mu(pr,T);
  fet->pair_mu(lam,T);
  fet->pair_mu(sigp,T);
  fet->pair_mu(sigz,T);
  fet->pair_mu(sigm,T);
  if (inc_cascade) {
    fet->pair_mu(casz,T);
    fet->pair_mu(casm,T);
  }

  sig2=sig*sig;
  sig4=sig2*sig2;
  ome2=ome*ome;
  ome4=ome2*ome2;
  rho2=lrho*lrho;

  duds=b*ne.m*gs2*gs*sig2+c*gs2*gs2*sig2*sig;
  
  fun=a1*sig+a2*sig2+a3*sig2*sig+a4*sig4+
    a5*sig4*sig+a6*sig4*sig2+b1*ome2+b2*ome4+b3*ome4*ome2;
  dfdw=2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome;
  
  us=b/3.0*ne.m*gs2*gs*sig2*sig+c/4.0*gs2*gs2*sig4;
  
  double ns_n=0.0;
  double ns_p=0.0;
  double ns_lam=0.0;
  double ns_sigp=0.0;
  double ns_sigz=0.0;
  double ns_sigm=0.0;
  double ns_casz=0.0;
  double ns_casm=0.0;
  if (ne.inc_rest_mass) {
    ns_n=1.0/ne.ms*(ne.ed-3.0*ne.pr);
  } else {
    ns_n=1.0/ne.ms*(ne.ed+ne.n*ne.m-3.0*ne.pr);
  }
  if (pr.inc_rest_mass) {
    ns_p=1.0/pr.ms*(pr.ed-3.0*pr.pr);
  } else {
    ns_p=1.0/pr.ms*(pr.ed+pr.n*pr.m-3.0*pr.pr);
  }
  if (lam.inc_rest_mass) {
    ns_lam=1.0/lam.ms*(lam.ed-3.0*lam.pr);
  } else {
    ns_lam=1.0/lam.ms*(lam.ed+lam.n*lam.m-3.0*lam.pr);
  }
  if (sigp.inc_rest_mass) {
    ns_sigp=1.0/sigp.ms*(sigp.ed-3.0*sigp.pr);
  } else {
    ns_sigp=1.0/sigp.ms*(sigp.ed+sigp.n*sigp.m-3.0*sigp.pr);
  }
  if (sigz.inc_rest_mass) {
    ns_sigz=1.0/sigz.ms*(sigz.ed-3.0*sigz.pr);
  } else {
    ns_sigz=1.0/sigz.ms*(sigz.ed+sigz.n*sigz.m-3.0*sigz.pr);
  }
  if (sigm.inc_rest_mass) {
    ns_sigm=1.0/sigm.ms*(sigm.ed-3.0*sigm.pr);
  } else {
    ns_sigm=1.0/sigm.ms*(sigm.ed+sigm.n*sigm.m-3.0*sigm.pr);
  }
  if (inc_cascade) {
    if (casz.inc_rest_mass) {
      ns_casz=1.0/casz.ms*(casz.ed-3.0*casz.pr);
    } else {
    ns_casz=1.0/casz.ms*(casz.ed+casz.n*casz.m-3.0*casz.pr);
    }
    if (casm.inc_rest_mass) {
      ns_casm=1.0/casm.ms*(casm.ed-3.0*casm.pr);
    } else {
      ns_casm=1.0/casm.ms*(casm.ed+casm.n*casm.m-3.0*casm.pr);
    }
  }

  if (inc_cascade) {
    
    f1=ms*ms*sig-gs*(ns_n+ns_p)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig)-
      gss*(ns_lam+ns_sigp+ns_sigz+ns_sigm+ns_casz+ns_casm);
    f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
      (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome)-
      gws*(lam.n+sigp.n+sigz.n+sigm.n+casz.n+casm.n);
    f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
      2.0*gr2*lrho*fun-grs*(sigp.n-sigm.n+0.5*(casz.n-casm.n));
    
    lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+gr2*rho2*fun+
      ne.pr+pr.pr+lam.pr+sigp.pr+sigz.pr+sigm.pr+casz.pr+casm.pr;
    
    // Isospins
    // Neutron -1/2
    // Proton +1/2
    // Lambda 0
    // Sigma- -1
    // Sigma0 0
    // Sigma+ +1
    // Xi- -1/2
    // Xi0 +1/2

    lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed+casz.ed+casm.ed+
      zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
      gr2*rho2*(fun+ome*dfdw);
          
  } else {

    f1=ms*ms*sig-gs*(ns_n+ns_p)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig)-
      gss*(ns_lam+ns_sigp+ns_sigz+ns_sigm);
    f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
      (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome)-
      gws*(lam.n+sigp.n+sigz.n+sigm.n);
    f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
      2.0*gr2*lrho*fun-grs*0.5*(sigp.n-sigm.n);

    lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+gr2*rho2*fun+
      ne.pr+pr.pr+lam.pr+sigp.pr+sigz.pr+sigm.pr;
    
    // Isospins
    // Neutron -1/2
    // Proton +1/2
    // Lambda 0
    // Sigma- -1
    // Sigma0 0
    // Sigma+ +1
    // Xi- -1/2
    // Xi0 +1/2

    lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed+
      zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
      gr2*rho2*(fun+ome*dfdw);
    
  }

  return 0;
}

void eos_had_rmf_hyp::calc_xs(double lam_be) {

  // Compute the fields at saturation
  saturation();

  // Now compute the proper value of xs
  double gs=ms*cs;
  double gw=mw*cw;
  xs=(-lam_be+xw*gw*omega)/gs/sigma;
      
  return;
}

void eos_had_rmf_hyp::calc_xw(double lam_be) {

  // Compute the fields at saturation
  saturation();

  // Now compute the proper value of xw
  double gs=ms*cs;
  double gw=mw*cw;
  xw=(xs*gs*sigma+lam_be)/gw/omega;
      
  return;
}

int eos_had_rmf_hyp::calc_e_solve_fun(size_t nv, const ubvector &ex, 
				      ubvector &ey) {
  
  double f1,f2,f3,sig,ome,lrho;

  neutron->mu=ex[0];
  proton->mu=ex[1];
  lambda->mu=neutron->mu;
  sigma_p->mu=proton->mu;
  sigma_z->mu=neutron->mu;
  sigma_m->mu=2.0*neutron->mu-proton->mu;
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];
  if (inc_cascade) {
    cascade_z->mu=neutron->mu;
    cascade_m->mu=2.0*neutron->mu-proton->mu;
  }
  
  calc_eq_hyp_p(*neutron,*proton,*lambda,*sigma_p,*sigma_z,
                *sigma_m,*cascade_z,*cascade_m,sig,ome,lrho,f1,f2,f3,
                *eos_thermo);
  
  // 11/5/08 - We don't want to call the error handler here, because
  // sometimes the solver may accidentally find a region where 
  // nu<ms, and can handle it automatically
  if (!ce_prot_matter && neutron->nu<neutron->ms) {
    return 1;
  }
  if (!ce_neut_matter && proton->nu<proton->ms) {
    return 2;
  }

  if (ce_neut_matter) {
    ey[0]=proton->nu-proton->ms;
    ey[1]=neutron->n-n_baryon;
  } else if (ce_prot_matter) {
    ey[0]=neutron->nu-neutron->ms;
    ey[1]=proton->n-n_baryon;
  } else {
    if (calc_e_relative) {
      ey[0]=(proton->n+neutron->n-n_baryon)/n_baryon;
      ey[1]=proton->n-n_charge;
    } else {
      ey[0]=proton->n+neutron->n-n_baryon;
      ey[1]=proton->n-n_charge;
    }
  }
  ey[2]=f1;
  ey[3]=f2;
  ey[4]=f3;

  for(int i=0;i<5;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      // 07/12/11 - We don't want to call the error handler here, because
      // sometimes the solver may be able to handle it automatically
      return 3;
    }
  }

  return 0;
}

int eos_had_rmf_hyp::calc_e_nobeta_fun(size_t nv, const ubvector &ex, 
                                       ubvector &ey, double nB,
                                       double Ye, double Ys) {
  
  double muu=ex[0];
  double mud=ex[1];
  double mus=ex[2];
  
  double sig=ex[3];
  double ome=ex[4];
  double lrho=ex[5];

  //cout << "1: " << ex[0] << " " << ex[1] << " " << ex[2] << " "
  //<< ex[3] << " " << ex[4] << " " << ex[5] << endl;

  neutron->mu=muu+mud+mud;
  proton->mu=muu+muu+mud;
  lambda->mu=muu+mud+mus;
  sigma_p->mu=muu+muu+mus;
  sigma_z->mu=muu+mud+mus;
  sigma_m->mu=mud+mud+mus;
  
  if (inc_cascade) {
    cascade_z->mu=muu+mus+mus;
    cascade_m->mu=mud+mus+mus;
  }

  double f1, f2, f3;
  calc_eq_hyp_p(*neutron,*proton,*lambda,*sigma_p,*sigma_z,
                *sigma_m,*cascade_z,*cascade_m,sig,ome,lrho,f1,f2,f3,
                *eos_thermo);

  double nB2=neutron->n+proton->n+lambda->n+sigma_p->n+sigma_z->n+
    sigma_m->n;
  if (inc_cascade) nB2+=cascade_z->n+cascade_m->n;

  double Ye2=proton->n+sigma_p->n-sigma_m->n;
  if (inc_cascade) Ye2-=cascade_m->n;
  Ye2/=nB2;
    
  double Ys2=lambda->n+sigma_p->n+sigma_z->n+sigma_m->n;
  if (inc_cascade) {
    Ys2+=2.0*(cascade_m->n+cascade_m->n);
  }
  Ys2/=nB2;

  ey[0]=(nB-nB2)/nB;
  ey[1]=Ye-Ye2;
  ey[2]=Ys-Ys2;
  ey[3]=f1;
  ey[4]=f2;
  ey[5]=f3;

  //cout << nB2 << " " << Ye2 << " " << Ys2 << endl;
  //cout << "2: " << ex[0] << " " << ex[1] << " " << ex[2] << " "
  //<< ex[3] << " " << ex[4] << " " << ex[5] << endl;
  
  for(int i=0;i<6;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      // 07/12/11 - We don't want to call the error handler here, because
      // sometimes the solver may be able to handle it automatically
      return 3;
    }
  }

  return 0;
}

int eos_had_rmf_hyp::calc_e_nobeta_np_fun(size_t nv, const ubvector &ex, 
                                          ubvector &ey, double nn,
                                          double np, double Ys) {
  
  double muu=ex[0];
  double mud=ex[1];
  double mus=ex[2];
  
  double sig=ex[3];
  double ome=ex[4];
  double lrho=ex[5];

  //cout << "1: " << ex[0] << " " << ex[1] << " " << ex[2] << " "
  //<< ex[3] << " " << ex[4] << " " << ex[5] << endl;

  neutron->mu=muu+mud+mud;
  proton->mu=muu+muu+mud;
  lambda->mu=muu+mud+mus;
  sigma_p->mu=muu+muu+mus;
  sigma_z->mu=muu+mud+mus;
  sigma_m->mu=mud+mud+mus;
  
  if (inc_cascade) {
    cascade_z->mu=muu+mus+mus;
    cascade_m->mu=mud+mus+mus;
  }

  double f1, f2, f3;
  calc_eq_hyp_p(*neutron,*proton,*lambda,*sigma_p,*sigma_z,
                *sigma_m,*cascade_z,*cascade_m,sig,ome,lrho,f1,f2,f3,
                *eos_thermo);
  
  double nB2=neutron->n+proton->n+lambda->n+sigma_p->n+sigma_z->n+
    sigma_m->n;
  if (inc_cascade) nB2+=cascade_z->n+cascade_m->n;

  double Ys2=lambda->n+sigma_p->n+sigma_z->n+sigma_m->n;
  if (inc_cascade) {
    Ys2+=2.0*(cascade_m->n+cascade_m->n);
  }
  Ys2/=nB2;

  ey[0]=(neutron->n-nn)/nn;
  ey[1]=(proton->n-np)/np;
  ey[2]=Ys-Ys2;
  ey[3]=f1;
  ey[4]=f2;
  ey[5]=f3;

  //cout << nB2 << " " << Ye2 << " " << Ys2 << endl;
  //cout << "2: " << ex[0] << " " << ex[1] << " " << ex[2] << " "
  //<< ex[3] << " " << ex[4] << " " << ex[5] << endl;
  
  for(int i=0;i<6;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      // 07/12/11 - We don't want to call the error handler here, because
      // sometimes the solver may be able to handle it automatically
      return 3;
    }
  }

  return 0;
}

int eos_had_rmf_hyp::calc_hyp_e(double n_baryon_loc, double n_charge_loc,
                                fermion &ne, fermion &pr,
				fermion &lam, fermion &sigp, fermion &sigz, 
				fermion &sigm, fermion &casz, fermion &casm,
				thermo &lth) {
  size_t nv=5;

  ubvector x(nv), y(nv);
  int ret;

  ne.non_interacting=false;
  pr.non_interacting=false;
  
  set_thermo(lth);
  set_n_and_p(ne,pr);
  lambda=&lam;
  sigma_p=&sigp;
  sigma_z=&sigz;
  sigma_m=&sigm;
  cascade_z=&casz;
  cascade_m=&casm;

  lambda->non_interacting=false;
  sigma_p->non_interacting=false;
  sigma_z->non_interacting=false;
  sigma_m->non_interacting=false;
  cascade_z->non_interacting=false;
  cascade_m->non_interacting=false;

  // If zero-density, then just return rest mass energy
  // Otherwise, set whether we are in neutron or proton matter mode
  if (ne.n<=0.0 && pr.n<=0.0) {
    ne.mu=ne.m;
    ne.ms=ne.m;
    ne.nu=ne.m;
    ne.ed=0.0;
    ne.pr=0.0;
    ne.en=0.0;
    pr.mu=pr.m;
    pr.ms=pr.m;
    pr.nu=pr.m;
    pr.ed=0.0;
    pr.pr=0.0;
    pr.en=0.0;
    lambda->mu=lambda->m;
    lambda->ms=lambda->m;
    lambda->nu=lambda->m;
    lambda->ed=0.0;
    lambda->pr=0.0;
    lambda->en=0.0;
    sigma_p->mu=sigma_p->m;
    sigma_p->ms=sigma_p->m;
    sigma_p->nu=sigma_p->m;
    sigma_p->ed=0.0;
    sigma_p->pr=0.0;
    sigma_p->en=0.0;
    sigma_z->mu=sigma_z->m;
    sigma_z->ms=sigma_z->m;
    sigma_z->nu=sigma_z->m;
    sigma_z->ed=0.0;
    sigma_z->pr=0.0;
    sigma_z->en=0.0;
    sigma_m->mu=sigma_m->m;
    sigma_m->ms=sigma_m->m;
    sigma_m->nu=sigma_m->m;
    sigma_m->ed=0.0;
    sigma_m->pr=0.0;
    sigma_m->en=0.0;
    cascade_z->mu=cascade_z->m;
    cascade_z->ms=cascade_z->m;
    cascade_z->nu=cascade_z->m;
    cascade_z->ed=0.0;
    cascade_z->pr=0.0;
    cascade_z->en=0.0;
    cascade_m->mu=cascade_m->m;
    cascade_m->ms=cascade_m->m;
    cascade_m->nu=cascade_m->m;
    cascade_m->ed=0.0;
    cascade_m->pr=0.0;
    cascade_m->en=0.0;
    return 0;
  } else if (ne.n<=0.0) {
    ne.n=0.0;
    ce_prot_matter=true;
  } else if (pr.n<=0.0) {
    pr.n=0.0;
    ce_neut_matter=true;
  } else {
    ce_neut_matter=false;
    ce_prot_matter=false;
  }
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf_hyp::calc_e_solve_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  if (guess_set) {
    
    // If an initial guess is given, then use it to directly compute
    // the EOS
    
    x[0]=ne.mu;
    x[1]=pr.mu;
    x[2]=sigma;
    x[3]=omega;
    x[4]=rho;
    guess_set=false;
    
    ret=eos_mroot->msolve(nv,x,fmf);
    
    int rt=calc_e_solve_fun(nv,x,y);
    if (rt!=0) {
      O2SCL_CONV2_RET("Final solution failed (user guess) in ",
		      "eos_had_rmf_hyp::calc_e().",exc_efailed,
		      this->err_nonconv);
    }
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    double nn=ne.n;
    double np=pr.n;
    
    x[0]=ne.m+0.05;
    x[1]=pr.m+0.01;
    x[2]=0.1;
    x[3]=0.07;
    x[4]=0.001;
    
    if (verbose>0) {
      cout << "Solving in eos_had_rmf_hyp::calc_e()." << endl;
      cout << "alpha      n_B        n_ch       mu_n       "
	   << "mu_p       sigma       omega      rho         ret" << endl;
      cout.precision(4);
    }

    for(double alpha=0.0;alpha<=1.0+1.0e-10;
	alpha+=1.0/((double)calc_e_steps)) {

      if (ce_prot_matter) {
	n_baryon_loc=0.12*(1.0-alpha)+np*alpha;
	n_charge_loc=n_baryon_loc;
      } else if (ce_neut_matter) {
	n_baryon_loc=0.12*(1.0-alpha)+nn*alpha;
	n_charge_loc=0.0;
      } else {
	n_baryon_loc=0.16*(1.0-alpha)+(nn+np)*alpha;
	n_charge_loc=0.08*(1.0-alpha)+np*alpha;
      }
    
      // 10/16/14: I think this was some previous debug code
      // if (fabs(alpha-0.1)<1.0e-8) {
      // x[0]*=1.0+1.0e-5;
      // x[4]=-1.0e-10;
      // }

      // If the chemical potentials are too small, shift them by
      // a little more than required to get positive densities. 
      int rt=calc_e_solve_fun(5,x,y);
      if (!ce_prot_matter && neutron->nu<neutron->ms) {
	neutron->mu+=(neutron->ms-neutron->mu)*1.01;
	rt=calc_e_solve_fun(5,x,y);
      }
      if (!ce_neut_matter && proton->nu<proton->ms) {
	proton->mu+=(proton->ms-proton->mu)*1.01;
	rt=calc_e_solve_fun(5,x,y);
      }

      // The initial point has n_n = n_p and thus rho=0, and the
      // solver has a little problem with the stepsize getting away
      // from the rho=0 point, so we give rho a small non-zero value
      if (fabs(x[4])<1.0e-8) {
	if (neutron->n>proton->n) {
	  x[4]=-1.0e-8;
	} else if (neutron->n<proton->n) {
	  x[4]=1.0e-8;
	}
      }
      
      // If the initial guess failed then we won't be able to solve
      if (rt!=0) {
	string s=((string)"Initial guess failed at (nn=")+
	  dtos(neutron->n)+" and np="+dtos(proton->n)+") in "+
	  "eos_had_rmf_hyp::calc_e().";
	O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
      }

      ret=eos_mroot->msolve(5,x,fmf);
      if (verbose>0.0) {
	cout << alpha << " " << n_baryon_loc << " " << n_charge_loc << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " " << x[4] << " " << ret << endl;
      }
    }
    if (verbose>0) {
      cout.precision(6);
      cout << endl;
    }
    
    int rt2=calc_e_solve_fun(5,x,y);
    if (rt2!=0) {
      O2SCL_CONV_RET("Final solution failed in eos_had_rmf_hyp::calc_e().",
		     exc_efailed,this->err_nonconv);
    }
    
  }

  sigma=x[2];
  omega=x[3];
  rho=x[4];

  // return neutron and proton densities to original values
  ne.n=n_baryon_loc-n_charge_loc;
  pr.n=n_charge_loc;
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf_hyp::calc_e",
		    "(fermion,fermion,thermo).",exc_efailed,this->err_nonconv);
  }

  return 0;
}

int eos_had_rmf_hyp::calc_temp_hyp_e
(double n_baryon_loc, double n_charge_loc,
 fermion &ne, fermion &pr,
 fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm,
 double T, thermo &lth) {
  
  size_t nv=5;

  ubvector x(nv), y(nv);
  int ret;

  ne.non_interacting=false;
  pr.non_interacting=false;
  
  set_thermo(lth);
  set_n_and_p(ne,pr);
  lambda=&lam;
  sigma_p=&sigp;
  sigma_z=&sigz;
  sigma_m=&sigm;
  cascade_z=&casz;
  cascade_m=&casm;

  lambda->non_interacting=false;
  sigma_p->non_interacting=false;
  sigma_z->non_interacting=false;
  sigma_m->non_interacting=false;
  cascade_z->non_interacting=false;
  cascade_m->non_interacting=false;

  mm_funct fmf;
  //=std::bind
  //(std::mem_fn<int(size_t,const ubvector &,ubvector &)>
  //(&eos_had_rmf_hyp::calc_temp_hyp_e_fun),
  //this,std::placeholders::_1,std::placeholders::_2,
  //std::placeholders::_3);

  if (guess_set) {
    
    // If an initial guess is given, then use it to directly compute
    // the EOS
    
    x[0]=ne.mu;
    x[1]=pr.mu;
    x[2]=sigma;
    x[3]=omega;
    x[4]=rho;
    guess_set=false;
    
    ret=eos_mroot->msolve(nv,x,fmf);
    
    int rt=calc_e_solve_fun(nv,x,y);
    if (rt!=0) {
      O2SCL_CONV2_RET("Final solution failed (user guess) in ",
		      "eos_had_rmf_hyp::calc_e().",exc_efailed,
		      this->err_nonconv);
    }
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    double nn=ne.n;
    double np=pr.n;
    
    x[0]=ne.m+0.05;
    x[1]=pr.m+0.01;
    x[2]=0.1;
    x[3]=0.07;
    x[4]=0.001;
    
    if (verbose>0) {
      cout << "Solving in eos_had_rmf_hyp::calc_temp_hyp_e()." << endl;
      cout << "alpha      n_B        n_ch       mu_n       "
	   << "mu_p       sigma       omega      rho         ret" << endl;
      cout.precision(4);
    }

    for(double alpha=0.0;alpha<=1.0+1.0e-10;
	alpha+=1.0/((double)calc_e_steps)) {

      if (ce_prot_matter) {
	n_baryon_loc=0.12*(1.0-alpha)+np*alpha;
	n_charge_loc=n_baryon_loc;
      } else if (ce_neut_matter) {
	n_baryon_loc=0.12*(1.0-alpha)+nn*alpha;
	n_charge_loc=0.0;
      } else {
	n_baryon_loc=0.16*(1.0-alpha)+(nn+np)*alpha;
	n_charge_loc=0.08*(1.0-alpha)+np*alpha;
      }
    
      // 10/16/14: I think this was some previous debug code
      // if (fabs(alpha-0.1)<1.0e-8) {
      // x[0]*=1.0+1.0e-5;
      // x[4]=-1.0e-10;
      // }

      // If the chemical potentials are too small, shift them by
      // a little more than required to get positive densities. 
      int rt=calc_e_solve_fun(5,x,y);
      if (!ce_prot_matter && neutron->nu<neutron->ms) {
	neutron->mu+=(neutron->ms-neutron->mu)*1.01;
	rt=calc_e_solve_fun(5,x,y);
      }
      if (!ce_neut_matter && proton->nu<proton->ms) {
	proton->mu+=(proton->ms-proton->mu)*1.01;
	rt=calc_e_solve_fun(5,x,y);
      }

      // The initial point has n_n = n_p and thus rho=0, and the
      // solver has a little problem with the stepsize getting away
      // from the rho=0 point, so we give rho a small non-zero value
      if (fabs(x[4])<1.0e-8) {
	if (neutron->n>proton->n) {
	  x[4]=-1.0e-8;
	} else if (neutron->n<proton->n) {
	  x[4]=1.0e-8;
	}
      }
      
      // If the initial guess failed then we won't be able to solve
      if (rt!=0) {
	string s=((string)"Initial guess failed at (nn=")+
	  dtos(neutron->n)+" and np="+dtos(proton->n)+") in "+
	  "eos_had_rmf_hyp::calc_e().";
	O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
      }

      ret=eos_mroot->msolve(5,x,fmf);
      if (verbose>0.0) {
	cout << alpha << " " << n_baryon_loc << " " << n_charge_loc << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " " << x[4] << " " << ret << endl;
      }
    }
    if (verbose>0) {
      cout.precision(6);
      cout << endl;
    }
    
    int rt2=calc_e_solve_fun(5,x,y);
    if (rt2!=0) {
      O2SCL_CONV_RET("Final solution failed in eos_had_rmf_hyp::calc_e().",
		     exc_efailed,this->err_nonconv);
    }
    
  }

  sigma=x[2];
  omega=x[3];
  rho=x[4];

  // return neutron and proton densities to original values
  ne.n=n_baryon_loc-n_charge_loc;
  pr.n=n_charge_loc;
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf_hyp::calc_e",
		    "(fermion,fermion,thermo).",exc_efailed,this->err_nonconv);
  }

  return 0;
}

int eos_had_rmf_hyp::calc_hyp_e_nobeta
(double nB, double Ye, double Ys, fermion &ne, fermion &pr,
 fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm, thermo &lth) {
  
  size_t nv=6;

  ubvector x(nv), y(nv);
  int ret;

  ne.non_interacting=false;
  pr.non_interacting=false;
  
  set_thermo(lth);
  set_n_and_p(ne,pr);
  lambda=&lam;
  sigma_p=&sigp;
  sigma_z=&sigz;
  sigma_m=&sigm;
  cascade_z=&casz;
  cascade_m=&casm;

  lambda->non_interacting=false;
  sigma_p->non_interacting=false;
  sigma_z->non_interacting=false;
  sigma_m->non_interacting=false;
  cascade_z->non_interacting=false;
  cascade_m->non_interacting=false;

  if (guess_set) {
    
    // If an initial guess is given, then use it to directly compute
    // the EOS
    
    x[0]=(2.0*proton->mu-neutron->mu)/3.0;
    x[1]=(2.0*neutron->mu-proton->mu)/3.0;
    x[2]=lambda->mu-x[0]-x[1];
    x[3]=sigma;
    x[4]=omega;
    x[5]=rho;
    guess_set=false;
    
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       double,double)>
       (&eos_had_rmf_hyp::calc_e_nobeta_fun),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,nB,Ye,Ys);
    
    ret=eos_mroot->msolve(nv,x,fmf);
    
    int rt=calc_e_solve_fun(nv,x,y);
    if (rt!=0) {
      O2SCL_CONV2_RET("Final solution failed (user guess) in ",
		      "eos_had_rmf_hyp::calc_hyp_e_nobeta().",exc_efailed,
		      this->err_nonconv);
    }
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    x[0]=(ne.m+0.05)/2.5;
    x[1]=(pr.m+0.01)/2.5;
    x[2]=x[1];
    x[3]=0.1;
    x[4]=0.07;
    x[5]=0.001;
    
    if (verbose>0) {
      cout << "Solving in eos_had_rmf_hyp::calc_e()." << endl;
      cout << "alpha      n_B        n_ch       mu_n       "
	   << "mu_p       sigma       omega      rho         ret" << endl;
      cout.precision(4);
    }

    for(double alpha=0.0;alpha<=1.0+1.0e-10;
	alpha+=1.0/((double)calc_e_steps)) {

      double nB2=0.60*(1.0-alpha)+nB*alpha;
      double Ye2=0.4*(1.0-alpha)+Ye*alpha;
      double Ys2=0.1*(1.0-alpha)+Ys*alpha;
    
      mm_funct fmf=std::bind
        (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                         double,double)>
         (&eos_had_rmf_hyp::calc_e_nobeta_fun),
         this,std::placeholders::_1,std::placeholders::_2,
         std::placeholders::_3,nB2,Ye2,Ys2);
      
      // If the chemical potentials are too small, shift them by
      // a little more than required to get positive densities. 
      int rt=calc_e_nobeta_fun(6,x,y,nB2,Ye2,Ys2);
      if (neutron->nu<neutron->ms) {
        cout << "adjusting." << endl;
	neutron->mu+=(neutron->ms-neutron->mu)*1.01;
	rt=calc_e_nobeta_fun(6,x,y,nB2,Ye2,Ys2);
      }
      if (proton->nu<proton->ms) {
        cout << "adjusting 2." << endl;
	proton->mu+=(proton->ms-proton->mu)*1.01;
	rt=calc_e_nobeta_fun(6,x,y,nB2,Ye2,Ys2);
      }

      // The initial point has n_n = n_p and thus rho=0, and the
      // solver has a little problem with the stepsize getting away
      // from the rho=0 point, so we give rho a small non-zero value
      if (fabs(x[5])<1.0e-8) {
	if (Ye2<0.5) {
          cout << "adjusting 3." << endl;
	  x[5]=-1.0e-8;
	} else {
          cout << "adjusting 4." << endl;
	  x[5]=1.0e-8;
	}
      }
      
      // If the initial guess failed then we won't be able to solve
      if (rt!=0) {
	string s=((string)"Initial guess failed at (nn=")+
	  dtos(neutron->n)+" and np="+dtos(proton->n)+") in "+
	  "eos_had_rmf_hyp::calc_e().";
	O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
      }

      ret=eos_mroot->msolve(6,x,fmf);
      if (verbose>0.0) {
	cout << alpha << " " << nB2 << " " << Ye2 << " " << Ys2 << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " " << x[4] << " " << x[5] << " " << ret << endl;
      }
    }
    if (verbose>0) {
      cout.precision(6);
      cout << endl;
    }
    
    int rt2=calc_e_nobeta_fun(6,x,y,nB,Ye,Ys);
    if (rt2!=0) {
      O2SCL_CONV_RET("Final solution failed in eos_had_rmf_hyp::calc_e().",
		     exc_efailed,this->err_nonconv);
    }
    
  }

  sigma=x[3];
  omega=x[4];
  rho=x[5];
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf_hyp::calc_e",
		    "(fermion,fermion,thermo).",exc_efailed,this->err_nonconv);
  }

  return 0;
}

int eos_had_rmf_hyp::calc_hyp_e_nobeta_np
(double Ys, fermion &ne, fermion &pr,
 fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm, thermo &lth) {

  double nn=ne.n;
  double np=pr.n;
  
  size_t nv=6;

  ubvector x(nv), y(nv);
  int ret;

  ne.non_interacting=false;
  pr.non_interacting=false;
  
  set_thermo(lth);
  set_n_and_p(ne,pr);
  lambda=&lam;
  sigma_p=&sigp;
  sigma_z=&sigz;
  sigma_m=&sigm;
  cascade_z=&casz;
  cascade_m=&casm;

  lambda->non_interacting=false;
  sigma_p->non_interacting=false;
  sigma_z->non_interacting=false;
  sigma_m->non_interacting=false;
  cascade_z->non_interacting=false;
  cascade_m->non_interacting=false;

  if (guess_set) {
    
    // If an initial guess is given, then use it to directly compute
    // the EOS
    
    x[0]=(2.0*proton->mu-neutron->mu)/3.0;
    x[1]=(2.0*neutron->mu-proton->mu)/3.0;
    x[2]=lambda->mu-x[0]-x[1];
    x[3]=sigma;
    x[4]=omega;
    x[5]=rho;
    guess_set=false;
    
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       double,double)>
       (&eos_had_rmf_hyp::calc_e_nobeta_np_fun),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,nn,np,Ys);
    
    ret=eos_mroot->msolve(nv,x,fmf);
    
    int rt=calc_e_solve_fun(nv,x,y);
    if (rt!=0) {
      O2SCL_CONV2_RET("Final solution failed (user guess) in ",
		      "eos_had_rmf_hyp::calc_hyp_e_nobeta_np().",exc_efailed,
		      this->err_nonconv);
    }
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    x[0]=(ne.m+0.05)/2.5;
    x[1]=(pr.m+0.01)/2.5;
    x[2]=x[1];
    x[3]=0.1;
    x[4]=0.07;
    x[5]=0.001;
    
    if (verbose>0) {
      cout << "Solving in eos_had_rmf_hyp::calc_e()." << endl;
      cout << "alpha      n_B        n_ch       mu_n       "
	   << "mu_p       sigma       omega      rho         ret" << endl;
      cout.precision(4);
    }

    for(double alpha=0.0;alpha<=1.0+1.0e-10;
	alpha+=1.0/((double)calc_e_steps)) {

      double nn2=0.3*(1.0-alpha)+nn*alpha;
      double np2=0.3*(1.0-alpha)+np*alpha;
      double Ys2=0.1*(1.0-alpha)+Ys*alpha;
    
      mm_funct fmf=std::bind
        (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                         double,double)>
         (&eos_had_rmf_hyp::calc_e_nobeta_np_fun),
         this,std::placeholders::_1,std::placeholders::_2,
         std::placeholders::_3,nn2,np2,Ys2);
      
      // If the chemical potentials are too small, shift them by
      // a little more than required to get positive densities. 
      int rt=calc_e_nobeta_fun(6,x,y,nn2,np2,Ys2);
      if (neutron->nu<neutron->ms) {
        cout << "adjusting." << endl;
	neutron->mu+=(neutron->ms-neutron->mu)*1.01;
	rt=calc_e_nobeta_fun(6,x,y,nn2,np2,Ys2);
      }
      if (proton->nu<proton->ms) {
        cout << "adjusting 2." << endl;
	proton->mu+=(proton->ms-proton->mu)*1.01;
	rt=calc_e_nobeta_fun(6,x,y,nn2,np2,Ys2);
      }

      // The initial point has n_n = n_p and thus rho=0, and the
      // solver has a little problem with the stepsize getting away
      // from the rho=0 point, so we give rho a small non-zero value
      if (fabs(x[5])<1.0e-8) {
	if (nn2-np2<0.5) {
          cout << "adjusting 3." << endl;
	  x[5]=-1.0e-8;
	} else {
          cout << "adjusting 4." << endl;
	  x[5]=1.0e-8;
	}
      }
      
      // If the initial guess failed then we won't be able to solve
      if (rt!=0) {
	string s=((string)"Initial guess failed at (nn=")+
	  dtos(neutron->n)+" and np="+dtos(proton->n)+") in "+
	  "eos_had_rmf_hyp::calc_e().";
	O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
      }

      ret=eos_mroot->msolve(6,x,fmf);
      if (verbose>0.0) {
	cout << alpha << " " << nn2 << " " << np2 << " " << Ys2 << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " " << x[4] << " " << x[5] << " " << ret << endl;
      }
    }
    if (verbose>0) {
      cout.precision(6);
      cout << endl;
    }
    
    int rt2=calc_e_nobeta_np_fun(6,x,y,nn,np,Ys);
    if (rt2!=0) {
      O2SCL_CONV_RET("Final solution failed in eos_had_rmf_hyp::calc_e().",
		     exc_efailed,this->err_nonconv);
    }
    
  }

  sigma=x[3];
  omega=x[4];
  rho=x[5];
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf_hyp::calc_e",
		    "(fermion,fermion,thermo).",exc_efailed,this->err_nonconv);
  }

  return 0;
}

void eos_had_rmf_hyp::set_hyp(fermion &lam, fermion &sigp, fermion &sigz, 
			      fermion &sigm, fermion &casz, fermion &casm) {
  lambda=&lam;
  sigma_p=&sigp;
  sigma_z=&sigz;
  sigma_m=&sigm;
  cascade_z=&casz;
  cascade_m=&casm;
  return;
}

int eos_had_rmf_hyp::beta_eq_T0(ubvector &nB_grid, ubvector &guess,
                                eos_leptons &elep,
				//fermion &e, bool include_muons,
				//fermion &mu, fermion_rel &frel,
				std::shared_ptr<table_units<> > results) {
  
  guess[0]=0.01;
  guess[1]=0.01;
  guess[2]=-0.001;
  guess[3]=5.0;
  guess[4]=5.0;

  double nB_temp;
      
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &, ubvector &, 
		     const double &, eos_leptons &)>
     (&eos_had_rmf_hyp::solve_beta_eq_T0),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::cref(nB_temp),std::ref(elep));
  //include_muons,std::ref(mu),std::ref(frel));
      
  results->clear();
  results->line_of_names(((std::string)"ed pr nb ne nmu nn np nlam ")+
			 "nsigp nsigz nsigm mun mup mulam musigp musigz "+
			 "musigm kfn kfp kflam kfsigp kfsigz kfsigm");
  results->line_of_units(((std::string)"1/fm^4 1/fm^4 1/fm^3 1/fm^3 ")+
			 "1/fm^3 1/fm^3 1/fm^3 1/fm^3 1/fm^3 "+
			 "1/fm^3 1/fm^3 1/fm 1/fm 1/fm 1/fm 1/fm 1/fm "+
			 "1/fm 1/fm 1/fm 1/fm 1/fm 1/fm");
  if (inc_cascade) {
    results->line_of_names("ncasz ncasm mucasz mucasm kfcasz kfcasm");
    results->set_unit("ncasz","1/fm^3");
    results->set_unit("ncasz","1/fm^3");
    results->set_unit("mucasz","1/fm");
    results->set_unit("mucasm","1/fm");
    results->set_unit("kfcasz","1/fm");
    results->set_unit("kfcasm","1/fm");
  }
      
  for(size_t i=0;i<nB_grid.size();i++) {
    nB_temp=nB_grid[i];

    beta_mroot.msolve(5,guess,fmf);
	
    // Final function evaluation to make sure, e.g.
    // eos_thermo object is correct
    ubvector y(5);
    fmf(5,guess,y);

    std::vector<double> line={eos_thermo->ed,eos_thermo->pr,nB_temp,
      elep.e.n,elep.mu.n,
      neutron->n,proton->n,lambda->n,
      sigma_p->n,sigma_z->n,sigma_m->n,
      neutron->mu,proton->mu,lambda->mu,
      sigma_p->mu,sigma_z->mu,sigma_m->mu,
      neutron->kf,proton->kf,lambda->kf,
      sigma_p->kf,sigma_z->kf,sigma_m->kf};
    results->line_of_data(line);
    if (inc_cascade) {
      size_t row=results->get_nlines()-1;
      results->set("ncasz",row,cascade_z->n);
      results->set("ncasm",row,cascade_z->n);
      results->set("mucasz",row,cascade_z->mu);
      results->set("mucasm",row,cascade_z->mu);
      results->set("kfcasz",row,cascade_z->kf);
      results->set("kfcasm",row,cascade_z->kf);
    }
	  
  }
      
  return 0;
}

int eos_had_rmf_hyp::solve_beta_eq_T0(size_t nv, const ubvector &x,
				      ubvector &y, const double &nB,
                                      eos_leptons &elep) {

  neutron->mu=x[3];
  proton->mu=x[4];
  lambda->mu=neutron->mu;
  sigma_p->mu=proton->mu;
  sigma_z->mu=neutron->mu;
  sigma_m->mu=2.0*neutron->mu-proton->mu;
  cascade_z->mu=neutron->mu;
  cascade_m->mu=2.0*neutron->mu-proton->mu;
  
  double f1, f2, f3;
  calc_eq_hyp_p(*neutron,*proton,*lambda,*sigma_p,*sigma_z,*sigma_m,
		*cascade_z,*cascade_m,x[0],x[1],x[2],f1,f2,f3,*eos_thermo);
  elep.e.mu=neutron->mu-proton->mu;
  elep.pair_mu(0.0);
  
  y[0]=proton->n+sigma_p->n-sigma_m->n-cascade_m->n-elep.e.n;
  if (elep.include_muons) {
    y[0]-=elep.mu.n;
  }
  y[1]=neutron->n+proton->n+lambda->n+sigma_p->n+sigma_z->n+sigma_m->n+
    cascade_z->n+cascade_m->n-nB;
  y[2]=f1;
  y[3]=f2;
  y[4]=f3;
  
  return 0;
}
