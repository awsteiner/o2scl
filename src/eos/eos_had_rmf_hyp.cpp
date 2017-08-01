/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2012, Andrew W. Steiner
  
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

#include <o2scl/eos_had_rmf_hyp.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_rmf_hyp::eos_had_rmf_hyp() {
  xs=1.0;
  xw=1.0;
  xr=1.0;
  inc_cascade=true;
}

int eos_had_rmf_hyp::calc_eq_p
(fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
 fermion &sigm, fermion &casz, fermion &casm, double sig, double ome, 
 double lrho, double &f1, double &f2, double &f3, thermo &lth) {
  
  double duds,us,fun,gs2,gr2,gw2,gs,gr,gw,sig2,ome2;
  double rho2,sig4,ome4,dfdw,fnn=0.0,fnp=0.0;

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
  sigz.ms=sigp.ms;
  sigm.ms=sigp.ms;
  if (inc_cascade) {
    casm.ms=casm.m-gss*sig;
    casz.ms=casm.ms;
  }

  if (ne.ms<0.0 || pr.ms<0.0 || lam.ms<0.0 || sigp.ms<0.0 || 
      (inc_cascade && casm.ms<0.0)) {
    O2SCL_ERR2("Hadron mass negative in ",
	       "eos_had_rmf_hyp::calc_eq_p().",o2scl::exc_efailed);
  }
  
  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;
  lam.nu=lam.mu-gws*ome;
  sigp.nu=sigp.mu-gws*ome-grs*rho;
  sigz.nu=sigz.nu-gws*ome;
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
    
    // The thermodynamic identity could be used instead of 
    // these explicit expressions, but it's nice to have them
    // available here.
    
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

    lth.ed=us+0.5*ms*ms*sig2-0.5*mw*mw*ome*ome-0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed+casz.ed+casm.ed-
      zeta/24.0*gw2*gw2*ome4-xi/24.0*gr2*gr2*rho2*rho2-fun*gr2*rho2+
      gw*ome*(ne.n+pr.n+lam.n+sigp.n+sigz.n+sigm.n+casz.n+casm.n)-
      0.5*gr*lrho*(ne.n-pr.n)-grs*lrho*(sigm.n-sigp.n-0.5*
					(casz.n-casm.n));
    
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

    // The thermodynamic identity could be used instead of 
    // these explicit expressions, but it's nice to have them
    // available here.
    
    lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
      zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+gr2*rho2*fun+
      ne.pr+pr.pr+lam.pr+sigp.pr+sigz.pr+sigm.pr;
    
    lth.ed=us+0.5*ms*ms*sig2-0.5*mw*mw*ome*ome-0.5*mr*mr*lrho*lrho+
      ne.ed+pr.ed+lam.ed+sigp.ed+sigz.ed+sigm.ed-
      zeta/24.0*gw2*gw2*ome4-xi/24.0*gr2*gr2*rho2*rho2-fun*gr2*rho2+
      gw*ome*(ne.n+pr.n+lam.n+sigp.n+sigz.n+sigm.n)-
      0.5*gr*lrho*(ne.n-pr.n)-grs*lrho*(sigm.n-sigp.n);
    
  }

  return 0;
}

void eos_had_rmf_hyp::calc_xs(double lam_be) {
      
  // Compute the saturation density
  saturation();

  // Nuclear matter
  def_neutron.n=n0/2.0;
  def_proton.n=n0/2.0;

  // Compute the fields at a fixed density
  thermo th;
  //calc_e_fields(def_neutron,def_proton,th,sigma,omega,rho);

  // Now compute the proper value of xs
  double gs=ms*cs;
  double gw=mw*cw;
  xs=(-lam_be+xw*gw*omega)/gs/sigma;
      
  return;
}

void eos_had_rmf_hyp::calc_xw(double lam_be) {
      
  // Compute the saturation density
  saturation();

  // Nuclear matter
  def_neutron.n=n0/2.0;
  def_proton.n=n0/2.0;

  // Compute the fields at a fixed density
  thermo th;
  //calc_e_fields(def_neutron,def_proton,th,sigma,omega,rho);

  // Now compute the proper value of xw
  double gs=ms*cs;
  double gw=mw*cw;
  xw=(xs*gs*sigma+lam_be)/gw/omega;
      
  return;
}
