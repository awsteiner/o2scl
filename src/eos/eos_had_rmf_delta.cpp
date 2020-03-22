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
#include <o2scl/eos_had_rmf_delta.h>

using namespace std;
using namespace o2scl;

int eos_had_rmf_delta::calc_e(fermion &ne, fermion &pr, thermo &lth) {
  ubvector x(6), y(6);
  int ret;
  
  set_thermo(lth);
  set_n_and_p(ne,pr);

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
  
  n_baryon=ne.n+pr.n;
  n_charge=pr.n;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf_delta::calc_e_solve_fun),
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
    x[5]=del;
    guess_set=false;

    ret=eos_mroot->msolve(6,x,fmf);
    calc_e_solve_fun(6,x,y);
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    double nn=ne.n;
    double np=pr.n;

    if (ce_neut_matter) {
      x[0]=4.9;
      x[1]=4.7;
      x[2]=0.06;
      x[3]=0.04;
      x[4]=-0.04;
      x[5]=0.01;
    } else if (ce_prot_matter) {
      x[0]=4.7;
      x[1]=4.9;
      x[2]=0.1;
      x[3]=0.07;
      x[4]=0.07;
      x[5]=0.01;
    } else {
      x[0]=4.8;
      x[1]=4.8;
      x[2]=0.1;
      x[3]=0.07;
      x[4]=-0.001;
      x[5]=0.01;
    }
    
    for(double alpha=0.0;alpha<=1.001;alpha+=0.1) {
      if (ce_prot_matter) {
	n_baryon=0.12*(1.0-alpha)+np*alpha;
	n_charge=n_baryon;
      } else if (ce_neut_matter) {
	n_baryon=0.12*(1.0-alpha)+nn*alpha;
	n_charge=0.0;
      } else {
	n_baryon=0.16*(1.0-alpha)+(nn+np)*alpha;
	n_charge=0.08*(1.0-alpha)+np*alpha;
      }
      
      ret=eos_mroot->msolve(6,x ,fmf);
    }
    calc_e_solve_fun(6,x,y );
    
  }
  
  sigma=x[2];
  omega=x[3];
  rho=x[4];
  del=x[5];
  
  if (ret!=0) {
    O2SCL_ERR2("msolve failed in eos_had_rmf_delta::calc_e",
	       "(fermion,fermion,thermo).",exc_efailed);
  }

  return 0;
}

int eos_had_rmf_delta::calc_temp_eqd_p
(fermion& ne, fermion& pr, double temper, double sig, double ome,
 double lrho, double delta, double &f1, double &f2, double &f3,
 double &f4, thermo& lth) {
  
  double gs2, gs, gr2, gr, gw2, gw, duds, gd, gd2;
  double us, fun, dfdw, fnn=0.0, fnp=0.0, sig2, sig4, ome2, ome4, rho2;
  
  if (temper<=0.0) {
    calc_eqd_p(ne,pr,sig,ome,lrho,delta,f1,f2,f3,f4,lth);
  }
  
  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gd=md*cd;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;
  gd2=gd*gd;
  
  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;
  
  ne.ms=ne.m-gs*sig+gd*delta;
  pr.ms=pr.m-gs*sig-gd*delta;

  if (ne.ms<=0.0) ne.ms=1.e-9;
  if (pr.ms<=0.0) pr.ms=1.e-9;
  
  fet->pair_mu(ne,temper);
  fet->pair_mu(pr,temper);
  
  sig2=sig*sig;
  sig4=sig2*sig2;
  ome2=ome*ome;
  ome4=ome2*ome2;
  rho2=lrho*lrho;

  duds=b*ne.m*gs2*gs*sig2+c*gs2*gs2*sig2*sig;

  fun=a1*sig+a2*sig2+a3*sig2*sig+a4*sig4+
    a5*sig4*sig+a6*sig4*sig2+b1*ome2+b2*ome4+b3*ome4*ome2;
  dfdw=2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome;

  // Scalar densities
  double nsn, nsp;
  if (ne.inc_rest_mass) {
    nsn=1.0/ne.ms*(ne.ed-3.0*ne.pr);
    nsp=1.0/pr.ms*(pr.ed-3.0*pr.pr);
  } else {
    nsn=1.0/ne.ms*(ne.ed+ne.n*ne.m-3.0*ne.pr);
    nsp=1.0/pr.ms*(pr.ed+pr.n*pr.m-3.0*pr.pr);
  }

  f1=ms*ms*sig-gs*(nsn+nsp)+duds-gr2*rho2*
    (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
     5.0*a5*sig4+6.0*a6*sig4*sig);
  f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
    (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome);
  f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
    2.0*gr2*lrho*fun;
  f4=md*md*delta-gd*(nsp-nsn);

  us=b/3.0*ne.m*gs2*gs*sig2*sig+c/4.0*gs2*gs2*sig4;
  
  lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
    zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+
    ne.pr+pr.pr+gr2*rho2*fun-0.5*md*md*delta*delta;
  lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
    ne.ed+pr.ed+zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
    gr2*rho2*(fun+ome*dfdw)+0.5*md*md*delta*delta;
  lth.en=(lth.ed+lth.pr-ne.n*ne.mu-pr.n*pr.mu)/temper;

  return success;
}

int eos_had_rmf_delta::calc_eqd_p(fermion &ne, fermion &pr, 
			  double sig, double ome, double rhof, double delta,
			  double &f1, double &f2, double &f3, double &f4,
			  thermo &th) {
  double duds, us, fun;
  double gs2,gr2,gw2,gs,gr,gw,gd,gd2,dfdw;
  char ch;

  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gd=md*cd;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;
  gd2=gd*gd;
  
  double sig2=sig*sig;
  double ome2=ome*ome;
  double ome4=ome2*ome2;
  double rho2=rhof*rhof;

  ne.ms=ne.m-gs*sig+gd*delta;
  pr.ms=pr.m-gs*sig-gd*delta;

  ne.nu=ne.mu-gw*ome+0.5*gr*rhof;
  pr.nu=pr.mu-gw*ome-0.5*gr*rhof;

  if (ne.nu<ne.ms || ne.ms<0.0) ne.kf=0.0;
  else ne.kf=sqrt(ne.nu*ne.nu-ne.ms*ne.ms);
  if (pr.nu<pr.ms || pr.ms<0.0) pr.kf=0.0;
  else pr.kf=sqrt(pr.nu*pr.nu-pr.ms*pr.ms);

  fet->calc_mu_zerot(ne);
  fet->calc_mu_zerot(pr);
  
  duds=b*ne.m*pow(gs,3.0)*pow(sig,2.0)+c*pow(gs,4.0)*pow(sig,3.0);

  fun=a1*sig+a2*sig*sig+a3*pow(sig,3.0)+a4*pow(sig,4.0)+
    a5*pow(sig,5.0)+a6*pow(sig,6.0)+b1*ome*ome+b2*pow(ome,4.0)+
    b3*pow(ome,6.0);
  dfdw=2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome;

  // Scalar densities
  double nsn, nsp;
  if (ne.inc_rest_mass) {
    nsn=1.0/ne.ms*(ne.ed-3.0*ne.pr);
    nsp=1.0/pr.ms*(pr.ed-3.0*pr.pr);
  } else {
    nsn=1.0/ne.ms*(ne.ed+ne.n*ne.m-3.0*ne.pr);
    nsp=1.0/pr.ms*(pr.ed+pr.n*pr.m-3.0*pr.pr);
  }

  f1=ms*ms*sig-gs*(nsn+nsp)+duds-gr2*rhof*rhof*
    (a1+2.0*a2*sig+3.0*a3*sig*sig+4.0*a4*pow(sig,3.0)+
     5.0*a5*pow(sig,4.0)+6.0*a6*pow(sig,5.0));
  f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*pow(ome,3.0)/6.0+gr2*rhof*rhof*
    (2.0*b1*ome+4.0*b2*pow(ome,3.0)+6.0*b3*pow(ome,5.0));
  f3=mr*mr*rhof-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*pow(rhof,3.0)/6.0+
    2.0*gr2*rhof*fun;
  f4=md*md*delta-gd*(nsp-nsn);

  us=b/3.0*ne.m*pow((gs*sig),3.0)+c/4.0*pow((gs*sig),4.0);
  
  th.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*rhof*rhof+
    zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+
    ne.pr+pr.pr+gr2*rho2*fun-0.5*md*md*delta*delta;
  th.ed=-th.pr+ne.mu*ne.n+pr.mu*pr.n;
  
  if (!std::isfinite(th.pr) || !std::isfinite(th.ed)) {
    O2SCL_ERR2("Pressure or energy not finite in ",
		   "eos_had_rmf_delta::calc_eqd_p().",exc_efailed);
  }
  
  return success;
}

int eos_had_rmf_delta::zero_pressure(size_t nv, const ubvector &ex, 
				 ubvector &ey) {
  double f1,f2,f3,f4,sig,ome,lrho,delta;
  fermion *n=neutron, *p=proton;
  int i;

  n->mu=ex[0];
  p->mu=ex[1];
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];
  delta=ex[5];

  for(i=0;i<6;i++) {
    if (ex[i]!=0.0 && !std::isfinite(ex[i])) {
      O2SCL_ERR("Variable not finite in zero_pressure()",exc_efailed);
    }
  }
  
  calc_eqd_p(*n,*p,sig,ome,lrho,delta,f1,f2,f3,f4,*eos_thermo);
  
  ey[0]=eos_thermo->pr;
  ey[1]=p->n/(n->n+p->n)-0.5;
  ey[2]=f1;
  ey[3]=f2;
  ey[4]=f3;
  ey[5]=f4;
  
  for(i=0;i<6;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      O2SCL_ERR((((string)"Eq. ")+itos(i)+
		   " not finite in zero_pressure()").c_str(),exc_efailed);
    }
  }
  return 0;
}

int eos_had_rmf_delta::calc_e_solve_fun(size_t nv, const ubvector &ex, 
				    ubvector &ey) {
  double f1,f2,f3,f4,sig,ome,lrho,delta;

  neutron->mu=ex[0];
  proton->mu=ex[1];
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];
  delta=ex[5];
  
  calc_eqd_p(*neutron,*proton,sig,ome,lrho,delta,f1,f2,f3,f4,*eos_thermo);
  
  if (!ce_prot_matter && neutron->nu<neutron->ms) {
    O2SCL_ERR("ne.nu<ne.ms in eos_had_rmf_delta::calc_e_solve_fun().",
	      exc_efailed);
  }
  if (!ce_neut_matter && proton->nu<proton->ms) {
    O2SCL_ERR("pr.nu<pr.ms in eos_had_rmf_delta::calc_e_solve_fun().",
	      exc_efailed);
  }

  if (ce_neut_matter) {
    ey[0]=proton->nu-proton->ms;
    ey[1]=neutron->n-n_baryon;
  } else if (ce_prot_matter) {
    ey[0]=neutron->nu-neutron->ms;
    ey[1]=proton->n-n_baryon;
  } else {
    ey[0]=proton->n+neutron->n-n_baryon;
    ey[1]=proton->n-n_charge;
  }
  ey[2]=f1;
  ey[3]=f2;
  ey[4]=f3;
  ey[5]=f4;

  for(int i=0;i<6;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      O2SCL_ERR((((string)"Eq. ")+itos(i)+
		 " not finite in eos_had_rmf_delta::"+
		 "calc_e_solve_fun().").c_str(),
		exc_efailed);
    }
  }

  return 0;
}

int eos_had_rmf_delta::saturation() {
  cout << "In saturation: " << endl;
  ubvector x(6);
  int test;

  x[0]=neutron->mu;
  x[1]=proton->mu;
  if (guess_set) {
    x[2]=sigma;
    x[3]=omega;
    x[4]=rho;
    x[5]=del;
  } else {
    x[2]=0.1;
    x[3]=0.07;
    x[4]=0.0;
    x[5]=0.0;
  }
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf_delta::zero_pressure),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  test=sat_mroot->msolve(6,x,fmf);
  
  sigma=x[2];
  omega=x[3];
  rho=x[4];
  del=x[5];

  if (test!=0) {
    O2SCL_CONV_RET("Solver failed in eos_had_rmf_delta::saturation().",
		   exc_efailed,this->err_nonconv);
  }
  
  n0=neutron->n+proton->n;
  msom=neutron->ms/neutron->m;
  eoa=(eos_thermo->ed/n0-neutron->m);
  esym=fesym(n0);
  
  return 0;
}

