/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

#include <o2scl/eos_had_rmf.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_rmf::eos_had_rmf() {

  mnuc=939.0/hc_mev_fm;
  ms=0.0;
  mr=0.0;
  mw=0.0;
  cs=0.0;
  cw=0.0;
  cr=0.0;
  b=0.0;
  c=0.0;

  zeta=0.0;
  xi=0.0;
  a1=0.0;
  a2=0.0;
  a3=0.0;
  a4=0.0;
  a5=0.0;
  a6=0.0;
  b1=0.0;
  b2=0.0;
  b3=0.0;

  zm_mode=false;
  guess_set=false;
  
  verbose=0;

  calc_e_relative=true;
  calc_e_steps=20;

  // AWS, 4/7/22: These settings have been moved from eos_had_rmf_ts
  // and are now the default, because they appear to improve results.
  // One problem was that the rho field is zero in nuclear matter and
  // so the step size has to be large enough to compute derivatives.
  def_mroot.def_jac.set_epsrel(1.0e-4);
  def_sat_mroot.def_jac.set_epsrel(1.0e-4);
  def_mroot.def_jac.set_epsmin(1.0e-15);
  
}

int eos_had_rmf::calc_eq_temp_p
(fermion &ne, fermion &pr, double temper, double sig, double ome, 
 double lrho, double &f1, double &f2, double &f3, thermo &lth) {

  double gs2, gs, gr2, gr, gw2, gw, duds;
  double us, fun, dfdw, fnn=0.0, fnp=0.0, sig2, sig4, ome2, ome4, rho2;

#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the chemical
  // potentials aren't finite, but I've found this extra checking of
  // the inputs useful for debugging.
  if (!std::isfinite(ne.mu) || !std::isfinite(ne.mu)) {
    O2SCL_ERR2("Chemical potentials not finite in ",
	       "eos_had_rmf::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_had_rmf::calc_eq_temp_p().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_rmf::calc_eq_temp_p().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true) {
    O2SCL_ERR2("Neutron or protons non-interacting in ",
	       "eos_had_rmf::calc_eq_temp_p().",exc_einval);
  }
#endif

  if (temper<=0.0) {
    calc_eq_p(ne,pr,sig,ome,rho,f1,f2,f3,lth);
  }
  
  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;
  
  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;
  
  if (zm_mode) {
    fnn=(1.0+gs*sig/ne.m);
    fnp=(1.0+gs*sig/pr.m);
    ne.ms=ne.m/fnn;
    pr.ms=pr.m/fnp;
  } else {
    ne.ms=ne.m-gs*sig;
    pr.ms=pr.m-gs*sig;
  }

  if (ne.ms<0.0 || pr.ms<0.0) {
    O2SCL_CONV2_RET("Neutron or proton mass negative in ",
		    "eos_had_rmf::calc_eq_temp_p().",exc_efailed,
		    this->err_nonconv);
  }
  
  ne.non_interacting=false;
  pr.non_interacting=false;

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
  } else {
    nsn=1.0/ne.ms*(ne.ed+ne.n*ne.m-3.0*ne.pr);
  }
  if (pr.inc_rest_mass) {
    nsp=1.0/pr.ms*(pr.ed-3.0*pr.pr);
  } else {
    nsp=1.0/pr.ms*(pr.ed+pr.n*pr.m-3.0*pr.pr);
  }

  if (zm_mode==true) {
    f1=ms*ms*sig-gs*(nsn/fnn/fnn+nsp/fnp/fnp)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig);
  } else {
    f1=ms*ms*sig-gs*(nsn+nsp)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig);
  }
  f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
    (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome);
  f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
    2.0*gr2*lrho*fun;

  us=b/3.0*ne.m*gs2*gs*sig2*sig+c/4.0*gs2*gs2*sig4;
  
  lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
    zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+
    ne.pr+pr.pr+gr2*rho2*fun;
  lth.ed=us+0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
    ne.ed+pr.ed+zeta/8.0*gw2*gw2*ome4+xi/8.0*gr2*gr2*rho2*rho2+
    gr2*rho2*(fun+ome*dfdw);
  lth.en=(lth.ed+lth.pr-ne.n*ne.mu-pr.n*pr.mu)/temper;

  return success;
}

int eos_had_rmf::check_derivs
(double &dPds, double &dPdw, double &dPdr, fermion &ne, fermion &pr,
 double sig, double ome, double lrho) {
  
  deriv_gsl<> dg;
  double feq1, feq2, feq3;

  funct f1=std::bind
    (std::mem_fn<int(fermion &,fermion &,double,
		     double,double,double &, double &, double &,
		     thermo &)>
     (&eos_had_rmf::calc_eq_p),this,std::ref(ne),std::ref(pr),
     std::placeholders::_1,ome,lrho,std::ref(feq1),std::ref(feq2),
     std::ref(feq3),std::ref(*eos_thermo));
  funct f2=std::bind
    (std::mem_fn<int(fermion &,fermion &,double,
		     double,double,double &, double &, double &,
		     thermo &)>
     (&eos_had_rmf::calc_eq_p),this,std::ref(ne),std::ref(pr),
     sig,std::placeholders::_1,lrho,std::ref(feq1),std::ref(feq2),
     std::ref(feq3),std::ref(*eos_thermo));
  funct f3=std::bind
    (std::mem_fn<int(fermion &,fermion &,double,
		     double,double,double &, double &, double &,
		     thermo &)>
     (&eos_had_rmf::calc_eq_p),this,std::ref(ne),std::ref(pr),
     sig,ome,std::placeholders::_1,std::ref(feq1),std::ref(feq2),
     std::ref(feq3),std::ref(*eos_thermo));

  double err;
  dg.deriv_err(sig,f1,dPds,err);
  dg.deriv_err(ome,f2,dPdw,err);
  dg.deriv_err(rho,f3,dPdr,err);
  
  return 0;
}

int eos_had_rmf::field_eqs(size_t nv, const ubvector &x, ubvector &y) {
#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the 
  // fields aren't finite, but I've found this extra checking of
  // the inputs useful for debugging.

  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) ||
      !std::isfinite(x[2])) {
    O2SCL_ERR("Fields not finite in eos_had_rmf::field_eqs().",
	      exc_efailed);
  }
#endif

  calc_eq_p(*neutron,*proton,x[0],x[1],x[2],y[0],y[1],y[2],
	    *eos_thermo);

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    return exc_ebadfunc;
  }
  return 0;
}

int eos_had_rmf::field_eqsT(size_t nv, const ubvector &x, ubvector &y) {
#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the 
  // fields aren't finite, but I've found this extra checking of
  // the inputs useful for debugging.
  if (!std::isfinite(x[0]) || !std::isfinite(x[1]) ||
      !std::isfinite(x[2])) {
    O2SCL_ERR("Fields not finite in eos_had_rmf::field_eqsT().",
	      exc_efailed);
  }
#endif

  double gs=ms*cs;
  if (zm_mode) {
    neutron->ms=neutron->m/(1.0+gs*x[0]/neutron->m);
    proton->ms=proton->m/(1.0+gs*x[0]/proton->m);
  } else {
    neutron->ms=neutron->m-gs*x[0];
    proton->ms=proton->m-gs*x[0];
  }
  if (neutron->ms<0.0 || proton->ms<0.0) return exc_ebadfunc;

  calc_eq_temp_p(*neutron,*proton,fe_temp,x[0],x[1],x[2],y[0],y[1],y[2],
		 *eos_thermo);

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    return exc_ebadfunc;
  }

  return 0;
}

int eos_had_rmf::calc_p(fermion &ne, fermion &pr, thermo &lth) {
  int ret;
  ubvector x(3);

  ne.non_interacting=false;
  pr.non_interacting=false;

  set_thermo(lth);
  set_n_and_p(ne,pr);
  
  if (guess_set) {
    x[0]=sigma;
    x[1]=omega;
    x[2]=rho;
    guess_set=false;
  } else {
    x[0]=0.3;
    x[1]=0.3;
    x[2]=-0.05;
  }
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf::field_eqs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  eos_mroot->verbose=2;
  ret=eos_mroot->msolve(3,x,fmf);
  if (ret!=0) {
    O2SCL_CONV_RET("Solver failed in eos_had_rmf::calc_p().",
		   exc_efailed,this->err_nonconv);
  }

  sigma=x[0];
  omega=x[1];
  rho=x[2];

  // 10/16/14: Final evaluation to store results in ne, pr, and lth
  // Use the x vector as a temporary. I'm not sure why this wasn't
  // here before, as it seems to be necessary in order to
  // properly report the correct pressure, energy density, etc.
  calc_eq_p(ne,pr,sigma,omega,rho,x[0],x[1],x[2],lth);
  
  return 0;
}

int eos_had_rmf::calc_temp_p(fermion &ne, fermion &pr, const double T,
			     thermo &lth) {
  int ret=0;
  ubvector x(3), y(3);

  if (!std::isfinite(ne.mu) || !std::isfinite(ne.mu)) {
    O2SCL_ERR2("Chemical potentials not finite in ",
	       "eos_had_rmf::calc_temp_p().",exc_efailed);
  }

  ne.non_interacting=false;
  pr.non_interacting=false;

  set_thermo(lth);
  set_n_and_p(ne,pr);
  fe_temp=T;
  
  if (guess_set) {
    x[0]=sigma;
    x[1]=omega;
    x[2]=rho;
    guess_set=false;
  } else {
    x[0]=0.3;
    x[1]=0.3;
    x[2]=-0.05;
  }
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf::field_eqsT),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  ret=eos_mroot->msolve(3,x,fmf);
  
  sigma=x[0];
  omega=x[1];
  rho=x[2];

  if (ret!=0) {
    O2SCL_CONV_RET("Solver failed in eos_had_rmf::calc_p().",
		   exc_efailed,this->err_nonconv);
  }

  // 10/16/14: Final evaluation to store results in ne, pr, and lth
  // Use the x vector as a temporary. I'm not sure why this wasn't
  // here before, as it seems to be necessary in order to
  // properly report the correct pressure, energy density, etc.
  calc_eq_temp_p(ne,pr,T,sigma,omega,rho,x[0],x[1],x[2],lth);
  
  return 0;
}

int eos_had_rmf::calc_e(fermion &ne, fermion &pr, thermo &lth) {
  ubvector x(5), y(5);
  int ret;

  ne.non_interacting=false;
  pr.non_interacting=false;

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
     (&eos_had_rmf::calc_e_solve_fun),
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

    ret=eos_mroot->msolve(5,x,fmf);

    int rt=calc_e_solve_fun(5,x,y);
    if (rt!=0) {
      O2SCL_CONV2_RET("Final solution failed (user guess) in ",
		      "eos_had_rmf::calc_e().",exc_efailed,this->err_nonconv);
    }
    
  } else {

    // If no initial guess is given, then create one by beginning
    // at saturated nuclear matter and proceeding incrementally.

    double nn=ne.n;
    double np=pr.n;
    
    x[0]=4.8;
    x[1]=4.8;
    x[2]=0.1;
    x[3]=0.07;
    x[4]=-0.001;

    if (verbose>0) {
      cout << "Solving in eos_had_rmf::calc_e(): Steps: "
           << calc_e_steps << endl;
      cout << "i  alpha      n_B        n_ch       mu_n       "
	   << "mu_p       sigma       omega      rho         ret" << endl;
    }

    int cnt=0;
    for(double alpha=0.0;alpha<=1.0+1.0e-10;
	alpha+=1.0/((double)calc_e_steps)) {

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

      // If the chemical potentials are too small, shift them to
      // get positive densities
      int rt=calc_e_solve_fun(5,x,y);
      
      if (!ce_prot_matter && neutron->nu<neutron->ms) {
	for(size_t j=0;j<5 && !ce_prot_matter && neutron->nu<neutron->ms;
	    j++) {
	  x[0]+=0.1;
	  rt=calc_e_solve_fun(5,x,y);
	}
      }
      
      if (!ce_neut_matter && proton->nu<proton->ms) {
	for(size_t j=0;j<5 && !ce_neut_matter && proton->nu<proton->ms;
	    j++) {
	  x[1]+=0.1;
	  rt=calc_e_solve_fun(5,x,y);
	}
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
	  "eos_had_rmf::calc_e().";
	O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
      }
      ret=eos_mroot->msolve(5,x,fmf);
      if (verbose>0) {
	cout.precision(4);
        cout.width(2);
	cout << cnt << " ";
        cout << alpha << " " << n_baryon << " " << n_charge << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " ";
        cout.setf(ios::showpos);
        cout << x[4] << " ";
        cout.unsetf(ios::showpos);
        cout << ret << endl;
	cout.precision(6);
      }
      cnt++;
    }
    if (verbose>0) {
      cout << endl;
    }
    
    int rt2=calc_e_solve_fun(5,x,y);
    if (rt2!=0) {
      O2SCL_CONV_RET("Final solution failed in eos_had_rmf::calc_e().",
		     exc_efailed,this->err_nonconv);
    }
    
  }

  sigma=x[2];
  omega=x[3];
  rho=x[4];

  // return neutron and proton densities to original values
  ne.n=n_baryon-n_charge;
  pr.n=n_charge;
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf::calc_e",
		    "(fermion,fermion,thermo).",exc_efailed,this->err_nonconv);
  }
  
  return 0;
}

int eos_had_rmf::calc_temp_e(fermion &ne, fermion &pr, const double T, 
			     thermo &lth) {
  
  if (T<=0.0) return calc_e(ne,pr,lth);

  ubvector x(5), y(5);
  
  ce_temp=T;

  ne.non_interacting=false;
  pr.non_interacting=false;

  set_thermo(lth);
  set_n_and_p(ne,pr);

  // If zero density, then just return rest mass energy
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
  }
  
  n_baryon=ne.n+pr.n;
  n_charge=pr.n;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf::calc_temp_e_solve_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  int ret=1;

  if (guess_set) {
    
    // If an initial guess is given, then use it to directly compute
    // the EOS

    x[0]=ne.mu;
    x[1]=pr.mu;
    x[2]=sigma;
    x[3]=omega;
    x[4]=rho;
    guess_set=false;

    if (verbose>0) {
      cout << "Solving in eos_had_rmf::calc_temp_e()." << endl;
      cout << " Init. guess: " << ne.mu << " " << pr.mu << " "
	   << sigma << " " << omega << " " << rho << endl;
    }

    bool ent=eos_mroot->err_nonconv;
    eos_mroot->err_nonconv=false;
    ret=eos_mroot->msolve(5,x,fmf);
    if (ret!=0) {
      cout << " Return value: " << ret << endl;
    }
    eos_mroot->err_nonconv=ent;
    calc_temp_e_solve_fun(5,x,y);
    
  } 

  // If no initial guess is given, or if the initial guess failed,
  // begin at saturated nuclear matter and proceeding incrementally.

  if (ret!=0) {

    double nn=ne.n;
    double np=pr.n;

    x[0]=4.8;
    x[1]=4.8;
    x[2]=0.1;
    x[3]=0.07;
    x[4]=-0.001;
    
    bool log_mode=false;
    if (fabs(log10(0.16/n_baryon))>6.0) log_mode=true;
    
    if (verbose>0) {
      cout << " Proceeding incrementally." << endl;
      cout << " alpha       n_B         n_ch        mu_n        "
	   << "mu_p        sigma       omega       rho         ret" << endl;
      cout.setf(ios::showpos);
      cout.precision(4);
    }
    
    for(size_t i=0;i<calc_e_steps;i++) {

      double alpha=((double)i)/((double)(calc_e_steps-1));

      // At the last point, there are finite precision problems when 
      // alpha=1, so we handle the last point separately
      if (i==calc_e_steps-1) {
	n_baryon=nn+np;
	n_charge=np;
      } else {
	if (log_mode) {
	  n_baryon=0.16*pow((nn+np)/0.16,alpha);
	  if (np>0.0) {
	    n_charge=0.08*pow(np/0.08,alpha);
	  } else {
	    n_charge=0.08*(1.0-alpha)+np*alpha;
	  }
	} else {
	  n_baryon=0.16*(1.0-alpha)+(nn+np)*alpha;
	  n_charge=0.08*(1.0-alpha)+np*alpha;
	}
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
      
      ret=eos_mroot->msolve(5,x,fmf);
      if (verbose>0) {
	cout << alpha << " " << n_baryon << " " << n_charge << " "
	     << x[0] << " " << x[1] << " " << x[2] << " " 
	     << x[3] << " " << x[4] << " " << ret << endl;
      }
    }
    if (verbose>0) {
      cout.unsetf(ios::showpos);
      cout.precision(6);
      cout << endl;
    }
    calc_temp_e_solve_fun(5,x,y);
    if (verbose>0) {
      cout.setf(ios::scientific);
      cout << "x: ";
      vector_out(cout,x,true);
      cout << "y: ";
      vector_out(cout,y,true);
      cout.unsetf(ios::scientific);
    }
  }
  
  sigma=x[2];
  omega=x[3];
  rho=x[4];
  
  if (ret!=0) {
    O2SCL_CONV2_RET("Solver failed in eos_had_rmf::calc_temp_e(fermion,",
		    "fermion,double,thermo).",exc_efailed,this->err_nonconv);
  }
  
  return 0;
}

int eos_had_rmf::calc_eq_p(fermion &ne, fermion &pr, double sig, double ome, 
			   double lrho, double &f1, double &f2, double &f3, 
			   thermo &lth) {

  ne.non_interacting=false;
  pr.non_interacting=false;
  
  double duds,us,fun,gs2,gr2,gw2,gs,gr,gw,sig2,ome2;
  double rho2,sig4,ome4,dfdw,fnn=0.0,fnp=0.0;

  gs=ms*cs;
  gw=mw*cw;
  gr=mr*cr;
  gs2=gs*gs;
  gw2=gw*gw;
  gr2=gr*gr;

  if (zm_mode) {
    fnn=(1.0+gs*sig/ne.m);
    fnp=(1.0+gs*sig/pr.m);
    ne.ms=ne.m/fnn;
    pr.ms=pr.m/fnp;
  } else {
    ne.ms=ne.m-gs*sig;
    pr.ms=pr.m-gs*sig;
  }
  
  if (ne.ms<0.0 || pr.ms<0.0) {
    O2SCL_CONV2_RET("Neutron or proton mass negative in ",
		    "eos_had_rmf::calc_eq_p().",exc_efailed,
		    this->err_nonconv);
  }

  ne.nu=ne.mu-gw*ome+0.5*gr*lrho;
  pr.nu=pr.mu-gw*ome-0.5*gr*lrho;

  if (ne.nu<ne.ms || ne.ms<0.0) ne.kf=0.0;
  else ne.kf=sqrt(ne.nu*ne.nu-ne.ms*ne.ms);
  if (pr.nu<pr.ms || pr.ms<0.0) pr.kf=0.0;
  else pr.kf=sqrt(pr.nu*pr.nu-pr.ms*pr.ms);

  // We don't record error values, since these functions usually
  // always succeed
  fet->calc_mu_zerot(ne);
  fet->calc_mu_zerot(pr);

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

  if (zm_mode==true) {
    f1=ms*ms*sig-gs*(nsn/fnn/fnn+nsp/fnp/fnp)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig);
  } else {
    f1=ms*ms*sig-gs*(nsn+nsp)+duds-gr2*rho2*
      (a1+2.0*a2*sig+3.0*a3*sig2+4.0*a4*sig2*sig+
       5.0*a5*sig4+6.0*a6*sig4*sig);
  }
  f2=mw*mw*ome-gw*(ne.n+pr.n)+zeta*gw2*gw2*ome2*ome/6.0+gr2*rho2*
    (2.0*b1*ome+4.0*b2*ome2*ome+6.0*b3*ome4*ome);
  f3=mr*mr*lrho-gr*(pr.n-ne.n)*0.5+xi*gr2*gr2*rho2*lrho/6.0+
    2.0*gr2*lrho*fun;

  us=b/3.0*ne.m*gs2*gs*sig2*sig+c/4.0*gs2*gs2*sig4;
  
  // The thermodynamic identity could be used instead of 
  // these explicit expressions, but it's nice to have them
  // available here.

  lth.pr=-us-0.5*ms*ms*sig2+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho+
    zeta/24.0*gw2*gw2*ome4+xi/24.0*gr2*gr2*rho2*rho2+
    ne.pr+pr.pr+gr2*rho2*fun;

  lth.ed=us+0.5*ms*ms*sig2-0.5*mw*mw*ome*ome-0.5*mr*mr*lrho*lrho+
    ne.ed+pr.ed-zeta/24.0*gw2*gw2*ome4-xi/24.0*gr2*gr2*rho2*rho2-
    fun*gr2*rho2+gw*ome*(ne.n+pr.n)-0.5*gr*lrho*(ne.n-pr.n);

  return success;
}

int eos_had_rmf::fix_saturation2_fun(size_t nv, const ubvector &x, 
				     ubvector &y, double fix_n0,
				     double fix_eoa, double fix_comp,
				     double fix_esym, double fix_msom) {
  
  cs=x[0];
  cw=x[1];
  cr=x[2];
  b=x[3];
  c=x[4];

  saturation();

  y[0]=(comp-fix_comp)/fix_comp;
  y[1]=(esym-fix_esym)/fix_esym;
  y[2]=(msom-fix_msom)/fix_msom;
  y[3]=(n0-fix_n0)/fix_n0;
  y[4]=(eoa-fix_eoa)/fix_comp;

  return 0;
}

int eos_had_rmf::fix_saturation_fun(size_t nv, const ubvector &x, 
				    ubvector &y) {
  
  double phi,power,ome,dome,one,two,tri;
  double cq,lr,lr13,sqt,fir,sec,kf,kf2,pr,ed,lfac;
  double efs,rhos,kbar,lbar,us,duds,uv,eps,denr,aknm,mns;
  
  cs=x[0];
  cw=x[1]; 
  b=x[2];
  c=x[3];

  // b and c can be negative but if cw or cs are negative
  // then omega or sigma are negative respectively
  if (cs<0.0 || cw<0.0) return 1;

  // Calculate phi from effective mass
  phi=mnuc*(1.0-msom);
  mns=mnuc*msom;

  // Calculate omega and domega We don't have to worry about the
  // non-linear rho interactions, since rho is zero for nuclear
  // matter.

  // This section is just the solution of the vector meson field
  // equation
  power=1.0/3.0;
  if (zeta>0.0) {
    cq=8.0/(9.0*pow(cw*cw,3.0)*n0*n0*zeta);
    lr=3*n0/zeta;
    lr13=pow(lr,power);
    sqt=sqrt(1.0+cq);
    fir=pow((1.0+sqt),power);
    sec=pow(fabs(1.0-sqt),power);
    ome=lr13*(fir-sec);
    dome=1.0/(1.0/(cw*cw)+zeta/2.0*ome*ome);
  } else {
    ome=cw*cw*n0;
    dome=cw*cw;
  }

  kf=pow(1.5*pi2*n0,power);
  kf2=kf*kf;
  efs=sqrt(kf*kf+mns*mns);

  // Scalar density
  rhos=1.0/pi2*mns*(kf*efs-mns*mns*log((kf+efs)/mns));

  kbar=2.0*b*mnuc;
  lbar=6.0*c;
  us=kbar/6.0*pow(phi,3.0)+lbar/24.0*pow(phi,4.0);
  duds=kbar/2.0*phi*phi+lbar/6.0*pow(phi,3.0);
  uv=-zeta/24.0*pow(ome,4.0);

  // Energy density and pressure
  lfac=(kf+efs)/mns;
  pr=1.0/24.0/pi2*(2.0*efs*pow(kf,3.0)-3.0*kf*efs*mns*mns
                   +3.0*pow(mns,4.0)*log(lfac));
  ed=1.0/8.0/pi2*(2.0*kf*pow(efs,3.0)-kf*efs*mns*mns
                  -pow(mns,4.0)*log(lfac));
  eps=phi*phi/(2.0*cs*cs)+ome*n0-ome*ome/(2.0*cw*cw)+
    us+uv+ed*2.0;

  // Compressibility - This expression for the compressibility
  // only works for nuclear matter density. 
  one=1.0/cs/cs+kbar*phi+lbar/2.0*phi*phi;
  two=3.0/mns*rhos;
  tri=3.0*n0/efs;
  denr=one+two-tri;
  aknm=9.0*n0*(dome+kf2/(3.0*efs*n0)-pow(mns/efs,2.0)/denr);

  // Scalar field equation
  y[0]=phi/cs/cs+duds-rhos;

  // Binding energy of nuclear matter
  y[1]=(eps/n0-mnuc)-eoa;

  // Pressure
  y[2]=-phi*phi/(2.0*cs*cs)+ome*ome/(2.0*cw*cw)-us-uv+pr*2.0;

  // Compressibilty
  y[3]=aknm-comp;

  if (!std::isfinite(y[1]) || !std::isfinite(y[2]) || 
      !std::isfinite(y[3]) || !std::isfinite(y[0])) {
    O2SCL_ERR2("Equation not finite in ",
	       "eos_had_rmf::fix_saturation_fun().",exc_efailed);
  }
  return 0;
}

int eos_had_rmf::fix_saturation2(double gcs, double gcw,
				 double gcr, double gb, double gc) {
  ubvector x(5);

  x[0]=gcs;
  x[1]=gcw;
  x[2]=gcr;
  x[3]=gb;
  x[4]=gc;

  double fix_n0=n0, fix_eoa=eoa, fix_comp=comp, fix_esym=esym;
  double fix_msom=msom;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
		     double, double, double, double)>
     (&eos_had_rmf::fix_saturation2_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,fix_n0,fix_eoa,fix_comp,fix_esym,
     fix_msom);

  int test=sat_mroot->msolve(5,x,fmf);
  if (test!=0) {
    O2SCL_ERR("Solve failed in fix_saturation2().",exc_efailed);
  }

  cs=x[0];
  cw=x[1];
  cr=x[2];
  b=x[3];
  c=x[4];

  saturation();
  
  return 0;
}

int eos_had_rmf::fix_saturation(double gcs, double gcw, double gb, double gc) {
  int nvar, test;
  ubvector x(4);
  double power,kf,kf2,kf3,efs;
  double mns,gw,ome,lr,lr13,cq,sqt,fir,sec,sig,gs;

  if (zm_mode==true) {
    O2SCL_ERR("Function fix_saturation() does not work for zm_mode.",
	      exc_eunimpl);
    return exc_eunimpl;
  }

  nvar=4;

  x[0]=gcs;
  x[1]=gcw;
  x[2]=gb;
  x[3]=gc;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf::fix_saturation_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  test=sat_mroot->msolve(4,x,fmf);
  if (test!=0) {
    O2SCL_ERR("Solve failed in fix_saturation().",exc_efailed);
  }

  cs=x[0];
  cw=x[1];
  b=x[2];
  c=x[3];

  power=1.0/3.0;
  kf=pow(1.5*pi2*n0,power);
  kf2=kf*kf;
  kf3=kf2*kf;
  mns=msom*mnuc;
  efs=sqrt(kf*kf+mns*mns);

  gw=mw*cw;
  if (zeta>0.0) {
    cq=8.0/(9.0*pow(cw*cw,3.0)*n0*n0*zeta);
    lr=3*n0/zeta;
    lr13=pow(lr,power);
    sqt=sqrt(1.0+cq);
    fir=pow((1.0+sqt),power);
    sec=pow(fabs(1.0-sqt),power);
    ome=lr13*(fir-sec);
  } else {
    ome=(cw*cw)*n0;
  }
  ome/=gw;
  
  gs=ms*cs;
  sig=mnuc*(1.0-msom)/gs;
  if (zm_mode) {
    O2SCL_ERR("Function fix_saturation() does not work with zm_mode=true.",
	      exc_efailed);
  }

  if (calc_cr(sig,ome,n0)==0) {
    return 0;
  } 
  
  O2SCL_ERR("calc_cr() failed in fix_saturation().",exc_efailed);
  return exc_efailed;
}

int eos_had_rmf::saturation() {

  ubvector x(5), y(5);
  int test;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_had_rmf::zero_pressure),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  if (guess_set) {
    x[0]=neutron->mu;
    x[1]=proton->mu;
    x[2]=sigma;
    x[3]=omega;
    x[4]=rho;
    guess_set=false;
  } else {
    x[0]=4.8;
    x[1]=4.8;
    x[2]=0.1;
    x[3]=0.07;
    x[4]=0.0;
  }

  // Ensure that the initial guess has a finite initial density
  int max_it=20;
  int it=0;
  if (x[0]<0.0) x[0]=0.0;
  if (x[1]<0.0) x[1]=0.0;
  zero_pressure(5,x,y);
  while (it<max_it && (neutron->n<=0.0 || proton->n<=0.0)) {
    x[0]+=neutron->m*0.05;
    x[1]+=proton->m*0.05;
    zero_pressure(5,x,y);
    it++;
    if (verbose>0) {
      cout << "eos_had_rmf::saturation() fixing density: " << it 
	   << "\n\tmun=" << x[0] << " mup=" << x[1] << " sig=" << x[2]
	   << " ome=" << x[3] << "\n\trho=" << x[4] << " y[0]=" << y[0]
	   << " y[1]=" << y[1] << "\n\ty[2]=" << y[2] << " y[3]="
	   << y[3] << " y[4]=" << y[4] << endl;
    }
  } 
  if (it==max_it) {
    O2SCL_CONV_RET("Failed to make initial density finite in saturation()",
                   exc_efailed,this->err_nonconv);
  }
  
  if (verbose>0) {
    cout << "eos_had_rmf::saturation() initial guess: \n\t" << x[0] << " "
	 << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << endl;
  }

  test=sat_mroot->msolve(5,x,fmf);
  
  sigma=x[2];
  omega=x[3];
  rho=x[4];

  if (test!=0) {
    O2SCL_CONV_RET("Solver failed in eos_had_rmf::saturation().",
		   exc_efailed,this->err_nonconv);
  }
  
  if (verbose>0) {
    cout << "eos_had_rmf::saturation() final value:\n\t" << x[0] << " "
	 << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << endl;
  }
  
  n0=neutron->n+proton->n;
  msom=neutron->ms/neutron->m;
  eoa=(eos_thermo->ed/n0-(neutron->m+proton->m)/2.0);
  
  comp=eos_had_base::fcomp(n0);
  kprime=eos_had_base::fkprime(n0);
  esym=eos_had_base::fesym(n0);

  // These can't be used because they don't work with unequal
  // neutron and proton masses. 
  // fkprime_fields(x[2],x[3],n0,comp,kprime);
  // esym=fesym_fields(x[2],x[3],n0);

  return 0;
}

int eos_had_rmf::calc_e_solve_fun(size_t nv, const ubvector &ex, 
				  ubvector &ey) {

  double f1,f2,f3,sig,ome,lrho;
  
  neutron->mu=ex[0];
  proton->mu=ex[1];
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];
  
  calc_eq_p(*neutron,*proton,sig,ome,lrho,f1,f2,f3,*eos_thermo);
  
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

int eos_had_rmf::calc_temp_e_solve_fun(size_t nv, const ubvector &ex, 
				       ubvector &ey) {
  double f1,f2,f3,sig,ome,lrho;

  neutron->mu=ex[0];
  proton->mu=ex[1];
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];
  
#if !O2SCL_NO_RANGE_CHECK
  // This may not be strictly necessary, because it should be clear
  // that this function will produce gibberish if the inputs aren't
  // finite, but I've found this extra checking of the inputs useful
  // for debugging.

  if (!std::isfinite(ex[0]) || !std::isfinite(ex[1])) {
    O2SCL_ERR2("Chemical potentials not finite in ",
	       "eos_had_rmf::calc_temp_e_solve_fun().",exc_efailed);
  }
  if (!std::isfinite(ex[2]) || !std::isfinite(ex[3]) ||
      !std::isfinite(ex[4])) {
    O2SCL_ERR2("Fields not finite in ",
	       "eos_had_rmf::calc_temp_e_solve_fun().",exc_efailed);
  }
#endif

  double gs=ms*cs;
  if (zm_mode) {
    neutron->ms=neutron->m/(1.0+gs*sig/neutron->m);
    proton->ms=proton->m/(1.0+gs*sig/proton->m);
  } else {
    neutron->ms=neutron->m-gs*sig;
    proton->ms=proton->m-gs*sig;
  }

  if (neutron->ms<0.0 || proton->ms<0.0) {
    return exc_ebadfunc;
  }

  calc_eq_temp_p(*neutron,*proton,
		 ce_temp,sig,ome,lrho,f1,f2,f3,*eos_thermo);

  if (calc_e_relative) {
    ey[0]=(proton->n+neutron->n-n_baryon)/n_baryon;
    if (n_charge>0.0) {
      ey[1]=(proton->n-n_charge)/n_charge;
    } else {
      ey[1]=proton->n-n_charge;
    }
  } else {
    ey[0]=proton->n+neutron->n-n_baryon;
    ey[1]=proton->n-n_charge;
  }
  ey[2]=f1;
  ey[3]=f2;
  ey[4]=f3;

  for(int i=0;i<5;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      O2SCL_ERR
	((((string)"Eq. ")+itos(i)+
	  " not finite in eos_had_rmf::calc_temp_e_solve_fun().").c_str(),
	 exc_efailed);
    }
  }

  return 0;
}

int eos_had_rmf::zero_pressure(size_t nv, const ubvector &ex, 
			       ubvector &ey) {

  double f1,f2,f3,sig,ome,lrho;
  fermion *n=neutron, *p=proton;
  int i;
  
  n->mu=ex[0];
  p->mu=ex[1];
  sig=ex[2];
  ome=ex[3];
  lrho=ex[4];

  for(i=0;i<5;i++) {
    if (!std::isfinite(ex[i])) {
      O2SCL_ERR("Variable not finite in zero_pressure()",exc_efailed);
    }
  }
  
  calc_eq_p(*n,*p,sig,ome,lrho,f1,f2,f3,*eos_thermo);

  ey[0]=eos_thermo->pr;
  ey[1]=p->n/(n->n+p->n)-0.5;
  ey[2]=f1;
  ey[3]=f2;
  ey[4]=f3;

  if (n->n<=0.0 || p->n<=0.0) {
    // 07/12/11 - We don't want to call the error handler here, because
    // sometimes the solver may be able to handle it automatically
    return 1;
  }

  for(i=0;i<5;i++) {
    if (!std::isfinite(ex[i]) || !std::isfinite(ey[i])) {
      // 07/12/11 - We don't want to call the error handler here, because
      // sometimes the solver may be able to handle it automatically
      return 2;
    }
  }

  return 0;
}

double eos_had_rmf::fesym_fields(double sig, double ome, double l_nb) {

  double kf, efs, mstar, ret, fun;
  
  kf=pow(1.5*l_nb*pi2,1.0/3.0);
  mstar=mnuc-cs*ms*sig;
  efs=sqrt(mstar*mstar+kf*kf);

  fun=a1*sig+a2*sig*sig+a3*pow(sig,3.0)+a4*pow(sig,4.0)+
    a5*pow(sig,5.0)+a6*pow(sig,6.0)+b1*ome*ome+b2*pow(ome,4.0)+
    b3*pow(ome,6.0);
  
  ret=kf*kf/6.0/efs+l_nb/8.0/(1.0/cr/cr+2.0*fun);
  
  return ret;
}

int eos_had_rmf::calc_cr(double sig, double ome, double l_nb) {

  double kf, efs, mstar, up, dn, fun;

  kf=pow(1.5*l_nb*pi2,1.0/3.0);
  mstar=mnuc-cs*ms*sig;
  efs=sqrt(mstar*mstar+kf*kf);

  fun=a1*sig+a2*sig*sig+a3*pow(sig,3.0)+a4*pow(sig,4.0)+
    a5*pow(sig,5.0)+a6*pow(sig,6.0)+b1*ome*ome+b2*pow(ome,4.0)+
    b3*pow(ome,6.0);

  up=4.0*kf*kf-24.0*esym*efs;
  dn=fun*(48.0*esym*efs-8.0*kf*kf)-3.0*l_nb*efs;

  cr=sqrt(up/dn);

  if (!std::isfinite(cr)) {
    O2SCL_ERR("Coupling not finite in eos_had_rmf::calc_cr()",exc_efailed);
  }

  return 0;
}

double eos_had_rmf::fcomp_fields(double sig, double ome, double l_nb) {

  double gs, gw, mstar, efs, d2u;
  double alpha, dsdn, dwdn, ret, rhos, kf;

  gs=ms*cs;
  gw=mw*cw;
  mstar=mnuc-gs*sig;
  kf=pow(1.5*pi2*l_nb,1.0/3.0);
  efs=sqrt(kf*kf+mstar*mstar);
  rhos=(mstar*kf*efs-pow(mstar,3.0)*log((kf+efs)/mstar))/pi2;

  d2u=pow(gs,3.0)*(2.0*b*mnuc*sig+3.0*c*gs*sig*sig);

  alpha=0.5/pi2/efs*(2.0*pow(kf,3.0)+6.0*mstar*mstar*kf-
		     6.0*mstar*mstar*efs*log((kf+efs)/mstar));
  dsdn=gs*mstar/efs/(ms*ms+d2u+gs*gs*alpha);
  dwdn=gw/(mw*mw+zeta/2.0*pow(gw*gw*ome,2.0));

  ret=dwdn*(mw*mw*ome+zeta/6.0*pow(gw*ome,3.0)*gw);
  ret+=kf*kf/3.0/efs;
  ret+=dsdn*(gs*mstar*alpha/3.0-gs*rhos);
  ret*=9.0;

  return ret;
}
  
void eos_had_rmf::fkprime_fields(double sig, double ome, double l_nb,
				 double &k, double &l_kprime) {

  double gs, gw, mstar, efs, d2u;
  double alpha, dsdn, dwdn, rhos, kf, lterm;
  double dpdn, d3u, dalphadms, d2sdn2, d2pdn2, d2wdn2;

  gs=ms*cs;
  gw=mw*cw;
  mstar=mnuc-gs*sig;
  kf=pow(1.5*pi2*l_nb,1.0/3.0);
  efs=sqrt(kf*kf+mstar*mstar);
  lterm=log((kf+efs)/mstar);
  rhos=(mstar*kf*efs-pow(mstar,3.0)*lterm)/pi2;

  // d2u is d^2(U(sigma))/d(sigma)
  d2u=pow(gs,3.0)*(2.0*b*mnuc*sig+3.0*c*gs*sig*sig);

  // alpha is d(scalar density)/d(mstar)
  alpha=1.0/pi2/efs*(pow(kf,3.0)+3.0*mstar*mstar*kf-
		     3.0*mstar*mstar*efs*lterm);
  // d(sigma)/dn
  dsdn=gs*mstar/efs/(ms*ms+d2u+gs*gs*alpha);
  // d(omega)/dn
  dwdn=gw/(mw*mw+zeta/2.0*pow(gw*gw*ome,2.0));

  // d(Pressure)/dn
  dpdn=dwdn*(mw*mw*ome+zeta/6.0*pow(gw*ome,3.0)*gw);
  dpdn+=kf*kf/3.0/efs;
  dpdn+=dsdn*(gs*mstar*alpha/3.0-gs*rhos);

  // d3u is d^3(U(sigma))/d(sigma)
  d3u=pow(gs,3.0)*(2.0*b*mnuc+6.0*c*gs*sig);

  dalphadms=-(8.0*pow(kf,3.0)+6.0*kf*mstar*mstar-6.0*pow(efs,3.0)*
	      lterm)*mstar/pi2/pow(efs,3.0);

  // d^2(sigma)/dn^2
  d2sdn2=-pi2/2.0/kf/efs/efs*dsdn+pow(dsdn,3.0)*
    (gs*efs/mstar*dalphadms-efs/mstar/gs*d3u);

  // d^2(omega)/dn^2
  d2wdn2=-pow(gw*dwdn,3.0)*zeta*ome;

  // d^2(Pressure)/dn^2
  d2pdn2=dsdn*(2.0*gs*kf*kf*mstar/3.0/pow(efs,3.0)-gs*mstar/efs)+
    dsdn*dsdn*(2*gs*gs*alpha/3.0-gs*gs*mstar/3.0*dalphadms)+
    d2sdn2*(gs*mstar*alpha/3.0-gs*rhos)+
    dwdn*dwdn*(mw*mw+zeta/2.0*pow(gw*gw*ome,2.0))+
    d2wdn2*(mw*mw*ome+zeta/6.0*gw*pow(gw*ome,3.0))+
    pi2/3.0/efs*(1.0/kf-kf/2.0/efs/efs);

  k=9.0*dpdn;
  l_kprime=27.0*l_nb*(d2pdn2-4.0*dpdn/l_nb);
  
  return;
}

