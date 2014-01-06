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

#include <o2scl/hadronic_eos.h>
// For unit conversions
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

hadronic_eos::hadronic_eos() {

  sat_deriv=&def_deriv;
  sat_deriv2=&def_deriv2;
  def_deriv.h=1.0e-3;
  def_deriv2.h=1.0e-3;

  def_neutron.init(o2scl_settings.get_convert_units().convert
		   ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  def_proton.init(o2scl_settings.get_convert_units().convert
		  ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  def_neutron.non_interacting=false;
  def_proton.non_interacting=false;
  neutron=&def_neutron;
  proton=&def_proton;
  
  eos_mroot=&def_mroot;
  sat_root=&def_sat_root;
}

double hadronic_eos::fcomp(double nb, const double &alpha) {
  double lcomp, err;
  
  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::calc_pressure_nb,alpha);
  lcomp=9.0*sat_deriv->deriv(nb,fmn);

  return lcomp;
}

double hadronic_eos::fcomp_err(double nb, double alpha, double &unc) {
  double lcomp;
  
  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::calc_pressure_nb,alpha);
  sat_deriv->deriv_err(nb,fmn,lcomp,unc);

  lcomp*=9.0;
  
  return lcomp;
}

double hadronic_eos::feoa(double nb, const double &alpha) {
  double leoa;

  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  leoa=(eos_thermo->ed-neutron->n*neutron->m-
	proton->n*proton->m)/nb;
  
  return leoa;
}

double hadronic_eos::fesym(double nb, const double &alpha) {

  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::calc_dmu_alpha,nb);
  return sat_deriv->deriv(alpha,fmn)/4.0;

  // * Old method using second derivative *
  //funct_mfptr_param<hadronic_eos,const double> 
  //fmn(this,&hadronic_eos::calc_edensity_alpha,nb);
  //return sat_deriv->calc2(alpha,fmn)/2.0/nb;
}

double hadronic_eos::fesym_err(double nb, double alpha, 
			       double &unc) {
  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::calc_dmu_alpha,nb);
  double val, err;
  sat_deriv->deriv_err(alpha,fmn,val,err);
  val/=4.0; 
  err/=4.0;
  return val;
}

double hadronic_eos::fesym_slope(double nb, const double &alpha) {
  
  if (false) {
    // The form below is effectively a second derivative since it must
    // compute the symmetry energy first. This form is also a second
    // derivative, and may or may not be less accurate. It might be
    // good to make this a separate function, to allow the user to
    // choose which way to evaluate L.
    funct_mfptr_param<hadronic_eos,const double> 
      fmn(this,&hadronic_eos::calc_musum_alpha,nb);
    return sat_deriv->deriv2(alpha,fmn)*0.75-3.0*fesym(nb,alpha);
  }

  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::fesym,alpha);
  return sat_deriv2->deriv(nb,fmn)*3.0*nb;
}

double hadronic_eos::fesym_curve(double nb, const double &alpha) {

  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::fesym,alpha);
  return sat_deriv2->deriv2(nb,fmn)*9.0*nb*nb;
}

double hadronic_eos::fesym_skew(double nb, const double &alpha) {

  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::fesym,alpha);
  return sat_deriv2->deriv3(nb,fmn)*27.0*nb*nb*nb;
}

double hadronic_eos::fesym_diff(double nb) {
  double eoa_neut, eoa_nuc;

  neutron->n=nb;
  proton->n=0.0;
  calc_e(*neutron,*proton,*eos_thermo);
  eoa_neut=eos_thermo->ed/nb-neutron->m;
  
  neutron->n=nb/2.0;
  proton->n=nb/2.0;
  calc_e(*neutron,*proton,*eos_thermo);
  eoa_nuc=eos_thermo->ed/nb-(neutron->m+proton->m)/2.0;
  
  return eoa_neut-eoa_nuc;
}

double hadronic_eos::feta(double nb) {
  double eoa_neut, eoa_nuc, eoa_mixed;

  neutron->n=nb;
  proton->n=0.0;
  calc_e(*neutron,*proton,*eos_thermo);
  eoa_neut=eos_thermo->ed/nb-neutron->m;
  
  neutron->n=nb/2.0;
  proton->n=nb/2.0;
  calc_e(*neutron,*proton,*eos_thermo);
  eoa_nuc=eos_thermo->ed/nb-(neutron->m+proton->m)/2.0;

  neutron->n=nb*0.75;
  proton->n=nb*0.25;
  calc_e(*neutron,*proton,*eos_thermo);
  eoa_mixed=eos_thermo->ed/nb-(neutron->m+proton->m)/2.0;
  
  return (eoa_neut-eoa_mixed)/3.0/(eoa_mixed-eoa_nuc);
}

double hadronic_eos::fkprime(double nb, const double &alpha) {
  double lkprime, err;
  int ret=0;
  
  funct_mfptr_param<hadronic_eos,const double> 
    fmn(this,&hadronic_eos::calc_press_over_den2,alpha);
  sat_deriv->deriv2_err(nb,fmn,lkprime,err);
  lkprime*=27.0*nb*nb*nb;
  
  return lkprime;
}

double hadronic_eos::fmsom(double nb, const double &alpha) {

  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->ms/neutron->m;
}

double hadronic_eos::fn0(double alpha, double &leoa) {
  double nb;
  int ret=0;
  
  // Initial guess
  nb=0.16;
  
  funct_mfptr_param<hadronic_eos,const double> 
    fmf(this,&hadronic_eos::calc_pressure_nb,alpha);
  
  sat_root->solve(nb,fmf);
  calc_pressure_nb(nb);
  leoa=eos_thermo->ed/nb-(neutron->m+proton->m)/2.0;
  
  return nb;
}

void hadronic_eos::saturation() {
  n0=fn0(0.0,eoa);
  comp=fcomp(n0);
  esym=fesym(n0);
  msom=fmsom(n0);
  kprime=fkprime(n0);

  return;
}

void hadronic_eos::gradient_qij(fermion &n, fermion &p, thermo &th,
			       double &qnn, double &qnp, double &qpp, 
			       double &dqnndnn, double &dqnndnp,
			       double &dqnpdnn, double &dqnpdnp,
			       double &dqppdnn, double &dqppdnp) {
  double nn=n.n, np=p.n, t1=0.0, t2=0.0, barn=nn+np, den;

  int vpx=0;

  t1=t1_fun(barn);
  t2=t2_fun(barn);
    
  funct_mfptr<hadronic_eos> t1fun(this,&hadronic_eos::t1_fun);
  funct_mfptr<hadronic_eos> t2fun(this,&hadronic_eos::t2_fun);

  set_n_and_p(n,p);
  set_thermo(th);
  double dt1=sat_deriv->deriv(barn,t1fun);
  double ddt1=sat_deriv->deriv2(barn,t1fun);
  double dt2=sat_deriv->deriv(barn,t2fun);
  
  qnn=0.0625*(3.0*t1-3.0*t2+2.0*(2.0*barn-nn)*dt1);
  qnp=0.0625*(6.0*t1-2.0*t2+3*barn*dt1);
  qpp=0.0625*(3.0*t1-3.0*t2+2.0*(2.0*barn-np)*dt1);
  
  dqnndnn=0.0625*(3.0*dt1-3.0*dt2+2.0*dt1+2.0*(2.0*barn-nn)*ddt1);
  dqnndnp=0.0625*(3.0*dt1-3.0*dt2+4.0*dt1+2.0*(2.0*barn-nn)*ddt1);
  dqnpdnn=0.0625*(6.0*dt1-2.0*dt2+3.0*barn*ddt1+3.0*dt1);
  dqnpdnp=0.0625*(6.0*dt1-2.0*dt2+3.0*barn*ddt1+3.0*dt1);
  dqppdnn=0.0625*(3.0*dt1-3.0*dt2+4.0*dt1+2.0*(2.0*barn-np)*ddt1);
  dqppdnp=0.0625*(3.0*dt1-3.0*dt2+2.0*dt1+2.0*(2.0*barn-np)*ddt1);

  return;
}

double hadronic_eos::t1_fun(double barn) {
  double xp=proton->n/(neutron->n+proton->n);
  neutron->n=(1.0-xp)*barn;
  proton->n=xp*barn;
  calc_e(*neutron,*proton,*eos_thermo);
  double den=neutron->m*proton->m*neutron->ms*proton->ms*barn*
    (neutron->n-proton->n);
  return ((proton->m-proton->ms)*neutron->m*neutron->ms*(neutron->n+2*barn)+
	  (neutron->ms-neutron->m)*proton->ms*proton->m*
	  (proton->n+2*barn))/den;
}

double hadronic_eos::t2_fun(double barn) {
  double xp=proton->n/(neutron->n+proton->n);
  neutron->n=(1.0-xp)*barn;
  proton->n=xp*barn;
  calc_e(*neutron,*proton,*eos_thermo);
  double den=neutron->m*proton->m*neutron->ms*proton->ms*barn*
    (neutron->n-proton->n);
  return ((proton->ms-proton->m)*neutron->m*neutron->ms*(2*barn-neutron->n)+
	  (neutron->m-neutron->ms)*proton->m*proton->ms*
	  (2*barn-proton->n))/den;
}

double hadronic_eos::calc_pressure_nb(double nb, const double &alpha) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->pr;
}

void hadronic_eos::const_pf_derivs(double nb, double pf, 
				   double &dednb_pf, double &dPdnb_pf) {

  // Take derivatives w.r.t. alpha and then multiply by -2 to get
  // derivatives w.r.t. x
  funct_mfptr_param<hadronic_eos,const double> 
    fmpp(this,&hadronic_eos::calc_pressure_nb,1.0-2.0*pf);
  funct_mfptr_param<hadronic_eos,const double> 
    fmpe(this,&hadronic_eos::calc_edensity_nb,1.0-2.0*pf);
  dPdnb_pf=-2.0*sat_deriv->deriv(nb,fmpp);
  dednb_pf=-2.0*sat_deriv->deriv(nb,fmpe);
  return;
}

double hadronic_eos::calc_press_over_den2(double nb, const double &alpha) {
  
  neutron->n=nb/2.0;
  proton->n=nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->pr/nb/nb;
}

double hadronic_eos::calc_edensity_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->ed;
}

double hadronic_eos::calc_dmu_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu-proton->mu;
}

double hadronic_eos::calc_musum_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu+proton->mu;
}

double hadronic_eos::calc_edensity_nb(double nb, const double &alpha) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->ed;
}

int hadronic_eos::nuc_matter_e(size_t nv, const ubvector &x, 
			       ubvector &y, double *&pa) {
  double mun=pa[0];
  double mup=pa[1];
  
  neutron->n=x[0];
  proton->n=x[1];
  
  calc_e(*neutron,*proton,*eos_thermo);

  y[0]=neutron->mu-mun;
  y[1]=proton->mu-mup;
  
  if (!o2scl::is_finite(neutron->mu) || !o2scl::is_finite(proton->mu)) {
    return exc_efailed;
  }

  return 0;
}

int hadronic_eos::nuc_matter_p(size_t nv, const ubvector &x, 
			       ubvector &y, double *&pa) {
  
  double nn=pa[0];
  double np=pa[1];
  
  neutron->mu=x[0];
  proton->mu=x[1];

  calc_p(*neutron,*proton,*eos_thermo);

  y[0]=neutron->n-nn;
  y[1]=proton->n-np;

  if (!o2scl::is_finite(neutron->n) || !o2scl::is_finite(proton->n)) {
    return exc_efailed;
  }
  
  return 0;
}

void hadronic_eos::set_n_and_p(fermion &n, fermion &p) {
  neutron=&n;
  proton=&p;
  return;
}

void hadronic_eos::set_mroot(mroot<mm_funct<>,
				  boost::numeric::ublas::vector<double>, 
				  jac_funct<> > &mr) {
  eos_mroot=&mr;
  return;
}

void hadronic_eos::set_sat_root(root<funct > &mr) {
  sat_root=&mr;
  return;
}

void hadronic_eos::set_sat_deriv(deriv_base<funct > &de) {
  sat_deriv=&de;
  return;
}

void hadronic_eos::set_sat_deriv2(deriv_base<funct > &de) {
  sat_deriv2=&de;
  return;
}

int hadronic_eos_eden::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;
    
  double pa[2]={n.mu,p.mu};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_eden,double *> 
    fmf(this,&hadronic_eos_eden::nuc_matter_e,pap);
  eos_mroot->msolve(2,x,fmf);
    
  th=*eos_thermo;

  return 0;
}

int hadronic_eos_pres::calc_e(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector mu(2);
  mu[0]=n.mu;
  mu[1]=p.mu;
    
  if (mu[1]<n.ms) mu[0]=n.ms+1.0e-4;
  if (mu[2]<p.ms) mu[1]=p.ms+1.0e-4;
    
  double pa[2]={n.n,p.n};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_pres,double *> 
    fmf(this,&hadronic_eos_pres::nuc_matter_p,pap);
  eos_mroot->msolve(2,mu,fmf);
    
  th=*eos_thermo;
    
  return 0;
}

int hadronic_eos_temp::nuc_matter_temp_e(size_t nv, const ubvector &x, 
					 ubvector &y, double *&pa) {
  neutron->n=x[0];
  proton->n=x[1];
  
  if (!o2scl::is_finite(neutron->n) || !o2scl::is_finite(proton->n)) {
    O2SCL_ERR2_RET("Density problem in ",
		   "hadronic_eos_temp::nuc_matter_e().",exc_esanity);
  }
  int ret=calc_temp_e(*neutron,*proton,lT,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR2("Function calc_e() failed in ",
	       "hadronic_eos_temp::nuc_matter_e().",exc_efailed);
  }

  y[0]=neutron->mu-pa[0];
  y[1]=proton->mu-pa[1];
  
  if (!o2scl::is_finite(neutron->mu) || !o2scl::is_finite(proton->mu)) {
    O2SCL_ERR2_RET("Chemical potential problem in ",
		   "hadronic_eos_temp::nuc_matter_e().",exc_esanity);
  }

  return ret;
}

int hadronic_eos_temp::nuc_matter_temp_p(size_t nv, const ubvector &x, 
					 ubvector &y, double *&pa) {
  
  neutron->mu=x[0];
  proton->mu=x[1];

  int ret=calc_temp_p(*neutron,*proton,lT,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR("calc_p() failed in hadronic_eos_temp::nuc_matter_p().",ret);
  }

  y[0]=neutron->n-pa[0];
  y[1]=proton->n-pa[1];
  
  return ret;
}


int hadronic_eos_temp_eden::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;

  double pa[2]={n.mu,p.mu};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_temp_eden,double *> 
    fmf(this,&hadronic_eos_temp_eden::nuc_matter_e,pap);
  ret=eos_mroot->msolve(2,x,fmf);
  if (ret!=0) {
    O2SCL_ERR_RET("Solver failed in hadronic_eos_temp_eden::calc_p().",ret);
  }
    
  th=*eos_thermo;
  
  return 0;
}

int hadronic_eos_temp_eden::calc_temp_p(fermion &n, fermion &p, 
					double T, thermo &th) {
  int ret;

  set_n_and_p(n,p);
  set_thermo(th);
  
  lT=T;
  
  ubvector den(2);
  den[0]=n.n;
  den[1]=p.n;

  double pa[2]={n.mu,p.mu};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_temp_eden,double *>
    fmf(this,&hadronic_eos_temp_eden::nuc_matter_temp_e,pap);
  ret=eos_mroot->msolve(2,den,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2_RET("Solver failed in ",
		   "hadronic_eos_temp_eden::calc_temp_p().",ret);
  }
  
  th=*eos_thermo;

  return 0;
}

int hadronic_eos_temp_pres::calc_e(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector mu(2);
  mu[0]=n.mu;
  mu[1]=p.mu;
    
  if (mu[1]<n.ms) mu[0]=n.ms+1.0e-4;
  if (mu[2]<p.ms) mu[1]=p.ms+1.0e-4;
    
  double pa[2]={n.n,p.n};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_temp_pres,double *> 
    fmf(this,&hadronic_eos_temp_pres::nuc_matter_p,pap);
  ret=eos_mroot->msolve(2,mu,fmf);
    
  if (ret!=0) {
    O2SCL_ERR_RET("Solver failed in hadronic_eos_temp_pres::calc_p().",ret);
  }
    
  th=*eos_thermo;
    
  return 0;
}

int hadronic_eos_temp_pres::calc_temp_e(fermion &n, fermion &p, 
					double T, thermo &th) {
  int ret;

  set_n_and_p(n,p);
  set_thermo(th);
  
  lT=T;

  ubvector mu(2);
  mu[0]=n.mu;
  mu[1]=p.mu;

  double pa[2]={n.n,p.n};
  double *pap=&(pa[0]);
  mm_funct_mfptr_param<hadronic_eos_temp_pres,double *>
    fmf(this,&hadronic_eos_temp_pres::nuc_matter_temp_p,pap);
  ret=eos_mroot->msolve(2,mu,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2_RET("Solver failed in ",
		   "hadronic_eos_temp_pres::calc_temp_e().",ret);
  }
  
  th=*eos_thermo;

  return 0;
}

