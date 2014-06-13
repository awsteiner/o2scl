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

#include <o2scl/eos_had_base.h>
// For unit conversions
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

eos_had_base::eos_had_base() {

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

double eos_had_base::fcomp(double nb, const double &alpha) {
  double lcomp, err;
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
		       (&eos_had_base::calc_pressure_nb),
		       this,std::placeholders::_1,alpha);

  lcomp=9.0*sat_deriv->deriv(nb,fmn);

  return lcomp;
}

double eos_had_base::fcomp_err(double nb, double alpha, double &unc) {
  double lcomp;
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
		       (&eos_had_base::calc_pressure_nb),
		       this,std::placeholders::_1,alpha);

  sat_deriv->deriv_err(nb,fmn,lcomp,unc);

  lcomp*=9.0;
  
  return lcomp;
}

double eos_had_base::feoa(double nb, const double &alpha) {
  double leoa;

  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  leoa=(eos_thermo->ed-neutron->n*neutron->m-
	proton->n*proton->m)/nb;
  
  return leoa;
}

double eos_had_base::fesym(double nb, const double &alpha) {
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::calc_dmu_alpha),
			this,std::placeholders::_1,nb);

  return sat_deriv->deriv(alpha,fmn)/4.0;

  // * Old method using second derivative *
  //funct_mfptr_param<eos_had_base,const double> 
  //fmn(this,&eos_had_base::calc_edensity_alpha,nb);
  //return sat_deriv->calc2(alpha,fmn)/2.0/nb;
}

double eos_had_base::fesym_err(double nb, double alpha, 
			       double &unc) {

  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::calc_dmu_alpha),
			this,std::placeholders::_1,nb);

  double val, err;
  sat_deriv->deriv_err(alpha,fmn,val,err);
  val/=4.0; 
  err/=4.0;
  return val;
}

double eos_had_base::fesym_slope(double nb, const double &alpha) {
  
  if (false) {
    // The form below is effectively a second derivative since it must
    // compute the symmetry energy first. This form is also a second
    // derivative, and may or may not be less accurate. It might be
    // good to make this a separate function, to allow the user to
    // choose which way to evaluate L.
    funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			  (&eos_had_base::calc_musum_alpha),
			  this,std::placeholders::_1,nb);
    
    return sat_deriv->deriv2(alpha,fmn)*0.75-3.0*fesym(nb,alpha);
  }
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::fesym),
			this,std::placeholders::_1,alpha);
  return sat_deriv2->deriv(nb,fmn)*3.0*nb;
}

double eos_had_base::fesym_curve(double nb, const double &alpha) {
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::fesym),
			this,std::placeholders::_1,alpha);
  
  return sat_deriv2->deriv2(nb,fmn)*9.0*nb*nb;
}

double eos_had_base::fesym_skew(double nb, const double &alpha) {

  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::fesym),
			this,std::placeholders::_1,alpha);

  return sat_deriv2->deriv3(nb,fmn)*27.0*nb*nb*nb;
}

double eos_had_base::fesym_diff(double nb) {
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

double eos_had_base::feta(double nb) {
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

double eos_had_base::fkprime(double nb, const double &alpha) {
  double lkprime, err;
  int ret=0;
  
  funct11 fmn=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::calc_press_over_den2),
			this,std::placeholders::_1,alpha);
  
  sat_deriv->deriv2_err(nb,fmn,lkprime,err);
  lkprime*=27.0*nb*nb*nb;
  
  return lkprime;
}

double eos_had_base::fmsom(double nb, const double &alpha) {

  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->ms/neutron->m;
}

double eos_had_base::fn0(double alpha, double &leoa) {
  double nb;
  int ret=0;
  
  // Initial guess
  nb=0.16;
  
  funct11 fmf=std::bind(std::mem_fn<double(double,const double &)>
			(&eos_had_base::calc_pressure_nb),
			this,std::placeholders::_1,alpha);
  
  sat_root->solve(nb,fmf);
  calc_pressure_nb(nb);
  leoa=eos_thermo->ed/nb-(neutron->m+proton->m)/2.0;
  
  return nb;
}

void eos_had_base::saturation() {
  n0=fn0(0.0,eoa);
  comp=fcomp(n0);
  esym=fesym(n0);
  msom=fmsom(n0);
  kprime=fkprime(n0);

  return;
}

void eos_had_base::gradient_qij(fermion &n, fermion &p, thermo &th,
			       double &qnn, double &qnp, double &qpp, 
			       double &dqnndnn, double &dqnndnp,
			       double &dqnpdnn, double &dqnpdnp,
			       double &dqppdnn, double &dqppdnp) {
  double nn=n.n, np=p.n, t1=0.0, t2=0.0, barn=nn+np, den;

  int vpx=0;

  t1=t1_fun(barn);
  t2=t2_fun(barn);
    
  funct11 t1fun=std::bind(std::mem_fn<double(double)>
			  (&eos_had_base::t1_fun),
			  this,std::placeholders::_1);
  funct11 t2fun=std::bind(std::mem_fn<double(double)>
			  (&eos_had_base::t2_fun),
			  this,std::placeholders::_1);

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

double eos_had_base::t1_fun(double barn) {
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

double eos_had_base::t2_fun(double barn) {
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

double eos_had_base::calc_pressure_nb(double nb, const double &alpha) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->pr;
}

void eos_had_base::const_pf_derivs(double nb, double pf, 
				   double &dednb_pf, double &dPdnb_pf) {

  // Take derivatives w.r.t. alpha and then multiply by -2 to get
  // derivatives w.r.t. x
  funct11 fmpp=std::bind(std::mem_fn<double(double,const double &)>
			 (&eos_had_base::calc_pressure_nb),
			 this,std::placeholders::_1,1.0-2.0*pf);
  funct11 fmpe=std::bind(std::mem_fn<double(double,const double &)>
			 (&eos_had_base::calc_edensity_nb),
			 this,std::placeholders::_1,1.0-2.0*pf);

  dPdnb_pf=-2.0*sat_deriv->deriv(nb,fmpp);
  dednb_pf=-2.0*sat_deriv->deriv(nb,fmpe);
  return;
}

double eos_had_base::calc_press_over_den2(double nb, const double &alpha) {
  
  neutron->n=nb/2.0;
  proton->n=nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->pr/nb/nb;
}

double eos_had_base::calc_edensity_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->ed;
}

double eos_had_base::calc_dmu_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu-proton->mu;
}

double eos_had_base::calc_musum_alpha(double alpha, const double &nb) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu+proton->mu;
}

double eos_had_base::calc_edensity_nb(double nb, const double &alpha) {
  
  neutron->n=(1.0+alpha)*nb/2.0;
  proton->n=(1.0-alpha)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->ed;
}

int eos_had_base::nuc_matter_e(size_t nv, const ubvector &x, 
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

int eos_had_base::nuc_matter_p(size_t nv, const ubvector &x, 
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

void eos_had_base::set_n_and_p(fermion &n, fermion &p) {
  neutron=&n;
  proton=&p;
  return;
}

void eos_had_base::set_mroot(mroot<mm_funct11,
				  boost::numeric::ublas::vector<double>, 
				  jac_funct<> > &mr) {
  eos_mroot=&mr;
  return;
}

void eos_had_base::set_sat_root(root<funct11> &mr) {
  sat_root=&mr;
  return;
}

void eos_had_base::set_sat_deriv(deriv_base<funct11> &de) {
  sat_deriv=&de;
  return;
}

void eos_had_base::set_sat_deriv2(deriv_base<funct11> &de) {
  sat_deriv2=&de;
  return;
}

int eos_had_eden_base::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;
    
  double pa[2]={n.mu,p.mu};
  double *pap=&(pa[0]);
  
  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>(&eos_had_eden_base::nuc_matter_e),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);
  eos_mroot->msolve(2,x,fmf);
    
  th=*eos_thermo;

  return 0;
}

int eos_had_pres_base::calc_e(fermion &n, fermion &p, thermo &th) {
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

  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>(&eos_had_eden_base::nuc_matter_p),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);
  eos_mroot->msolve(2,mu,fmf);
    
  th=*eos_thermo;
    
  return 0;
}

int eos_had_temp_base::nuc_matter_temp_e(size_t nv, const ubvector &x, 
					 ubvector &y, double *&pa) {
  neutron->n=x[0];
  proton->n=x[1];
  
  if (!o2scl::is_finite(neutron->n) || !o2scl::is_finite(proton->n)) {
    O2SCL_ERR2("Density problem in ",
		   "eos_had_temp_base::nuc_matter_e().",exc_esanity);
  }
  int ret=calc_temp_e(*neutron,*proton,lT,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR2("Function calc_e() failed in ",
	       "eos_had_temp_base::nuc_matter_e().",exc_efailed);
  }

  y[0]=neutron->mu-pa[0];
  y[1]=proton->mu-pa[1];
  
  if (!o2scl::is_finite(neutron->mu) || !o2scl::is_finite(proton->mu)) {
    O2SCL_ERR2("Chemical potential problem in ",
		   "eos_had_temp_base::nuc_matter_e().",exc_esanity);
  }

  return ret;
}

int eos_had_temp_base::nuc_matter_temp_p(size_t nv, const ubvector &x, 
					 ubvector &y, double *&pa) {
  
  neutron->mu=x[0];
  proton->mu=x[1];

  int ret=calc_temp_p(*neutron,*proton,lT,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR("calc_p() failed in eos_had_temp_base::nuc_matter_p().",ret);
  }

  y[0]=neutron->n-pa[0];
  y[1]=proton->n-pa[1];
  
  return ret;
}


int eos_had_temp_eden_base::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;

  double pa[2]={n.mu,p.mu};
  double *pap=&(pa[0]);

  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>(&eos_had_eden_base::nuc_matter_e),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);

  ret=eos_mroot->msolve(2,x,fmf);
  if (ret!=0) {
    O2SCL_ERR("Solver failed in eos_had_temp_eden_base::calc_p().",ret);
  }
    
  th=*eos_thermo;
  
  return 0;
}

int eos_had_temp_eden_base::calc_temp_p(fermion &n, fermion &p, 
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
  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>
    (&eos_had_temp_eden_base::nuc_matter_temp_e),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);

  ret=eos_mroot->msolve(2,den,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2("Solver failed in ",
		   "eos_had_temp_eden_base::calc_temp_p().",ret);
  }
  
  th=*eos_thermo;

  return 0;
}

int eos_had_temp_pres_base::calc_e(fermion &n, fermion &p, thermo &th) {
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

  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>
    (&eos_had_eden_base::nuc_matter_p),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);

  ret=eos_mroot->msolve(2,mu,fmf);
    
  if (ret!=0) {
    O2SCL_ERR("Solver failed in eos_had_temp_pres_base::calc_p().",ret);
  }
    
  th=*eos_thermo;
    
  return 0;
}

int eos_had_temp_pres_base::calc_temp_e(fermion &n, fermion &p, 
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

  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double *&)>
    (&eos_had_temp_pres_base::nuc_matter_temp_p),
    this,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,pap);

  ret=eos_mroot->msolve(2,mu,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2("Solver failed in ",
		   "eos_had_temp_pres_base::calc_temp_e().",ret);
  }
  
  th=*eos_thermo;

  return 0;
}

