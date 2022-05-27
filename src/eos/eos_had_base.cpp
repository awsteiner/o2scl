/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
  sat_mroot=&def_sat_mroot;

  err_nonconv=true;
}

double eos_had_base::fcomp(double nb, double delta) {
  double lcomp, err;
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::calc_pressure_nb),
		      this,std::placeholders::_1,delta);

  lcomp=9.0*sat_deriv->deriv(nb,fmn);

  return lcomp;
}

double eos_had_base::fcomp_err(double nb, double delta, double &unc) {
  double lcomp;
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::calc_pressure_nb),
		      this,std::placeholders::_1,delta);

  sat_deriv->deriv_err(nb,fmn,lcomp,unc);

  lcomp*=9.0;
  
  return lcomp;
}

double eos_had_base::feoa(double nb, double delta) {
  double leoa;

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  leoa=(eos_thermo->ed-neutron->n*neutron->m-proton->n*proton->m)/nb;
  
  return leoa;
}

double eos_had_base::fesym(double nb, double delta) {
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::calc_dmu_delta),
		      this,std::placeholders::_1,nb);

  return sat_deriv->deriv(delta,fmn)/4.0;

  // * Old method using second derivative *
  //funct_mfptr_param<eos_had_base,const double> 
  //fmn(this,&eos_had_base::calc_edensity_delta,nb);
  //return sat_deriv->calc2(delta,fmn)/2.0/nb;
}

double eos_had_base::fesym_err(double nb, double delta, 
			       double &unc) {

  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::calc_dmu_delta),
		      this,std::placeholders::_1,nb);

  double val, err;
  sat_deriv->deriv_err(delta,fmn,val,err);
  val/=4.0; 
  err/=4.0;
  return val;
}

double eos_had_base::fesym_slope(double nb, double delta) {
  
  if (false) {
    // The form below is effectively a second derivative since it must
    // compute the symmetry energy first. This form is also a second
    // derivative, and may or may not be less accurate. It might be
    // good to make this a separate function, to allow the user to
    // choose which way to evaluate L.
    funct fmn=std::bind(std::mem_fn<double(double,double)>
			(&eos_had_base::calc_musum_delta),
			this,std::placeholders::_1,nb);
    
    return sat_deriv->deriv2(delta,fmn)*0.75-3.0*fesym(nb,delta);
  }
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::fesym),
		      this,std::placeholders::_1,delta);
  return sat_deriv2->deriv(nb,fmn)*3.0*nb;
}

double eos_had_base::fesym_curve(double nb, double delta) {
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::fesym),
		      this,std::placeholders::_1,delta);
  
  return sat_deriv2->deriv2(nb,fmn)*9.0*nb*nb;
}

double eos_had_base::fesym_skew(double nb, double delta) {

  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::fesym),
		      this,std::placeholders::_1,delta);

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

double eos_had_base::feta_prime(double nb) {
  funct fmn=std::bind(std::mem_fn<double(double)>
		      (&eos_had_base::feta),
		      this,std::placeholders::_1);
  
  return sat_deriv->deriv(nb,fmn);
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
  eoa_mixed=eos_thermo->ed/nb-(3.0*neutron->m+proton->m)/4.0;
  
  return (eoa_neut-eoa_mixed)/3.0/(eoa_mixed-eoa_nuc);
}

double eos_had_base::fkprime(double nb, double delta) {
  double lkprime, err;
  int ret=0;
  
  funct fmn=std::bind(std::mem_fn<double(double,double)>
		      (&eos_had_base::calc_press_over_den2),
		      this,std::placeholders::_1,delta);
  
  sat_deriv->deriv2_err(nb,fmn,lkprime,err);
  lkprime*=27.0*nb*nb*nb;
  
  return lkprime;
}

double eos_had_base::fmsom(double nb, double delta) {

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->ms/neutron->m;
}

double eos_had_base::f_effm_neut(double nb, double delta) {

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->ms/neutron->m;
}

double eos_had_base::f_effm_prot(double nb, double delta) {

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return proton->ms/proton->m;
}

double eos_had_base::f_effm_scalar(double nb, double delta) {

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  double mn=neutron->ms/neutron->m;
  double mp=proton->ms/proton->m;
  
  return 2.0*mn*mp/(mn+mp);
}

double eos_had_base::f_effm_vector(double nb, double delta) {

  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  double mn=neutron->ms/neutron->m;
  double mp=proton->ms/proton->m;
  
  return 2.0*mn*mp*delta/((mn+mp)*delta+mn-mp);
}

double eos_had_base::fn0(double delta, double &leoa) {
  double nb;
  int ret=0;
  
  // Initial guess
  nb=0.16;
  
  mm_funct fmf=std::bind(std::mem_fn<int(size_t,const ubvector &,
                                         ubvector &,double)>
                         (&eos_had_base::calc_pressure_nb_mroot),
                         this,std::placeholders::_1,
                         std::placeholders::_2,std::placeholders::_3,
                         delta);

  ubvector x(1);
  x[0]=0.16;
  
  int mret=sat_mroot->msolve(1,x,fmf);
  if (mret!=0) {
    O2SCL_ERR2("Solver failed in ",
	       "eos_had_base::fn0().",exc_efailed);
  }
  
  nb=x[0];
  
  calc_pressure_nb(nb);
  leoa=eos_thermo->ed/nb;
  
  if (neutron->inc_rest_mass) {
    leoa-=neutron->m/2.0;
  }
  if (proton->inc_rest_mass) {
    leoa-=proton->m/2.0;
  }
  
  return nb;
}

int eos_had_base::saturation() {
  std::cout << "sat 1" << std::endl;
  n0=fn0(0.0,eoa);
  if (n0<0.08 || n0>0.24) {
    O2SCL_CONV_RET("Function eos_had_base::saturation() found an ",
                   "unphysical saturation density.",
                   o2scl::exc_efailed,err_nonconv);
  }
  std::cout << "sat 2" << std::endl;
  comp=fcomp(n0);
  std::cout << "sat 3" << std::endl;
  esym=fesym(n0);
  std::cout << "sat 4" << std::endl;
  msom=fmsom(n0);
  std::cout << "sat 5" << std::endl;
  kprime=fkprime(n0);
  std::cout << "sat 6" << std::endl;

  return 0;
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
    
  funct t1fun=std::bind(std::mem_fn<double(double)>
			(&eos_had_base::t1_fun),
			this,std::placeholders::_1);
  funct t2fun=std::bind(std::mem_fn<double(double)>
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

double eos_had_base::calc_pressure_nb(double nb, double delta) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  std::cout << "H: " << neutron->n << " " << proton->n << std::endl;
  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->pr;
}

int eos_had_base::calc_pressure_nb_mroot(size_t nv, const ubvector &x,
                                         ubvector &y, double delta) {

  if ((1.0+delta)*x[0]/2.0<0.0 || (1.0-delta)*x[0]/2.0<0.0) {
    // Avoid negative densities
    return 1;
  }
  y[0]=calc_pressure_nb(x[0],delta);
  
  return 0;
}

void eos_had_base::const_pf_derivs(double nb, double pf, 
				   double &dednb_pf, double &dPdnb_pf) {

  // Take derivatives w.r.t. delta and then multiply by -2 to get
  // derivatives w.r.t. x
  funct fmpp=std::bind(std::mem_fn<double(double,double)>
		       (&eos_had_base::calc_pressure_nb),
		       this,std::placeholders::_1,1.0-2.0*pf);
  
  dPdnb_pf=-2.0*sat_deriv->deriv(nb,fmpp);

  // The other derivative can be obtained simply from the
  // chemical potentials
  neutron->n=nb*(1.0-pf);
  proton->n=nb*pf;
  calc_e(*neutron,*proton,*eos_thermo);
  dednb_pf=-2.0*(neutron->mu*(1.0-pf)+proton->mu*pf);
  
  return;
}

double eos_had_base::calc_press_over_den2(double nb, double delta) {
  
  neutron->n=nb/2.0;
  proton->n=nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->pr/nb/nb;
}

double eos_had_base::calc_edensity_delta(double delta, double nb) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->ed;
}

void eos_had_base::f_number_suscept(double mun, double mup, double &dPdnn, 
				    double &dPdnp, double &dPdpp) {

  // For (d^2 P)/(d mun d mun)
  funct fnn=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_nn_p),
     this,std::placeholders::_1,mup);
  dPdnn=sat_deriv->deriv(mun,fnn);
  
  // For (d^2 P)/(d mun d mup)
  funct fnp=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_nn_p),
     this,mun,std::placeholders::_1);
  dPdnp=sat_deriv->deriv(mup,fnn);

  // For (d^2 P)/(d mup d mup)
  funct fpp=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_np_p),
     this,mun,std::placeholders::_1);
  dPdpp=sat_deriv->deriv(mup,fpp);
  
  return;
}

void eos_had_base::f_inv_number_suscept(double nn, double np, double &dednn, 
					double &dednp, double &dedpp) {

  // For (d^2 ed)/(d mun d mun)
  funct fnn=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_mun_e),
     this,std::placeholders::_1,np);
  dednn=sat_deriv->deriv(nn,fnn);
  
  // For (d^2 ed)/(d mun d mup)
  funct fnp=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_mun_e),
     this,nn,std::placeholders::_1);
  dednp=sat_deriv->deriv(np,fnn);

  // For (d^2 ed)/(d mup d mup)
  funct fpp=std::bind
    (std::mem_fn<double(double,double)>(&eos_had_base::calc_mup_e),
     this,nn,std::placeholders::_1);
  dedpp=sat_deriv->deriv(np,fpp);
  
  return;
}

double eos_had_base::calc_mun_e(double nn, double np) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu;
}

double eos_had_base::calc_ed(double nn, double np) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return eos_thermo->ed;
}

double eos_had_base::calc_pr(double mun, double mup) {
  
  neutron->mu=mun;  
  proton->mu=mup;
  
  calc_p(*neutron,*proton,*eos_thermo);
  cout << "calc_pr(): " << mun << " " << mup << " " 
       << eos_thermo->pr << endl;

  return eos_thermo->pr;
}

double eos_had_base::calc_mup_e(double nn, double np) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu;
}

double eos_had_base::calc_nn_p(double mun, double mup) {
  
  neutron->mu=mun;  
  proton->mu=mup;
  
  calc_p(*neutron,*proton,*eos_thermo);

  return neutron->n;
}

double eos_had_base::calc_np_p(double mun, double mup) {
  
  neutron->n=mun;  
  proton->n=mup;
  
  calc_p(*neutron,*proton,*eos_thermo);

  return neutron->n;
}

double eos_had_base::calc_dmu_delta(double delta, double nb) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu-proton->mu;
}

double eos_had_base::calc_musum_delta(double delta, double nb) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;
  
  calc_e(*neutron,*proton,*eos_thermo);

  return neutron->mu+proton->mu;
}

double eos_had_base::calc_edensity_nb(double nb, double delta) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;

  calc_e(*neutron,*proton,*eos_thermo);
  
  return eos_thermo->ed;
}

int eos_had_base::nuc_matter_e(size_t nv, const ubvector &x, 
			       ubvector &y, double mun0, double mup0) {

  neutron->n=x[0];
  proton->n=x[1];
  
  calc_e(*neutron,*proton,*eos_thermo);

  y[0]=neutron->mu-mun0;
  y[1]=proton->mu-mup0;
  
  if (!std::isfinite(neutron->mu) || !std::isfinite(proton->mu)) {
    return exc_efailed;
  }

  return 0;
}

int eos_had_base::nuc_matter_p(size_t nv, const ubvector &x, 
			       ubvector &y, double nn0, double np0) {
  
  neutron->mu=x[0];
  proton->mu=x[1];

  calc_p(*neutron,*proton,*eos_thermo);

  y[0]=neutron->n-nn0;
  y[1]=proton->n-np0;

  if (!std::isfinite(neutron->n) || !std::isfinite(proton->n)) {
    return exc_efailed;
  }
  
  return 0;
}

void eos_had_base::set_n_and_p(fermion &n, fermion &p) {
  neutron=&n;
  proton=&p;
  return;
}

void eos_had_base::set_mroot(mroot<mm_funct,
			     boost::numeric::ublas::vector<double>, 
			     jac_funct> &mr) {
  eos_mroot=&mr;
  return;
}

void eos_had_base::set_sat_mroot(mroot<> &mr) {
  sat_mroot=&mr;
  return;
}

void eos_had_base::set_sat_deriv(deriv_base<funct> &de) {
  sat_deriv=&de;
  return;
}

void eos_had_base::set_sat_deriv2(deriv_base<funct> &de) {
  sat_deriv2=&de;
  return;
}

void eos_had_base::check_mu(fermion &n, fermion &p, thermo &th,
			    double &mun_deriv, double &mup_deriv,
			    double &mun_err, double &mup_err) {

  set_n_and_p(n,p);
  set_thermo(th);
  double nn=n.n;
  double np=p.n;

  funct fn=std::bind
    (std::mem_fn<double(double,double)>
     (&eos_had_base::calc_ed),this,std::placeholders::_1,p.n);
  sat_deriv->deriv_err(n.n,fn,mun_deriv,mun_err);

  n.n=nn;
  p.n=np;

  funct fp=std::bind
    (std::mem_fn<double(double,double)>
     (&eos_had_base::calc_ed),this,n.n,std::placeholders::_1);
  sat_deriv->deriv_err(p.n,fp,mup_deriv,mup_err);

  n.n=nn;
  p.n=np;

  calc_e(n,p,th);

  return;
}

void eos_had_base::check_den(fermion &n, fermion &p, thermo &th,
			     double &nn_deriv, double &np_deriv,
			     double &nn_err, double &np_err) {

  set_n_and_p(n,p);
  set_thermo(th);
  double mun=n.mu;
  double mup=p.mu;

  funct fn=std::bind
    (std::mem_fn<double(double,double)>
     (&eos_had_base::calc_pr),this,std::placeholders::_1,p.mu);
  sat_deriv->deriv_err(n.mu,fn,nn_deriv,nn_err);

  n.mu=mun;
  p.mu=mup;

  funct fp=std::bind
    (std::mem_fn<double(double,double)>
     (&eos_had_base::calc_pr),this,n.mu,std::placeholders::_1);
  sat_deriv->deriv_err(p.mu,fp,np_deriv,np_err);

  n.mu=mun;
  p.mu=mup;

  calc_p(n,p,th);

  return;
}

void eos_had_temp_base::check_en(fermion &n, fermion &p, double T,
				 thermo &th,
				 double &en_deriv, double &en_err) {

  set_n_and_p(n,p);
  set_thermo(th);

  funct fn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_fr),this,n.n,p.n,std::placeholders::_1);
  sat_deriv->deriv_err(T,fn,en_deriv,en_err);
  en_deriv*=-1.0;
  
  calc_temp_e(n,p,T,th);

  return;
}

void eos_had_temp_base::check_mu_T(fermion &n, fermion &p, double T, thermo &th,
				   double &mun_deriv, double &mup_deriv,
				   double &mun_err, double &mup_err) {

  set_n_and_p(n,p);
  set_thermo(th);
  double nn=n.n;
  double np=p.n;

  // In order to help ensure the derivatives don't cause
  // negative densities
  def_deriv.h=n.n/10.0;
  
  funct fn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_fr),this,std::placeholders::_1,p.n,T);
  sat_deriv->deriv_err(n.n,fn,mun_deriv,mun_err);

  n.n=nn;
  p.n=np;

  // In order to help ensure the derivatives don't cause
  // negative densities
  def_deriv.h=p.n/10.0;
  
  funct fp=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_fr),this,n.n,std::placeholders::_1,T);
  sat_deriv->deriv_err(p.n,fp,mup_deriv,mup_err);
  
  n.n=nn;
  p.n=np;

  calc_temp_e(n,p,T,th);

  return;
}

int eos_had_eden_base::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;
    
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &, ubvector &, double, double)>
     (&eos_had_eden_base::nuc_matter_e),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.mu,p.mu);
#ifdef O2SCL_NEVER_DEFINED
}{
#endif
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
    
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &, ubvector &, double, double)>
     (&eos_had_eden_base::nuc_matter_p),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.n,p.n);
  
#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  
  eos_mroot->msolve(2,mu,fmf);

  th=*eos_thermo;

  return 0;
}

double eos_had_temp_base::calc_fr(double nn, double np, double T) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_temp_e(*neutron,*proton,T,*eos_thermo);
  
  return eos_thermo->ed-T*eos_thermo->en;
}

int eos_had_temp_base::calc_liqgas_dens_temp_e
(fermion &n1, fermion &p1, fermion &n2, fermion &p2,
 double T, thermo &th1, thermo &th2) {
  
  ubvector x(3);
  x[0]=p1.n;
  x[1]=n2.n;
  x[2]=p2.n;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t, const ubvector &, ubvector &, fermion &,
		     fermion &, fermion &, fermion &, 
		     double, thermo &, thermo &)>
     (&eos_had_temp_base::liqgas_dens_solve),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(n1),std::ref(p1),
     std::ref(n2),std::ref(p2),T,std::ref(th1),std::ref(th2));

#ifdef O2SCL_NEVER_DEFINED
}{
#endif

  int ret=eos_mroot->msolve(3,x,fmf);

  p1.n=x[0];
  n2.n=x[1];
  p2.n=x[2];
      
  return ret;
}

int eos_had_temp_base::liqgas_dens_solve(size_t nv, const ubvector &x, 
					 ubvector &y, fermion &n1, fermion &p1,
					 fermion &n2, fermion &p2, double T,
					 thermo &th1, thermo &th2) {
      
  p1.n=x[0];
  calc_temp_e(n1,p1,T,th1);

  n2.n=x[1];
  p2.n=x[2];
  calc_temp_e(n2,p2,T,th2);

  y[0]=n1.mu-n2.mu;
  y[1]=p1.mu-p2.mu;
  y[2]=th1.pr-th2.pr;

  return 0;
}

int eos_had_temp_base::liqgas_solve(size_t nv, const ubvector &x, 
				    ubvector &y, fermion &n1, fermion &p1,
				    fermion &n2, fermion &p2, double nB0,
				    double Ye0, double T, 
				    thermo &th1, thermo &th2) {
      
  n1.n=x[0];
  p1.n=x[1];
  calc_temp_e(n1,p1,T,th1);

  n2.n=x[2];
  p2.n=x[3];
  calc_temp_e(n2,p2,T,th2);

  double chi=x[4];

  y[0]=n1.mu-n2.mu;
  y[1]=p1.mu-p2.mu;
  y[2]=th1.pr-th2.pr;
  y[3]=(n1.n+p1.n)*chi+(n2.n+p2.n)*(1.0-chi)-nB0;
  y[4]=p1.n*chi+p2.n*(1.0-chi)-Ye0*nB0;

  return 0;
}

int eos_had_temp_base::liqgas_beta_solve(size_t nv, const ubvector &x, 
					 ubvector &y, fermion &n1, fermion &p1,
					 fermion &n2, fermion &p2, 
					 double nB0, double T, 
					 thermo &th1, thermo &th2, fermion &e) {
      
  n1.n=x[0];
  p1.n=x[1];
  calc_temp_e(n1,p1,T,th1);

  n2.n=x[2];
  p2.n=x[3];
  calc_temp_e(n2,p2,T,th2);

  double chi=x[4];

  e.n=p1.n*chi+p2.n*(1.0-chi);
  fet->calc_density(e,T);
      
  y[0]=n1.mu-n2.mu;
  y[1]=p1.mu-p2.mu;
  y[2]=th1.pr-th2.pr;
  y[3]=(n1.n+p1.n)*chi+(n2.n+p2.n)*(1.0-chi)-nB0;
  y[4]=n1.mu-p1.mu-e.mu;

  return 0;
}

int eos_had_temp_base::calc_liqgas_temp_e
(fermion &n1, fermion &p1, fermion &n2, fermion &p2,
 double nB, double Ye, double T, thermo &th1, thermo &th2,
 double &chi) {

  ubvector x(5);
  x[0]=n1.n;
  x[1]=p1.n;
  x[2]=n2.n;
  x[3]=p2.n;
  x[4]=chi;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t, const ubvector &, ubvector &, fermion &,
		     fermion &, fermion &, fermion &, 
		     double, double, double, thermo &, thermo &)>
     (&eos_had_temp_base::liqgas_solve),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(n1),std::ref(p1),
     std::ref(n2),std::ref(p2),nB,Ye,T,std::ref(th1),std::ref(th2));
#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  int ret=eos_mroot->msolve(5,x,fmf);

  n1.n=x[0];
  p1.n=x[1];
  n2.n=x[2];
  p2.n=x[3];
  chi=x[4];
      
  return ret;
}

int eos_had_temp_base::calc_liqgas_beta_temp_e
(fermion &n1, fermion &p1, fermion &n2, fermion &p2,
 double nB, double T, thermo &th1, thermo &th2,
 double &Ye, double &chi) {
      
  fermion electron(o2scl_settings.get_convert_units().convert
		   ("kg","1/fm",o2scl_mks::mass_electron),2.0);

  ubvector x(5);
  x[0]=n1.n;
  x[1]=p1.n;
  x[2]=n2.n;
  x[3]=p2.n;
  x[4]=chi;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t, const ubvector &, ubvector &, fermion &,
		     fermion &, fermion &, fermion &, 
		     double, double, thermo &, thermo &, fermion &)>
     (&eos_had_temp_base::liqgas_beta_solve),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(n1),std::ref(p1),
     std::ref(n2),std::ref(p2),nB,T,std::ref(th1),std::ref(th2),
     std::ref(electron));
#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  int ret=eos_mroot->msolve(5,x,fmf);

  n1.n=x[0];
  p1.n=x[1];
  n2.n=x[2];
  p2.n=x[3];
  chi=x[4];
      
  Ye=(p1.n*chi+p2.n*(1.0-chi))/nB;
      
  return ret;
}

double eos_had_temp_base::calc_temp_mun_e(double nn, double np, double T) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_temp_e(*neutron,*proton,T,*eos_thermo);

  return neutron->mu;
}

double eos_had_temp_base::calc_temp_mup_e(double nn, double np, double T) {
  
  neutron->n=nn;  
  proton->n=np;
  
  calc_temp_e(*neutron,*proton,T,*eos_thermo);

  return neutron->mu;
}

double eos_had_temp_base::calc_temp_nn_p(double mun, double mup, double T) {
  
  neutron->mu=mun;  
  proton->mu=mup;
  
  calc_temp_p(*neutron,*proton,T,*eos_thermo);

  return neutron->n;
}

double eos_had_temp_base::calc_temp_np_p(double mun, double mup, double T) {
  
  neutron->n=mun;  
  proton->n=mup;
  
  calc_temp_p(*neutron,*proton,T,*eos_thermo);

  return neutron->n;
}

void eos_had_temp_base::f_number_suscept_T
(double mun, double mup, double T, double &dPdnn, double &dPdnp, 
 double &dPdpp) {

  // For (d^2 P)/(d mun d mun)
  funct fnn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_nn_p),
     this,std::placeholders::_1,mup,T);
  dPdnn=sat_deriv->deriv(mun,fnn);
  
  // For (d^2 P)/(d mun d mup)
  funct fnp=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_nn_p),
     this,mun,std::placeholders::_1,T);
  dPdnp=sat_deriv->deriv(mup,fnn);
  
  // For (d^2 P)/(d mup d mup)
  funct fpp=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_np_p),
     this,mun,std::placeholders::_1,T);
  dPdpp=sat_deriv->deriv(mup,fpp);
  
  return;
}

void eos_had_temp_base::f_inv_number_suscept_T
(double nn, double np, double T, double &dednn, double &dednp, double &dedpp) {

  // For (d^2 ed)/(d mun d mun)
  funct fnn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_mun_e),
     this,std::placeholders::_1,np,T);
  dednn=sat_deriv->deriv(nn,fnn);
  
  // For (d^2 ed)/(d mun d mup)
  funct fnp=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_mun_e),
     this,nn,std::placeholders::_1,T);
  dednp=sat_deriv->deriv(np,fnn);
  
  // For (d^2 ed)/(d mup d mup)
  funct fpp=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_temp_mup_e),
     this,nn,std::placeholders::_1,T);
  dedpp=sat_deriv->deriv(np,fpp);
  
  return;
}

int eos_had_temp_base::nuc_matter_temp_e(size_t nv, const ubvector &x, 
					 ubvector &y, double mun0, double mup0,
					 double T) {
  neutron->n=x[0];
  proton->n=x[1];
  
  if (x[0]<0.0 || x[1]<0.0) {
    return exc_ebadfunc;
  }

  if (!std::isfinite(neutron->n) || !std::isfinite(proton->n)) {
    O2SCL_ERR2("Density problem in ",
	       "eos_had_temp_base::nuc_matter_temp_e().",exc_esanity);
  }
  int ret=calc_temp_e(*neutron,*proton,T,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR2("Function calc_e() failed in ",
	       "eos_had_temp_base::nuc_matter_temp_e().",exc_efailed);
  }

  y[0]=(neutron->mu-mun0)/mun0;
  y[1]=(proton->mu-mup0)/mup0;
  
  if (!std::isfinite(neutron->mu) || !std::isfinite(proton->mu)) {
    O2SCL_ERR2("Chemical potential problem in ",
	       "eos_had_temp_base::nuc_matter_temp_e().",exc_esanity);
  }

  return ret;
}

int eos_had_temp_base::nuc_matter_temp_p(size_t nv, const ubvector &x, 
					 ubvector &y, double nn0,
					 double np0, double T) {
  
  neutron->mu=x[0];
  proton->mu=x[1];

  int ret=calc_temp_p(*neutron,*proton,T,*eos_thermo);
  if (ret!=0) {
    O2SCL_ERR("calc_p() failed in eos_had_temp_base::nuc_matter_p().",ret);
  }

  y[0]=neutron->n-nn0;
  y[1]=proton->n-np0;
  
  return ret;
}

double eos_had_temp_base::calc_entropy_delta(double delta, double nb, 
					     double T) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;
  
  calc_temp_e(*neutron,*proton,T,*eos_thermo);
  
  return eos_thermo->en;
}    

double eos_had_temp_base::calc_dmu_delta_T(double delta, double nb, 
					   double T) {
  
  neutron->n=(1.0+delta)*nb/2.0;
  proton->n=(1.0-delta)*nb/2.0;
  
  calc_temp_e(*neutron,*proton,T,*eos_thermo);
  
  return neutron->mu-proton->mu;
}

double eos_had_temp_base::fesym_T(double nb, double T, double delta) {
  
  funct fmn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_dmu_delta_T),this,std::placeholders::_1,
     nb,T);
  
  return sat_deriv->deriv(delta,fmn)/4.0;
}

double eos_had_temp_base::fsyment_T(double nb, double T, double delta) {
  
  funct fmn=std::bind
    (std::mem_fn<double(double,double,double)>
     (&eos_had_temp_base::calc_entropy_delta),this,std::placeholders::_1,
     nb,T);
  
  return sat_deriv->deriv2(delta,fmn)/2.0/nb;
}

int eos_had_temp_eden_base::calc_p(fermion &n, fermion &p, thermo &th) {
  int ret;
  
  set_n_and_p(n,p);
  set_thermo(th);
    
  ubvector x(2);
  x[0]=n.n;
  x[1]=p.n;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &, double, double)>
     (&eos_had_eden_base::nuc_matter_e),
     
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.mu,p.mu);
  
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
  
  // Replace a bad guess if necessary
  if (n.n<=0.0) n.n=0.08;
  if (p.n<=0.0) p.n=0.08;

  ubvector den(2);
  den[0]=n.n;
  den[1]=p.n;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double, double, double)>
     (&eos_had_temp_eden_base::nuc_matter_temp_e),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.mu,p.mu,T);

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  
  ret=eos_mroot->msolve(2,den,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2("Solver failed in ",
	       "eos_had_temp_eden_base::calc_temp_p().",ret);
  }
  
  th=*eos_thermo;
  
  return ret;
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
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &, double, double)>
     (&eos_had_eden_base::nuc_matter_p),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.n,p.n);
  
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
  
  ubvector mu(2);
  mu[0]=n.mu;
  mu[1]=p.mu;

  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,
		     ubvector &, double, double, double)>
     (&eos_had_temp_pres_base::nuc_matter_temp_p),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,n.n,p.n,T);

#ifdef O2SCL_NEVER_DEFINED
}{
#endif

  ret=eos_mroot->msolve(2,mu,fmf);
  
  if (ret!=0) {
    O2SCL_ERR2("Solver failed in ",
	       "eos_had_temp_pres_base::calc_temp_e().",ret);
  }
  
  th=*eos_thermo;

  return 0;
}

int eos_had_base::beta_eq_T0(ubvector &nB_grid, ubvector &guess,
                             eos_leptons &elep,
			     //fermion &e, bool include_muons,
			     //fermion &mu, fermion_rel &frel,
			     std::shared_ptr<table_units<> > results) {
  
  // Ensure initial guess is valid
  if (guess[0]<=0.0 || guess[0]>=nB_grid[0]) guess[0]=nB_grid[0]/2.0;

  results->clear();
  results->line_of_names("ed pr nb nn np mun mup kfn kfp");
  results->line_of_units(((std::string)"1/fm^4 1/fm^4 1/fm^3 1/fm^3 ")+
			 "1/fm^3 1/fm 1/fm 1/fm 1/fm");
  
  for(size_t i=0;i<nB_grid.size();i++) {
    
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &, ubvector &, 
                       const double &, eos_leptons &)>
       (&eos_had_base::solve_beta_eq_T0),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,std::cref(nB_grid[i]),std::ref(elep));
    
    int iret=beta_mroot.msolve(1,guess,fmf);
    if (iret!=0) return iret;
	
    // Final function evaluation to make sure, e.g.
    // eos_thermo object is correct
    ubvector y(1);
    fmf(1,guess,y);
	
    std::vector<double> line={eos_thermo->ed,eos_thermo->pr,nB_grid[i],
			      neutron->n,proton->n,
			      neutron->mu,proton->mu,
			      neutron->kf,proton->kf};
    results->line_of_data(line);
  }

  // Store the initial guess for the first density in the guess vector
  // since we used the guess vector as storage space above
  guess[0]=results->get("np",0);
      
  return 0;
}

int eos_had_base::solve_beta_eq_T0(size_t nv, const ubvector &x,
				   ubvector &y, const double &nB,
                                   eos_leptons &elep) {
  
  if (x[0]<0.0) return 1;
  double n_charge=x[0];
  proton->n=n_charge;
  neutron->n=nB-n_charge;
  if (neutron->n<0.0) return 2;
  this->calc_e(*neutron,*proton,*eos_thermo);
  elep.e.mu=neutron->mu-proton->mu;
  elep.pair_mu(0.0);
  y[0]=n_charge-elep.e.n;
  if (elep.include_muons) {
    y[0]=n_charge-elep.e.n-elep.mu.n;
  }
  return 0;
}
