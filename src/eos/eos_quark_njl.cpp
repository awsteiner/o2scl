/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#include <o2scl/eos_quark_njl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_quark_njl::eos_quark_njl() {

  up_default_mass=5.5/hc_mev_fm;
  down_default_mass=5.5/hc_mev_fm;
  strange_default_mass=140.7/hc_mev_fm;

  def_up.non_interacting=true;
  def_down.non_interacting=true;
  def_strange.non_interacting=true;
  def_up.init(up_default_mass,6.0);
  def_down.init(down_default_mass,6.0);
  def_strange.init(strange_default_mass,6.0);
  
  up=&def_up;
  down=&def_down;
  strange=&def_strange;

  B0=0.0;

  L=602.3/hc_mev_fm;
  G=1.835/L/L;
  K=12.36/pow(L,5.0);
  
  solver=&def_solver;
  it=&def_it;

  // We don't call set_parameters() here, because we can't call a
  // virtual member function from the constructor
  limit=20.0;

  fromqq=true;
}

int eos_quark_njl::set_quarks(quark &u, quark &d, quark &s) {
  up=&u;
  down=&d;
  strange=&s;

  return 0;
}

int eos_quark_njl::set_parameters(double lambda, double fourferm, 
				  double sixferm) {
  
  if (lambda!=0.0) L=lambda;
  else L=602.3/hc_mev_fm;
  if (fourferm!=0.0) G=fourferm;
  else G=1.835/L/L;
  if (sixferm!=0.0) K=sixferm;
  else K=12.36/pow(L,5.0);

  ubvector bx(3);
  bool fromqqold=fromqq;
  fromqq=false;

  // Solve zero-density gap-equations:
  bx[0]=1.0;
  bx[1]=1.0;
  bx[2]=2.0;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_quark_njl::B0fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  solver->msolve(3,bx,fmf);
  
  // Make the appropriate correction
  B0+=eos_thermo->pr;

  // Return fromqq to its original value
  fromqq=fromqqold;
  return 0;
}

int eos_quark_njl::calc_p(quark &u, quark &d, quark &s, thermo &th) {
  ubvector x(3);
  int ret;

  up=&u;
  down=&d;
  strange=&s;
  eos_thermo=&th;
  
  if (fromqq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::gapfunqq),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    ret=solver->msolve(3,x,fmf);
    
  } else {

    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::gapfunms),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    ret=solver->msolve(3,x,fmf);

  }

  return ret;
}

int eos_quark_njl::calc_temp_p(quark &u, quark &d, quark &s, 
			       double T, thermo &th) {
  ubvector x(3);
  int vp=0;
  int ret;

  up=&u;
  down=&d;
  strange=&s;
  eos_thermo=&th;
  cp_temp=T;
  
  if (fromqq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::gapfunqqT),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    ret=solver->msolve(3,x,fmf);
    
  } else {

    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::gapfunmsT),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    ret=solver->msolve(3,x,fmf);

  }

  double gap1, gap2, gap3;
  
  return ret;
}

int eos_quark_njl::calc_eq_p(quark &u, quark &d, quark &s, double &gap1,
			     double &gap2, double &gap3, thermo &th) {
  
  if (fromqq==true) {
    
    //---------------------------------------------------------------
    // Calculate everything from the quark condensates and the
    // chemical potentials, then gap1, gap2, and gap3 contain
    // Eq. 3 from Buballa99

    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
    
    if (u.mu>u.ms) u.n=1.0/pi2*pow(u.mu*u.mu-u.ms*u.ms,1.5);
    else u.n=0.0;
    if (d.mu>d.ms) d.n=1.0/pi2*pow(d.mu*d.mu-d.ms*d.ms,1.5);
    else d.n=0.0;
    if (s.mu>s.ms) s.n=1.0/pi2*pow(s.mu*s.mu-s.ms*s.ms,1.5);
    else s.n=0.0;

    fet.kf_from_density(u);
    fet.kf_from_density(d);
    fet.kf_from_density(s);
  
    gap1=3.0/pi2*(-u.ms/2.0*L*sqrt(u.ms*u.ms+L*L)+pow(u.ms,3.0)/2.0*
		  log(L+sqrt(L*L+u.ms*u.ms)));
    gap1+=3.0/pi2*(u.ms/2.0*u.kf*sqrt(u.ms*u.ms+u.kf*u.kf)-
		   pow(u.ms,3.0)/2.0*
		   log(u.kf+sqrt(u.kf*u.kf+u.ms*u.ms)));
    gap1-=u.qq;
    gap2=3.0/pi2*(-d.ms/2.0*L*sqrt(d.ms*d.ms+L*L)+pow(d.ms,3.0)/2.0*
		  log(L+sqrt(L*L+d.ms*d.ms)));
    gap2+=3.0/pi2*(d.ms/2.0*d.kf*sqrt(d.ms*d.ms+d.kf*d.kf)-
		   pow(d.ms,3.0)/2.0*
		   log(d.kf+sqrt(d.kf*d.kf+d.ms*d.ms)));
    gap2-=d.qq;
    gap3=3.0/pi2*(-s.ms/2.0*L*sqrt(s.ms*s.ms+L*L)+pow(s.ms,3.0)/2.0*
		  log(L+sqrt(L*L+s.ms*s.ms)));
    gap3+=3.0/pi2*(s.ms/2.0*s.kf*sqrt(s.ms*s.ms+s.kf*s.kf)-
		   pow(s.ms,3.0)/2.0*
		   log(s.kf+sqrt(s.kf*s.kf+s.ms*s.ms)));
    gap3-=s.qq;

  } else {

    //---------------------------------------------------------------
    // Calculate everything from the dynamical masses and the chemical
    // potentials, and then gap1, gap2, and gap3 contain Eq. 2 from
    // Buballa99.

    if (u.mu>u.ms) u.kf=sqrt(u.mu*u.mu-u.ms*u.ms);
    else u.kf=0.0;
    if (d.mu>d.ms) d.kf=sqrt(d.mu*d.mu-d.ms*d.ms);
    else d.kf=0.0;
    if (s.mu>s.ms) s.kf=sqrt(s.mu*s.mu-s.ms*s.ms);
    else s.kf=0.0;

    u.qq=3.0/pi2*(-u.ms/2.0*L*sqrt(u.ms*u.ms+L*L)+pow(u.ms,3.0)/2.0*
		  log(L+sqrt(L*L+u.ms*u.ms)));
    u.qq+=3.0/pi2*(u.ms/2.0*u.kf*sqrt(u.ms*u.ms+u.kf*u.kf)-
		   pow(u.ms,3.0)/2.0*
		   log(u.kf+sqrt(u.kf*u.kf+u.ms*u.ms)));
    d.qq=3.0/pi2*(-d.ms/2.0*L*sqrt(d.ms*d.ms+L*L)+pow(d.ms,3.0)/2.0*
		  log(L+sqrt(L*L+d.ms*d.ms)));
    d.qq+=3.0/pi2*(d.ms/2.0*d.kf*sqrt(d.ms*d.ms+d.kf*d.kf)-
		   pow(d.ms,3.0)/2.0*
		   log(d.kf+sqrt(d.kf*d.kf+d.ms*d.ms)));
    s.qq=3.0/pi2*(-s.ms/2.0*L*sqrt(s.ms*s.ms+L*L)+pow(s.ms,3.0)/2.0*
		  log(L+sqrt(L*L+s.ms*s.ms)));
    s.qq+=3.0/pi2*(s.ms/2.0*s.kf*sqrt(s.ms*s.ms+s.kf*s.kf)-
		   pow(s.ms,3.0)/2.0*
		   log(s.kf+sqrt(s.kf*s.kf+s.ms*s.ms)));

    if (u.mu>u.ms) u.n=1.0/pi2*pow(u.mu*u.mu-u.ms*u.ms,1.5);
    else u.n=0.0;
    if (d.mu>d.ms) d.n=1.0/pi2*pow(d.mu*d.mu-d.ms*d.ms,1.5);
    else d.n=0.0;
    if (s.mu>s.ms) s.n=1.0/pi2*pow(s.mu*s.mu-s.ms*s.ms,1.5);
    else s.n=0.0;

    gap1=-u.ms+u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    gap2=-d.ms+d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    gap3=-s.ms+s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;

  }    
  
  fet.energy_density_zerot(u);
  fet.energy_density_zerot(d);
  fet.energy_density_zerot(s);
  
  njbag(u);
  njbag(d);
  njbag(s);
  
  u.ed-=u.B;
  d.ed-=d.B;
  s.ed-=s.B;
  
  u.pr=-u.ed+u.n*u.mu;
  d.pr=-d.ed+d.n*d.mu;
  s.pr=-s.ed+s.n*s.mu;
  
  th.ed=u.ed+d.ed+s.ed+
    2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)+B0-4.0*K*u.qq*d.qq*s.qq;
  th.pr=-th.ed+u.n*u.mu+d.n*d.mu+s.n*s.mu;
  th.en=0.0;

  return 0;
}

int eos_quark_njl::calc_eq_e(quark &u, quark &d, quark &s, double &gap1,
			     double &gap2, double &gap3, thermo &th) {
  
  if (fromqq==true) {

    //--------------------------------------------
    // Calculate everything from the quark 
    // condensates and the chemical potentials

    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;

    fet.kf_from_density(u);
    fet.kf_from_density(d);
    fet.kf_from_density(s);
  
    u.mu=sqrt(u.kf*u.kf+u.ms*u.ms);
    d.mu=sqrt(d.kf*d.kf+d.ms*d.ms);
    s.mu=sqrt(s.kf*s.kf+s.ms*s.ms);
    
    gap1=3.0/pi2*(-u.ms/2.0*L*sqrt(u.ms*u.ms+L*L)+pow(u.ms,3.0)/2.0*
		  log(L+sqrt(L*L+u.ms*u.ms)));
    gap1+=3.0/pi2*(u.ms/2.0*u.kf*sqrt(u.ms*u.ms+u.kf*u.kf)-
		   pow(u.ms,3.0)/2.0*
		   log(u.kf+sqrt(u.kf*u.kf+u.ms*u.ms)));
    gap1-=u.qq;
    gap2=3.0/pi2*(-d.ms/2.0*L*sqrt(d.ms*d.ms+L*L)+pow(d.ms,3.0)/2.0*
		  log(L+sqrt(L*L+d.ms*d.ms)));
    gap2+=3.0/pi2*(d.ms/2.0*d.kf*sqrt(d.ms*d.ms+d.kf*d.kf)-
		   pow(d.ms,3.0)/2.0*
		   log(d.kf+sqrt(d.kf*d.kf+d.ms*d.ms)));
    gap2-=d.qq;
    gap3=3.0/pi2*(-s.ms/2.0*L*sqrt(s.ms*s.ms+L*L)+pow(s.ms,3.0)/2.0*
		  log(L+sqrt(L*L+s.ms*s.ms)));
    gap3+=3.0/pi2*(s.ms/2.0*s.kf*sqrt(s.ms*s.ms+s.kf*s.kf)-
		   pow(s.ms,3.0)/2.0*
		   log(s.kf+sqrt(s.kf*s.kf+s.ms*s.ms)));
    gap3-=s.qq;

  } else {

    //--------------------------------------------
    // Calculate everything from the dynamical 
    // masses and the chemical potentials

    fet.kf_from_density(u);
    fet.kf_from_density(d);
    fet.kf_from_density(s);
  
    u.mu=sqrt(u.kf*u.kf+u.ms*u.ms);
    d.mu=sqrt(d.kf*d.kf+d.ms*d.ms);
    s.mu=sqrt(s.kf*s.kf+s.ms*s.ms);

    u.qq=3.0/pi2*(-u.ms/2.0*L*sqrt(u.ms*u.ms+L*L)+pow(u.ms,3.0)/2.0*
		  log(L+sqrt(L*L+u.ms*u.ms)));
    u.qq+=3.0/pi2*(u.ms/2.0*u.kf*sqrt(u.ms*u.ms+u.kf*u.kf)-
		   pow(u.ms,3.0)/2.0*
		   log(u.kf+sqrt(u.kf*u.kf+u.ms*u.ms)));
    d.qq=3.0/pi2*(-d.ms/2.0*L*sqrt(d.ms*d.ms+L*L)+pow(d.ms,3.0)/2.0*
		  log(L+sqrt(L*L+d.ms*d.ms)));
    d.qq+=3.0/pi2*(d.ms/2.0*d.kf*sqrt(d.ms*d.ms+d.kf*d.kf)-
		   pow(d.ms,3.0)/2.0*
		   log(d.kf+sqrt(d.kf*d.kf+d.ms*d.ms)));
    s.qq=3.0/pi2*(-s.ms/2.0*L*sqrt(s.ms*s.ms+L*L)+pow(s.ms,3.0)/2.0*
		  log(L+sqrt(L*L+s.ms*s.ms)));
    s.qq+=3.0/pi2*(s.ms/2.0*s.kf*sqrt(s.ms*s.ms+s.kf*s.kf)-
		   pow(s.ms,3.0)/2.0*
		   log(s.kf+sqrt(s.kf*s.kf+s.ms*s.ms)));

    gap1=-u.ms+u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    gap2=-d.ms+d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    gap3=-s.ms+s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;

  }    

  fet.energy_density_zerot(u);
  fet.energy_density_zerot(d);
  fet.energy_density_zerot(s);

  njbag(u);
  njbag(d);
  njbag(s);

  u.ed-=u.B;
  d.ed-=d.B;
  s.ed-=s.B;
  
  u.pr=-u.ed+u.n*u.mu;
  d.pr=-d.ed+d.n*d.mu;
  s.pr=-s.ed+s.n*s.mu;
  
  th.ed=u.ed+d.ed+s.ed+
    2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)+B0-4.0*K*u.qq*d.qq*s.qq;
  th.pr=-th.ed+u.n*u.mu+d.n*d.mu+s.n*s.mu;
  th.en=0.0;

  return 0;
}

void eos_quark_njl::njbag(quark &pp) {

  //if (fabs(pp.ms)<1.0e-6) {
  //pp.B=0.0;
  //} else {
  //if (fabs(pp.ms)<1.0e-6) {
  pp.B=3.0/pi2*(1.0/4.0*L*pow(pp.ms*pp.ms+L*L,1.5)
                -1.0/8.0*pp.ms*pp.ms*L*sqrt(pp.ms*pp.ms+L*L)
                -1.0/8.0*pow(pp.ms,4.0)*log(L+sqrt(pp.ms*pp.ms+L*L))
                +1.0/16.0*pow(pp.ms,4.0)*log(pp.ms*pp.ms));
  //} else {
  //pp.B=3.0/pi2*(1.0/4.0*L*pow(pp.ms*pp.ms+L*L,1.5)
  //-1.0/8.0*pp.ms*pp.ms*L*sqrt(pp.ms*pp.ms+L*L)
  //-1.0/8.0*pow(pp.ms,4.0)*log(L+sqrt(pp.ms*pp.ms+L*L))
  //+1.0/16.0*pow(pp.ms,4.0)*log(pp.ms*pp.ms));
  //}
  //}

  return;
}

int eos_quark_njl::gapfunqq(size_t nv, const ubvector &x, ubvector &y) {

  double gap1,gap2,gap3;
  
  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  
  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) return 1;
  
  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,*eos_thermo);
  
  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;
  
  return 0;
}

int eos_quark_njl::gapfunms(size_t nv, const ubvector &x, ubvector &y) {

  double gap1,gap2,gap3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];

  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,*eos_thermo);

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;

  return 0;
}

int eos_quark_njl::gapfunqqT(size_t nv, const ubvector &x, ubvector &y) {

  double gap1,gap2,gap3;
  
  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  
  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) return 1;

  calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,*eos_thermo,cp_temp);
  
  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;
  
  return 0;
}

int eos_quark_njl::gapfunmsT(size_t nv, const ubvector &x, ubvector &y) {

  double gap1,gap2,gap3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];

  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,*eos_thermo,cp_temp);

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;

  return 0;
}

int eos_quark_njl::B0fun(size_t nv, const ubvector &x, ubvector &y) {
		      
  double gap1,gap2,gap3;
  
  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];
  
  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  up->mu=up->ms;
  down->mu=down->ms;
  strange->mu=strange->ms;
  
  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,*eos_thermo);

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;

  return 0;
}

int eos_quark_njl::calc_eq_temp_p(quark &u, quark &d, quark &s, 
				  double &gap1, double &gap2, double &gap3, 
				  thermo &qb, double temper) {
  double ierr;
  int iret1, iret2, iret3, iret4;

  if (temper<=0.0) {
    calc_eq_p(u,d,s,gap1,gap2,gap3,qb);
    return 0;
  }

  if (u.ms<=0.0) u.ms=1.e-9;
  if (d.ms<=0.0) d.ms=1.e-9;
  if (s.ms<=0.0) s.ms=1.e-9;

  njtp pa;
  pa.temper=temper;
  pa.limit=limit;

  // -----------------------------------------------------------------
  // Some of these integrals (iqq, and ide) converge better when they
  // are non-zero, so we add 1 to the integrand and subtract off the
  // contribution after the fact. This may cause inaccuracies when the
  // densities are small.

  if (fromqq==true) {
    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    pa.ms=u.ms;
    pa.mu=u.mu;
    pa.m=u.m;
    iret1=0;
  } else {
    pa.ms=u.ms;
    pa.mu=u.mu;
    pa.m=u.m;
    funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::iqq),
			this,std::placeholders::_1,pa);
    iret1=it->integ_err(fqq,0.0,L,u.qq,ierr);
    u.qq-=L;
  }

  {
    funct fde=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ide),
			this,std::placeholders::_1,pa);
    funct fed=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ied),
			this,std::placeholders::_1,pa);
    funct fpr=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ipr),
			this,std::placeholders::_1,pa);
    iret2=it->integ_err(fde,0.0,L,u.n,ierr);
    u.n-=L;
    fet.kf_from_density(u);
    iret3=it->integ_err(fpr,0.0,L,u.pr,ierr);
    iret4=it->integ_err(fed,0.0,L,u.ed,ierr);
    if (iret1!=0 || iret2!=0 || iret3!=0 || iret4!=0) {
      O2SCL_ERR("Up quark failed in eos_quark_njl::calc_eq_temp_p().",
		exc_efailed);
    }
  }
  
  if (fromqq==true) {
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    pa.ms=d.ms;
    pa.mu=d.mu;
    pa.m=d.m;
    iret1=0;
  } else {
    pa.ms=d.ms;
    pa.mu=d.mu;
    pa.m=d.m;
    funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::iqq),
			this,std::placeholders::_1,pa);
    iret1=it->integ_err(fqq,0.0,L,d.qq,ierr);
    d.qq-=L;
  }
  {
    funct fde=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ide),
			this,std::placeholders::_1,pa);
    funct fed=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ied),
			this,std::placeholders::_1,pa);
    funct fpr=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ipr),
			this,std::placeholders::_1,pa);

    iret2=it->integ_err(fde,0.0,L,d.n,ierr);
    d.n-=L;
    fet.kf_from_density(d);
    iret3=it->integ_err(fpr,0.0,L,d.pr,ierr);
    iret4=it->integ_err(fed,0.0,L,d.ed,ierr);
    if (iret1!=0 || iret2!=0 || iret3!=0 || iret4!=0) {
      O2SCL_ERR("Down quark failed in eos_quark_njl::calc_eq_temp_p().",
		exc_efailed);
    }
  }
  
  if (fromqq==true) {
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
    pa.ms=s.ms;
    pa.mu=s.mu;
    pa.m=s.m;
    iret1=0;
  } else {
    pa.ms=s.ms;
    pa.mu=s.mu;
    pa.m=s.m;
    funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::iqq),
			this,std::placeholders::_1,pa);
    iret1=it->integ_err(fqq,0.0,L,s.qq,ierr);
    s.qq-=L;
  }
  {
    funct fde=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ide),
			this,std::placeholders::_1,pa);
    funct fed=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ied),
			this,std::placeholders::_1,pa);
    funct fpr=std::bind(std::mem_fn<double(double,const njtp &)>
			(&eos_quark_njl::ipr),
			this,std::placeholders::_1,pa);
    iret2=it->integ_err(fde,0.0,L,s.n,ierr);
    s.n-=L;
    fet.kf_from_density(s);
    iret3=it->integ_err(fpr,0.0,L,s.pr,ierr);
    iret4=it->integ_err(fed,0.0,L,s.ed,ierr);
    if (iret1!=0 || iret2!=0 || iret3!=0 || iret4!=0) {
      O2SCL_ERR("Strange quark failed in eos_quark_njl::calc_eq_temp_p().",
		exc_efailed);
    }
  }

  qb.ed=u.ed+d.ed+s.ed+2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)+
    B0-4.0*K*u.qq*d.qq*s.qq;
  qb.pr=u.pr+d.pr+s.pr-2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)-
    B0+4.0*K*u.qq*d.qq*s.qq;
  qb.en=(qb.ed+qb.pr-u.n*u.mu-d.n*d.mu-s.n*s.mu)/
    temper;
  
  if (fromqq==true) {
    double qqt;

    pa.ms=u.ms;
    pa.mu=u.mu;
    pa.m=u.m;
    {
      funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			  (&eos_quark_njl::iqq),
			  this,std::placeholders::_1,pa);
      iret1=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap1=-u.qq+qqt-L;
    
    pa.ms=d.ms;
    pa.mu=d.mu;
    pa.m=d.m;
    {
      funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			  (&eos_quark_njl::iqq),
			  this,std::placeholders::_1,pa);
      iret2=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap2=-d.qq+qqt-L;
    
    pa.ms=s.ms;
    pa.mu=s.mu;
    pa.m=s.m;
    {
      funct fqq=std::bind(std::mem_fn<double(double,const njtp &)>
			  (&eos_quark_njl::iqq),
			  this,std::placeholders::_1,pa);
      iret3=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap3=-s.qq+qqt-L;

    if (iret1!=0 || iret2!=0 || iret3!=0) {
      O2SCL_ERR("Strange quark failed in eos_quark_njl::calc_eq_temp_p().",
		exc_efailed);
    }
  } else {
    gap1=-u.ms+u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    gap2=-d.ms+d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    gap3=-s.ms+s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
  }

  return 0;
}

double eos_quark_njl::iqq(double x, const njtp &p) {
  double en, ret;

  en=sqrt(x*x+p.ms*p.ms);
  ret=1.0;
  ret-=fermi_function(en,p.mu,p.temper,p.limit)+
    fermi_function(en,-p.mu,p.temper,p.limit);
  ret*=-3.0*x*x/en/pi2*p.ms;
  ret+=1.0;
  return ret;
}

double eos_quark_njl::ide(double x, const njtp &p) {
  double en, ret;


  en=sqrt(x*x+p.ms*p.ms);
  ret=fermi_function(en,p.mu,p.temper,p.limit)-
    fermi_function(en,-p.mu,p.temper,p.limit);
  ret*=3.0*x*x/pi2;
  ret+=1.0;
  return ret;
}

double eos_quark_njl::ied(double x, const njtp &p) {
  double en, ret;

  en=sqrt(x*x+p.ms*p.ms);
  ret=-en;
  ret+=en*fermi_function(en,p.mu,p.temper,p.limit)+
    en*fermi_function(en,-p.mu,p.temper,p.limit);
  ret*=3.0/pi2*x*x;
  return ret;
}

double eos_quark_njl::ipr(double x, const njtp &p) {
  double en, ret;

  en=sqrt(x*x+p.ms*p.ms);
  ret=en;
  if ((p.mu-en)/p.temper>p.limit) {
    ret+=p.mu-en;
  } else if ((p.mu-en)/p.temper<-p.limit) {
    ret+=0.0;
  } else {
    ret+=p.temper*log(1+exp((p.mu-en)/p.temper));
  }
  if ((-en-p.mu)/p.temper>p.limit) {
    ret+=-en-p.mu;
  } else if ((-en-p.mu)/p.temper<-p.limit) {
    ret+=0.0;
  } else {
    ret+=p.temper*log(1+exp((-en-p.mu)/p.temper));
  }
  ret*=3.0/pi2*x*x;
  return ret;
}

double eos_quark_njl::f_therm_pot(double qqu, double qqd, double qqs,
                                  double msu, double msd, double mss,
                                  bool vac_terms) {
  
  double g1, g2, g3;
  
  up->qq=qqu;
  down->qq=qqd;
  strange->qq=qqs;
  up->ms=msu;
  down->ms=msd;
  strange->ms=mss;
  
  calc_eq_p(*up,*down,*strange,g1,g2,g3,*eos_thermo);
  if (vac_terms==false) {
    double ret=-eos_thermo->pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                  strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq;
    return ret;
  }
  return -eos_thermo->pr;
}

double eos_quark_njl::f_therm_pot_T(double qqu, double qqd, double qqs,
                                    double msu, double msd, double mss,
                                    double temper, bool vac_terms) {
  
  double g1, g2, g3;
  
  up->qq=qqu;
  down->qq=qqd;
  strange->qq=qqs;
  up->ms=msu;
  down->ms=msd;
  strange->ms=mss;
  
  calc_eq_temp_p(*up,*down,*strange,g1,g2,g3,*eos_thermo,temper);
  if (vac_terms==false) {
    double ret=-eos_thermo->pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                  strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq;
    return ret;
  }
  return -eos_thermo->pr;
}

int eos_quark_njl_vec::calc_eq_p(quark &tu, quark &td, quark &ts,
                                 double &gap1, double &gap2, double &gap3,
                                 double &vec1, double &vec2, double &vec3,
                                 thermo &th) {

  if (fromqq==true) {
    O2SCL_ERR("Unimplemented.",o2scl::exc_eunimpl);
  }
  
  /*
    if (tu.m<0.0) tu.m*=-1.0;
    if (td.m<0.0) td.m*=-1.0;
    if (ts.m<0.0) ts.m*=-1.0;
  */

  // Following Klimt 90 for now. We use the "nu" field to store what
  // is referred to in that paper as {\bar{\mu}}. Note that a
  // different notation for G and K is used in Buballa et al. 1999,
  // and we use that same notation for the non-vector terms as Buballa
  // et al. 1999 in order to be consistent with the eos_quark_njl
  // class. Then Klimt's "delta mu" is mu-nu

  // This function presumes that initial guesses have been given
  // for both the effective masses (ms) and the effective chemical
  // potentials (nu)

  if (tu.nu>tu.ms) {
    tu.n=1.0/pi2*pow(tu.nu*tu.nu-tu.ms*tu.ms,1.5);
  } else {
    tu.n=0.0;
  }
  if (td.nu>td.ms) {
    td.n=1.0/pi2*pow(td.nu*td.nu-td.ms*td.ms,1.5);
  } else {
    td.n=0.0;
  }
  if (ts.nu>ts.ms) {
    ts.n=1.0/pi2*pow(ts.nu*ts.nu-ts.ms*ts.ms,1.5);
  } else {
    ts.n=0.0;
  }
  
  fet.kf_from_density(tu);
  fet.kf_from_density(td);
  fet.kf_from_density(ts);
  
  tu.qq=3.0/pi2*(-tu.ms/2.0*L*sqrt(tu.ms*tu.ms+L*L)+pow(tu.ms,3.0)/2.0*
                 log(L+sqrt(L*L+tu.ms*tu.ms)));
  tu.qq+=3.0/pi2*(tu.ms/2.0*tu.kf*sqrt(tu.ms*tu.ms+tu.kf*tu.kf)-
                  pow(tu.ms,3.0)/2.0*
                  log(tu.kf+sqrt(tu.kf*tu.kf+tu.ms*tu.ms)));
  td.qq=3.0/pi2*(-td.ms/2.0*L*sqrt(td.ms*td.ms+L*L)+pow(td.ms,3.0)/2.0*
                 log(L+sqrt(L*L+td.ms*td.ms)));
  td.qq+=3.0/pi2*(td.ms/2.0*td.kf*sqrt(td.ms*td.ms+td.kf*td.kf)-
                  pow(td.ms,3.0)/2.0*
                  log(td.kf+sqrt(td.kf*td.kf+td.ms*td.ms)));
  ts.qq=3.0/pi2*(-ts.ms/2.0*L*sqrt(ts.ms*ts.ms+L*L)+pow(ts.ms,3.0)/2.0*
                 log(L+sqrt(L*L+ts.ms*ts.ms)));
  ts.qq+=3.0/pi2*(ts.ms/2.0*ts.kf*sqrt(ts.ms*ts.ms+ts.kf*ts.kf)-
                  pow(ts.ms,3.0)/2.0*
                  log(ts.kf+sqrt(ts.kf*ts.kf+ts.ms*ts.ms)));

  fet.energy_density_zerot(tu);
  fet.energy_density_zerot(td);
  fet.energy_density_zerot(ts);
  
  njvecbag(tu);
  njvecbag(td);
  njvecbag(ts);
  
  tu.ed-=tu.B;
  td.ed-=td.B;
  ts.ed-=ts.B;

  tu.pr=-tu.ed+tu.n*tu.mu;
  td.pr=-td.ed+td.n*td.mu;
  ts.pr=-ts.ed+ts.n*ts.mu;
  
  th.ed=tu.ed+td.ed+ts.ed+
    2.0*G*(tu.qq*tu.qq+td.qq*td.qq+ts.qq*ts.qq)+B0-
    4.0*K*tu.qq*td.qq*ts.qq-2.0*GV*tu.n*tu.n-2.0*GV*td.n*td.n-
    2.0*GV*ts.n*ts.n;
  th.pr=-th.ed+tu.n*tu.mu+td.n*td.mu+ts.n*ts.mu;
  th.en=0.0;
  
  gap1=-tu.ms+tu.m-4.0*G*tu.qq+2.0*K*td.qq*ts.qq;
  gap2=-td.ms+td.m-4.0*G*td.qq+2.0*K*tu.qq*ts.qq;
  gap3=-ts.ms+ts.m-4.0*G*ts.qq+2.0*K*td.qq*tu.qq;

  // Fix the relationship between \mu_R and \mu, following eq 8b in
  // Klimt 90
  vec1=tu.mu-4.0*GV*tu.n-tu.nu;
  vec2=td.mu-4.0*GV*td.n-td.nu;
  vec3=ts.mu-4.0*GV*ts.n-ts.nu;
  
  return 0;
}

void eos_quark_njl_vec::njvecbag(quark &q) {
  //static const double zero=1.0e-9;
  //if (fabs(q.ms)<zero) {
  //q.B=0.0;
  //} else {
  //if (fabs(q.ms)<zero) {
  q.B=3.0/pi2*(1.0/4.0*L*pow(q.ms*q.ms+L*L,1.5)
               -1.0/8.0*q.ms*q.ms*L*sqrt(q.ms*q.ms+L*L)
               -1.0/8.0*pow(q.ms,4.0)*log(L+sqrt(q.ms*q.ms+L*L))
               +1.0/16.0*pow(q.ms,4.0)*log(q.ms*q.ms));
  //} else {
  //q.B=3.0/pi2*(1.0/4.0*L*pow(q.ms*q.ms+L*L,1.5)
  //-1.0/8.0*q.ms*q.ms*L*sqrt(q.ms*q.ms+L*L)
  //-1.0/8.0*pow(q.ms,4.0)*log(L+sqrt(q.ms*q.ms+L*L))
  //+1.0/16.0*pow(q.ms,4.0)*log(q.ms*q.ms));
  //}
  //}
  return;
}

/*
  #define LIMIT 20
  #define DEBUG 0
  
  double iqq(double x, int np, double *param);
  double ide(double x, int np, double *param);
  double ied(double x, int np, double *param);
  double ipr(double x, int np, double *param);
  //void njbag3(fermion &pp, double L, double G);
  */

int eos_quark_njl_vec::calc_eq_temp_p(quark &tu, quark &td, quark &ts,
                                      double &gap1, double &gap2,
                                      double &gap3, thermo &th,
                                      double temper) {
  /*
    double B0, double temper, int *check,
    double *energy, double *press, double *entropy,
    double *charge, double *barn, double *gap1, double *gap2,
    double *gap3, double *vec1, double *vec2, double *vec3) {
  */

  double *iparam, test, test2;
  int inp;

  //inp=4;
  //iparam=dvector(1,inp);

  tu.m=fabs(tu.m);
  td.m=fabs(td.m);
  ts.m=fabs(ts.m);

  /*
    iparam[3]=temper;
    
    iparam[1]=tu.m;
    iparam[2]=tu.nu;
    iparam[4]=tu.m;
  */
  
  //if (DEBUG) printf("uqq\n");
  //tu.qq=qrombaws(iqq,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //if (DEBUG) printf("un\n");
  //tu.n=qrombaws(ide,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //tu.kf=fermimom(tu.n,3.0);
  //if (DEBUG) printf("upr\n");
  //tu.pr=qrombaws(ipr,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  //  if (DEBUG) printf("ued\n");
  //tu.ed=qrombaws(ied,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  tu.ed+=2.0*G*tu.qq*tu.qq-2.0*GV*tu.n*tu.n;
  tu.pr-=2.0*G*tu.qq*tu.qq-2.0*GV*tu.n*tu.n;

  /*
  iparam[1]=td.m;
  iparam[2]=td.nu;
  iparam[4]=td.m;
  */
  //if (DEBUG) printf("dqq\n");
  //td.qq=qrombaws(iqq,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //if (DEBUG) printf("dn\n");
  //td.n=qrombaws(ide,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //td.kf=fermimom(td.n,3.0);
  //if (DEBUG) printf("dpr\n");
  //td.pr=qrombaws(ipr,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  //if (DEBUG) printf("ded\n");
  //td.ed=qrombaws(ied,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  td.ed+=2.0*G*td.qq*td.qq+2.0*GV*td.n*td.n;
  td.pr-=2.0*G*td.qq*td.qq+2.0*GV*td.n*td.n;

  /*
  iparam[1]=ts.m;
  iparam[2]=ts.nu;
  iparam[4]=ts.m;
  */
  //if (DEBUG) printf("sqq\n");
  //ts.qq=qrombaws(iqq,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //if (DEBUG) printf("sn\n");
  //ts.n=qrombaws(ide,0.0,L,inp,iparam,1.0e-7,400,-1.0,0)-L;
  //ts.kf=fermimom(ts.n,3.0);
  //if (DEBUG) printf("spr\n");
  //ts.pr=qrombaws(ipr,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  //if (DEBUG) printf("sed\n");
  //ts.ed=qrombaws(ied,0.0,L,inp,iparam,1.0e-7,400,-1.0,0);
  ts.ed+=2.0*G*ts.qq*ts.qq-2.0*GV*ts.n*ts.n;
  ts.pr-=2.0*G*ts.qq*ts.qq-2.0*GV*ts.n*ts.n;

  /*
    njbag(tu,L,G);
    njbag(td,L,G);
    njbag(ts,L,G);
  */

  th.ed=tu.ed+td.ed+ts.ed+B0-4.0*K*tu.qq*td.qq*ts.qq;
  th.pr=tu.pr+td.pr+ts.pr-B0-10.0*K*tu.qq*td.qq*ts.qq;
  th.en=(th.ed+th.pr-tu.n*tu.mu-td.n*td.mu-ts.n*ts.mu)/temper;
  
  //*charge=(tu.n*2.0-td.n-ts.n)/3.0;
  //*barn=(tu.n+td.n+ts.n)/3.0;
  
  gap1=-1.0*tu.m+tu.m-4.0*G*tu.qq-4.0*K*td.qq*ts.qq;
  gap2=-1.0*td.m+td.m-4.0*G*td.qq-4.0*K*tu.qq*ts.qq;
  gap3=-1.0*ts.m+ts.m-4.0*G*ts.qq-4.0*K*td.qq*tu.qq;

  double vec1, vec2, vec3;
  vec1=-1.0*tu.nu+tu.mu-4.0*GV*tu.n;
  vec2=-1.0*td.nu+td.mu-4.0*GV*td.n;
  vec3=-1.0*ts.nu+ts.mu-4.0*GV*ts.n;

  //free_dvector(iparam,1,inp);
  return 0;
}

double eos_quark_njl_vec::iqq(double x, int np, double *param) {
  double m, mu, t, en, ret;
  char ch;
  static const double limit=20.0;
  m=param[1];
  mu=param[2];
  t=param[3];
  en=sqrt(x*x+m*m);
  ret=1.0;
  if ((en-mu)/t>limit) {
    ret-=0.0;
  } else if ((en-mu)/t<-limit) {
    ret-=1.0;
  } else {
    ret-=1.0/(1+exp((en-mu)/t));
  }
  if ((en+mu)/t>limit) {
    ret-=0.0;
  } else if ((en+mu)/t<-limit) {
    ret-=1.0;
  } else {
    ret-=1.0/(1+exp((en+mu)/t));
  }
  ret*=-3.0*x*x/en/pi2*m;
  ret+=1.0;
  return ret;
}

double eos_quark_njl_vec::ide(double x, int np, double *param) {
  double m, mu, t, en, ret;
  char ch;
  m=param[1];
  mu=param[2];
  t=param[3];
  en=sqrt(x*x+m*m);
  ret=0.0;
  static const double limit=20.0;
  if ((en-mu)/t>limit) {
    ret+=0.0;
  } else if ((en-mu)/t<-limit) {
    ret+=1.0;
  } else {
    ret+=1.0/(1+exp((en-mu)/t));
  }
  if ((en+mu)/t>limit) {
    ret-=0.0;
  } else if ((en+mu)/t<-limit) {
    ret-=1.0;
  } else {
    ret-=1.0/(1+exp((en+mu)/t));
  }
  ret*=3.0*x*x/pi2;
  ret+=1.0;
  return ret;
}

double eos_quark_njl_vec::ied(double x, int np, double *param) {
  double m, mu, t, en, ret, en0, m0;
  m=param[1];
  mu=param[2];
  t=param[3];
  m0=param[4];
  en=sqrt(x*x+m*m);
  //  en0=sqrt(x*x+m0*m0);
  ret=-en;
  static const double limit=20.0;
  if ((en-mu)/t>limit) {
    ret+=0.0;
  } else if ((en-mu)/t<-limit) {
    ret+=en;
  } else {
    ret+=en/(1+exp((en-mu)/t));
  }
  if ((en+mu)/t>limit) {
    ret+=0.0;
  } else if ((en+mu)/t<-limit) {
    ret+=en;
  } else {
    ret+=en/(1+exp((en+mu)/t));
  }
  ret*=3.0/pi2*x*x;
  return ret;
}

double eos_quark_njl_vec::ipr(double x, int np, double *param) {
  double m, mu, t, en, ret, m0, en0;
  m=param[1];
  mu=param[2];
  t=param[3];
  m0=param[4];
  en=sqrt(x*x+m*m);
  // en0=sqrt(x*x+m0*m0);
  ret=en;
  if ((mu-en)/t>limit) {
    ret+=mu-en;
  } else if ((mu-en)/t<-limit) {
    ret+=0.0;
  } else {
    ret+=t*log(1+exp((mu-en)/t));
  }
  if ((-en-mu)/t>limit) {
    ret+=-en-mu;
  } else if ((-en-mu)/t<-limit) {
    ret+=0.0;
  } else {
    ret+=t*log(1+exp((-en-mu)/t));
  }
  ret*=3.0/pi2*x*x;
  return ret;
}

int eos_quark_njl_vec::gapfunmsvec(size_t nv, const ubvector &x,
                                   ubvector &y) {

  double gap1, gap2, gap3, vec1, vec2, vec3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];
  up->nu=x[3];
  down->nu=x[4];
  strange->nu=x[5];

  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,vec1,vec2,vec3,*eos_thermo);

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;
  y[3]=vec1;
  y[4]=vec2;
  y[5]=vec3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2]) || !std::isfinite(y[3]) ||
      !std::isfinite(y[4]) || !std::isfinite(y[5])) return 2;

  return 0;
}

int eos_quark_njl_vec::calc_p(quark &u, quark &d, quark &s, thermo &th) {
  ubvector x(6);
  int ret;

  if (fromqq==true) {
    O2SCL_ERR("Unimplemented.",o2scl::exc_eunimpl);
  }
  
  up=&u;
  down=&d;
  strange=&s;
  eos_thermo=&th;
  
  x[0]=u.ms;
  x[1]=d.ms;
  x[2]=s.ms;
  x[3]=u.nu;
  x[4]=d.nu;
  x[5]=s.nu;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_quark_njl_vec::gapfunmsvec),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  
  ret=solver->msolve(6,x,fmf);
  
  return ret;
}
