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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_had_apr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_apr::eos_had_apr() {
  pion=best; 

  select(a18_uix_deltav);
  
  neutron->init(939.0/hc_mev_fm,2.0);
  proton->init(939.0/hc_mev_fm,2.0);
  
  lp=0;

  parent_method=false;

  fet=&nrf;
}

eos_had_apr::~eos_had_apr() {
}

int eos_had_apr::gradient_qij2(double nn, double np, 
			   double &qnn, double &qnp, double &qpp, 
			   double &dqnndnn, double &dqnndnp,
			   double &dqnpdnn, double &dqnpdnp,
			   double &dqppdnn, double &dqppdnp) {

  double rho=nn+np, p3=par[3], p4=par[4], p5=par[5], a1=p3-2.0*p5, 
    ex=exp(-p4*rho);
  
  qnn=0.25*ex*(-6.0*p5-p4*a1*nn-2.0*p4*a1*np);
  qnp=0.125*ex*(4.0*(p3-4.0*p5)-3.0*p4*a1*nn-3.0*p4*a1*np);
  qpp=0.25*ex*(-6.0*p5-2.0*p4*a1*nn-p4*a1*np);

  dqnndnn=0.25*ex*p4*(-p3+8.0*p5+p4*a1*nn+2.0*p4*a1*np);
  dqnndnp=0.25*ex*p4*(p4*a1*nn+2.0*(-p3+5.0*p5+p4*a1*np));
  dqnpdnn=0.125*ex*p4*(-7.0*p3+22.0*p5+3.0*p4*a1*(nn+np));
  dqnpdnp=dqnpdnn;
  dqppdnn=0.25*ex*p4*(-2.0*(p3-5.0*p5)+2.0*p4*a1*nn+p4*a1*np);
  dqppdnp=0.25*ex*p4*(-p3+8.0*p5+2.0*p4*a1*nn+p4*a1*np);
  
  return 0;
}

int eos_had_apr::calc_e(fermion &ne, fermion &pr, thermo &lth) {

  double barn, xp, t1, t2, t3, t4, t5, t6, t7, t8, t9;
  double dt4, nb2, kin, low, high, gl1, gl2;
  double dmsndnn, dmspdnp, dkindnn, dkindnp, dmsndnp, dmspdnn;
  double gh1, gh2, dgl1, dgl2, dgh1, dgh2;

#if !O2SCL_NO_RANGE_CHECK
  if (!std::isfinite(ne.n) || !std::isfinite(ne.n)) {
    O2SCL_ERR2("Nucleon densities not finite in ",
	       "eos_had_apr::calc_e().",exc_einval);
  }
  if (ne.n<0.0 || pr.n<0.0) {
    O2SCL_ERR2("Nucleon densities negative in ",
	       "eos_had_apr::calc_e().",exc_einval);
  }
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_had_apr::calc_e().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_apr::calc_e().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true) {
    O2SCL_ERR2("Neutron or protons non-interacting in ",
	       "eos_had_apr::calc_e().",exc_einval);
  }
#endif
  
  //---------------------------------------
  // Some local variables of interest:
  //
  // kin is just the kinetic part of the hamiltonian hbar^2 tau / (2 mstar)
  // low and high are the energy density independent parts of the 
  //   hamiltonian;
  // Note: Remember that the chemical potential is exactly the single 
  //   particle energy evaluated at the Fermi momentum

  ne.non_interacting=false;
  pr.non_interacting=false;

  if (ne.n<=0.0 && pr.n<=0.0) {
    ne.ms=ne.m;
    pr.ms=pr.m;
    ne.mu=ne.m;
    pr.mu=pr.m;
    ne.pr=0.0;
    pr.pr=0.0;
    ne.ed=0.0;
    pr.ed=0.0;
    lth.pr=0.0;
    lth.ed=0.0;
    lth.en=0.0;
    return success;
  } else if (ne.n<=0.0) {
    ne.n=0.0;
  } else if (pr.n<=0.0) {
    pr.n=0.0;
  }

  barn=ne.n+pr.n;
  xp=pr.n/barn;

  // Partial derivatives of x wrt nn and np. 
  // Note that d(barn)/dnn=d(barn)dnp=1.0
  double dxdnn=-xp/barn, dxdnp=(1.0-xp)/barn;

  // Landau effective masses
  t1=barn*exp(-par[4]*barn);
  t2=par[3]+par[5]*(1-xp);
  t3=par[3]+xp*par[5];

  ne.ms=1.0/(1.0/ne.m+2.0*t1*t2);
  pr.ms=1.0/(1.0/pr.m+2.0*t1*t3);

  // We don't record error values, since these functions usually
  // always succeed
  nrf.calc_density_zerot(ne);
  nrf.calc_density_zerot(pr);
  
  if (ne.n==0 && pr.n==0) {
    lth.ed=0.0;
    ne.mu=ne.m; 
    pr.mu=pr.m;
  } else {
    kin=ne.ed+pr.ed;

    // Derivative of kinetic energy term:
    dmsndnn=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-t1*par[5]*dxdnn);
    dmsndnp=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-t1*par[5]*dxdnp);
    dmspdnn=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+t1*par[5]*dxdnn);
    dmspdnp=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+t1*par[5]*dxdnp);
    
    dkindnn=ne.nu-(ne.ed-ne.n*ne.m)/ne.ms*dmsndnn-
      (pr.ed-pr.n*pr.m)/pr.ms*dmspdnn;
    dkindnp=pr.nu-(pr.ed-pr.n*pr.m)/pr.ms*dmspdnp-
      (ne.ed-ne.n*ne.m)/ne.ms*dmsndnp;
    
    nb2=barn*barn;
    t4=exp(-par[9]*par[9]*nb2);
    dt4=-par[9]*par[9]*2.0*barn*t4;
    t5=barn-par[19];
    t6=barn-par[20];
    t9=(1.0-2.0*xp);
    t8=t9*t9;
    t7=1.0-t8;

    gl1=-nb2*(par[1]+par[2]*barn+par[6]*nb2+(par[10]+par[11]*barn)*t4);
    gl2=-nb2*(par[12]/barn+par[7]+par[8]*barn+par[13]*t4);
    gh1=gl1-nb2*(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5);
    gh2=gl2-nb2*(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6);

    low=gl1*t7+gl2*t8;
    high=gh1*t7+gh2*t8;

    dgl1=gl1*2.0/barn-nb2*(par[2]+2.0*par[6]*barn+t4*par[11]+
			   (par[10]+par[11]*barn)*dt4);
    dgl2=gl2*2.0/barn-nb2*(-par[12]/nb2+par[8]+par[13]*dt4);
    
    if ((low<=high && pion==best) || pion==ldp || barn<0.16) {
      
      // low-density phase
      lth.ed=kin+low;
      
      // This is not the fastest way, but it makes the contributions
      // to the derivatives from the various terms clear.
      ne.mu=dkindnn+dgl1*t7+dgl2*t8+dxdnn*4.0*t9*(gl1-gl2);
      pr.mu=dkindnp+dgl1*t7+dgl2*t8+dxdnp*4.0*t9*(gl1-gl2);

      lp=ldp;

    } else {

      // high-density phase
      lth.ed=kin+high;
      
      dgh1=dgl1-barn*(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5)*
	(2.0+par[18]*barn)-nb2*(par[17]+2.0*t5*par[21])*exp(par[18]*t5);
      dgh2=dgl2-barn*(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6)*
	(2.0+par[16]*barn)-nb2*(par[15]+2.0*t6*par[14])*exp(par[16]*t6);
      
      // This is not the fastest way, but it makes the contributions
      // to the derivatives from the various terms clear.
      ne.mu=dkindnn+dgh1*t7+dgh2*t8+dxdnn*4.0*t9*(gh1-gh2);
      pr.mu=dkindnp+dgh1*t7+dgh2*t8+dxdnp*4.0*t9*(gh1-gh2);

      lp=hdp;
      
    }
  }

  // Thermodynamics
  lth.pr=-lth.ed+ne.mu*ne.n+pr.mu*pr.n;
  lth.en=0.0;

  if (!std::isfinite(lth.pr)) {
    cout << ne.n << " " << pr.n << endl;
    O2SCL_ERR("Pressure not finite in calc_e().",exc_efailed);
    return exc_efailed;
  }

  return 0;
}

double eos_had_apr::fesym_diff(double nb) {

  double ret, t1, t2_neut, t2_nuc, nb2, t4, t5, t6;
  double gl1, gl2, gh1, gh2;
  
  if (parent_method) {
    return eos_had_base::fesym_diff(nb);
  }

  // Landau effective masses
  t1=nb*exp(-par[4]*nb);
  t2_neut=par[3]+par[5];
  t2_nuc=par[3]+par[5]*0.5;
  
  double msn_neut=1.0/(1.0/neutron->m+2.0*t1*t2_neut);
  double msn_nuc=1.0/(1.0/neutron->m+2.0*t1*t2_nuc);
  
  double kfn_neut=pow(3.0*pi2*nb,1.0/3.0);
  double kfn_nuc=pow(1.5*pi2*nb,1.0/3.0);
  
  double ed_neut=pow(kfn_neut,5.0)/10.0/pi2/msn_neut+nb*neutron->m;
  double ed_nuc=2.0*pow(kfn_nuc,5.0)/10.0/pi2/msn_nuc+
    nb/2.0*(neutron->m+proton->m);
  
  nb2=nb*nb;
  t4=exp(-par[9]*par[9]*nb2);
  t5=nb-par[19];
  t6=nb-par[20];

  gl1=-nb2*(par[1]+par[2]*nb+par[6]*nb2+(par[10]+par[11]*nb)*t4);
  gl2=-nb2*(par[12]/nb+par[7]+par[8]*nb+par[13]*t4);
  gh1=gl1-nb2*(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5);
  gh2=gl2-nb2*(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6);

  if (pion==hdp) {
    
    ret=((ed_neut+gh2)/nb-neutron->m)-
      ((ed_nuc+gh1)/nb-(neutron->m+proton->m)/2.0);

    lp=hdp;
    
  } else {
    
    ret=((ed_neut+gl2)/nb-neutron->m)-
      ((ed_nuc+gl1)/nb-(neutron->m+proton->m)/2.0);
    
    lp=ldp;
    
  }

  return ret;
}

double eos_had_apr::fcomp(double nb) {

  if (parent_method) {
    return eos_had_base::fcomp(nb);
  }
  
  /// Compute the compressibility directly with 9*nb*d^2(epsilon)/d(nb^2)

  double t1, t2, t3, t4, t5, t6, t7, t8, t9;
  double dt4, nb2, low, high, gl1, gl2, msn, msp;
  double gh1, gh2, ret, x=0.5;
  
  // Landau effective masses

  t1=nb*exp(-par[4]*nb);
  t2=par[3]+par[5]*(1-x);
  t3=par[3]+x*par[5];
  
  msn=1.0/(1.0/neutron->m+2.0*t1*t2);
  msp=1.0/(1.0/proton->m+2.0*t1*t3);

  double dt1=-par[4]*t1+t1/nb;
  double ddt1=-par[4]*dt1+dt1/nb-t1/nb/nb;
  
  double dmn=-msn*msn*2.0*t2*dt1;
  double dmp=-msp*msp*2.0*t3*dt1;
  double ddmn=pow(2.0*msn,3.0)*t2*t2*dt1*dt1-2.0*t2*ddt1*msn*msn;
  double ddmp=pow(2.0*msp,3.0)*t3*t3*dt1*dt1-2.0*t3*ddt1*msp*msp;
  
  /// Kinetic energy part

  double kfn=pow(3.0*pi2*nb*(1.0-x),1.0/3.0);
  double kfp=pow(3.0*pi2*nb*x,1.0/3.0);
  
  double ednpp=pi2/msn/kfn/4.0-kfn*kfn*dmn/msn/msn/2.0+
    pow(kfn,5.0)/10.0/pi2*(2.0*dmn*dmn/pow(msn,3.0)-ddmn/msn/msn);
  double edppp=pi2/msp/kfp/4.0-kfp*kfp*dmp/msp/msp/2.0+
    pow(kfp,5.0)/10.0/pi2*(2.0*dmp*dmp/pow(msp,3.0)-ddmp/msp/msp);
  
  /// Potential energy part

  nb2=nb*nb;
  t4=exp(-par[9]*par[9]*nb2);
  t5=nb-par[19];
  t6=nb-par[20];
  t9=(1.0-2.0*x);
  t8=t9*t9;
  t7=1.0-t8;
  
  gl1=-nb2*(par[1]+par[2]*nb+par[6]*nb2+(par[10]+par[11]*nb)*t4);
  gl2=-nb2*(par[12]/nb+par[7]+par[8]*nb+par[13]*t4);

  double fh1=(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5);
  double fh2=(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6);

  gh1=gl1-nb2*fh1;
  gh2=gl2-nb2*fh2;
  
  low=gl1*t7+gl2*t8;
  high=gh1*t7+gh2*t8;

  /// First and second derivatives of the potential energy part

  dt4=-par[9]*par[9]*2.0*nb*t4;
  double ddt4=-par[9]*par[9]*2.0*(t4+nb*dt4);
  
  double dfl1=(par[2]+2.0*par[6]*nb+t4*par[11]+(par[10]+par[11]*nb)*dt4);
  double ddfl1=(2.0*par[6]+2.0*dt4*par[11]+(par[10]+par[11]*nb)*ddt4);
  
  double dfl2=(-par[12]/nb2+par[8]+par[13]*dt4);
  double ddfl2=(2.0*par[12]/nb2/nb+par[13]*ddt4);
  
  double dfh1=(par[17]*(1.0+par[18]*t5)+par[21]*t5*(2.0+t5*par[18]))*
    exp(par[18]*t5);
  double dfh2=(par[15]*(1.0+par[16]*t6)+par[14]*t6*(2.0+t6*par[16]))*
    exp(par[16]*t6);

  double ddfh1=(par[17]*(2.0*par[18]+t5*par[18]*par[18])+
		par[21]*(2.0+4.0*par[18]*t5+t5*t5*par[18]*par[18]))*
    exp(par[18]*t5);
  double ddfh2=(par[15]*(2.0*par[16]+t6*par[16]*par[16])+
		par[14]*(2.0+4.0*par[16]*t6+t6*t6*par[16]*par[16]))*
    exp(par[16]*t6);
  
  
  double ddgl1=2.0*gl1/nb/nb-4.0*nb*dfl1-nb*nb*ddfl1;
  double ddgl2=2.0*gl2/nb/nb-4.0*nb*dfl2-nb*nb*ddfl2;
  
  double ddgh1=-2.0*fh1-4.0*nb*dfh1-nb*nb*ddfh1+ddgl1;
  double ddgh2=-2.0*fh2-4.0*nb*dfh2-nb*nb*ddfh2+ddgl2;
  
  double ddlow=ddgl1*t7+ddgl2*t8;
  double ddhigh=ddgh1*t7+ddgh2*t8;

  /// Return the correct value for the LDP or HDP

  if ((low<=high && pion==best) || pion==ldp) {
    
    ret=9.0*nb*(ednpp+edppp+ddlow);
    
  } else {
    
    ret=9.0*nb*(ednpp+edppp+ddhigh);

  }
  
  return ret;
}

void eos_had_apr::select(int model_index) {
  
  choice=model_index;
  par[3]=89.8/hc_mev_fm;
  par[4]=0.457;
  par[5]=-59.0/hc_mev_fm;

  if (choice==a18_uix) {
    /*A18+UIX*/
    par[1]=328.8;
    par[2]=-404.6;
    par[6]=-34.0;
    par[7]=217.5;
    par[8]=-385.6;
    par[9]=6.35;
    par[10]=25.4;
    par[11]=0.0;
    par[12]=0.47;
    par[13]=-0.9;
    par[14]=-452;
    par[15]=217.1;
    par[16]=-1.0;
    par[17]=100.3;
    par[18]=-1.19;
    par[19]=0.32;
    par[20]=0.2;
    par[21]=-275;
  } else if (choice==a18_deltav) {
    /* A18+deltav */
    par[1]=281.0;
    par[2]=-151.1;
    par[6]=-10.6;
    par[7]=210.1;
    par[8]=-158.0;
    par[9]=5.88;
    par[10]=58.8;
    par[11]=-15.0;
    par[12]=-.2;
    par[13]=-0.9;
    par[14]=0.0;
    par[15]=0.0;
    par[16]=0.0;
    par[17]=0.0;
    par[18]=0.0;
    par[19]=0.0;
    par[20]=0.0;
    par[21]=0.0;
  } else if (choice==a18) {
    /* A18 */
    par[1]=297.6;
    par[2]=-134.6;
    par[6]=-15.9;
    par[7]=215.0;
    par[8]=-116.5;
    par[9]=6.42;
    par[10]=51.0;
    par[11]=-35.0;
    par[12]=-.2;
    par[13]=0.0;
    par[14]=0.0;
    par[15]=0.0;
    par[16]=0.0;
    par[17]=0.0;
    par[18]=0.0;
    par[19]=0.0;
    par[20]=0.0;
    par[21]=0.0;
  } else if (choice==a18_uix_deltav) {
    /* A18+UIX*+deltav */
    par[1]=337.2;
    par[2]=-382.0;
    par[6]=-19.1;
    par[7]=214.6;
    par[8]=-384.0;
    par[9]=6.4;
    par[10]=69.0;
    par[11]=-33.0;
    par[12]=0.35;
    par[13]=0.0;
    par[14]=0.0;
    par[15]=287.0;
    par[16]=-1.54;
    par[17]=175.0;
    par[18]=-1.45;
    par[19]=0.32;
    par[20]=0.195;
    par[21]=0.0;
  } else {
    O2SCL_ERR("Improper value in eos_had_apr::select().",exc_einval);
  }
  
  par[1]/=hc_mev_fm;
  par[2]/=hc_mev_fm;
  par[6]/=hc_mev_fm;
  par[7]/=hc_mev_fm;
  par[8]/=hc_mev_fm;
  par[10]/=hc_mev_fm;
  par[11]/=hc_mev_fm;
  par[12]/=hc_mev_fm;
  par[13]/=hc_mev_fm;
  par[14]/=hc_mev_fm;
  par[15]/=hc_mev_fm;
  par[17]/=hc_mev_fm;
  par[21]/=hc_mev_fm;

  return;
}

int eos_had_apr::calc_temp_e(fermion &ne, fermion &pr, const double temper, 
			 thermo &lth) {

#if !O2SCL_NO_RANGE_CHECK
  if (!std::isfinite(ne.n) || !std::isfinite(ne.n) ||
      !std::isfinite(temper)) {
    O2SCL_ERR2("Nucleon densities or temperature not finite in ",
	       "eos_had_apr::calc_temp_e().",exc_einval);
  }
  if (ne.n<0.0 || pr.n<0.0) {
    O2SCL_ERR2("Nucleon densities negative in ",
	       "eos_had_apr::calc_temp_e().",exc_einval);
  }
  if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
    O2SCL_ERR2("Neutron or proton spin degeneracies wrong in ",
	       "eos_had_apr::calc_temp_e().",exc_einval);
  }
  if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
    O2SCL_ERR2("Neutron or proton masses wrong in ",
	       "eos_had_apr::calc_temp_e().",exc_einval);
  }
  if (ne.non_interacting==true || pr.non_interacting==true) {
    O2SCL_ERR2("Neutron or protons non-interacting in ",
	       "eos_had_apr::calc_temp_e().",exc_einval);
  }
#endif

  double barn, xp, t1, t2, t3, t4, t5, t6, t7, t8, t9;
  double dt4, nb2, kin, low, high, gl1, gl2;
  double dmsndnn, dmspdnp, dkindnn, dkindnp, dmsndnp, dmspdnn;
  double gh1, gh2, dgl1, dgl2, dgh1, dgh2;
  
  //---------------------------------------
  // Some local variables of interest:
  //
  // kin is just the kinetic part of the hamiltonian hbar^2 tau / (2 mstar)
  // low and high are the energy density independent parts of the 
  //   hamiltonian;
  // Note: Remember that the chemical potential is exactly the single 
  //   particle energy evaluated at the Fermi momentum
  
  ne.non_interacting=false;
  pr.non_interacting=false;

  //---------------------------------------
  // If the temperature is too small, just use 
  // the zero-temperature   

  if (temper<=0.0) {
    calc_e(ne,pr,lth);
    return success;
  }
  
  if (ne.n==0.0 && pr.n==0.0) {
    ne.ms=ne.m;
    pr.ms=pr.m;
    ne.mu=ne.m;
    pr.mu=pr.m;
    ne.pr=0.0;
    pr.pr=0.0;
    ne.ed=0.0;
    pr.ed=0.0;
    lth.pr=0.0;
    lth.ed=0.0;
    lth.en=0.0;
    return success;
  }

  barn=ne.n+pr.n;
  xp=pr.n/barn;

  // Partial derivatives of x wrt nn and np. 
  // Note that d(barn)/dnn=d(barn)dnp=1.0
  double dxdnn=-xp/barn, dxdnp=(1.0-xp)/barn;

  // Landau effective masses
  t1=barn*exp(-par[4]*barn);
  t2=par[3]+par[5]*(1-xp);
  t3=par[3]+xp*par[5];

  ne.ms=1.0/(1.0/ne.m+2.0*t1*t2);
  pr.ms=1.0/(1.0/pr.m+2.0*t1*t3);

  if (ne.n>0.0) {
    nrf.calc_density(ne,temper);
  } else {
    // If the neutron density is zero, we just assume we're 
    // computing pure proton matter
    ne.ed=0.0;
    ne.en=0.0;
    ne.pr=0.0;
    ne.nu=0.0;
  }
  if (pr.n>0.0) {
    nrf.calc_density(pr,temper);
  } else {
    // If the proton density is zero, we just assume we're 
    // computing pure neutron matter
    pr.ed=0.0;
    pr.en=0.0;
    pr.pr=0.0;
    pr.nu=0.0;
  }

  if (ne.n==0 && pr.n==0) {
    lth.ed=0.0;
    ne.mu=ne.m; 
    pr.mu=pr.m;
  } else {
    kin=ne.ed+pr.ed;

    // Derivative of kinetic energy term:
    dmsndnn=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-t1*par[5]*dxdnn);
    dmsndnp=-2.0*ne.ms*ne.ms*(t2*t1*(1.0/barn-par[4])-t1*par[5]*dxdnp);
    dmspdnn=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+t1*par[5]*dxdnn);
    dmspdnp=-2.0*pr.ms*pr.ms*(t3*t1*(1.0/barn-par[4])+t1*par[5]*dxdnp);
    dkindnn=ne.nu-(ne.ed-ne.n*ne.m)/ne.ms*dmsndnn-
      (pr.ed-pr.n*pr.m)/pr.ms*dmspdnn;
    dkindnp=pr.nu-(pr.ed-pr.n*pr.m)/pr.ms*dmspdnp-
      (ne.ed-ne.n*ne.m)/ne.ms*dmsndnp;
    
    nb2=barn*barn;
    t4=exp(-par[9]*par[9]*nb2);
    dt4=-par[9]*par[9]*2.0*barn*t4;
    t5=barn-par[19];
    t6=barn-par[20];
    t9=(1.0-2.0*xp);
    t8=t9*t9;
    t7=1.0-t8;

    gl1=-nb2*(par[1]+par[2]*barn+par[6]*nb2+(par[10]+par[11]*barn)*t4);
    gl2=-nb2*(par[12]/barn+par[7]+par[8]*barn+par[13]*t4);
    gh1=gl1-nb2*(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5);
    gh2=gl2-nb2*(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6);

    low=gl1*t7+gl2*t8;
    high=gh1*t7+gh2*t8;

    dgl1=gl1*2.0/barn-nb2*(par[2]+2.0*par[6]*barn+t4*par[11]+
			   (par[10]+par[11]*barn)*dt4);
    dgl2=gl2*2.0/barn-nb2*(-par[12]/nb2+par[8]+par[13]*dt4);
    
    if ((low<=high && pion==best) || pion==ldp || barn<0.16) {
      // low-density phase
      lth.ed=kin+low;
      
      // This is not the fastest way, but it makes the contributions
      // to the derivatives from the various terms clear.
      ne.mu=dkindnn+dgl1*t7+dgl2*t8+dxdnn*4.0*t9*(gl1-gl2);
      pr.mu=dkindnp+dgl1*t7+dgl2*t8+dxdnp*4.0*t9*(gl1-gl2);
      
    } else if ((high<low && pion==best) || pion==hdp) {
      // high-density phase
      lth.ed=kin+high;
      
      dgh1=dgl1-barn*(par[17]*t5+par[21]*t5*t5)*exp(par[18]*t5)*
	(2.0+par[18]*barn)-nb2*(par[17]+2.0*t5*par[21])*exp(par[18]*t5);
      dgh2=dgl2-barn*(par[15]*t6+par[14]*t6*t6)*exp(par[16]*t6)*
	(2.0+par[16]*barn)-nb2*(par[15]+2.0*t6*par[14])*exp(par[16]*t6);
      
      // This is not the fastest way, but it makes the contributions
      // to the derivatives from the various terms clear.
      ne.mu=dkindnn+dgh1*t7+dgh2*t8+dxdnn*4.0*t9*(gh1-gh2);
      pr.mu=dkindnp+dgh1*t7+dgh2*t8+dxdnp*4.0*t9*(gh1-gh2);
      
    } else {
      O2SCL_ERR("Bad value for pion in calc_temp_e",exc_efailed);
    }
  }
  
  // Thermodynamics
  lth.en=ne.en+pr.en;
  lth.pr=-lth.ed+ne.mu*ne.n+pr.mu*pr.n+temper*lth.en;
  
  if (!std::isfinite(lth.pr)) {
    O2SCL_ERR("Pressure not finite in calc_e().",exc_efailed);
  }

  return success;
}

