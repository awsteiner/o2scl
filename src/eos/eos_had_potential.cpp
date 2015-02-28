/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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

#include <o2scl/eos_had_potential.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_potential::eos_had_potential() {
  def_mu_deriv.h=0.01;
  mu_deriv_set=false;
  mu_deriv_ptr=&def_mu_deriv;

  neutron->init(939.0/hc_mev_fm,2.0);
  proton->init(939.0/hc_mev_fm,2.0);
  
  bpal_esym=30.0/hc_mev_fm;
}

double eos_had_potential::mom_integral(double pft, double pftp) {
  double qf=(pft+pftp)/2.0, result;
  double K=16.0/3.0/pow(2.0*pi,6.0)*pi2*Lambda*Lambda;
  double tx=1.0+4.0*qf*qf/Lambda/Lambda;
  result=K*((qf-Lambda/2.0*atan(2.0*qf/Lambda))*
	    4.0*(pow(pft,3.0)+pow(pftp,3.0))
	    -(3.0*(pft*pft+pftp*pftp)+Lambda*Lambda/2.0)*qf*qf+pow(qf,4.0)+
	    (0.75*Lambda*Lambda*(pft*pft+pftp*pftp)+pow(Lambda,4.0)/8.0-
	     0.375*pow(pft*pft-pftp*pftp,2.0))*
	    log(tx));
  return result;
}

double eos_had_potential::energy(double var) {
  double n, hamk, ham, ham1, ham2, ham3=0.0, xp;

  //---------------------------------------
  // Some local variables of interest:
  //
  // hamk is just the kinetic part of the hamiltonian 
  //   hbar^2 tau / (2 m^{star})
  //
  // ham{1-3} are remaining parts of the hamiltonian
  //

  if (mode==nmode) neutron->n=var;
  else if (mode==pmode) proton->n=var;

  if (neutron->n<0.0) neutron->n=0.01;
  if (proton->n<0.0) proton->n=0.01;
  
  n=neutron->n+proton->n;
  xp=proton->n/n;
  double delta=1.0-2.0*xp;
  
  neutron->ms=neutron->m;
  proton->ms=proton->m;

  nrf.calc_density_zerot(*neutron);
  nrf.calc_density_zerot(*proton);

  hamk=neutron->ed+proton->ed;

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
    
    ham=hamk+ham1+ham2+ham3;

  } else if (form==bpal_form) {
    
    double u=n/rho0;
    double cu=cbrt(u);
    double pf=cbrt(1.5*pi2*n);
    double pf0=cbrt(1.5*pi2*rho0);
    
    double ef0=pf0*pf0/2.0/neutron->m;
    // this "hamk" is a little different than the one above because
    // it doesn't include the symmetry energy, which is included 
    // separately below
    hamk=(0.6*ef0*cu*cu+neutron->m)*n;
    
    ham1=0.5*A/rho0*n*n;
    
    ham2=B*n*pow(u,sigma)/(1.0+Bp*pow(u,sigma-1.0));
    
    double g1=n*3.0*C1*pow(Lambda/pf0,3.0)*(pf/Lambda-atan(pf/Lambda));
    double g2=n*3.0*C2*pow(Lambda2/pf0,3.0)*(pf/Lambda2-atan(pf/Lambda2));
    
    ham3=g1+g2;

    double c2=cbrt(2.0);
    double fu;
    if (sym_index==2) fu=2.0*u*u/(1.0+u);
    else if (sym_index==3) fu=sqrt(u);
    else fu=u;
    double hsym=(c2*c2-1.0)*0.6*pf*pf/(neutron->m+proton->m)*
      (cu*cu-fu)+bpal_esym*fu;
    hsym*=delta*delta*n;
    
    ham=hamk+ham1+ham2+ham3+hsym;

  } else {

    ham1=2.0/3.0*A/rho0*((1.0+0.5*x0)*n*n-(0.5+x0)*
			 (neutron->n*neutron->n+proton->n*proton->n));
    double term=((1.0+0.5*x3)*n*n-
		 (0.5+x3)*(neutron->n*neutron->n+proton->n*proton->n))*
      pow(n,sigma-1.0);
    ham2=4.0/3.0*B/pow(rho0,sigma)*term/
      (1.0+4.0/3.0*Bp/pow(rho0,sigma-1.0)*term/n/n);
    double u=n/rho0;

    if (form==bgbd_form) {

      double gn, gp;
      gn=Lambda*Lambda/pi2*(neutron->kf-Lambda*atan(neutron->kf/Lambda));
      gp=Lambda*Lambda/pi2*(proton->kf-Lambda*atan(proton->kf/Lambda));
      ham3=0.8/rho0*(C1+2.0*z1)*n*(gn+gp)+0.4/rho0*(C1-8.0*z1)*
	(neutron->n*gn+proton->n*gp);

    } else if (form==bpalb_form) {

      double gn1, gp1, gn2, gp2;
      gn1=Lambda*Lambda/pi2*(neutron->kf-Lambda*atan(neutron->kf/Lambda));
      gp1=Lambda*Lambda/pi2*(proton->kf-Lambda*atan(proton->kf/Lambda));
      gn2=Lambda2*Lambda2/pi2*(neutron->kf-Lambda2*atan(neutron->kf/Lambda2));
      gp2=Lambda2*Lambda2/pi2*(proton->kf-Lambda2*atan(proton->kf/Lambda2));
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

    ham=hamk+ham1+ham2+ham3;

  }
  
  return ham;
}

int eos_had_potential::calc_e(fermion &ne, fermion &pr, 
			      thermo &locth) {
  double xp, n, hamk, ham, ham1, ham2, ham3;
  double dhdnn, dhdnp, na, npa, nna, term, term2, gn, gp;

  if (ne.n<=0.0 && pr.n<=0.0) {
    ne.ms=ne.m;
    pr.ms=pr.m;
    ne.mu=ne.m;
    pr.mu=pr.m;
    ne.pr=0.0;
    pr.pr=0.0;
    ne.ed=0.0;
    pr.ed=0.0;
    locth.pr=0.0;
    locth.ed=0.0;
    locth.en=0.0;
    return success;
  } else if (ne.n<=0.0) {
    ne.n=0.0;
  } else if (pr.n<=0.0) {
    pr.n=0.0;
  }

  ne.non_interacting=false;
  pr.non_interacting=false;

  n=ne.n+pr.n;

  set_n_and_p(ne,pr);
  
  double tmp;
  funct11 df=std::bind(std::mem_fn<double(double)>
		       (&eos_had_potential::energy),
		       this,std::placeholders::_1);

  mode=nmode;
  tmp=ne.n;
  ne.mu=mu_deriv_ptr->deriv(ne.n,df);
  ne.n=tmp;

  mode=pmode;
  tmp=pr.n;
  pr.mu=mu_deriv_ptr->deriv(pr.n,df);
  pr.n=tmp;

  mode=normal;
  ham=energy(pr.n);
  
  if (form==mdi_form) {
    double termn=-(Cl+Cu)*ne.m*Lambda*Lambda/4.0/pow(ne.kf,3.0)/rho0/
      pi2*(-4.0*ne.kf*ne.kf+(2.0*ne.kf*ne.kf+Lambda*Lambda)*
	   log(1.0+4.0*ne.kf*ne.kf/Lambda/Lambda));
    double termp=-(Cl+Cu)*pr.m*Lambda*Lambda/4.0/pow(pr.kf,3.0)/rho0/
      pi2*(-4.0*pr.kf*pr.kf+(2.0*pr.kf*pr.kf+Lambda*Lambda)*
	   log(1.0+4.0*pr.kf*pr.kf/Lambda/Lambda));
    ne.ms=ne.m/(1.0+termn);
    pr.ms=pr.m/(1.0+termp);
  } else if (form==gbd_form) {
    double termn=-2.0*ne.m*Lambda*Lambda/rho0/
      pow(ne.kf*ne.kf+Lambda*Lambda,2.0)*(Cl*ne.n+Cu*pr.n);
    double termp=-2.0*pr.m*Lambda*Lambda/rho0/
      pow(pr.kf*pr.kf+Lambda*Lambda,2.0)*(Cl*pr.n+Cu*ne.n);
    ne.ms=ne.m/(1.0+termn);
    pr.ms=pr.m/(1.0+termp);
  } else if (form==bgbd_form) {
    double termn=-0.8*ne.m/rho0*Lambda*Lambda/
      pow(ne.kf*ne.kf+Lambda*Lambda,2.0)*
      (3.0*C1*ne.n+2.0*C1*pr.n-4.0*z1*ne.n+4.0*z1*pr.n);
    double termp=-0.8*pr.m/rho0*Lambda*Lambda/
      pow(pr.kf*pr.kf+Lambda*Lambda,2.0)*
      (3.0*C1*pr.n+2.0*C1*ne.n-4.0*z1*pr.n+4.0*z1*ne.n);
    ne.ms=ne.m/(1.0+termn);
    pr.ms=pr.m/(1.0+termp);
  } else if (form==bpalb_form) {
    double termn=(-ne.n*(C1-8.0*z1)*Lambda*Lambda-
		  2.0*n*(C1+2.0*z1)*Lambda*Lambda)/
      pow(ne.kf*ne.kf+Lambda*Lambda,2.0);
    termn+=(-ne.n*(C2-8.0*z2)*Lambda2*Lambda2-
	    2.0*n*(C2+2.0*z2)*Lambda2*Lambda2)/
      pow(ne.kf*ne.kf+Lambda2*Lambda2,2.0);
    termn*=0.8*ne.m/rho0;
    double termp=(-pr.n*(C1-8.0*z1)*Lambda*Lambda-
		  2.0*n*(C1+2.0*z1)*Lambda*Lambda)/
      pow(pr.kf*pr.kf+Lambda*Lambda,2.0);
    termp+=(-pr.n*(C2-8.0*z2)*Lambda2*Lambda2-
	    2.0*n*(C2+2.0*z2)*Lambda2*Lambda2)/
      pow(pr.kf*pr.kf+Lambda2*Lambda2,2.0);
    termp*=0.8*pr.m/rho0;
    ne.ms=ne.m/(1.0+termn);
    pr.ms=pr.m/(1.0+termp);
  } else if (form==sl_form) {
    double termn=(-4.0*ne.kf*ne.n*(C1-8.0*z1)-8.0*ne.kf*n*(C1+2.0*z1))/
      (5.0*rho0*Lambda*Lambda);
    termn+=(-4.0*ne.kf*ne.n*(C2-8.0*z2)-8.0*ne.kf*n*(C2+2.0*z2))/
      (5.0*rho0*Lambda2*Lambda2*pow(1.0+ne.kf*ne.kf/Lambda2/Lambda2,2.0));
    termn*=ne.m/ne.kf;
    double termp=(-4.0*pr.kf*pr.n*(C1-8.0*z1)-8.0*pr.kf*n*(C1+2.0*z1))/
      (5.0*rho0*Lambda*Lambda);
    termp+=(-4.0*pr.kf*pr.n*(C2-8.0*z2)-8.0*pr.kf*n*(C2+2.0*z2))/
      (5.0*rho0*Lambda2*Lambda2*pow(1.0+pr.kf*pr.kf/Lambda2/Lambda2,2.0));
    termp*=pr.m/pr.kf;
    ne.ms=ne.m/(1.0+termn);
    pr.ms=pr.m/(1.0+termp);
  } else {
    ne.ms=ne.m;
    pr.ms=pr.m;
  }

  if (ne.kf<1.0e-12) ne.ms=ne.m;
  if (pr.kf<1.0e-12) pr.ms=pr.m;
  
  locth.ed=ham;
  locth.pr=-locth.ed+ne.mu*ne.n+pr.mu*pr.n;
  locth.en=0.0;
  
  return success;
}

