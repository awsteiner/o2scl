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

#include <o2scl/part_deriv.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

bool fermion_deriv_thermo::calc_mu_deg
(fermion_deriv &f, double temper, double prec) {

  // Double check to ensure T and mass are positive
  if (temper<0.0 || f.ms<0.0) {
    O2SCL_ERR2("Temperature or mass negative in fermion_deriv_thermo",
	       "::calc_mu_deg().",exc_einval);
  }
      
  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
      
  // Compute psi and tt
  double psi;
  if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
  else psi=(f.nu+f.m-f.ms)/temper;
  double tt=temper/f.ms;
      
  // Return false immediately psi<0 where the expressions below
  // don't work because of the square roots
  if (psi<0.0) return false;
      
  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);
      
  // Define x = psi * t = (mu/m - 1) and related values
  double x=psi*tt;
  double sx=sqrt(x);
  double s2x=sqrt(2.0+x);
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  
  // Evaluate the first and last term for the pressure
  double pterm1;
  if (x>1.0e-5) {
    pterm1=(x*(1.0+x)*(2.0+x)*(-3.0+2.0*x*(2.0+x))+6.0*sx*s2x*
	    log((sx+s2x)/sqrt(2.0)))/24.0/sx/s2x;
  } else {
    pterm1=x2*sx*(29568.0+15840.0*x+1540.0*x2-105.0*x3)/55440.0/sqrt(2.0);
  }
  double pterm4=-31.0*pow(pi*tt,6.0)/1008.0*(1.0+x)*
    sx*s2x/pow(x*(2.0+x),4.0);

  // Check if we're going to succeed
  if (fabs(pterm4)/fabs(pterm1)>prec) {
    return false;
  }
  
  // First order density term (first order entropy term is zero)
  double nterm1=sx*s2x*x*(2.0+x)/3.0/f.ms;
  double dndmu_term1=sx*s2x*(1.0+x)/f.ms/f.ms;
  
  // Second order terms
  double pterm2=tt*tt*pi2/6.0*(1.0+x)*sx*s2x;
  double nterm2=tt*tt*pi2/6.0*(1.0+4.0*x+2.0*x2)/
    f.ms/sx/s2x;
  double enterm2=tt*pi2/3.0*(1.0+x)*sx*s2x/f.ms;
  double dndmu_term2=tt*tt*pi2/6.0*(1.0+x)*(-1.0+2.0*x*(2.0+x))/
    f.ms/f.ms/sx/s2x/x/(2.0+x);
  double dndT_term2=tt*pi2/3.0*(1.0+2.0*x*(2.0+x))/
    f.ms/f.ms/sx/s2x;
  double dsdT_term2=pi2/3.0*(1.0+x)*sx*s2x/
    f.ms/f.ms;

  // Third order terms
  double pterm3=7.0*pow(pi*tt,4.0)/360.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/pow(x*(2.0+x),1.5);
  double nterm3=7.0*pow(pi*tt,4.0)/120.0/sx/s2x/
    x2/(x+2.0)/(x+2.0)/f.ms;
  double enterm3=7.0*pow(pi*tt,4.0)/tt/90.0*(1.0+x)*
    (-1.0+4.0*x+2.0*x2)/f.ms/sx/s2x/x/(x+2.0);
  double dndmu_term3=-7.0*pow(pi*tt,4.0)/24.0*(1.0+x)/sx/s2x/
    x3/(x+2.0)/(x+2.0)/(x+2.0)/f.ms/f.ms;
  double dndT_term3=7.0*pow(pi*tt,4.0)/tt/30.0/
    f.ms/f.ms/pow(x*(2.0+x),2.5);
  double dsdT_term3=7.0*pow(pi*tt,2.0)*pi2/30.0/
    f.ms/f.ms*(1.0+x)*(-1.0+2.0*x*(2.0+x))/x/(2.0+x)/sx/s2x;

  // Fourth order terms for density and entropy
  double nterm4=31.0*pow(pi*tt,6.0)/1008.0*sx*s2x*
    (7.0+12.0*x+6.0*x2)/f.ms/pow(x*(2.0+x),5.0);
  double enterm4=-31.0*pow(pi*tt,6.0)/tt/168.0*sx*s2x*
    (1.0+x)/pow(x*(2.0+x),4.0);
  double dndmu_term4=-31.0*pow(pi*tt,6.0)/48.0*(1.0+x)*
    (3.0+2.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),5.5);
  double dndT_term4=31.0*pow(pi*tt,6.0)/tt/168.0*
    (7.0+6.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),4.5);
  double dsdT_term4=-155.0*pow(pi*tt,4.0)*pi2/168.0*
    (1.0+x)/f.ms/f.ms/pow(x*(2.0+x),3.5);

  // Add up all the terms
  f.pr=prefac*(pterm1+pterm2+pterm3+pterm4);
  f.n=prefac*(nterm1+nterm2+nterm3+nterm4);
  f.en=prefac*(enterm2+enterm3+enterm4);
  f.ed=-f.pr+f.nu*f.n+temper*f.en;
  f.dndmu=prefac*(dndmu_term1+dndmu_term2+dndmu_term3+dndmu_term4);
  f.dndT=prefac*(dndT_term2+dndT_term3+dndT_term4);
  f.dsdT=prefac*(dsdT_term2+dsdT_term3+dsdT_term4);

  return true;
}

bool fermion_deriv_thermo::calc_mu_ndeg
(fermion_deriv &f, double temper, double prec, bool inc_antip) {

  if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

  // Compute psi and tt
  double psi, psi_num;
  if (f.inc_rest_mass) {
    psi_num=f.nu-f.ms;
  } else {
    psi_num=f.nu+f.m-f.ms;
  }
  psi=psi_num/temper;
  double tt=temper/f.ms;
  double xx=psi*tt;

  // Return false immediately if we're degenerate
  if (inc_antip==false && psi>-1.0) return false;

  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);

  // One term is always used, so only values of max_term greater than
  // 0 are useful.
  static const size_t max_term=200;
  
  // Maximum argument for exponential
  // double log_dbl_max=709.78;

  // Return zero if psi+1/t is too small
  if (psi+1.0/tt<-700.0) {
    f.n=0.0;
    f.ed=0.0;
    f.pr=0.0;
    f.en=0.0;
    return true;
  }

  // -----------------------------------------------------
  // Return early if the last term is going to be too large.
  
  // Ratio of last term to first term in the pressure expansion
  double rat;
  double dj1=((double)max_term), jot1=max_term/tt;
  double dj2=1.0, jot2=1.0/tt;
  if (inc_antip==false) {
    rat=exp(dj1*psi)/jot1/jot1*gsl_sf_bessel_Kn_scaled(2.0,jot1);
    rat/=exp(dj2*psi)/jot2/jot2*gsl_sf_bessel_Kn_scaled(2.0,jot2);
  } else {
    if (f.inc_rest_mass) {
      rat=exp(-jot1)*2.0*cosh(dj1*f.nu/temper)/jot1/jot1*
	gsl_sf_bessel_Kn_scaled(2.0,jot1);
      rat/=exp(-jot2)*2.0*cosh(dj2*f.nu/temper)/jot2/jot2*
	gsl_sf_bessel_Kn_scaled(2.0,jot2);
    } else {
      rat=exp(-jot1)*2.0*cosh(dj1*(f.nu+f.m)/temper)/jot1/jot1*
	gsl_sf_bessel_Kn_scaled(2.0,jot1);
      rat/=exp(-jot2)*2.0*cosh(dj2*(f.nu+f.m)/temper)/jot2/jot2*
	gsl_sf_bessel_Kn_scaled(2.0,jot2);
    }
  }

  // If the ratio between the last term and the first term is 
  // not small enough, return false
  if (std::isfinite(rat) && rat>prec) {
    return false;
  }
  
  double first_term=0.0;
  f.pr=0.0;
  f.n=0.0;
  f.en=0.0;

  double nu2=f.nu;
  if (f.inc_rest_mass==false) nu2+=f.m;
  
  for(size_t j=1;j<=max_term;j++) {
    
    double dj=((double)j);
    double jot=dj/tt;

    double pterm, nterm, enterm;
    double dndmu_term, dndT_term, dsdT_term;

    if (inc_antip==false) {
      pterm=exp(dj*psi)/jot/jot*gsl_sf_bessel_Kn_scaled(2.0,jot);
      if (j%2==0) {
	pterm*=-1.0;
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj-
		exp(dj*psi)/2.0/dj*(gsl_sf_bessel_Kn_scaled(1.0,jot)+
				    gsl_sf_bessel_Kn_scaled(3.0,jot)))/f.ms-
	  pterm*dj*psi_num/temper/temper;
      } else {
	enterm=(pterm*2.0/tt-pterm/tt/tt*dj+
		exp(dj*psi)/2.0/dj*(gsl_sf_bessel_Kn_scaled(1.0,jot)+
				    gsl_sf_bessel_Kn_scaled(3.0,jot)))/f.ms-
	  pterm*dj*psi_num/temper/temper;
      }
      nterm=pterm*dj/temper;
      dndmu_term=nterm*dj/temper;
      dndT_term=-dj*pterm/temper/temper+dj/temper*enterm;
      if (j%2==0) {
	dsdT_term=((xx+1.0)/2.0/tt/tt*exp(dj*xx/tt)*
		   (gsl_sf_bessel_Kn_scaled(1.0,jot)+
		    gsl_sf_bessel_Kn_scaled(3.0,jot))-
		   1.0/4.0/tt/tt*exp(dj*xx/tt)*
		   (gsl_sf_bessel_Kn_scaled(0.0,jot)+
		    2.0*gsl_sf_bessel_Kn_scaled(2.0,jot)+
		    gsl_sf_bessel_Kn_scaled(4.0,jot)))/f.ms/f.ms;
      } else {
	dsdT_term=(-(xx+1.0)/2.0/tt/tt*exp(dj*xx/tt)*
		   (gsl_sf_bessel_Kn_scaled(1.0,jot)+
		    gsl_sf_bessel_Kn_scaled(3.0,jot))+
		   1.0/4.0/tt/tt*exp(dj*xx/tt)*
		   (gsl_sf_bessel_Kn_scaled(0.0,jot)+
		    2.0*gsl_sf_bessel_Kn_scaled(2.0,jot)+
		    gsl_sf_bessel_Kn_scaled(4.0,jot)))/f.ms/f.ms;
      }
      dsdT_term+=(2.0*dj*(xx+1.0)-2.0*tt)/tt/temper/temper*pterm+
	(2.0*tt-dj*(xx+1.0))/tt/temper*enterm;
    } else {
      if (f.inc_rest_mass) {
	pterm=exp(-jot)*2.0*cosh(dj*f.nu/temper)/jot/jot*
	  gsl_sf_bessel_Kn_scaled(2.0,jot);
	if (j%2==0) {
	  pterm*=-1.0;
	}
	nterm=pterm*tanh(dj*f.nu/temper)*dj/temper;
      } else {
	pterm=exp(-jot)*2.0*cosh(dj*(f.nu+f.m)/temper)/jot/jot*
	  gsl_sf_bessel_Kn_scaled(2.0,jot);
	if (j%2==0) {
	  pterm*=-1.0;
	}
	nterm=pterm*tanh(dj*(f.nu+f.m)/temper)*dj/temper;
      }
      if (j%2==0) {
	enterm=(pterm*2.0/tt-cosh(dj*nu2/temper)/dj*exp(-jot)*
		(gsl_sf_bessel_Kn_scaled(1.0,jot)+
		 gsl_sf_bessel_Kn_scaled(3.0,jot))+2.0*pterm*nu2*dj/tt/tt*
		tanh(dj*nu2/temper)/f.ms)/f.ms;
      } else {
	enterm=(pterm*2.0/tt+cosh(dj*nu2/temper)/dj*exp(-jot)*
		(gsl_sf_bessel_Kn_scaled(1.0,jot)+
		 gsl_sf_bessel_Kn_scaled(3.0,jot))+2.0*pterm*nu2*dj/tt/tt*
		tanh(dj*nu2/temper)/f.ms)/f.ms;
      }
      dndmu_term=pterm*dj*dj/temper/temper;
      dndT_term=(dj/temper*enterm-dj/temper/temper*pterm)*
	tanh(dj*(xx+1.0)/tt)-dj*dj*(xx+1.0)/temper/temper*pterm/
	pow(cosh(dj*(xx+1.0)/tt),2.0);
      if (j%2==0) {
	dsdT_term=-(xx+1.0)/2.0/tt/tt*exp(-dj/tt)*sinh(dj*(xx+1.0)/tt)*
	  (gsl_sf_bessel_Kn_scaled(1.0,jot)+
	   gsl_sf_bessel_Kn_scaled(3.0,jot))-
	  1.0/2.0/dj*exp(-dj/tt)*cos(dj*(xx+1.0)/tt)*
	  (gsl_sf_bessel_Kn_scaled(0.0,jot)+
	   2.0*gsl_sf_bessel_Kn_scaled(2.0,jot)+
	   gsl_sf_bessel_Kn_scaled(4.0,jot));
      } else {
	dsdT_term=(xx+1.0)/2.0/tt/tt*exp(-dj/tt)*sinh(dj*(xx+1.0)/tt)*
	  (gsl_sf_bessel_Kn_scaled(1.0,jot)+
	   gsl_sf_bessel_Kn_scaled(3.0,jot))+
	  1.0/2.0/dj*exp(-dj/tt)*cos(dj*(xx+1.0)/tt)*
	  (gsl_sf_bessel_Kn_scaled(0.0,jot)+
	   2.0*gsl_sf_bessel_Kn_scaled(2.0,jot)+
	   gsl_sf_bessel_Kn_scaled(4.0,jot));
      }
      dsdT_term+=(-2.0*tt+dj)/tt/tt/tt*pterm+(2.0*tt+dj)/tt/tt*enterm-
	(dj*(xx+1.0)*enterm*tanh(dj*(xx+1.0)/tt))/tt/tt+
	(dj*dj*(xx+1.0)*(xx+1.0)*pterm/pow(cosh(dj*(xx+1.0)/tt),2.0))/
	pow(tt,4.0);
      dsdT_term/=f.ms;
    }
    
    if (j==1) first_term=pterm;
    f.pr+=pterm;
    f.n+=nterm;
    f.en+=enterm;
    f.dndmu+=dndmu_term;
    f.dndT+=dndT_term;
    f.dsdT+=dsdT_term;

    // If the first term is zero, then the rest of the terms
    // will be zero so just return early
    if (first_term==0.0) {
      f.pr=0.0;
      f.n=0.0;
      f.ed=0.0;
      f.en=0.0;
      f.dndmu=0.0;
      f.dndT=0.0;
      f.dsdT=0.0;
      return true;
    }

    // Stop if the last term is sufficiently small compared to
    // the first term
    if (j>1 && fabs(pterm)<prec*fabs(first_term)) {
      f.pr*=prefac;
      f.n*=prefac;
      f.en*=prefac;
      f.ed=-f.pr+f.nu*f.n+temper*f.en;
      f.dndT*=prefac;
      f.dndmu*=prefac;
      f.dsdT*=prefac;
      return true;
    }

    // End of 'for(size_t j=1;j<=max_term;j++)'
  }

  // We failed to add enough terms, so return false
  return false;
}
