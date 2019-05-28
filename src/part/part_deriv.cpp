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
  
  fermion_rel fr;
  if (fr.calc_mu_deg_tlate<fermion_deriv>(f,temper,prec)==false) {
    return false;
  }
  
  // Compute psi and tt
  double psi;
  if (f.inc_rest_mass) psi=(f.nu-f.ms)/temper;
  else psi=(f.nu+f.m-f.ms)/temper;
  double tt=temper/f.ms;
      
  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);
      
  // Define x = psi * t = (mu/m - 1) and related values
  double x=psi*tt;
  double sx=sqrt(x);
  double s2x=sqrt(2.0+x);
  double x2=x*x;
  double x3=x2*x;
  double x4=x2*x2;
  
  // First order density term (first order entropy term is zero)
  double dndmu_term1=sx*s2x*(1.0+x)/f.ms/f.ms;
  
  // Second order terms
  double dndmu_term2=tt*tt*pi2/6.0*(1.0+x)*(-1.0+2.0*x*(2.0+x))/
    f.ms/f.ms/sx/s2x/x/(2.0+x);
  double dndT_term2=tt*pi2/3.0*(1.0+2.0*x*(2.0+x))/
    f.ms/f.ms/sx/s2x;
  double dsdT_term2=pi2/3.0*(1.0+x)*sx*s2x/
    f.ms/f.ms;

  // Third order terms
  double dndmu_term3=-7.0*pow(pi*tt,4.0)/24.0*(1.0+x)/sx/s2x/
    x3/(x+2.0)/(x+2.0)/(x+2.0)/f.ms/f.ms;
  double dndT_term3=7.0*pow(pi*tt,4.0)/tt/30.0/
    f.ms/f.ms/pow(x*(2.0+x),2.5);
  double dsdT_term3=7.0*pow(pi*tt,2.0)*pi2/30.0/
    f.ms/f.ms*(1.0+x)*(-1.0+2.0*x*(2.0+x))/x/(2.0+x)/sx/s2x;

  // Fourth order terms for density and entropy
  double dndmu_term4=-31.0*pow(pi*tt,6.0)/48.0*(1.0+x)*
    (3.0+2.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),5.5);
  double dndT_term4=31.0*pow(pi*tt,6.0)/tt/168.0*
    (7.0+6.0*x*(2.0+x))/f.ms/f.ms/pow(x*(2.0+x),4.5);
  double dsdT_term4=-155.0*pow(pi*tt,4.0)*pi2/168.0*
    (1.0+x)/f.ms/f.ms/pow(x*(2.0+x),3.5);

  // Add up all the terms
  f.dndmu=prefac*(dndmu_term1+dndmu_term2+dndmu_term3+dndmu_term4);
  f.dndT=prefac*(dndT_term2+dndT_term3+dndT_term4);
  f.dsdT=prefac*(dsdT_term2+dsdT_term3+dsdT_term4);

  return true;
}

bool fermion_deriv_thermo::calc_mu_ndeg
(fermion_deriv &f, double temper, double prec, bool inc_antip) {

  fermion_rel fr;
  if (fr.calc_mu_ndeg_tlate<fermion_deriv>(f,temper,prec,
					   inc_antip)==false) {
    return false;
  }

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

  // Prefactor 'd' in Johns96
  double prefac=f.g/2.0/pi2*pow(f.ms,4.0);

  // One term is always used, so only values of max_term greater than
  // 0 are useful.
  static const size_t max_term=200;
  
  double first_dndT=0.0;
  double first_dsdT=0.0;
  double first_dndmu=0.0;

  double nu2=f.nu;
  if (f.inc_rest_mass==false) nu2+=f.m;

  f.dndmu=0.0;
  f.dndT=0.0;
  f.dsdT=0.0;
  
  for(size_t j=1;j<=max_term;j++) {
    
    double dj=((double)j);
    double jot=dj/tt;

    // Here, we are only computing the derivatives, but we need to
    // compute the terms in the pressure, density, and entropy because
    // they are used in computing the terms for the derivatives.
    double pterm, nterm, enterm;
    double dndmu_term, dndT_term, dsdT_term;

    fr.ndeg_terms(j,tt,psi*tt,f.ms,f.inc_rest_mass,inc_antip,
		  pterm,nterm,enterm);
    
    if (inc_antip==false) {
      dndmu_term=nterm*jot;
      dndT_term=jot*enterm-nterm/tt;
      dsdT_term=(3.0*tt-2.0*dj*xx-2.0*dj)/tt/tt*enterm+
	(5.0*dj*tt-2.0*dj*dj*xx+5.0*dj*tt*xx-dj*dj*xx*xx)/
	dj/tt/tt/tt*nterm;
    } else {
      dndmu_term=nterm*jot;
      dndT_term=jot*enterm*tanh(jot*(xx+1.0))-
	(tt+2.0*dj*(1.0+xx))/sinh(jot*(xx+1.0))*nterm/tt/tt;
      dsdT_term=(2.0*dj*(1.0+xx)*tanh(jot*(xx+1.0))-3.0*tt)*enterm/tt/tt+
	(2.0*pow(dj*1.0+xx,2.0)*tanh(jot*(xx+1.0))-
	 dj*dj*(2.0+2.0*xx+xx*xx)*cosh(jot*(xx+1.0))-
	 5.0*dj*(1.0+xx)*tt)*nterm/dj/tt/tt/tt;
    }
    dndmu_term/=f.ms;
    dndT_term/=f.ms;
    dsdT_term/=f.ms;

    if (j==1) {
      first_dndT=dndT_term;
      first_dsdT=dsdT_term;
      first_dndmu=dndmu_term;
    }
    f.dndmu+=dndmu_term;
    f.dndT+=dndT_term;
    f.dsdT+=dsdT_term;

    /*
      cout << j << " " << dj << " " << tt << " " << xx << " "
      << f.ms << " " << pterm << " " << nterm << " "
      << enterm << " "
      << dndmu_term << " " << dndT_term << endl;
    */

    // If the first terms are zero, then the rest of the terms
    // will be zero so just return early
    if (first_dndT==0.0 && first_dndmu==0.0 && first_dsdT==0.0) {
      f.dndmu=0.0;
      f.dndT=0.0;
      f.dsdT=0.0;
      return true;
    }
    
    // Stop if the last term is sufficiently small compared to
    // the first term
    if (j>1 &&
	fabs(dndT_term)<prec*fabs(first_dndT) &&
	fabs(dndmu_term)<prec*fabs(first_dndmu) &&
	fabs(dsdT_term)<prec*fabs(first_dsdT)) {
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
