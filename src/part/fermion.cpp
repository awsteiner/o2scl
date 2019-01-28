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

#include <o2scl/fermion.h>
// For gsl_hypot()
#include <gsl/gsl_sys.h>

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

fermion::fermion(double mass, double dof) : part(mass,dof) {
  kf=0.0;
  del=0.0;
}

void fermion_zerot::kf_from_density(fermion &f) {
  f.kf=cbrt(6.0*pi2/f.g*f.n);
  return;
}

void fermion_zerot::energy_density_zerot(fermion &f) {
  double r,efs;
  if (f.kf>0.0) {
    if (f.ms<=0.0) {
      f.ed=f.g*(pow(f.kf,4.0)/8.0/pi2);
    } else {
      efs=gsl_hypot(f.kf,f.ms);
      r=(f.kf+efs)/f.ms;
      f.ed=f.g/16.0/pi2*(2.0*f.kf*pow(efs,3.0)-f.kf*efs*f.ms*f.ms
			 -pow(f.ms,4.0)*log(r));
    }
  } else {
    f.ed=0.0;
  }
  return;
}

void fermion_zerot::pressure_zerot(fermion &f) {
  double r,efs;
  if (f.kf>0.0) {
    if (f.ms<=0.0) {
      f.pr=f.g*(pow(f.kf,4.0)/24.0/pi2);
    } else {
      efs=gsl_hypot(f.kf,f.ms);
      r=(f.kf+efs)/f.ms;
      f.pr=f.g/48.0/pi2*(2.0*efs*pow(f.kf,3.0)-3.0*f.kf*efs*f.ms*f.ms
			 +3.0*pow(f.ms,4.0)*log(r));
    }
  } else {
    f.pr=0.0;
  }
  return;
}

fermion_thermo::fermion_thermo() {
  massless_root=&def_massless_root;
}

void fermion_zerot::calc_mu_zerot(fermion &f) {
  bool nulessthan0;
  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  
  if (f.inc_rest_mass) {
    if (f.nu>f.ms) {
      nulessthan0=false;
      f.kf=sqrt(f.nu*f.nu-f.ms*f.ms);
    } else {
      nulessthan0=false;
      f.kf=0.0;
    }
  } else {
    double nupm=f.nu+f.m;
    if ((nupm)>f.ms) {
      nulessthan0=false;
      f.kf=sqrt(nupm*nupm-f.ms*f.ms);
    } else {
      nulessthan0=false;
      f.kf=0.0;
    }
  }
  f.n=f.g/6.0/pi2*pow(f.kf,3.0);
  energy_density_zerot(f);
  pressure_zerot(f);
  f.en=0.0;
  if (!f.inc_rest_mass) f.ed-=f.n*f.m;
  if (nulessthan0==true) {
    f.n*=-1.0;
    f.kf*=-1.0;
  }

  return;
}

void fermion_zerot::calc_density_zerot(fermion &f) {
  if (f.non_interacting) { f.ms=f.m; }

  f.kf=cbrt(6.0*pi2/f.g*f.n);
  f.nu=gsl_hypot(f.kf,f.ms);
  energy_density_zerot(f);
  pressure_zerot(f);
  f.en=0.0;

  if (!f.inc_rest_mass) {
    f.nu-=f.m;
    f.ed-=f.n*f.m;
  }

  if (f.non_interacting) { f.mu=f.nu; }

  return;
}

void fermion_thermo::massless_calc_mu(fermion &f, double temper) {
  
  double fm2, fm3;

  if (f.non_interacting) { f.nu=f.mu; }

  fm2=gsl_sf_fermi_dirac_int(2,f.nu/temper);
  fm3=gsl_sf_fermi_dirac_int(3,f.nu/temper);
  
  f.n=f.g/pi2*pow(temper,3.0)*fm2;
  f.ed=f.g*3.0/pi2*pow(temper,4.0)*fm3;
  f.pr=f.ed/3.0;
  f.en=(f.ed+f.pr-f.n*f.nu)/temper;

  return;
}

double fermion_thermo::massless_solve_fun
(double x, fermion &f, double temper) {
  double fm2=gsl_sf_fermi_dirac_int(2,x/(temper));
  return f.g*pow(temper,3.0)*fm2/pi2/f.n-1.0;
}

void fermion_thermo::massless_calc_density(fermion &f, double temper) {
  double x, T=temper;
  
  x=f.ms+temper;
  funct mf2=std::bind(std::mem_fn<double(double,fermion &,double)>
		      (&fermion_thermo::massless_solve_fun),
		      this,std::placeholders::_1,std::ref(f),temper);
  massless_root->solve(x,mf2);
  f.nu=x;

  massless_calc_mu(f,temper);

  // If the particle is non-interacting, then need to set
  // mu=nu to get the entropy right
  if (f.non_interacting) { f.mu=f.nu; }

  return;
}

void fermion_thermo::massless_pair_mu(fermion &f, double temper) {
  double pitmu, pitmu2, nu2;

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
  if (f.nu==0.0) {
    f.n=0.0;
    f.ed=f.g/8.0/pi2*7.0/15.0*pow(pi*temper,4.0);
    f.pr=f.ed/3.0;
    f.en=(f.ed+f.pr-f.n*f.mu)/temper;
  } else {
    nu2=f.nu*f.nu;
    pitmu=pi*temper/f.nu;
    pitmu2=pitmu*pitmu;
    f.n=f.g*f.nu*nu2/6.0/pi2*(1.0+pitmu2);
    f.ed=f.g*nu2*nu2/8.0/pi2*(1.0+2.0*pitmu2+7.0/15.0*pitmu2*pitmu2);
    f.pr=f.ed/3.0;
    f.en=(f.ed+f.pr-f.n*f.mu)/temper;
    
    // Might the following work better for the energy density?
    // pit=pi*temper;
    // pit2=pit*pit;
    // ed=g/8.0/pi2*(nu2*nu2+2.0*pit2*nu2+7.0/15.0*pit2*pit2);
    
  }

  return;
}

void fermion_thermo::massless_pair_density(fermion &f, double temper) {

  double t2=temper*temper,pitmu,pitmu2,nu2;
  double cbt, alpha, two13, alpha16;

  if (f.non_interacting) { f.ms=f.m; }
  if (f.n<=0.0) {
    f.nu=0.0;
    f.ed=f.g/8.0/pi2*7.0/15.0*pow(pi*temper,4.0);
    f.pr=f.ed/3.0;
  } else {
    alpha=f.g*f.g*pi2*t2*t2*t2/243.0/f.n/f.n;
    if (alpha>1.0e4) {
      f.nu=(2.0/3.0/sqrt(alpha)-8.0/81.0/pow(alpha,1.5)+
	    32.0/729.0/pow(alpha,2.5))*pi*temper/sqrt(3.0);
    } else if (alpha<3.0e-4) {
      two13=cbrt(2.0);
      alpha16=pow(alpha,1.0/6.0);
      f.nu=(two13/alpha16-alpha16/two13+alpha/alpha16/6.0/two13/two13
	    +alpha*alpha16/12.0/two13-alpha*alpha/alpha16/18.0/two13/two13-
	    5.0*alpha*alpha*alpha16/144.0/two13+
	    77.0/2592.0*alpha*alpha*alpha/alpha16/two13/two13)*
	pi*temper/sqrt(3.0);
    } else {
      cbt=pow(-1.0+sqrt(1.0+alpha),1.0/3.0)/pow(alpha,1.0/6.0);
      f.nu=pi*temper/sqrt(3.0)*(1.0/cbt-cbt);
    }
    pitmu=pi*temper/f.nu;
    pitmu2=pitmu*pitmu;
    nu2=f.nu*f.nu;
    f.ed=f.g*nu2*nu2/8.0/pi2*(1.0+2.0*pitmu2+7.0/15.0*pitmu2*pitmu2);
    f.pr=f.ed/3.0;

    if (!std::isfinite(f.nu)) {
      string str="Chemical potential not finite ("+dtos(f.nu)+
	") in fermion::massless_pair_density().";
      O2SCL_ERR(str.c_str(),exc_efailed);
    }
  }

  if (f.non_interacting) { f.mu=f.nu; }
  f.en=(f.ed+f.pr-f.n*f.nu)/temper;

  return;
}

bool fermion_thermo::calc_mu_deg(fermion &f, double temper, 
				 double prec) {
  return calc_mu_deg_tlate<fermion>(f,temper,prec);
}

bool fermion_thermo::calc_mu_ndeg(fermion &f, double temper, 
				  double prec, bool inc_antip) {
  return calc_mu_ndeg_tlate<fermion>(f,temper,prec,inc_antip);
}
