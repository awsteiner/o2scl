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
#include <o2scl/fermion_mag_zerot.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

void fermion_mag_zerot::calc_mu_zerot_mag(fermion &f, double qB, 
					  double kappa) {

  if (qB==0.0) {
    calc_mu_zerot(f);
    return;
  }

  if (f.g!=2.0) {
    O2SCL_ERR("Class fermion_mag_zerot only works for g=2.",
	      exc_eunimpl);
  }

  double sign=1.0;
  // If qB<0, then q<0, so set sign to -1 and use qB=|qB|
  if (qB<0.0) {
    qB=-qB;
    sign=-1.0;
  }

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
    
  if (f.inc_rest_mass) {
    
    double nmax_upd=(f.nu*f.nu-f.ms*f.ms-qB+qB*sign)/2.0/qB;
    if (nmax_upd>0.0) {
      nmax_up=(int)(nmax_upd+1.0e-12);
    } else {
      nmax_up=-1;
    }
    
    double nmax_dnd=(f.nu*f.nu-f.ms*f.ms-qB-qB*sign)/2.0/qB;
    if (nmax_dnd>0.0) {
      nmax_dn=(int)(nmax_dnd+1.0e-12);
    } else {
      nmax_dn=-1;
    }
    
  } else {

    // If we're not including the rest mass, handle the m=m* case
    // separately, so we can take advantage of the cancelation

    if (f.ms==f.m) {
      
      double nmax_upd=(f.nu*f.nu+2.0*f.m*f.nu-qB+qB*sign)/2.0/qB;
      if (nmax_upd>0.0) {
	nmax_up=(int)(nmax_upd+1.0e-12);
      } else {
	nmax_up=-1;
      }
      
      double nmax_dnd=(f.nu*f.nu+2.0*f.m*f.nu-qB-qB*sign)/2.0/qB;
      if (nmax_dnd>0.0) {
	nmax_dn=(int)(nmax_dnd+1.0e-12);
      } else {
	nmax_dn=-1;
      }
    
    } else {
      
      double nmax_upd=((f.nu+f.m)*(f.nu+f.m)-f.ms*f.ms-qB+qB*sign)/2.0/qB;
      if (nmax_upd>0.0) {
	nmax_up=(int)(nmax_upd+1.0e-12);
      } else {
	nmax_up=-1;
      }
      
      double nmax_dnd=((f.nu+f.m)*(f.nu+f.m)-f.ms*f.ms-qB-qB*sign)/2.0/qB;
      if (nmax_dnd>0.0) {
	nmax_dn=(int)(nmax_dnd+1.0e-12);
      } else {
	nmax_dn=-1;
      }
    
    }

  }
  
  // Doesn't matter really if we choose nmax_up or nmax_dn here.
  if (nmax_up<sum_limit) {

    f.n=0.0;
    f.ed=0.0;

    // Spin up part
    for(int i=0;i<=nmax_up;i++) {
      double mt;

      if (f.inc_rest_mass) {

	double mt2=f.ms*f.ms+2.0*(i+0.5-0.5*sign)*qB;
	mt=sqrt(mt2);
	f.kf=sqrt(f.nu*f.nu-mt2);
	f.n+=f.kf;
	f.ed+=f.nu*f.kf+mt2*log(fabs((f.nu+f.kf)/mt));
	
      } else {

	double X=2.0*(i+0.5-0.5*sign)*qB;
	double mt2=f.ms*f.ms+X;
	mt=sqrt(mt2);
	if (f.ms==f.m) {
	  f.kf=sqrt(f.nu*f.nu+2.0*f.m*f.nu-X);
	} else {
	  f.kf=sqrt((f.nu+f.m)*(f.nu+f.m)-mt2);
	}
	f.n+=f.kf;
	f.ed+=(f.nu+f.m)*f.kf+mt2*log(fabs((f.nu+f.m+f.kf)/mt));
	
      }

    }
  
    // Spin down part
    for(int i=0;i<=nmax_dn;i++) {

      double mt;

      if (f.inc_rest_mass) {
	
	double mt2=f.ms*f.ms+2.0*(i+0.5+0.5*sign)*qB;
	mt=sqrt(mt2);
	f.kf=sqrt(f.nu*f.nu-mt2);
	f.n+=f.kf;
	f.ed+=f.nu*f.kf+mt2*log(fabs((f.nu+f.kf)/mt));
	
      } else {
	
	double X=2.0*(i+0.5+0.5*sign)*qB;
	double mt2=f.ms*f.ms+X;
	mt=sqrt(mt2);
	if (f.ms==f.m) {
	  f.kf=sqrt(f.nu*f.nu+2.0*f.m*f.nu-X);
	} else {
	  f.kf=sqrt((f.nu+f.m)*(f.nu+f.m)-mt2);
	}
	f.n+=f.kf;
	f.ed+=(f.nu+f.m)*f.kf+mt2*log(fabs((f.nu+f.m+f.kf)/mt));

      }

    }
    
    f.n*=qB/2.0/pi2;
    f.ed*=qB/4.0/pi2;
    if (f.inc_rest_mass==false) f.ed-=f.n*f.m;
    f.pr=f.n*f.nu-f.ed;

  } else {
    
    calc_mu_zerot(f);
    
  }

  return;
}

void fermion_mag_zerot::calc_density_zerot_mag
(fermion &f, double qB, double kappa) {

  if (qB==0.0) {
    calc_density_zerot(f);
    return;
  }

  if (f.g!=2.0) {
    O2SCL_ERR("Class fermion_mag_zerot only works for g=2.",
	      exc_eunimpl);
  }

  // Handle low density limit gracefully.
  // Can probably improve this section by directly computing
  // things instead of calling calc_mu_zerot_mag().
  if (f.n<pow(fabs(qB),1.5)*sqrt(2.0)/2.0/pi2) {
    double dentest=f.n;
    f.kf=2.0*pi2*f.n/fabs(qB);
    if (f.non_interacting) {
      if (f.inc_rest_mass) {
	f.mu=gsl_hypot(f.kf,f.m);
      } else {
	f.mu=gsl_hypot(f.kf,f.m)-f.m;
      }
    } else {
      if (f.inc_rest_mass) {
	f.nu=gsl_hypot(f.kf,f.ms);
      } else {
	f.nu=gsl_hypot(f.kf,f.ms)-f.m;
      }
    }
    if (f.non_interacting && f.inc_rest_mass) {
      double term;
      if (f.kf/f.m<1.0e-4) {
	double kf3=f.kf*f.kf*f.kf;
	double kf5=kf3*f.kf*f.kf;
	term=f.m*f.kf-kf3/6.0/f.m+3.0*kf5/40.0/f.m/f.m/f.m;
      } else {
	term=f.m*f.m*log((f.mu+f.kf)/f.m);
      }
      f.ed=f.mu*f.n/2+fabs(qB)/4.0/pi2*term;
      f.pr=f.mu*f.n/2-fabs(qB)/4.0/pi2*term;
    } else {
      calc_mu_zerot_mag(f,qB,kappa);
    }
    if (fabs(f.n-dentest)/dentest>1.0e-6) {
      O2SCL_ERR("Low density one-level failure in calc_density().",
		exc_esanity);
    }
    return;
  }

  qBt=qB;
  kt=kappa;
  dent=f.n;

  ubvector x(1), y(1);

  if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

  mm_funct11 mf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,fermion &f)>
     (&fermion_mag_zerot::solve_fun),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(f));

  // Construct an initial guess from the B=0 result
  calc_density_zerot(f);
  
  // Ensure that f.nu>f.ms
  if (f.inc_rest_mass && f.nu<f.ms) f.nu=f.ms*1.01;

  x[0]=f.nu;

  // Check initial guess 
  int ret2=solve_fun(1,x,y,f);
  if (ret2!=0) {
    O2SCL_ERR("Initial guess failed in calc_density_zerot_mag().",
	      exc_efailed);
  }

  density_root->err_nonconv=false;
  int ret=density_root->msolve(1,x,mf);
  if (ret!=0) {
    O2SCL_ERR("Solver failed in calc_density_zerot_mag().",
	      exc_efailed);
  }
  f.nu=x[0];
   
  // Compute final values
  solve_fun(1,x,y,f);

  if (f.non_interacting) { f.mu=f.nu; }

  return;
}

int fermion_mag_zerot::solve_fun(size_t nv, const ubvector &x,
				 ubvector &y, fermion &f) {
  
  if (!o2scl::is_finite(x[0])) {
    return 3;
  }

  if (f.non_interacting) {
    f.mu=x[0];
    if (f.inc_rest_mass && f.mu<f.m) return 1;
  } else {
    f.nu=x[0];
    if (f.inc_rest_mass && f.nu<f.ms) return 2;
  }
  
  calc_mu_zerot_mag(f,qBt,kt);
  y[0]=(f.n-dent)/dent;

  if (!o2scl::is_finite(y[0])) {
    return 4;
  }
  return success;
}

