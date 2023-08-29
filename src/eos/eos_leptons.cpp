/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_leptons.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_leptons::eos_leptons() {
  include_muons=true;
  include_photons=false;
  include_deriv=false;
  err_nonconv=false;
      
  convert_units<double> &cu=o2scl_settings.get_convert_units();
  cu_ld.default_conversions();
  cu_cdf25.default_conversions();
  
  e.init(cu.convert("kg","1/fm",mass_electron_f<double>()),2);
  mu.init(cu.convert("kg","1/fm",mass_muon_f<double>()),2);
  eld.init(cu_ld.convert("kg","1/fm",mass_electron_f<long double>()),2);
  muld.init(cu_ld.convert("kg","1/fm",mass_muon_f<long double>()),2);
  ecdf25.init(cu_cdf25.convert("kg","1/fm",
                               mass_electron_f<cpp_dec_float_25>()),2);
  mucdf25.init(cu_cdf25.convert("kg","1/fm",
                                mass_muon_f<cpp_dec_float_25>()),2);

  
  ph.init(0.0,2.0);

  pde_from_density=true;
  verbose=0;
  accuracy=acc_default;

}

int eos_leptons::electron_density(double T) {

  int retx;

  bool inc_rest_mass=false;
  if (e.inc_rest_mass) {
    
    // I find that the calculation without the rest mass is a bit more
    // stable, so we use that method and add the rest mass back in
    // after the fact.
    inc_rest_mass=true;
    e.inc_rest_mass=false;
    e.mu-=e.m;

  }

  if (accuracy==acc_fp_25) {
    ecdf25.n=e.n;
    ecdf25.mu=e.mu;
    ecdf25.inc_rest_mass=e.inc_rest_mass;
    retx=frel_cdf25.pair_density(ecdf25,T);
    e.mu=static_cast<double>(ecdf25.mu);
    e.ed=static_cast<double>(ecdf25.ed);
    e.pr=static_cast<double>(ecdf25.pr);
    e.en=static_cast<double>(ecdf25.en);
  } else if (accuracy==acc_ld) {
    eld.n=e.n;
    eld.mu=e.mu;
    eld.inc_rest_mass=e.inc_rest_mass;
    retx=frel_ld.pair_density(eld,T);
    e.mu=static_cast<double>(eld.mu);
    e.ed=static_cast<double>(eld.ed);
    e.pr=static_cast<double>(eld.pr);
    e.en=static_cast<double>(eld.en);
  } else {
    retx=frel.pair_density(e,T);
  }
      
  // Sometimes the solver fails, but we can recover by adjusting
  // the upper limit for degenerate fermions and tightening the
  // integration tolerances
  if (retx!=0 && accuracy==acc_default) {
        
    frel.upper_limit_fac=40.0;
    frel.fri.dit.tol_rel=1.0e-10;
    frel.fri.dit.tol_abs=1.0e-10;
    frel.fri.nit.tol_rel=1.0e-10;
    frel.fri.nit.tol_abs=1.0e-10;
        
    retx=frel.pair_density(e,T);

    // If it still fails, then we don't call the error handler here
    // because this function is used in pair_density_eq_fun().
        
    frel.upper_limit_fac=20.0;
    frel.fri.dit.tol_rel=1.0e-8;
    frel.fri.dit.tol_abs=1.0e-8;
    frel.fri.nit.tol_rel=1.0e-8;
    frel.fri.nit.tol_abs=1.0e-8;
        
  }

  if (inc_rest_mass) {
    e.inc_rest_mass=true;
    e.mu+=e.m;
    e.ed+=e.m*e.n;
  }

  return retx;
}

int eos_leptons::pair_density_eq_fun(size_t nv, const ubvector &x,
                                     ubvector &y, double T, double nq) {

  if (pde_from_density) {

    if (accuracy==acc_ld) {
      eld.n=x[0]*nq;
      int retx=electron_density(T);
      if (retx!=0) return retx;
    } else if (accuracy==acc_fp_25) {
      ecdf25.n=x[0]*nq;
      int retx=electron_density(T);
      if (retx!=0) return retx;
    } else {
      e.n=x[0]*nq;
      int retx=electron_density(T);
      if (retx!=0) return retx;
    }
    
  } else {
    
    e.mu=x[0];

    bool inc_rest_mass=false;
    if (e.inc_rest_mass) {
      inc_rest_mass=true;
      e.inc_rest_mass=false;
      e.mu-=e.m;
    }

    if (accuracy==acc_ld) {
      eld.mu=e.mu;
      frel_ld.pair_mu(eld,T);
      e.n=eld.n;
      e.ed=eld.ed;
      e.pr=eld.pr;
      e.en=eld.en;
    } else if (accuracy==acc_fp_25) {
      ecdf25.mu=e.mu;
      frel_cdf25.pair_mu(ecdf25,T);
      e.n=static_cast<double>(ecdf25.n);
      e.ed=static_cast<double>(ecdf25.ed);
      e.pr=static_cast<double>(ecdf25.pr);
      e.en=static_cast<double>(ecdf25.en);
    } else {
      frel.pair_mu(e,T);
    }

    if (inc_rest_mass) {
      e.inc_rest_mass=true;
      e.mu+=e.m;
      e.ed+=e.n*e.m;
    }
    
  }

  if (e.inc_rest_mass) {
    if (mu.inc_rest_mass) {
      mu.mu=e.mu;
    } else {
      mu.mu=e.mu-mu.m;
    }
  } else {
    if (mu.inc_rest_mass) {
      mu.mu=e.mu+e.m;
    } else {
      mu.mu=e.mu+e.m-mu.m;
    }
  }
      
  if (mu.inc_rest_mass) {
    mu.inc_rest_mass=false;
    mu.mu-=mu.m;
    frel.pair_mu(mu,T);
    mu.inc_rest_mass=true;
    mu.mu+=mu.m;
    mu.ed+=mu.m*mu.n;
  } else {
    frel.pair_mu(mu,T);
  }

  y[0]=(e.n+mu.n-nq)/fabs(nq);

  return 0;
}

int eos_leptons::pair_mu(double T) {

  // Electron section
  
  if (include_deriv) {
    
    fermion_deriv fd=e;
    
    if (accuracy==acc_ld || accuracy==acc_fp_25) {
      fdrel.multip=true;
    } else {
      fdrel.multip=false;
    }
    
    fdrel.pair_mu(fd,T);

    // Copy results from the fermion_deriv object
    e.n=fd.n;
    e.ed=fd.ed;
    e.pr=fd.pr;
    e.en=fd.en;
    
    ed.dndT=fd.dndT;
    ed.dndmu=fd.dndmu;
    ed.dsdT=fd.dsdT;
    
    // Collect the total derivative quantities in the thd object, the
    // totals for the non-derivative quantities are computed below.
    thd.dndT=fd.dndT;
    thd.dndmu=fd.dndmu;
    thd.dsdT=fd.dsdT;
    
  } else {
    
    frel.pair_mu(e,T);
    
  }

  // Collect the total for the non-derivative quantities
  th.ed=e.ed;
  th.pr=e.pr;
  th.en=e.en;

  // Muon section
  
  if (include_muons) {
    
    if (include_deriv) {
      
      fermion_deriv fd=mu;
      
      if (accuracy==acc_ld || accuracy==acc_fp_25) {
        fdrel.multip=true;
      } else {
        fdrel.multip=false;
      }
      
      fdrel.pair_mu(fd,T);
      
      mu.n=fd.n;
      mu.ed=fd.ed;
      mu.pr=fd.pr;
      mu.en=fd.en;
      
      mud.dndT=fd.dndT;
      mud.dndmu=fd.dndmu;
      mud.dsdT=fd.dsdT;
      
      thd.dndT+=fd.dndT;
      thd.dndmu+=fd.dndmu;
      thd.dsdT+=fd.dsdT;
      
    } else {
      
      frel.pair_mu(mu,T);
      
    }
    
    th.ed+=mu.ed;
    th.pr+=mu.pr;
    th.en+=mu.en;
    
  }

  // Photon section
  
  if (include_photons) {
    
    ph.massless_calc(T);
    
    if (include_deriv) {
      
      phd.dsdT=ph.g*pi2*3.0*T*T/22.5;
      phd.dndT=ph.g*zeta3_f<double>()/pi2*3.0*T*T;
      phd.dndmu=0.0;
      
      thd.dndT+=phd.dndT;
      thd.dndmu+=phd.dndmu;
      thd.dsdT+=phd.dsdT;
      
    }
    
    th.ed+=ph.ed;
    th.pr+=ph.pr;
    th.en+=ph.en;
  }

  return 0;
}

int eos_leptons::pair_mu_eq(double T) {
  // Set muon chemical potential from the electron chemical
  // potential
  if (include_muons) {
    if (e.inc_rest_mass) {
      if (mu.inc_rest_mass) {
        mu.mu=e.mu;
      } else {
        mu.mu=e.mu-mu.m;
      }
    } else {
      if (mu.inc_rest_mass) {
        mu.mu=e.mu+e.m;
      } else {
        mu.mu=e.mu+e.m-mu.m;
      }
    }
  }
  pair_mu(T);
  return 0;
}

int eos_leptons::pair_density(double T) {
      
  bool fr_en=frel.err_nonconv;
  frel.err_nonconv=false;

  int retx=electron_density(T);
      
  th.ed=e.ed;
  th.pr=e.pr;
  th.en=e.en;

  if (include_deriv) {
    fermion_deriv fd;
    fd=e;
    if (accuracy==acc_ld || accuracy==acc_fp_25) {
      fdrel.multip=true;
    } else {
      fdrel.multip=false;
    }
    fdrel.pair_mu(fd,T);
    ed.dndmu=fd.dndmu;
    ed.dndT=fd.dndT;
    ed.dsdT=fd.dsdT;
  }
  
  if (include_muons) {

    if (mu.inc_rest_mass) {
      mu.inc_rest_mass=false;
      mu.mu-=mu.m;
      retx=frel.pair_density(mu,T);
      mu.inc_rest_mass=true;
      mu.mu+=mu.m;
      mu.ed+=mu.m*mu.n;
    } else {
      retx=frel.pair_density(mu,T);
    }
        
    // Sometimes the solver fails, but we can recover by adjusting
    // the upper limit for degenerate fermions and tightening the
    // integration tolerances
    if (retx!=0) {
          
      frel.upper_limit_fac=40.0;
      frel.fri.dit.tol_rel=1.0e-10;
      frel.fri.dit.tol_abs=1.0e-10;
      frel.fri.nit.tol_rel=1.0e-10;
      frel.fri.nit.tol_abs=1.0e-10;
          
      if (mu.inc_rest_mass) {
        mu.inc_rest_mass=false;
        mu.mu-=mu.m;
        retx=frel.pair_density(mu,T);
        mu.inc_rest_mass=true;
        mu.mu+=mu.m;
        mu.ed+=mu.m*mu.n;
      } else {
        retx=frel.pair_density(mu,T);
      }
          
      if (retx!=0) {
        O2SCL_ERR2("Function pair_density() for muons failed in ",
                   "class eos_leptons().",o2scl::exc_efailed);
      }
        
      frel.upper_limit_fac=20.0;
      frel.fri.dit.tol_rel=1.0e-8;
      frel.fri.dit.tol_abs=1.0e-8;
      frel.fri.nit.tol_rel=1.0e-8;
      frel.fri.nit.tol_abs=1.0e-8;
          
    }
    
    if (include_deriv) {
      fermion_deriv fd;
      fd=mu;
    if (accuracy==acc_ld || accuracy==acc_fp_25) {
      fdrel.multip=true;
    } else {
      fdrel.multip=false;
    }
      fdrel.pair_mu(fd,T);
      mud.dndmu=fd.dndmu;
      mud.dndT=fd.dndT;
      mud.dsdT=fd.dsdT;
    }
    
    th.ed+=mu.ed;
    th.pr+=mu.pr;
    th.en+=mu.en;
  }
      
  if (include_photons) {
    ph.massless_calc(T);
    th.ed+=ph.ed;
    th.pr+=ph.pr;
    th.en+=ph.en;
  }

  frel.err_nonconv=fr_en;
      
  return 0;
}

int eos_leptons::pair_density_eq(double nq, double T) {
      
  bool fr_en=frel.err_nonconv;
  frel.err_nonconv=false;

  int retx;
  if (include_muons) {
    if (verbose>1) {
      std::cout << "pair_density_eq() with muons, pde_from_density="
                << pde_from_density << std::endl;
    }

    ubvector x(1), y(1);
    if (pde_from_density) {
      x[0]=e.n/nq;
    } else {
      x[0]=e.mu;
    }

    mm_funct mf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,double)>
       (&eos_leptons::pair_density_eq_fun),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,T,nq);
    mh.err_nonconv=false;
    mh.def_jac.err_nonconv=false;
    mh.tol_rel=1.0e-6;
    size_t maxj=10;
    if (verbose>1) {
      cout << "Initial guess: " << x[0] << endl;
    }
    int mret=mh.msolve(1,x,mf);
    for(size_t j=0;j<maxj && mret!=0;j++) {
      if (verbose>1) {
        cout << "Attempt " << j+2 << " with guess " << x[0]
             << " and tolerance: " << mh.tol_rel << std::endl;
      }
      mret=mh.msolve(1,x,mf);
      mh.tol_rel*=pow(10.0,1.0/2.0);
    }
    if (mret!=0) {
      std::cout << "nq,T,T_MeV: " << nq << " " << T << " "
                << T*o2scl_const::hc_mev_fm << std::endl;
      O2SCL_ERR2("Failed to compute muons in ",
                 "eos_leptons::pair_density_eq()",o2scl::exc_einval);
    }
    if (verbose>1) {
      cout << "Solution: " << x[0] << endl;
    }

    mf(1,x,y);
    e.n=x[0]*nq;

    if (include_deriv) {
      
    if (accuracy==acc_ld || accuracy==acc_fp_25) {
      fdrel.multip=true;
    } else {
      fdrel.multip=false;
    }

    fermion_deriv fd;
      fd=e;
      fdrel.pair_mu(fd,T);
      ed.dndmu=fd.dndmu;
      ed.dndT=fd.dndT;
      ed.dsdT=fd.dsdT;
      fd=mu;
      fdrel.pair_mu(fd,T);
      mud.dndmu=fd.dndmu;
      mud.dndT=fd.dndT;
      mud.dsdT=fd.dsdT;
    }
        
  } else {
    if (verbose>1) {
      std::cout << "pair_density_eq() just electrons." << std::endl;
    }
        
    e.n=nq;
    mu.n=0.0;
        
    retx=electron_density(T);
        
    if (include_deriv) {
      fermion_deriv fd;
      fd=e;
      if (accuracy==acc_ld || accuracy==acc_fp_25) {
	fdrel.multip=true;
      } else {
	fdrel.multip=false;
      }
      fdrel.pair_mu(fd,T);
      ed.dndmu=fd.dndmu;
      ed.dndT=fd.dndT;
      ed.dsdT=fd.dsdT;
    }
        
  }
      
  th.ed=e.ed;
  th.pr=e.pr;
  th.en=e.en;

  if (include_muons) {
    th.ed+=mu.ed;
    th.pr+=mu.pr;
    th.en+=mu.en;
  }
      
  if (include_photons) {
    ph.massless_calc(T);
    th.ed+=ph.ed;
    th.pr+=ph.pr;
    th.en+=ph.en;
    if (include_deriv) {
      phd.dsdT=ph.g*pi2*3.0*T*T/22.5;
      phd.dndT=ph.g*zeta3_f<double>()/pi2*3.0*T*T;
      phd.dndmu=0.0;
    }
  }

  if (include_deriv) {
    thd.dndmu=ed.dndmu;
    thd.dndT=ed.dndT;
    thd.dsdT=ed.dsdT;
    if (include_muons) {
      thd.dndmu+=mud.dndmu;
      thd.dndT+=mud.dndT;
      thd.dsdT+=mud.dsdT;
    }
    if (include_photons) {
      thd.dndmu+=phd.dndmu;
      thd.dndT+=phd.dndT;
      thd.dsdT+=phd.dsdT;
    }
  }
  
  frel.err_nonconv=fr_en;
      
  return 0;
}

