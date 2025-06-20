/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
  
  e.init(cu.convert("kg","1/fm",mass_electron_f<double>()),2);
  mu.init(cu.convert("kg","1/fm",mass_muon_f<double>()),2);
  tau.init(cu.convert("kg","1/fm",mass_tau_f<double>()),2);

  nu_e.init(0,2);
  nu_mu.init(0,2);
  nu_tau.init(0,2);
  
  ph.init(0.0,2.0);

  pde_from_density=true;
  verbose=0;
  accuracy=acc_default;

}

int eos_leptons::electron_density(double T) {

  bool inc_rest_mass=false;
  if (e.inc_rest_mass) {
    
    // I find that the calculation without the rest mass is a bit more
    // stable, so we use that method and add the rest mass back in
    // after the fact.
    inc_rest_mass=true;
    e.inc_rest_mass=false;
    e.mu-=e.m;

  }

  int retx=1;
  
  if (accuracy==acc_fp_25) {
    O2SCL_ERR("The eos_leptons object doesn't support multiprecision.",
              o2scl::exc_eunimpl);
  } else if (accuracy==acc_ld) {
    O2SCL_ERR("The eos_leptons object doesn't support multiprecision.",
              o2scl::exc_eunimpl);
  } else {
    retx=frel.pair_density(e,T);
  }
      
  // Sometimes the solver fails, but we can recover by adjusting
  // the upper limit for degenerate fermions and tightening the
  // integration tolerances
  if (retx!=0 && accuracy==acc_default) {
        
    frel.upper_limit_fac=40.0;
    frel.dit.tol_rel=1.0e-10;
    frel.dit.tol_abs=1.0e-10;
    frel.nit.tol_rel=1.0e-10;
    frel.nit.tol_abs=1.0e-10;
        
    retx=frel.pair_density(e,T);

    // If it still fails, then we don't call the error handler here
    // because this function is used in pair_density_eq_fun().
        
    frel.upper_limit_fac=20.0;
    frel.dit.tol_rel=1.0e-8;
    frel.dit.tol_abs=1.0e-8;
    frel.nit.tol_rel=1.0e-8;
    frel.nit.tol_abs=1.0e-8;
        
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
      O2SCL_ERR("This object doesn't do multip.",o2scl::exc_einval);
    } else if (accuracy==acc_fp_25) {
      O2SCL_ERR("This object doesn't do multip.",o2scl::exc_einval);
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
      O2SCL_ERR("This object doesn't do multip.",o2scl::exc_einval);
    } else if (accuracy==acc_fp_25) {
      O2SCL_ERR("This object doesn't do multip.",o2scl::exc_einval);
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

    mu.mu=e.mu;
    
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
      frel.dit.tol_rel=1.0e-10;
      frel.dit.tol_abs=1.0e-10;
      frel.nit.tol_rel=1.0e-10;
      frel.nit.tol_abs=1.0e-10;
          
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
      frel.dit.tol_rel=1.0e-8;
      frel.dit.tol_abs=1.0e-8;
      frel.nit.tol_rel=1.0e-8;
      frel.nit.tol_abs=1.0e-8;
          
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
      std::cout << "eos_leptons::pair_density_eq(): No muons."
                << std::endl;
    }
    mu.n=0.0;
    
    e.n=nq;
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

#ifdef O2SCL_SET_MULTIP

eos_leptons_multip::eos_leptons_multip() {

  cu_ld.default_conversions();
  cu_cdf25.default_conversions();
  
  eld.init(cu_ld.convert("kg","1/fm",mass_electron_f<long double>()),2);
  muld.init(cu_ld.convert("kg","1/fm",mass_muon_f<long double>()),2);
  tauld.init(cu_ld.convert("kg","1/fm",mass_tau_f<long double>()),2);
  
  ecdf25.init(cu_cdf25.convert("kg","1/fm",
                               mass_electron_f<cpp_dec_float_25>()),2);
  mucdf25.init(cu_cdf25.convert("kg","1/fm",
                                mass_muon_f<cpp_dec_float_25>()),2);
  taucdf25.init(cu_cdf25.convert("kg","1/fm",
                                 mass_tau_f<cpp_dec_float_25>()),2);

  ph_ld.init(0,2);
  ph_cdf25.init(0,2);
  
}

int eos_leptons_multip::electron_density(double T) {

  int retx;

  // If we're just using double precision, then use the parent
  // function
  if (accuracy!=acc_fp_25 && accuracy!=acc_ld) {
    return eos_leptons::electron_density(T);
  }
  
  // I find that the calculation without the rest mass is a bit more
  // stable, so we use that method and add the rest mass back in
  // later if necessary.
  bool inc_rest_mass=false;
  if (e.inc_rest_mass) {
    
    inc_rest_mass=true;
    e.inc_rest_mass=false;
    e.mu-=e.m;
  }

  if (accuracy==acc_fp_25) {
    O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
               "electron_density().",o2scl::exc_einval);
  } else if (accuracy==acc_ld) {
    O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
               "electron_density().",o2scl::exc_einval);
  }
      
  if (inc_rest_mass) {
    e.inc_rest_mass=true;
    eld.inc_rest_mass=true;
    ecdf25.inc_rest_mass=true;
    e.mu+=e.m;
    e.ed+=e.m*e.n;
  }

  return retx;
}

int eos_leptons_multip::electron_density_ld(long double T) {

  int retx;

  // I find that the calculation without the rest mass is a bit more
  // stable, so we use that method and add the rest mass back in
  // later if necessary.
  bool inc_rest_mass=false;
  if (eld.inc_rest_mass) {
    
    inc_rest_mass=true;
    eld.inc_rest_mass=false;
    eld.mu-=eld.m;
  }

  if (accuracy==acc_fp_25) {
    O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
               "electron_density_ld().",o2scl::exc_einval);
  } else {
    retx=frel_ld.pair_density(eld,T);
  }
      
  if (inc_rest_mass) {
    eld.inc_rest_mass=true;
    ecdf25.inc_rest_mass=true;
    eld.mu+=eld.m;
    eld.ed+=eld.m*eld.n;
  }

  return retx;
}

int eos_leptons_multip::electron_density_cdf25(cpp_dec_float_25 T) {

  int retx;

  // I find that the calculation without the rest mass is a bit more
  // stable, so we use that method and add the rest mass back in
  // later if necessary.
  bool inc_rest_mass=false;
  if (ecdf25.inc_rest_mass) {
    inc_rest_mass=true;
    ecdf25.inc_rest_mass=false;
    ecdf25.mu-=ecdf25.m;
  }
  
  retx=frel_cdf25.pair_density(ecdf25,T);
  
  if (inc_rest_mass) {
    ecdf25.inc_rest_mass=true;
    ecdf25.mu+=ecdf25.m;
    ecdf25.ed+=ecdf25.m*ecdf25.n;
  }

  return retx;
}

int eos_leptons_multip::pair_density_eq(double nq, double T) {
      
  bool fr_en=frel.err_nonconv;
  frel.err_nonconv=false;
  
  int retx=-99;
  
  if (include_muons) {
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq(): "
                << "with muons, pde_from_density="
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
       (&eos_leptons_multip::pair_density_eq_fun),
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
        O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
                   "pair_density_eq().",o2scl::exc_einval);
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
      std::cout << "eos_leptons_multip::pair_density_eq(): No muons."
                << std::endl;
    }
    mu.n=0.0;
    
    if (accuracy==acc_ld) {
      O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
                 "pair_density_eq().",o2scl::exc_einval);
    } else if (accuracy==acc_fp_25) {
      O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
                 "pair_density_eq().",o2scl::exc_einval);
    } else {
      e.n=nq;
      retx=electron_density(T);
    }
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq(): "
                << "Return value " << retx << std::endl;
    }
    
    if (include_deriv) {
      if (verbose>1) {
        std::cout << "eos_leptons_multip::pair_density_eq(): "
                  << "Including derivatives." << std::endl;
      }
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

int eos_leptons_multip::pair_density_eq_ld(long double nq, long double T) {
      
  bool fr_en=frel_ld.err_nonconv;
  frel_ld.err_nonconv=false;
  
  int retx=-99;
  if (include_muons) {
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_ld(): "
                << "with muons, pde_from_density="
                << pde_from_density << std::endl;
    }

    long double x, y;
    if (pde_from_density) {
      x=eld.n/nq;
    } else {
      x=eld.mu;
    }
    
    funct_ld mf=std::bind
      (std::mem_fn<long double(long double,long double,long double)>
       (&eos_leptons_multip::pair_density_eq_ld_fun),
       this,std::placeholders::_1,T,nq);
    size_t maxj=10;
    if (verbose>1) {
      cout << "Initial guess: " << x << endl;
    }
    int mret=rc_ld.solve(x,mf);
    for(size_t j=0;j<maxj && mret!=0;j++) {
      if (verbose>1) {
        cout << "Attempt " << j+2 << " with guess " << x
             << " and tolerance: " << rc_ld.tol_rel << std::endl;
      }
      mret=rc_ld.solve(x,mf);
      rc_ld.tol_rel*=pow(10.0,1.0/2.0);
    }
    if (mret!=0) {
      std::cout << "nq,T,T_MeV: " << nq << " " << T << " "
                << T*o2scl_const::hc_mev_fm << std::endl;
      O2SCL_ERR2("Failed to compute muons in ",
                 "eos_leptons::pair_density_eq()",o2scl::exc_einval);
    }
    if (verbose>1) {
      cout << "Solution: " << x << endl;
    }
    
    y=mf(x);
    eld.n=x*nq;
    
    if (include_deriv) {
      if (accuracy==acc_fp_25) {
        O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
                   "pair_density_eq_ld().",o2scl::exc_einval);
      } else {
        fdrel_ld.multip=true;
      }
      fermion_deriv_ld fd;
      fd=eld;
      fdrel_ld.pair_mu(fd,T);
      ed_ld.dndmu=fd.dndmu;
      ed_ld.dndT=fd.dndT;
      ed_ld.dsdT=fd.dsdT;
      fd=muld;
      fdrel_ld.pair_mu(fd,T);
      mud_ld.dndmu=fd.dndmu;
      mud_ld.dndT=fd.dndT;
      mud_ld.dsdT=fd.dsdT;
    }
    
  } else {
    
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_ld(): No muons."
                << std::endl;
    }
    muld.n=0.0;
    
    if (accuracy==acc_fp_25) {
      O2SCL_ERR2("This accuracy not supported in eos_leptons_multip::",
                 "pair_density_eq_ld().",o2scl::exc_einval);
    } else {
      eld.n=nq;
      retx=electron_density_ld(T);
    }
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_ld(): "
                << "Return value " << retx << std::endl;
    }
    
    if (include_deriv) {
      if (verbose>1) {
        std::cout << "eos_leptons_multip::pair_density_eq_ld(): "
                  << "Including derivatives." << std::endl;
      }
      fermion_deriv_ld fd;
      fd=eld;
      if (accuracy==acc_fp_25) {
	fdrel_cdf25.multip=true;
      } else {
	fdrel_ld.multip=true;
      }
      fdrel_ld.pair_mu(fd,T);
      ed_ld.dndmu=fd.dndmu;
      ed_ld.dndT=fd.dndT;
      ed_ld.dsdT=fd.dsdT;

    }
        
  }
      
  th_ld.ed=eld.ed;
  th_ld.pr=eld.pr;
  th_ld.en=eld.en;

  if (include_muons) {
    th_ld.ed+=muld.ed;
    th_ld.pr+=muld.pr;
    th_ld.en+=muld.en;
  }
      
  if (include_photons) {
    ph_ld.massless_calc(T);
    th_ld.ed+=ph_ld.ed;
    th_ld.pr+=ph_ld.pr;
    th_ld.en+=ph_ld.en;
    if (include_deriv) {
      phd_ld.dsdT=ph.g*pi2*3.0*T*T/22.5;
      phd_ld.dndT=ph.g*zeta3_f<double>()/pi2*3.0*T*T;
      phd_ld.dndmu=0.0;
    }
  }

  if (include_deriv) {
    thd_ld.dndmu=ed_ld.dndmu;
    thd_ld.dndT=ed_ld.dndT;
    thd_ld.dsdT=ed_ld.dsdT;
    if (include_muons) {
      thd_ld.dndmu+=mud_ld.dndmu;
      thd_ld.dndT+=mud_ld.dndT;
      thd_ld.dsdT+=mud_ld.dsdT;
    }
    if (include_photons) {
      thd_ld.dndmu+=phd_ld.dndmu;
      thd_ld.dndT+=phd_ld.dndT;
      thd_ld.dsdT+=phd_ld.dsdT;
    }
  }
  
  frel_ld.err_nonconv=fr_en;
      
  return 0;
}

int eos_leptons_multip::pair_density_eq_cdf25
(cpp_dec_float_25 nq, cpp_dec_float_25 T) {
      
  bool fr_en=frel_cdf25.err_nonconv;
  frel_cdf25.err_nonconv=false;
  
  int retx;
  if (include_muons) {
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_cdf25(): "
                << "with muons, pde_from_density="
                << pde_from_density << std::endl;
    }

    cpp_dec_float_25 x, y;
    if (pde_from_density) {
      x=ecdf25.n/nq;
    } else {
      x=ecdf25.mu;
    }

    funct_cdf25 mf=std::bind
      (std::mem_fn<cpp_dec_float_25(cpp_dec_float_25,
                                    cpp_dec_float_25,cpp_dec_float_25)>
       (&eos_leptons_multip::pair_density_eq_cdf25_fun),
       this,std::placeholders::_1,T,nq);
    size_t maxj=10;
    if (verbose>1) {
      cout << "Initial guess: " << x << endl;
    }
    int mret=rc_cdf25.solve(x,mf);
    for(size_t j=0;j<maxj && mret!=0;j++) {
      if (verbose>1) {
        cout << "Attempt " << j+2 << " with guess " << x
             << " and tolerance: " << rc_cdf25.tol_rel << std::endl;
      }
      mret=rc_cdf25.solve(x,mf);
      rc_cdf25.tol_rel*=pow(10.0,1.0/2.0);
    }
    if (mret!=0) {
      std::cout << "nq,T,T_MeV: " << nq << " " << T << " "
                << T*o2scl_const::hc_mev_fm << std::endl;
      O2SCL_ERR2("Failed to compute muons in ",
                 "eos_leptons::pair_density_eq()",o2scl::exc_einval);
    }
    if (verbose>1) {
      cout << "Solution: " << x << endl;
    }
    
    y=mf(x);
    ecdf25.n=x*nq;
    
    if (include_deriv) {
      fdrel_cdf25.multip=true;
      
      fermion_deriv_cdf25 fd;
      fd=ecdf25;
      fdrel_cdf25.pair_mu(fd,T);
      ed_cdf25.dndmu=fd.dndmu;
      ed_cdf25.dndT=fd.dndT;
      ed_cdf25.dsdT=fd.dsdT;
      fd=mucdf25;
      fdrel_cdf25.pair_mu(fd,T);
      mud_cdf25.dndmu=fd.dndmu;
      mud_cdf25.dndT=fd.dndT;
      mud_cdf25.dsdT=fd.dsdT;
    }
    
  } else {
    
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_cdf25(): No muons."
                << std::endl;
    }
    mucdf25.n=0.0;
    
    ecdf25.n=nq;
    retx=electron_density_cdf25(T);

    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_cdf25(): "
                << "Return value " << retx << std::endl;
    }
    
    if (include_deriv) {
      if (verbose>1) {
        std::cout << "eos_leptons_multip::pair_density_eq_cdf25(): "
                  << "Including derivatives." << std::endl;
      }
      fermion_deriv_cdf25 fd;
      fd=ecdf25;
      fdrel_cdf25.multip=true;
      fdrel_cdf25.pair_mu(fd,T);
      ed_cdf25.dndmu=fd.dndmu;
      ed_cdf25.dndT=fd.dndT;
      ed_cdf25.dsdT=fd.dsdT;

    }
        
  }
      
  th_cdf25.ed=ecdf25.ed;
  th_cdf25.pr=ecdf25.pr;
  th_cdf25.en=ecdf25.en;

  if (include_muons) {
    th_cdf25.ed+=mucdf25.ed;
    th_cdf25.pr+=mucdf25.pr;
    th_cdf25.en+=mucdf25.en;
  }
      
  if (include_photons) {
    ph_cdf25.massless_calc(T);
    th_cdf25.ed+=ph_cdf25.ed;
    th_cdf25.pr+=ph_cdf25.pr;
    th_cdf25.en+=ph_cdf25.en;
    if (include_deriv) {
      phd_cdf25.dsdT=ph.g*pi2*3.0*T*T/22.5;
      phd_cdf25.dndT=ph.g*zeta3_f<double>()/pi2*3.0*T*T;
      phd_cdf25.dndmu=0.0;
    }
  }

  if (include_deriv) {
    thd_cdf25.dndmu=ed_cdf25.dndmu;
    thd_cdf25.dndT=ed_cdf25.dndT;
    thd_cdf25.dsdT=ed_cdf25.dsdT;
    if (include_muons) {
      thd_cdf25.dndmu+=mud_cdf25.dndmu;
      thd_cdf25.dndT+=mud_cdf25.dndT;
      thd_cdf25.dsdT+=mud_cdf25.dsdT;
    }
    if (include_photons) {
      thd_cdf25.dndmu+=phd_cdf25.dndmu;
      thd_cdf25.dndT+=phd_cdf25.dndT;
      thd_cdf25.dsdT+=phd_cdf25.dsdT;
    }
  }
  
  frel_cdf25.err_nonconv=fr_en;
      
  return 0;
}

int eos_leptons_multip::pair_density_eq_fun
(size_t nv, const ubvector &x, ubvector &y, double T, double nq) {

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

long double eos_leptons_multip::pair_density_eq_ld_fun
(long double x, long double T, long double nq) {

  if (pde_from_density) {

    if (accuracy==acc_fp_25) {
      ecdf25.n=x*nq;
      int retx=electron_density_cdf25(T);
      if (retx!=0) return retx;
    } else {
      eld.n=x*nq;
      int retx=electron_density_ld(T);
      if (retx!=0) return retx;
    }
    
  } else {
    
    eld.mu=x;

    bool inc_rest_mass=false;
    if (eld.inc_rest_mass) {
      inc_rest_mass=true;
      eld.inc_rest_mass=false;
      eld.mu-=eld.m;
    }

    if (accuracy==acc_fp_25) {
      ecdf25.mu=eld.mu;
      frel_cdf25.pair_mu(ecdf25,T);
      eld.n=static_cast<double>(ecdf25.n);
      eld.ed=static_cast<double>(ecdf25.ed);
      eld.pr=static_cast<double>(ecdf25.pr);
      eld.en=static_cast<double>(ecdf25.en);
    } else {
      frel_ld.pair_mu(eld,T);
    }

    if (inc_rest_mass) {
      eld.inc_rest_mass=true;
      eld.mu+=eld.m;
      eld.ed+=eld.n*eld.m;
    }
    
  }

  if (eld.inc_rest_mass) {
    if (muld.inc_rest_mass) {
      muld.mu=eld.mu;
    } else {
      muld.mu=eld.mu-muld.m;
    }
  } else {
    if (muld.inc_rest_mass) {
      muld.mu=eld.mu+eld.m;
    } else {
      muld.mu=eld.mu+eld.m-muld.m;
    }
  }
      
  if (muld.inc_rest_mass) {
    muld.inc_rest_mass=false;
    muld.mu-=muld.m;
    frel_ld.pair_mu(muld,T);
    muld.inc_rest_mass=true;
    muld.mu+=muld.m;
    muld.ed+=muld.m*muld.n;
  } else {
    frel_ld.pair_mu(muld,T);
  }

  long double y=(eld.n+muld.n-nq)/fabs(nq);

  return y;
}

cpp_dec_float_25 eos_leptons_multip::pair_density_eq_cdf25_fun
(cpp_dec_float_25 x, cpp_dec_float_25 T, cpp_dec_float_25 nq) {

  if (pde_from_density) {

    ecdf25.n=x*nq;
    int retx=electron_density_cdf25(T);
    if (retx!=0) return retx;
    
  } else {
    
    ecdf25.mu=x;

    bool inc_rest_mass=false;
    if (ecdf25.inc_rest_mass) {
      inc_rest_mass=true;
      ecdf25.inc_rest_mass=false;
      ecdf25.mu-=ecdf25.m;
    }

    frel_cdf25.pair_mu(ecdf25,T);

    if (inc_rest_mass) {
      ecdf25.inc_rest_mass=true;
      ecdf25.mu+=ecdf25.m;
      ecdf25.ed+=ecdf25.n*ecdf25.m;
    }
    
  }

  if (ecdf25.inc_rest_mass) {
    if (mucdf25.inc_rest_mass) {
      mucdf25.mu=ecdf25.mu;
    } else {
      mucdf25.mu=ecdf25.mu-mucdf25.m;
    }
  } else {
    if (mucdf25.inc_rest_mass) {
      mucdf25.mu=ecdf25.mu+ecdf25.m;
    } else {
      mucdf25.mu=ecdf25.mu+ecdf25.m-mucdf25.m;
    }
  }
      
  if (mucdf25.inc_rest_mass) {
    mucdf25.inc_rest_mass=false;
    mucdf25.mu-=mucdf25.m;
    frel_cdf25.pair_mu(mucdf25,T);
    mucdf25.inc_rest_mass=true;
    mucdf25.mu+=mucdf25.m;
    mucdf25.ed+=mucdf25.m*mucdf25.n;
  } else {
    frel_cdf25.pair_mu(mucdf25,T);
  }

  cpp_dec_float_25 y=(ecdf25.n+mucdf25.n-nq)/fabs(nq);

  return y;
}

#endif
