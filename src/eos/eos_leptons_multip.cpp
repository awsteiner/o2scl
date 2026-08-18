/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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

#include <o2scl/eos_leptons_multip.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

#ifdef O2SCL_SET_MULTIP

eos_leptons_multip::eos_leptons_multip() {

  cu_ld.default_conversions();
  cu_fp25.default_conversions();
  
  eld.init(cu_ld.convert("kg","1/fm",mass_electron_f<long double>()),2);
  muld.init(cu_ld.convert("kg","1/fm",mass_muon_f<long double>()),2);
  tauld.init(cu_ld.convert("kg","1/fm",mass_tau_f<long double>()),2);
  
  efp25.init(cu_fp25.convert("kg","1/fm",
                               mass_electron_f<o2fp_25>()),2);
  mufp25.init(cu_fp25.convert("kg","1/fm",
                                mass_muon_f<o2fp_25>()),2);
  taufp25.init(cu_fp25.convert("kg","1/fm",
                                 mass_tau_f<o2fp_25>()),2);

  ph_ld.init(0,2);
  ph_fp25.init(0,2);
  
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
    efp25.inc_rest_mass=true;
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
    efp25.inc_rest_mass=true;
    eld.mu+=eld.m;
    eld.ed+=eld.m*eld.n;
  }

  return retx;
}

int eos_leptons_multip::electron_density_fp25(o2fp_25 T) {

  int retx;

  // I find that the calculation without the rest mass is a bit more
  // stable, so we use that method and add the rest mass back in
  // later if necessary.
  bool inc_rest_mass=false;
  if (efp25.inc_rest_mass) {
    inc_rest_mass=true;
    efp25.inc_rest_mass=false;
    efp25.mu-=efp25.m;
  }
  
  retx=frel_fp25.pair_density(efp25,T);
  
  if (inc_rest_mass) {
    efp25.inc_rest_mass=true;
    efp25.mu+=efp25.m;
    efp25.ed+=efp25.m*efp25.n;
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
	fdrel_fp25.multip=true;
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

int eos_leptons_multip::pair_density_eq_fp25
(o2fp_25 nq, o2fp_25 T) {
      
  bool fr_en=frel_fp25.err_nonconv;
  frel_fp25.err_nonconv=false;
  
  int retx;
  if (include_muons) {
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_fp25(): "
                << "with muons, pde_from_density="
                << pde_from_density << std::endl;
    }

    o2fp_25 x, y;
    if (pde_from_density) {
      x=efp25.n/nq;
    } else {
      x=efp25.mu;
    }

    funct_fp25 mf=std::bind
      (std::mem_fn<o2fp_25(o2fp_25,
                                    o2fp_25,o2fp_25)>
       (&eos_leptons_multip::pair_density_eq_fp25_fun),
       this,std::placeholders::_1,T,nq);
    size_t maxj=10;
    if (verbose>1) {
      cout << "Initial guess: " << x << endl;
    }
    int mret=rc_fp25.solve(x,mf);
    for(size_t j=0;j<maxj && mret!=0;j++) {
      if (verbose>1) {
        cout << "Attempt " << j+2 << " with guess " << x
             << " and tolerance: " << rc_fp25.tol_rel << std::endl;
      }
      mret=rc_fp25.solve(x,mf);
      rc_fp25.tol_rel*=pow(10.0,1.0/2.0);
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
    efp25.n=x*nq;
    
    if (include_deriv) {
      fdrel_fp25.multip=true;
      
      fermion_deriv_fp25 fd;
      fd=efp25;
      fdrel_fp25.pair_mu(fd,T);
      ed_fp25.dndmu=fd.dndmu;
      ed_fp25.dndT=fd.dndT;
      ed_fp25.dsdT=fd.dsdT;
      fd=mufp25;
      fdrel_fp25.pair_mu(fd,T);
      mud_fp25.dndmu=fd.dndmu;
      mud_fp25.dndT=fd.dndT;
      mud_fp25.dsdT=fd.dsdT;
    }
    
  } else {
    
    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_fp25(): No muons."
                << std::endl;
    }
    mufp25.n=0.0;
    
    efp25.n=nq;
    retx=electron_density_fp25(T);

    if (verbose>1) {
      std::cout << "eos_leptons_multip::pair_density_eq_fp25(): "
                << "Return value " << retx << std::endl;
    }
    
    if (include_deriv) {
      if (verbose>1) {
        std::cout << "eos_leptons_multip::pair_density_eq_fp25(): "
                  << "Including derivatives." << std::endl;
      }
      fermion_deriv_fp25 fd;
      fd=efp25;
      fdrel_fp25.multip=true;
      fdrel_fp25.pair_mu(fd,T);
      ed_fp25.dndmu=fd.dndmu;
      ed_fp25.dndT=fd.dndT;
      ed_fp25.dsdT=fd.dsdT;

    }
        
  }
      
  th_fp25.ed=efp25.ed;
  th_fp25.pr=efp25.pr;
  th_fp25.en=efp25.en;

  if (include_muons) {
    th_fp25.ed+=mufp25.ed;
    th_fp25.pr+=mufp25.pr;
    th_fp25.en+=mufp25.en;
  }
      
  if (include_photons) {
    ph_fp25.massless_calc(T);
    th_fp25.ed+=ph_fp25.ed;
    th_fp25.pr+=ph_fp25.pr;
    th_fp25.en+=ph_fp25.en;
    if (include_deriv) {
      phd_fp25.dsdT=ph.g*pi2*3.0*T*T/22.5;
      phd_fp25.dndT=ph.g*zeta3_f<double>()/pi2*3.0*T*T;
      phd_fp25.dndmu=0.0;
    }
  }

  if (include_deriv) {
    thd_fp25.dndmu=ed_fp25.dndmu;
    thd_fp25.dndT=ed_fp25.dndT;
    thd_fp25.dsdT=ed_fp25.dsdT;
    if (include_muons) {
      thd_fp25.dndmu+=mud_fp25.dndmu;
      thd_fp25.dndT+=mud_fp25.dndT;
      thd_fp25.dsdT+=mud_fp25.dsdT;
    }
    if (include_photons) {
      thd_fp25.dndmu+=phd_fp25.dndmu;
      thd_fp25.dndT+=phd_fp25.dndT;
      thd_fp25.dsdT+=phd_fp25.dsdT;
    }
  }
  
  frel_fp25.err_nonconv=fr_en;
      
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
      efp25.n=x[0]*nq;
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
      efp25.mu=e.mu;
      frel_fp25.pair_mu(efp25,T);
      e.n=static_cast<double>(efp25.n);
      e.ed=static_cast<double>(efp25.ed);
      e.pr=static_cast<double>(efp25.pr);
      e.en=static_cast<double>(efp25.en);
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
      efp25.n=x*nq;
      int retx=electron_density_fp25(T);
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
      efp25.mu=eld.mu;
      frel_fp25.pair_mu(efp25,T);
      eld.n=static_cast<double>(efp25.n);
      eld.ed=static_cast<double>(efp25.ed);
      eld.pr=static_cast<double>(efp25.pr);
      eld.en=static_cast<double>(efp25.en);
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

o2fp_25 eos_leptons_multip::pair_density_eq_fp25_fun
(o2fp_25 x, o2fp_25 T, o2fp_25 nq) {

  if (pde_from_density) {

    efp25.n=x*nq;
    int retx=electron_density_fp25(T);
    if (retx!=0) return retx;
    
  } else {
    
    efp25.mu=x;

    bool inc_rest_mass=false;
    if (efp25.inc_rest_mass) {
      inc_rest_mass=true;
      efp25.inc_rest_mass=false;
      efp25.mu-=efp25.m;
    }

    frel_fp25.pair_mu(efp25,T);

    if (inc_rest_mass) {
      efp25.inc_rest_mass=true;
      efp25.mu+=efp25.m;
      efp25.ed+=efp25.n*efp25.m;
    }
    
  }

  if (efp25.inc_rest_mass) {
    if (mufp25.inc_rest_mass) {
      mufp25.mu=efp25.mu;
    } else {
      mufp25.mu=efp25.mu-mufp25.m;
    }
  } else {
    if (mufp25.inc_rest_mass) {
      mufp25.mu=efp25.mu+efp25.m;
    } else {
      mufp25.mu=efp25.mu+efp25.m-mufp25.m;
    }
  }
      
  if (mufp25.inc_rest_mass) {
    mufp25.inc_rest_mass=false;
    mufp25.mu-=mufp25.m;
    frel_fp25.pair_mu(mufp25,T);
    mufp25.inc_rest_mass=true;
    mufp25.mu+=mufp25.m;
    mufp25.ed+=mufp25.m*mufp25.n;
  } else {
    frel_fp25.pair_mu(mufp25,T);
  }

  o2fp_25 y=(efp25.n+mufp25.n-nq)/fabs(nq);

  return y;
}

#endif
