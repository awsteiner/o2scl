/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  from_qq=true;

  verbose=0;
  err_on_fail=true;
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
  bool from_qq_old=from_qq;
  from_qq=false;

  // Solve zero-density gap-equations:
  bx[0]=1.0;
  bx[1]=1.0;
  bx[2]=2.0;

  thermo th;
  // Fix uninit'ed var warnings
  th.pr=0.0;
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,
                     thermo &)>(&eos_quark_njl::B0_func),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::ref(th));

  solver->msolve(3,bx,fmf);
  
  // Make the appropriate correction
  B0+=th.pr;

  // Return from_qq to its original value
  from_qq=from_qq_old;
  return 0;
}

int eos_quark_njl::calc_p(quark &u, quark &d, quark &s, thermo &th) {
  ubvector x(3);
  int ret;

  up=&u;
  down=&d;
  strange=&s;

  if (err_on_fail==false) {
    solver->err_nonconv=false;
    def_solver.def_jac.err_nonconv=false;
  } else {
    solver->err_nonconv=true;
    def_solver.def_jac.err_nonconv=true;
  }
  
  if (from_qq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
                       thermo &)>
       (&eos_quark_njl::gap_func_qq),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,th);

    ret=solver->msolve(3,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=true) in ",
                   "eos_quark_njl::calc_p().",o2scl::exc_efailed);
      }
      return ret;
    }
    
  } else {

    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
                       thermo &)>
       (&eos_quark_njl::gap_func_ms),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,th);

    ret=solver->msolve(3,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=false) in ",
                   "eos_quark_njl::calc_p().",o2scl::exc_efailed);
      }
      return ret;
    }

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
  
  if (verbose>0) {
    cout << "eos_quark_njl::calc_temp_p() from_qq=" << from_qq << endl;
  }
  
  if (err_on_fail==false) {
    solver->err_nonconv=false;
  } else {
    solver->err_nonconv=true;
  }
  
  if (from_qq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       thermo &)>
       (&eos_quark_njl::gap_func_qq_T),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,T,th);

    ret=solver->msolve(3,x,fmf);
    if (ret!=0) {
      // AWS, 5/25/22: I've found that the solver sometimes fails
      // and we can recover just by calling it a second time
      ret=solver->msolve(3,x,fmf);
      if (ret!=0) {
        if (err_on_fail) {
          O2SCL_ERR2("Solver failed (from_qq=true) in ",
                     "eos_quark_njl::calc_temp_p().",o2scl::exc_efailed);
        }
        return ret;
      }
    }
    
  } else {

    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       thermo &)>
       (&eos_quark_njl::gap_func_ms_T),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,T,th);

    ret=solver->msolve(3,x,fmf);
    if (ret!=0) {
      // AWS, 5/25/22: I've found that the solver sometimes fails
      // and we can recover just by calling it a second time.
      ret=solver->msolve(3,x,fmf);
      if (ret!=0) {
        if (err_on_fail) {
          O2SCL_ERR2("Solver failed (from_qq=false) in ",
                     "eos_quark_njl::calc_temp_p().",o2scl::exc_efailed);
        }
        return ret;
      }
    }

  }

  return ret;
}

int eos_quark_njl::calc_eq_p(quark &u, quark &d, quark &s, double &gap1,
			     double &gap2, double &gap3, thermo &th) {
  
  if (from_qq==true) {
    
    //---------------------------------------------------------------
    // Calculate everything from the quark condensates and the
    // chemical potentials, then gap1, gap2, and gap3 contain
    // Eq. 3 from Buballa99

    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;

  }
    
  if (u.mu>u.ms) {
    u.kf=sqrt(u.mu*u.mu-u.ms*u.ms);
    u.n=u.kf*u.kf*u.kf/pi2;
  } else {
    u.kf=0.0;
    u.n=0.0;
  }
  if (d.mu>d.ms) {
    d.kf=sqrt(d.mu*d.mu-d.ms*d.ms);
    d.n=d.kf*d.kf*d.kf/pi2;
  } else {
    d.kf=0.0;
    d.n=0.0;
  }
  if (s.mu>s.ms) {
    s.kf=sqrt(s.mu*s.mu-s.ms*s.ms);
    s.n=s.kf*s.kf*s.kf/pi2;
  } else {
    s.kf=0.0;
    s.n=0.0;
  }

  if (from_qq==true) {
    
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
  
  njl_bag(u);
  njl_bag(d);
  njl_bag(s);
  
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
  
  if (from_qq==true) {

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

  njl_bag(u);
  njl_bag(d);
  njl_bag(s);

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

void eos_quark_njl::njl_bag(quark &pp) {

  pp.B=3.0/pi2*(1.0/4.0*L*pow(pp.ms*pp.ms+L*L,1.5)
                -1.0/8.0*pp.ms*pp.ms*L*sqrt(pp.ms*pp.ms+L*L)
                -1.0/8.0*pow(pp.ms,4.0)*log(L+sqrt(pp.ms*pp.ms+L*L))
                +1.0/16.0*pow(pp.ms,4.0)*log(pp.ms*pp.ms));
  
  return;
}

int eos_quark_njl::gap_func_qq(size_t nv, const ubvector &x, ubvector &y,
                               thermo &th) {

  double gap1,gap2,gap3;
  
  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  
  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) return 1;
  
  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,th);
  
  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;
  
  return 0;
}

int eos_quark_njl::gap_func_ms(size_t nv, const ubvector &x, ubvector &y,
                               thermo &th) {

  double gap1,gap2,gap3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];

  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,th);

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) return 2;

  return 0;
}

int eos_quark_njl::gap_func_qq_T(size_t nv, const ubvector &x,
                                 ubvector &y, double T,
                                 thermo &th) {

  double gap1, gap2, gap3;
  
  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  
  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) {
    if (verbose>0) {
      cout << "Quark condensate positive." << endl;
    }
    return 1;
  }

  int ret=calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,
                         th,T);
  if (ret!=0) {
    if (verbose>0) {
      cout << "cetp returned non-zero." << endl;
    }
    return 3;
  }
  
  if (up->ms<=0.0 || down->ms<=0.0 || strange->ms<=0.0) {
    if (verbose>0) {
      cout << "Quark mass negative." << endl;
    }
    return 4;
  }

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    cout << "Gap equation not finite (gap_func_qq_T)." << endl;
    exit(-1);
    return 2;
  }
  
  return 0;
}

int eos_quark_njl::gap_func_ms_T(size_t nv, const ubvector &x,
                                 ubvector &y, double T,
                                 thermo &th) {

  double gap1,gap2,gap3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];

  if (x[0]<=0.0 || x[1]<=0.0 || x[2]<=0.0) {
    if (verbose>0) {
      cout << "Quark mass negative." << endl;
    }
    return 4;
  }

  int ret=calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,th,T);
  if (ret!=0) {
    if (verbose>0) {
      cout << "cetp returned non-zero." << endl;
    }
    return 3;
  }

  if (up->qq>0.0 || down->qq>0.0 || strange->qq>0.0) {
    if (verbose>0) {
      cout << "Quark condensate positive." << endl;
    }
    return 1;
  }

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    cout << "Gap equation not finite (gap_func_ms_T)." << endl;
    exit(-1);
    return 2;
  }

  return 0;
}

int eos_quark_njl::B0_func(size_t nv, const ubvector &x, ubvector &y,
                           thermo &th) {
  
  double gap1,gap2,gap3;
  
  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];
  
  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  up->mu=up->ms;
  down->mu=down->ms;
  strange->mu=strange->ms;
  
  calc_eq_p(*up,*down,*strange,gap1,gap2,gap3,th);
  
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

  // AWS, 2/17/22, I think this is commented out in order to
  // allow the effective masses to be negative
  // which may happen while solving the gap equations.
  
  /*
    if (u.ms<=0.0 || d.ms<=0.0 || s.ms<=0.0) {
    if (err_on_fail) {
    O2SCL_ERR("Effective masses negative.",o2scl::exc_einval);
    }
    return exc_efailed;
    }
  */
  
  // -----------------------------------------------------------------
  // Some of these integrals (integ_qq, and integ_density) converge
  // better when they are non-zero, so we add 1 to the integrand and
  // subtract off the contribution after the fact. This may cause
  // inaccuracies when the densities are small.

  if (from_qq==true) {
    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_qq),
			this,std::placeholders::_1,temper,u.mu,u.m,u.ms);
    iret1=it->integ_err(fqq,0.0,L,u.qq,ierr);
    u.qq-=L;
  }

  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_density),
			this,std::placeholders::_1,temper,u.mu,u.m,u.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_edensity),
			this,std::placeholders::_1,temper,u.mu,u.m,u.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_pressure),
			this,std::placeholders::_1,temper,u.mu,u.m,u.ms);
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

  if (from_qq==true) {
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_qq),
			this,std::placeholders::_1,temper,d.mu,d.m,d.ms);
    iret1=it->integ_err(fqq,0.0,L,d.qq,ierr);
    d.qq-=L;
  }
  
  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_density),
			this,std::placeholders::_1,temper,d.mu,d.m,d.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_edensity),
			this,std::placeholders::_1,temper,d.mu,d.m,d.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_pressure),
			this,std::placeholders::_1,temper,d.mu,d.m,d.ms);

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
  
  if (from_qq==true) {
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_qq),
			this,std::placeholders::_1,temper,s.mu,s.m,s.ms);
    iret1=it->integ_err(fqq,0.0,L,s.qq,ierr);
    s.qq-=L;
  }
  
  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_density),
			this,std::placeholders::_1,temper,s.mu,s.m,s.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_edensity),
			this,std::placeholders::_1,temper,s.mu,s.m,s.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl::integ_pressure),
			this,std::placeholders::_1,temper,s.mu,s.m,s.ms);
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
  
  if (from_qq==true) {
    double qqt;

    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl::integ_qq),
			  this,std::placeholders::_1,temper,u.mu,u.m,u.ms);
      iret1=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap1=-u.qq+qqt-L;
    
    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl::integ_qq),
			  this,std::placeholders::_1,temper,d.mu,d.m,d.ms);
      iret2=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap2=-d.qq+qqt-L;
    
    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl::integ_qq),
			  this,std::placeholders::_1,temper,s.mu,s.m,s.ms);
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

double eos_quark_njl::integ_qq(double x, double temper, double mu,
                               double m, double ms) {

  double en=sqrt(x*x+ms*ms);
  double ret=1.0;
  ret-=fermi_function((en-mu)/temper)+
    fermi_function((en+mu)/temper);
  ret*=-3.0*x*x/en/pi2*ms;
  ret+=1.0;
  return ret;
}

double eos_quark_njl::integ_density(double x, double temper, double mu,
                                    double m, double ms) {

  double en=sqrt(x*x+ms*ms);
  double ret=fermi_function((en-mu)/temper)-
    fermi_function((en+mu)/temper);
  ret*=3.0*x*x/pi2;
  ret+=1.0;
  return ret;
}

double eos_quark_njl::integ_edensity(double x, double temper, double mu,
                                     double m, double ms) {

  double en=sqrt(x*x+ms*ms);
  double ret=-en;
  ret+=en*fermi_function((en-mu)/temper)+
    en*fermi_function((en+mu)/temper);
  ret*=3.0/pi2*x*x;
  return ret;
}

double eos_quark_njl::integ_pressure(double x, double temper, double mu,
                                     double m, double ms) {
  
  double en=sqrt(x*x+ms*ms);
  double ret=en;
  if ((mu-en)/temper>limit) {
    ret+=mu-en;
  } else if ((mu-en)/temper<-limit) {
    ret+=0.0;
  } else {
    ret+=temper*log(1+exp((mu-en)/temper));
  }
  if ((-en-mu)/temper>limit) {
    ret+=-en-mu;
  } else if ((-en-mu)/temper<-limit) {
    ret+=0.0;
  } else {
    ret+=temper*log(1+exp((-en-mu)/temper));
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

  thermo th;
  calc_eq_p(*up,*down,*strange,g1,g2,g3,th);
  if (vac_terms==false) {
    double ret=-th.pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                      strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq;
    return ret;
  }
  return -th.pr;
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

  thermo th;
  calc_eq_temp_p(*up,*down,*strange,g1,g2,g3,th,temper);
  if (vac_terms==false) {
    double ret=-th.pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                      strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq;
    return ret;
  }
  return -th.pr;
}

int eos_quark_njl_vec::calc_eq_p_vec(quark &tu, quark &td, quark &ts,
                                     double &gap1, double &gap2, double &gap3,
                                     double &vec1, double &vec2, double &vec3,
                                     thermo &th) {

  if (from_qq==true) {

    //--------------------------------------------
    // Calculate everything from the quark 
    // condensates and the chemical potentials

    tu.ms=tu.m-4.0*G*tu.qq+2.0*K*td.qq*ts.qq;
    td.ms=td.m-4.0*G*td.qq+2.0*K*tu.qq*ts.qq;
    ts.ms=ts.m-4.0*G*ts.qq+2.0*K*td.qq*tu.qq;
    
  }
  
  // This function presumes that initial guesses have been given
  // for both the effective masses (ms) and the effective chemical
  // potentials (nu)
  
  if (tu.nu>tu.ms) {
    tu.kf=sqrt(tu.nu*tu.nu-tu.ms*tu.ms);
    tu.n=tu.kf*tu.kf*tu.kf/pi2;
  } else {
    tu.kf=0.0;
    tu.n=0.0;
  }
  if (td.nu>td.ms && td.nu<L) {
    td.kf=sqrt(td.nu*td.nu-td.ms*td.ms);
    td.n=td.kf*td.kf*td.kf/pi2;
  } else {
    td.kf=0.0;
    td.n=0.0;
  }
  if (ts.nu>ts.ms) {
    ts.kf=sqrt(ts.nu*ts.nu-ts.ms*ts.ms);
    ts.n=ts.kf*ts.kf*ts.kf/pi2;
  } else {
    ts.kf=0.0;
    ts.n=0.0;
  }
  
  if (from_qq==true) {

    gap1=3.0/pi2*(-tu.ms/2.0*L*sqrt(tu.ms*tu.ms+L*L)+pow(tu.ms,3.0)/2.0*
                  log(L+sqrt(L*L+tu.ms*tu.ms)));
    gap1+=3.0/pi2*(tu.ms/2.0*tu.kf*sqrt(tu.ms*tu.ms+tu.kf*tu.kf)-
                   pow(tu.ms,3.0)/2.0*
                   log(tu.kf+sqrt(tu.kf*tu.kf+tu.ms*tu.ms)));
    gap1-=tu.qq;
    gap2=3.0/pi2*(-td.ms/2.0*L*sqrt(td.ms*td.ms+L*L)+pow(td.ms,3.0)/2.0*
                  log(L+sqrt(L*L+td.ms*td.ms)));
    gap2+=3.0/pi2*(td.ms/2.0*td.kf*sqrt(td.ms*td.ms+td.kf*td.kf)-
                   pow(td.ms,3.0)/2.0*
                   log(td.kf+sqrt(td.kf*td.kf+td.ms*td.ms)));
    gap2-=td.qq;
    gap3=3.0/pi2*(-ts.ms/2.0*L*sqrt(ts.ms*ts.ms+L*L)+pow(ts.ms,3.0)/2.0*
                  log(L+sqrt(L*L+ts.ms*ts.ms)));
    gap3+=3.0/pi2*(ts.ms/2.0*ts.kf*sqrt(ts.ms*ts.ms+ts.kf*ts.kf)-
                   pow(ts.ms,3.0)/2.0*
                   log(ts.kf+sqrt(ts.kf*ts.kf+ts.ms*ts.ms)));
    gap3-=ts.qq;
    
  } else {
    
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
    
    gap1=-tu.ms+tu.m-4.0*G*tu.qq+2.0*K*td.qq*ts.qq;
    gap2=-td.ms+td.m-4.0*G*td.qq+2.0*K*tu.qq*ts.qq;
    gap3=-ts.ms+ts.m-4.0*G*ts.qq+2.0*K*td.qq*tu.qq;
    
  }

  fet.energy_density_zerot(tu);
  fet.energy_density_zerot(td);
  fet.energy_density_zerot(ts);
  
  njl_bag(tu);
  njl_bag(td);
  njl_bag(ts);
  
  tu.ed-=tu.B;
  td.ed-=td.B;
  ts.ed-=ts.B;

  tu.pr=-tu.ed+tu.n*tu.mu;
  td.pr=-td.ed+td.n*td.mu;
  ts.pr=-ts.ed+ts.n*ts.mu;
  
  // Fix the relationship between mu and nu
  vec1=tu.mu-4.0*GV*tu.n-tu.nu;
  vec2=td.mu-4.0*GV*td.n-td.nu;
  vec3=ts.mu-4.0*GV*ts.n-ts.nu;

  th.ed=tu.ed+td.ed+ts.ed+
    2.0*G*(tu.qq*tu.qq+td.qq*td.qq+ts.qq*ts.qq)+B0-
    4.0*K*tu.qq*td.qq*ts.qq+2.0*GV*tu.n*tu.n+2.0*GV*td.n*td.n+
    2.0*GV*ts.n*ts.n;
  th.pr=-th.ed+tu.n*tu.mu+td.n*td.mu+ts.n*ts.mu;
  th.en=0.0;
  
  return 0;
}

int eos_quark_njl_vec::calc_temp_p(quark &u, quark &d, quark &s, 
                                   double T, thermo &th) {
  ubvector x(6);
  int ret;

  up=&u;
  down=&d;
  strange=&s;
  
  if (verbose>0) {
    cout << "eos_quark_njl::calc_temp_p() from_qq=" << from_qq << endl;
  }
  
  if (err_on_fail==false) {
    solver->err_nonconv=false;
  } else {
    solver->err_nonconv=true;
  }
  
  if (from_qq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;
    x[3]=u.nu;
    x[4]=d.nu;
    x[5]=s.nu;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       thermo &)>
       (&eos_quark_njl_vec::gap_func_qq_T),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,T,th);

    ret=solver->msolve(6,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=true) in ",
                   "eos_quark_njl::calc_temp_p().",o2scl::exc_efailed);
      }
      return ret;
    }
    
  } else {

    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;
    x[3]=u.nu;
    x[4]=d.nu;
    x[5]=s.nu;

    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
                       thermo &th)>
       (&eos_quark_njl_vec::gap_func_ms_T),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,T,th);

    ret=solver->msolve(6,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=false) in ",
                   "eos_quark_njl::calc_temp_p().",o2scl::exc_efailed);
      }
      return ret;
    }

  }

  return ret;
}

int eos_quark_njl_vec::calc_eq_temp_p(quark &u, quark &d, quark &s,
                                      double &gap1, double &gap2,
                                      double &gap3, double &vec1,
                                      double &vec2, double &vec3,
                                      thermo &th, double temper) {

  double ierr;
  int iret1, iret2, iret3, iret4, iret5;
  
  if (temper<=0.0) {
    calc_eq_p_vec(u,d,s,gap1,gap2,gap3,vec1,vec2,vec3,th);
    return 0;
  }

  // AWS, 2/17/22, I think this is commented out in order to
  // allow the effective masses to be negative
  // which may happen while solving the gap equations.
  
  /*
    if (u.ms<=0.0 || d.ms<=0.0 || s.ms<=0.0) {
    if (err_on_fail) {
    O2SCL_ERR("Effective masses negative.",o2scl::exc_einval);
    }
    return exc_efailed;
    }
  */
  
  // -----------------------------------------------------------------
  // Some of these integrals (integ_qq, and integ_density) converge
  // better when they are non-zero, so we add 1 to the integrand and
  // subtract off the contribution after the fact. This may cause
  // inaccuracies when the densities are small.

  if (from_qq==true) {
    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_qq),
			this,std::placeholders::_1,temper,u.nu,u.m,u.ms);
    iret1=it->integ_err(fqq,0.0,L,u.qq,ierr);
    u.qq-=L;
  }

  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_density),
			this,std::placeholders::_1,temper,u.nu,u.m,u.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_edensity),
			this,std::placeholders::_1,temper,u.nu,u.m,u.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_pressure),
			this,std::placeholders::_1,temper,u.nu,u.m,u.ms);
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

  if (from_qq==true) {
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_qq),
			this,std::placeholders::_1,temper,d.nu,d.m,d.ms);
    iret1=it->integ_err(fqq,0.0,L,d.qq,ierr);
    d.qq-=L;
  }
  
  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_density),
			this,std::placeholders::_1,temper,d.nu,d.m,d.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_edensity),
			this,std::placeholders::_1,temper,d.nu,d.m,d.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_pressure),
			this,std::placeholders::_1,temper,d.nu,d.m,d.ms);
    iret2=it->integ_err(fde,0.0,L,d.n,ierr);
    d.n-=L;
    fet.kf_from_density(d);
    iret3=it->integ_err(fpr,0.0,L,d.pr,ierr);
    iret4=it->integ_err(fed,0.0,L,d.ed,ierr);
    if (iret1!=0 || iret2!=0 || iret3!=0 || iret4!=0) {
      O2SCL_ERR("Down quark failed in eos_quark_njl_vec::calc_eq_temp_p().",
		exc_efailed);
    }
  }
  
  if (from_qq==true) {
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
    iret1=0;
  } else {
    funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_qq),
			this,std::placeholders::_1,temper,s.nu,s.m,s.ms);
    iret1=it->integ_err(fqq,0.0,L,s.qq,ierr);
    s.qq-=L;
  }
  
  {
    funct fde=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_density),
			this,std::placeholders::_1,temper,s.nu,s.m,s.ms);
    funct fed=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_edensity),
			this,std::placeholders::_1,temper,s.nu,s.m,s.ms);
    funct fpr=std::bind(std::mem_fn<double(double,double,double,
                                           double,double)>
			(&eos_quark_njl_vec::integ_pressure),
			this,std::placeholders::_1,temper,s.nu,s.m,s.ms);
    iret2=it->integ_err(fde,0.0,L,s.n,ierr);
    s.n-=L;
    fet.kf_from_density(s);
    iret3=it->integ_err(fpr,0.0,L,s.pr,ierr);
    iret4=it->integ_err(fed,0.0,L,s.ed,ierr);
    if (iret1!=0 || iret2!=0 || iret3!=0 || iret4!=0) {
      O2SCL_ERR("Strange quark failed in eos_quark_njl_vec::calc_eq_temp_p().",
		exc_efailed);
    }
  }

  u.en=(u.ed+u.pr-u.nu*u.n)/temper;
  d.en=(d.ed+d.pr-d.nu*d.n)/temper;
  s.en=(s.ed+s.pr-s.nu*s.n)/temper;

  th.ed=u.ed+d.ed+s.ed+2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)+
    B0-4.0*K*u.qq*d.qq*s.qq+2.0*GV*u.n*u.n+2.0*GV*d.n*d.n+
    2.0*GV*s.n*s.n;
  th.pr=u.pr+d.pr+s.pr-2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)-
    B0+4.0*K*u.qq*d.qq*s.qq+2.0*GV*u.n*u.n+2.0*GV*d.n*d.n+
    2.0*GV*s.n*s.n;
  th.en=(th.ed+th.pr-u.n*u.mu-d.n*d.mu-s.n*s.mu)/temper;
  
  if (from_qq==true) {
    double qqt;

    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl_vec::integ_qq),
			  this,std::placeholders::_1,temper,u.nu,u.m,u.ms);
      iret1=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap1=-u.qq+qqt-L;
    
    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl_vec::integ_qq),
			  this,std::placeholders::_1,temper,d.nu,d.m,d.ms);
      iret2=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap2=-d.qq+qqt-L;
    
    {
      funct fqq=std::bind(std::mem_fn<double(double,double,double,
                                             double,double)>
			  (&eos_quark_njl_vec::integ_qq),
			  this,std::placeholders::_1,temper,s.nu,s.m,s.ms);
      iret3=it->integ_err(fqq,0.0,L,qqt,ierr);
    }
    gap3=-s.qq+qqt-L;

    if (iret1!=0 || iret2!=0 || iret3!=0) {
      O2SCL_ERR("Strange quark failed in eos_quark_njl_vec::calc_eq_temp_p().",
		exc_efailed);
    }
  } else {
    gap1=-u.ms+u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    gap2=-d.ms+d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    gap3=-s.ms+s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
  }

  vec1=-u.nu+u.mu-4.0*GV*u.n;
  vec2=-d.nu+d.mu-4.0*GV*d.n;
  vec3=-s.nu+s.mu-4.0*GV*s.n;

  return 0;
}

int eos_quark_njl_vec::gap_func_ms_vec(size_t nv, const ubvector &x,
                                       ubvector &y, thermo &th) {

  double gap1, gap2, gap3, vec1, vec2, vec3;

  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];
  up->nu=x[3];
  down->nu=x[4];
  strange->nu=x[5];

  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) return 1;

  calc_eq_p_vec(*up,*down,*strange,gap1,gap2,gap3,
                vec1,vec2,vec3,th);

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

int eos_quark_njl_vec::gap_func_qq_vec(size_t nv, const ubvector &x,
                                       ubvector &y,
                                       thermo &th) {
  
  double gap1, gap2, gap3, vec1, vec2, vec3;

  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  up->nu=x[3];
  down->nu=x[4];
  strange->nu=x[5];

  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) return 1;

  calc_eq_p_vec(*up,*down,*strange,gap1,gap2,gap3,
                vec1,vec2,vec3,th);

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

  up=&u;
  down=&d;
  strange=&s;
  
  if (err_on_fail==false) {
    solver->err_nonconv=false;
    def_solver.def_jac.err_nonconv=false;
  } else {
    solver->err_nonconv=true;
    def_solver.def_jac.err_nonconv=true;
  }

  if (from_qq==true) {

    x[0]=u.qq;
    x[1]=d.qq;
    x[2]=s.qq;
    x[3]=u.nu;
    x[4]=d.nu;
    x[5]=s.nu;
    
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,thermo &)>
       (&eos_quark_njl_vec::gap_func_qq_vec),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,th);
    
    ret=solver->msolve(6,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=true) in ",
                   "eos_quark_njl_vec::calc_p().",o2scl::exc_efailed);
      }
      return ret;
    }
    
  } else {
    
    x[0]=u.ms;
    x[1]=d.ms;
    x[2]=s.ms;
    x[3]=u.nu;
    x[4]=d.nu;
    x[5]=s.nu;
    
    mm_funct fmf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &, thermo &)>
       (&eos_quark_njl_vec::gap_func_ms_vec),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,th);
    
    ret=solver->msolve(6,x,fmf);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Solver failed (from_qq=false) in ",
                   "eos_quark_njl_vec::calc_p().",o2scl::exc_efailed);
      }
      return ret;
    }
  }
  
  return ret;
}

double eos_quark_njl_vec::f_therm_pot_vec(double qqu, double qqd, double qqs,
                                          double msu, double msd, double mss,
                                          double nuu, double nud, double nus,
                                          bool vac_terms) {
  
  double g1, g2, g3, h1, h2, h3;
  
  up->qq=qqu;
  down->qq=qqd;
  strange->qq=qqs;
  up->ms=msu;
  down->ms=msd;
  strange->ms=mss;
  up->nu=nuu;
  down->nu=nud;
  strange->nu=nus;

  thermo th;
  calc_eq_p_vec(*up,*down,*strange,g1,g2,g3,h1,h2,h3,th);
  
  if (vac_terms==false) {
    double ret=-th.pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                      strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq+
      2.0*GV*up->n*up->n+2.0*GV*down->n*down->n+2.0*GV*strange->n*strange->n;
    return ret;
  }
  return -th.pr;
}

double eos_quark_njl_vec::f_therm_pot_T_vec(double qqu, double qqd, double qqs,
                                            double msu, double msd, double mss,
                                            double nuu, double nud, double nus,
                                            double T, bool vac_terms) {
  
  double g1, g2, g3, h1, h2, h3;
  
  up->qq=qqu;
  down->qq=qqd;
  strange->qq=qqs;
  up->ms=msu;
  down->ms=msd;
  strange->ms=mss;
  up->nu=nuu;
  down->nu=nud;
  strange->nu=nus;

  thermo th;
  calc_eq_temp_p(*up,*down,*strange,g1,g2,g3,h1,h2,h3,th,T);
  
  if (vac_terms==false) {
    double ret=-th.pr-2.0*G*(up->qq*up->qq+down->qq*down->qq+
                                      strange->qq*strange->qq)+
      4.0*K*up->qq*down->qq*strange->qq+
      2.0*GV*up->n*up->n+2.0*GV*down->n*down->n+2.0*GV*strange->n*strange->n;
    return ret;
  }
  return -th.pr;
}

int eos_quark_njl_vec::gap_func_qq_T(size_t nv, const ubvector &x,
                                     ubvector &y, double T,
                                     thermo &th) {

  double gap1, gap2, gap3, vec1, vec2, vec3;
  
  up->qq=x[0];
  down->qq=x[1];
  strange->qq=x[2];
  up->nu=x[3];
  down->nu=x[4];
  strange->nu=x[5];
  
  if (x[0]>0.0 || x[1]>0.0 || x[2]>0.0) {
    if (verbose>0) {
      cout << "Quark condensate positive." << endl;
    }
    return 1;
  }

  int ret=calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,
                         vec1,vec2,vec3,th,T);
  
  if (ret!=0) {
    if (verbose>0) {
      cout << "cetp returned non-zero." << endl;
    }
    return 3;
  }
  
  if (up->ms<=0.0 || down->ms<=0.0 || strange->ms<=0.0) {
    if (verbose>0) {
      cout << "Quark mass negative." << endl;
    }
    return 4;
  }

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;
  y[3]=vec1;
  y[4]=vec2;
  y[5]=vec3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    cout << "Gap equation not finite (gap_func_qq_T)." << endl;
    exit(-1);
    return 2;
  }
  
  return 0;
}

int eos_quark_njl_vec::gap_func_ms_T(size_t nv, const ubvector &x,
                                     ubvector &y, double T,
                                     thermo &th) {

  double gap1, gap2, gap3, vec1, vec2, vec3;
  
  up->ms=x[0];
  down->ms=x[1];
  strange->ms=x[2];
  up->nu=x[3];
  down->nu=x[4];
  strange->nu=x[5];
  
  if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0) {
    if (verbose>0) {
      cout << "Masses negative." << endl;
    }
    return 1;
  }

  int ret=calc_eq_temp_p(*up,*down,*strange,gap1,gap2,gap3,
                         vec1,vec2,vec3,th,T);
  
  if (ret!=0) {
    if (verbose>0) {
      cout << "cetp returned non-zero." << endl;
    }
    return 3;
  }
  
  if (up->qq>0.0 || down->qq>0.0 || strange->qq>0.0) {
    if (verbose>0) {
      cout << "Quark condensate positive." << endl;
    }
    return 1;
  }

  y[0]=gap1;
  y[1]=gap2;
  y[2]=gap3;
  y[3]=vec1;
  y[4]=vec2;
  y[5]=vec3;

  if (!std::isfinite(y[0]) || !std::isfinite(y[1]) ||
      !std::isfinite(y[2])) {
    cout << "Gap equation not finite (gap_func_qq_T)." << endl;
    exit(-1);
    return 2;
  }
  
  return 0;
}

