/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>
#include <o2scl/eos_had_sym4.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
 
int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  thermo th;

  fermion n, p;

  bool optv;

  eos_had_apr apr;
  eos_had_sym4_apr apr4;
  eos_had_sym4 s4e;
  eos_had_rmf rmf;
  eos_had_rmf rmf_b;
  eos_had_rmf rmf_c;
  eos_had_sym4_rmf rmf4;
  eos_had_sym4_rmf rmf4_b;
  eos_had_sym4_rmf rmf4_c;
  eos_had_sym4_rmf rnew;
  eos_had_sym4_rmf rhi;
  eos_had_sym4_rmf rlo;
  eos_had_skyrme sk;
  eos_had_sym4_skyrme sk4;
  eos_had_potential mdi;
  eos_had_sym4_mdi mdi4;

  nstar_cold cnx;
  
  /// ----------------------------------------------------------------
  /// Initialize particles

  n.init(939.0/hc_mev_fm,2.0);
  p.init(939.0/hc_mev_fm,2.0);
  n.non_interacting=false;
  p.non_interacting=false;

  o2scl_hdf::skyrme_load(sk,"../../data/o2scl/skdata/SLy230a.o2",1);
  o2scl_hdf::skyrme_load(sk4,"../../data/o2scl/skdata/SLy230a.o2",1);

  /// Initialize

  o2scl_hdf::rmf_load(rmf,"../../data/o2scl/rmfdata/NL4.o2",1);
  o2scl_hdf::rmf_load(rmf4,"../../data/o2scl/rmfdata/NL4.o2",1);
  rmf.mnuc=939.0/197.33;
  rmf4.mnuc=939.0/197.33;
  o2scl_hdf::rmf_load(rmf_b,"../../data/o2scl/rmfdata/NL4.o2",1);
  o2scl_hdf::rmf_load(rmf4_b,"../../data/o2scl/rmfdata/NL4.o2",1);
  rmf_b.mnuc=939.0/197.33;
  rmf4_b.mnuc=939.0/197.33;
  o2scl_hdf::rmf_load(rmf_c,"../../data/o2scl/rmfdata/NL4.o2",1);
  o2scl_hdf::rmf_load(rmf4_c,"../../data/o2scl/rmfdata/NL4.o2",1);
  rmf_c.mnuc=939.0/197.33;
  rmf4_c.mnuc=939.0/197.33;

  rmf.esym=34.0/hc_mev_fm;
  rmf.n0=0.16;
  rmf.eoa=-16.0/hc_mev_fm;
  rmf.comp=230.0/hc_mev_fm;
  rmf.msom=0.8;
  rmf.xi=0.0;
  rmf.zeta=0.0;
  rmf4.esym=34.0/hc_mev_fm;
  rmf4.n0=0.16;
  rmf4.eoa=-16.0/hc_mev_fm;
  rmf4.comp=230.0/hc_mev_fm;
  rmf4.msom=0.8;
  rmf4.xi=0.0;
  rmf4.zeta=0.0;
  rmf_b.esym=34.0/hc_mev_fm;
  rmf_b.n0=0.16;
  rmf_b.eoa=-16.0/hc_mev_fm;
  rmf_b.comp=230.0/hc_mev_fm;
  rmf_b.msom=0.8;
  rmf_b.xi=1.5;
  rmf_b.zeta=0.0;
  rmf4_b.esym=34.0/hc_mev_fm;
  rmf4_b.n0=0.16;
  rmf4_b.eoa=-16.0/hc_mev_fm;
  rmf4_b.comp=230.0/hc_mev_fm;
  rmf4_b.msom=0.8;
  rmf4_b.xi=1.5;
  rmf4_b.zeta=0.0;
  rmf_c.esym=34.0/hc_mev_fm;
  rmf_c.n0=0.16;
  rmf_c.eoa=-16.0/hc_mev_fm;
  rmf_c.comp=230.0/hc_mev_fm;
  rmf_c.msom=0.8;
  rmf_c.xi=0.0;
  rmf_c.zeta=0.0;
  rmf_c.a2=3.0;
  rmf4_c.esym=34.0/hc_mev_fm;
  rmf4_c.n0=0.16;
  rmf4_c.eoa=-16.0/hc_mev_fm;
  rmf4_c.comp=230.0/hc_mev_fm;
  rmf4_c.msom=0.8;
  rmf4_c.xi=0.0;
  rmf4_c.zeta=0.0;
  rmf4_c.a2=3.0;
  
  rmf.def_mroot.def_jac.set_epsrel(1.0e-4);
  rmf4.def_mroot.def_jac.set_epsrel(1.0e-4);
  rmf_b.def_mroot.def_jac.set_epsrel(1.0e-4);
  rmf_c.def_mroot.def_jac.set_epsrel(1.0e-4);
  rmf4_b.def_mroot.def_jac.set_epsrel(1.0e-4);
  rmf4_c.def_mroot.def_jac.set_epsrel(1.0e-4);
  rnew.def_mroot.def_jac.set_epsrel(1.0e-4);
  rlo.def_mroot.def_jac.set_epsrel(1.0e-4);
  rhi.def_mroot.def_jac.set_epsrel(1.0e-4);

  double dx, eoat;

  n.mu=4.8;
  p.mu=4.8;
  rmf.set_n_and_p(n,p);
  rmf.set_fields(0.1,0.07,-0.001);
  rmf.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf4.set_n_and_p(n,p);
  rmf4.set_fields(0.1,0.07,-0.001);
  rmf4.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf.set_fields(0.1,0.07,-0.001);
  rmf.saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf4.set_fields(0.1,0.07,-0.001);
  rmf4.saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf_b.set_n_and_p(n,p);
  rmf_b.set_fields(0.1,0.07,-0.001);
  rmf_b.fix_saturation();
  n.mu=4.8;
  p.mu=4.8;
  rmf4_b.set_n_and_p(n,p);
  rmf4_b.set_fields(0.1,0.07,-0.001);
  rmf4_b.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf_b.set_fields(0.1,0.07,-0.001);
  rmf_b.saturation();
  n.mu=4.8;
  p.mu=4.8;
  rmf4_b.set_fields(0.1,0.07,-0.001);
  rmf4_b.saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf_c.set_n_and_p(n,p);
  rmf_c.set_fields(0.1,0.07,-0.001);
  rmf_c.fix_saturation();
  n.mu=4.8;
  p.mu=4.8;
  rmf4_c.set_n_and_p(n,p);
  rmf4_c.set_fields(0.1,0.07,-0.001);
  rmf4_c.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rmf_c.set_fields(0.1,0.07,-0.001);
  rmf_c.saturation();
  n.mu=4.8;
  p.mu=4.8;
  rmf4_c.set_fields(0.1,0.07,-0.001);
  rmf4_c.saturation();

  mdi.Au=-95.98/hc_mev_fm;
  mdi.Al=-120.57/hc_mev_fm;
  mdi.B=106.35/hc_mev_fm;
  mdi.Cu=-103.40/hc_mev_fm;
  mdi.Cl=-11.70/hc_mev_fm;
  mdi.sigma=4.0/3.0;
  mdi.x=0.0;
  mdi.rho0=0.16;
  mdi.Lambda=cbrt(1.5*pi2*mdi.rho0);
  mdi.form=mdi.mdi_form;

  mdi4.Au=-95.98/hc_mev_fm;
  mdi4.Al=-120.57/hc_mev_fm;
  mdi4.B=106.35/hc_mev_fm;
  mdi4.Cu=-103.40/hc_mev_fm;
  mdi4.Cl=-11.70/hc_mev_fm;
  mdi4.sigma=4.0/3.0;
  mdi4.x=0.0;
  mdi4.rho0=0.16;
  mdi4.Lambda=cbrt(1.5*pi2*mdi4.rho0);
  mdi4.form=mdi4.mdi_form;
  
  o2scl_hdf::rmf_load(rnew,"../../data/o2scl/rmfdata/NL3.o2",1);

  rnew.mnuc=939.0/197.33;
  rnew.esym=34.0/hc_mev_fm;
  rnew.n0=0.16;
  rnew.eoa=-16.0/hc_mev_fm;
  rnew.comp=230.0/hc_mev_fm;
  rnew.msom=0.8;
  rnew.xi=0.0;
  rnew.zeta=0.02;
  rnew.b1=1.707452;
  rnew.a1=-0.1826542;
  rnew.a2=0.171308;

  n.mu=4.8;
  p.mu=4.8;
  rnew.set_n_and_p(n,p);
  rnew.set_fields(0.1,0.07,-0.001);
  rnew.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rnew.set_fields(0.1,0.07,-0.001);
  rnew.saturation();

  o2scl_hdf::rmf_load(rlo,"../../data/o2scl/rmfdata/NL3.o2",1);

  rlo.mnuc=939.0/197.33;
  rlo.esym=34.0/hc_mev_fm;
  rlo.n0=0.16;
  rlo.eoa=-16.0/hc_mev_fm;
  rlo.comp=230.0/hc_mev_fm;
  rlo.msom=0.8;
  rlo.zeta=3.800923e-02;
  rlo.xi=1.498865e+00;
  rlo.b1=5.251033e-02;
  rlo.a1=3.525083e-01;
  rlo.a2=1.051629e+00;
  rlo.b2=-3.148076e-01;
  rlo.a3=1.666687e+00;
  rlo.a4=-8.033841e-01;

  n.mu=4.8;
  p.mu=4.8;
  rlo.set_n_and_p(n,p);
  rlo.set_fields(0.1,0.07,-0.001);
  rlo.fix_saturation();

  n.mu=4.8;
  p.mu=4.8;
  rlo.set_fields(0.1,0.07,-0.001);
  rlo.saturation();

  o2scl_hdf::rmf_load(rhi,"../../data/o2scl/rmfdata/NL3.o2",1);

  rhi.mnuc=939.0/197.33;
  rhi.esym=34.0/hc_mev_fm;
  rhi.n0=0.16;
  rhi.eoa=-16.0/hc_mev_fm;
  rhi.comp=230.0/hc_mev_fm;
  rhi.msom=0.8;
  rhi.zeta=5.879708e-02;
  rhi.xi=2.869631e-02;
  rhi.b1=2.719245e-01;
  rhi.a1=-1.790506e-01;
  rhi.a2=-6.660515e-01;
  rhi.b2=4.390330e-01;
  rhi.a3=2.527117e-01;
  rhi.a4=8.027970e-01;

  // For some reason this next call requires
  // a couple more iterations than the default

  rhi.def_sat_mroot.ntrial*=10;
  n.mu=4.8;
  p.mu=4.8;
  rhi.set_n_and_p(n,p);
  rhi.set_fields(0.1,0.07,-0.001);
  rhi.fix_saturation();
  rhi.def_sat_mroot.ntrial/=10;

  n.mu=4.8;
  p.mu=4.8;
  rhi.set_fields(0.1,0.07,-0.001);
  rhi.saturation();

  double t1, t2, t3;

  cout << endl;

  /// ----------------------------------------------------------------
  /// Test RMF EOS

  cout << rmf.n0 << " " << rmf.eoa*hc_mev_fm << endl;
  cout << rmf4.n0 << " " << rmf4.eoa*hc_mev_fm << '\n' << endl;
  cout << rmf_b.n0 << " " << rmf_b.eoa*hc_mev_fm << endl;
  cout << rmf4_b.n0 << " " << rmf4_b.eoa*hc_mev_fm << '\n' << endl;
  cout << rmf_c.n0 << " " << rmf_c.eoa*hc_mev_fm << endl;
  cout << rmf4_c.n0 << " " << rmf4_c.eoa*hc_mev_fm << '\n' << endl;

  t.test_rel(rmf.n0,0.16,1.0e-4,"");
  t.test_rel(rmf4.n0,0.16,1.0e-4,"");
  t.test_rel(rmf.eoa*hc_mev_fm,-16.0,1.0e-4,"");
  t.test_rel(rmf4.eoa*hc_mev_fm,-16.0,1.0e-4,"");

  t.test_rel(rmf_b.n0,0.16,1.0e-4,"");
  t.test_rel(rmf4_b.n0,0.16,1.0e-4,"");
  t.test_rel(rmf_b.eoa*hc_mev_fm,-16.0,1.0e-4,"");
  t.test_rel(rmf4_b.eoa*hc_mev_fm,-16.0,1.0e-4,"");

  t.test_rel(rmf_c.n0,0.16,1.0e-4,"");
  t.test_rel(rmf4_c.n0,0.16,1.0e-4,"");
  t.test_rel(rmf_c.eoa*hc_mev_fm,-16.0,1.0e-4,"");
  t.test_rel(rmf4_c.eoa*hc_mev_fm,-16.0,1.0e-4,"");

  /// ----------------------------------------------------------------
  /// Test MDI EOS

  n.n=0.08;
  p.n=0.08;
  //mdi4.test_separation(n,p,t);

  n.n=0.12;
  p.n=0.04;
  //mdi4.test_separation(n,p,t);

  n.n=0.16;
  p.n=0.0;
  //mdi4.test_separation(n,p,t);
  cout << endl;

  /// ----------------------------------------------------------------
  /// Try out APR EOS

  s4e.set_base_eos(apr4);
  s4e.alpha=3.0;
  
  cout << "APR:" << endl;

  /// Check nuclear matter

  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  apr.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;  
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  apr4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"apr nuc");
  t.test_rel(n.mu,t2,1.0e-6,"apr nuc");
  t.test_rel(p.mu,t3,1.0e-6,"apr nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"apr nuc");
  t.test_rel(n.mu,t2,1.0e-6,"apr nuc");
  t.test_rel(p.mu,t3,1.0e-6,"apr nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  apr.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  apr4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"apr neut");
  t.test_rel(n.mu,t2,1.0e-5,"apr neut");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"apr neut");
  t.test_rel(n.mu,t2,1.0e-5,"apr neut");
  
  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  apr.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  apr4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"apr nmix");
  t.test_rel(n.mu,t2,1.0e-5,"apr nmix");
  t.test_rel(p.mu,t3,1.0e-5,"apr nmix");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"apr nmix");
  t.test_rel(n.mu,t2,1.0e-5,"apr nmix");
  t.test_rel(p.mu,t3,1.0e-5,"apr nmix");

  /// Check neutron-rich matter with non-trivial alpha

  cout << "Neutron-rich matter, alpha!=3." << endl;
  n.n=0.12;
  p.n=0.04;
  s4e.alpha=2.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  cout << endl;

  /// ----------------------------------------------------------------
  /// A Skyrme EOS

  cout << "SLy230a:" << endl;
  s4e.set_base_eos(sk4);
  s4e.alpha=3.0;

  /// Check nuclear matter

  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  sk.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  sk4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk nuc");
  t.test_rel(n.mu,t2,1.0e-5,"sk nuc");
  t.test_rel(p.mu,t3,1.0e-5,"sk nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk nuc");
  t.test_rel(n.mu,t2,1.0e-5,"sk nuc");
  t.test_rel(p.mu,t3,1.0e-5,"sk nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  sk.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  sk4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk neut");
  t.test_rel(n.mu,t2,1.0e-5,"sk neut");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk neut");
  t.test_rel(n.mu,t2,1.0e-5,"sk neut");

  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  sk.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  sk4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk nrich");
  t.test_rel(n.mu,t2,1.0e-5,"sk nrich");
  t.test_rel(p.mu,t3,1.0e-5,"sk nrich");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"sk nrich");
  t.test_rel(n.mu,t2,1.0e-5,"sk nrich");
  t.test_rel(p.mu,t3,1.0e-5,"sk nrich");

  /// Check neutron-rich matter with non-trivial alpha

  cout << "Neutron-rich matter, alpha!=3." << endl;
  n.n=0.12;
  p.n=0.04;
  s4e.alpha=2.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " " 
       << n.mu << " " << p.mu << endl;
  cout << endl;

  /// ----------------------------------------------------------------
  /// RMF EOS

  s4e.set_base_eos(rmf4);
  s4e.alpha=3.0;

  /// Check nuclear matter
  
  cout << "RMF:" << endl;
  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  rmf.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  rmf4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  rmf.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.16;
  p.n=0.0;
  rmf4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf neut");
  n.n=0.16;
  p.n=0.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf neut");

  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  rmf.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.12;
  p.n=0.04;
  rmf4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf nrich");
  t.test_rel(n.mu,t2,1.0e-6,"rmf nrich");
  t.test_rel(p.mu,t3,1.0e-6,"rmf nrich");
  
  n.n=0.12;
  p.n=0.04;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,4.0e-2,"rmf nrich");
  t.test_rel(n.mu,t2,4.0e-2,"rmf nrich");
  t.test_rel(p.mu,t3,4.0e-2,"rmf nrich");
  
  /// Check neutron-rich matter with non-trivial alpha

  cout << "Neutron-rich matter, alpha!=3." << endl;
  n.n=0.12;
  p.n=0.04;
  s4e.alpha=2.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  cout << endl;
  
  /// ----------------------------------------------------------------
  /// RMF (B) EOS

  s4e.set_base_eos(rmf4_b);
  s4e.alpha=3.0;

  /// Check nuclear matter
  
  cout << "RMF_B:" << endl;
  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  rmf_b.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  rmf4_b.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_b nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_b nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf_b nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_b nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_b nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf_b nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  rmf_b.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.16;
  p.n=0.0;
  rmf4_b.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_b neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_b neut");
  n.n=0.16;
  p.n=0.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_b neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_b neut");

  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  rmf_b.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.12;
  p.n=0.04;
  rmf4_b.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_b nrich");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_b nrich");
  t.test_rel(p.mu,t3,1.0e-6,"rmf_b nrich");
  n.n=0.12;
  p.n=0.04;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,4.0e-2,"rmf_b nrich");
  t.test_rel(n.mu,t2,4.0e-2,"rmf_b nrich");
  t.test_rel(p.mu,t3,4.0e-2,"rmf_b nrich");

  /// Check neutron-rich matter with non-trivial alpha

  cout << "Neutron-rich matter, alpha!=3." << endl;
  n.n=0.12;
  p.n=0.04;
  s4e.alpha=2.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  cout << endl;

  /// ----------------------------------------------------------------
  /// RMF (C) EOS

  s4e.set_base_eos(rmf4_c);
  s4e.alpha=3.0;

  /// Check nuclear matter
  
  cout << "RMF_C:" << endl;
  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  rmf_c.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  rmf4_c.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_c nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_c nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf_c nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_c nuc");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_c nuc");
  t.test_rel(p.mu,t3,1.0e-6,"rmf_c nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  rmf_c.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.16;
  p.n=0.0;
  rmf4_c.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_c neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_c neut");
  n.n=0.16;
  p.n=0.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_c neut");
  t.test_rel(n.mu,t2,1.0e-6,"rmf_c neut");

  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  rmf_c.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.12;
  p.n=0.04;
  rmf4_c.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-6,"rmf_c nrich");
  //t.test_rel(n.mu,t2,1.0e-6,"rmf_c nrich");
  //t.test_rel(p.mu,t3,1.0e-6,"rmf_c nrich");
  n.n=0.12;
  p.n=0.04;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,4.0e-1,"rmf_c nrich");
  //t.test_rel(n.mu,t2,4.0e-1,"rmf_c nrich");
  //t.test_rel(p.mu,t3,4.0e-1,"rmf_c nrich");

  /// Check neutron-rich matter with non-trivial alpha

  cout << "Neutron-rich matter, alpha!=3." << endl;
  n.n=0.12;
  p.n=0.04;
  s4e.alpha=2.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  cout << endl;
  
  /// ----------------------------------------------------------------
  /// MDI EOS

  s4e.set_base_eos(mdi4);
  s4e.alpha=3.0;

  /// Check nuclear matter
  
  cout << "MDI:" << endl;
  cout << "Nuclear matter." << endl;
  n.n=0.08;
  p.n=0.08;
  mdi.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  mdi4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"mdi nuc");
  //t.test_rel(n.mu,t2,1.0e-5,"mdi nuc");
  //t.test_rel(p.mu,t3,1.0e-5,"mdi nuc");
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"mdi nuc");
  //t.test_rel(n.mu,t2,1.0e-5,"mdi nuc");
  //t.test_rel(p.mu,t3,1.0e-5,"mdi nuc");

  /// Check neutron matter

  cout << "Neutron matter." << endl;
  n.n=0.16;
  p.n=0.0;
  mdi.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.16;
  p.n=0.0;
  mdi4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"mdi neut");
  //t.test_rel(n.mu,t2,1.0e-5,"mdi neut");
  n.n=0.16;
  p.n=0.0;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"mdi neut");
  //t.test_rel(n.mu,t2,1.0e-5,"mdi neut");

  /// Check neutron-rich matter

  cout << "Neutron-rich matter." << endl;
  n.n=0.12;
  p.n=0.04;
  mdi.calc_e(n,p,th);
  t1=(th.ed/0.16-n.m)*hc_mev_fm;
  t2=n.mu;
  t3=p.mu;
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  n.n=0.12;
  p.n=0.04;
  mdi4.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,1.0e-5,"mdi nrich");
  t.test_rel(n.mu,t2,1.0e-5,"mdi nrich");
  t.test_rel(p.mu,t3,1.0e-5,"mdi nrich");
  n.n=0.12;
  p.n=0.04;
  s4e.calc_e(n,p,th);
  cout << n.n << " " << p.n << " " 
       << (th.ed/0.16-n.m)*hc_mev_fm << " "
       << n.mu << " " << p.mu << endl;
  //t.test_rel((th.ed/0.16-n.m)*hc_mev_fm,t1,4.0e-2,"mdi nrich");
  t.test_rel(n.mu,t2,4.0e-2,"mdi nrich");
  t.test_rel(p.mu,t3,4.0e-2,"mdi nrich");

  /// Check neutron-rich matter with non-trivial alpha

  /*
    cout << "Neutron-rich matter, alpha!=3." << endl;
    n.n=0.12;
    p.n=0.04;
    s4e.alpha=2.0;
    s4e.calc_e(n,p,th);
    cout << n.n << " " << p.n << " " 
    << (th.ed/0.16-n.m)*hc_mev_fm << " "
    << n.mu << " " << p.mu << endl;
    cout << endl;
  */

  /// ----------------------------------------------------------------
  /// Check that we get the same threshold for URCA in both 
  /// APR and the new version with a sym4 object
  
  /*
  cnx.set_n_and_p(n,p);
  cnx.set_eos(apr);

  double uden, uden2, nnu, npu;
  bool turning_on;
  nnu=0.6;
  npu=0.4;
  cnx.calc_eos();
  uden=cnx.allow_urca;
  
  s4e.alpha=3.0;
  s4e.set_base_eos(apr4);
  cnx.set_eos(s4e);
  nnu=0.6;
  npu=0.4;
  cnx.calc_eos();
  uden2=cnx.allow_urca;
  t.test_rel(uden,uden2,1.0e-6,"urca density comparison");
  cout << "APR urca density (no muons): " << uden << endl;
  */

  /// ----------------------------------------------------------------

  t.report();
  
  return 0;

}

