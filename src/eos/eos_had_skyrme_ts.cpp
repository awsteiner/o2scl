/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#include <o2scl/eos_had_skyrme.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int load_sly4(eos_had_skyrme &sk) {
  sk.t0=-2488.913/hc_mev_fm;
  sk.t1=486.818/hc_mev_fm;
  sk.t2=-546.395/hc_mev_fm;
  sk.t3=13777.0/hc_mev_fm;
  sk.x0=0.8340;
  sk.x1=-0.3438;
  sk.x2=-1.0; 
  sk.x3=1.3540;
  sk.a=0.0;
  sk.b=1.0;
  sk.alpha=0.1666666666667;
  sk.W0=123/hc_mev_fm;
  return 0;
}

int load_skms(eos_had_skyrme &sk) {
  sk.t0=-2645.0/hc_mev_fm;
  sk.t1=410.0/hc_mev_fm;
  sk.t2=-135.0/hc_mev_fm;
  sk.t3=15595.0/hc_mev_fm;
  sk.x0=0.090;
  sk.x1=0.0;
  sk.x2=0.0;
  sk.x3=0.0; 
  sk.a=0.0;
  sk.b=1.0;
  sk.alpha=0.1666666666667;
  sk.W0=130.0/hc_mev_fm;
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  double n0, eoa2;
  eos_had_skyrme sk;

  fermion n(939.0/197.33,2.0), p(939.0/197.33,2.0);
  n.non_interacting=false;
  p.non_interacting=false;
  thermo th;

  load_sly4(sk);

  // ------------------------------------------------------------
  // Use check_mu() to check chemical potentials
  // ------------------------------------------------------------

  n.n=0.08;
  p.n=0.08;
  double mund, mupd, munde, mupde;
  sk.check_mu(n,p,th,mund,mupd,munde,mupde);
  t.test_abs(n.mu,mund,fabs(munde),"neutron chem pot.");
  t.test_abs(p.mu,mupd,fabs(mupde),"proton chem pot.");

  sk.check_mu_T(n,p,1.0/hc_mev_fm,th,mund,mupd,munde,mupde);
  t.test_abs(n.mu,mund,fabs(munde),"neutron chem pot. T=1");
  t.test_abs(p.mu,mupd,fabs(mupde),"proton chem pot. T=1");

  sk.check_mu_T(n,p,10.0/hc_mev_fm,th,mund,mupd,munde,mupde);
  t.test_abs(n.mu,mund,fabs(munde),"neutron chem pot. T=10");
  t.test_abs(p.mu,mupd,fabs(mupde),"proton chem pot. T=10");

  double end, ende;
  sk.check_en(n,p,5.0/hc_mev_fm,th,end,ende);
  t.test_abs(th.en,end,ende,"entropy");
  
  // ------------------------------------------------------------
  // Test that calc_p() and calc_temp_p() are working
  // ------------------------------------------------------------

  n.n=0.02;
  p.n=0.02;

  sk.calc_e(n,p,th);
  cout << n.mu << " " << p.mu << endl;
  sk.calc_p(n,p,th);
  cout << n.n << " " << p.n << endl;
  t.test_rel(n.n,0.02,1.0e-6,"nn");
  t.test_rel(p.n,0.02,1.0e-6,"np");

  sk.calc_temp_e(n,p,0.0,th);
  cout << n.mu << " " << p.mu << endl;
  sk.calc_temp_p(n,p,0.0,th);
  cout << n.n << " " << p.n << endl;
  t.test_rel(n.n,0.02,1.0e-6,"nn");
  t.test_rel(p.n,0.02,1.0e-6,"np");

  sk.calc_temp_e(n,p,1.0/hc_mev_fm,th);
  cout << n.mu << " " << p.mu << endl;
  sk.calc_temp_p(n,p,1.0/hc_mev_fm,th);
  cout << n.n << " " << p.n << endl;
  t.test_rel(n.n,0.02,1.0e-6,"nn");
  t.test_rel(p.n,0.02,1.0e-6,"np");
  
  // ------------------------------------------------------------

  load_skms(sk);

  cout << endl;
  cout << "Compare saturation properties of SkMs with known values:"
       << endl;
  n0=sk.fn0(0.0,eoa2);
  cout << n0 << " " << 0.1603 << endl;
  cout << -sk.feoa(n0)*hc_mev_fm << " " << 15.78 << endl;
  cout << sk.fmsom(n0) << " " << 0.79 << endl;
  cout << sk.fcomp(n0)*hc_mev_fm << " " << 216.7 << endl;
  cout << sk.fesym(n0)*hc_mev_fm << " " << 30.03 << endl;
  cout << sk.fkprime(n0)*hc_mev_fm << endl;
  t.test_rel(-sk.feoa(n0)*hc_mev_fm,-eoa2*hc_mev_fm,1.0e-4,"eoa1");
  t.test_rel(-eoa2*hc_mev_fm,15.78,1.0e-3,"eoa2");
  t.test_rel(sk.fmsom(n0),0.79,4.0e-2,"msom");
  t.test_rel(sk.fcomp(n0)*hc_mev_fm,216.7,1.0e-3,"comp");
  t.test_rel(sk.fesym(n0)*hc_mev_fm,30.03,1.0e-4,"esym");
  cout << endl;

  cout << "Testing new fractional power of alpha:" << endl;
  t.test_rel(sk.alpha,1.0/6.0,1.0e-12,"frac. alpha");
  cout << endl;

  cout << "Compare with results from eos_had_base():" << endl;
  cout << "With a non-zero value of a" << endl;
  sk.a=0.1;

  sk.n0=sk.fn0(0.0,eoa2);
  n0=sk.n0;
  sk.eoa=eoa2;
  sk.msom=sk.fmsom(n0);
  sk.comp=sk.fcomp(n0);
  sk.esym=sk.fesym(n0);
  sk.kprime=sk.fkprime(n0);

  sk.parent_method=true;
  n0=sk.fn0(0.0,eoa2);
  cout << sk.fesym(n0)*hc_mev_fm << endl;
  /// FIXME n0 doesn't look right here?
  cout << "n0=" << n0 << endl;
  t.test_rel(n0,sk.n0,1.0e-6,"n0 pm");
  t.test_rel(sk.feoa(n0),sk.eoa,1.0e-6,"eoa1 pm");
  t.test_rel(sk.fmsom(n0),sk.msom,1.0e-6,"msom pm");
  t.test_rel(sk.fcomp(n0),sk.comp,1.0e-6,"comp pm");
  t.test_rel(sk.fesym(n0),sk.esym,2.0e-6,"esym pm");
  t.test_rel(sk.fkprime(n0),sk.kprime,1.0e-3,"kprime pm");
  sk.parent_method=false;
  cout << endl;
  //exit(-1);

  cout << "See if the symmetry energy works at all densities:" << endl;

  n0=0.08;
  sk.esym=sk.fesym(n0);
  sk.parent_method=true;
  t.test_rel(sk.fesym(n0),sk.esym,1.0e-5,"esym density 1");
  sk.parent_method=false;

  n0=0.32;
  sk.esym=sk.fesym(n0);
  sk.parent_method=true;
  t.test_rel(sk.fesym(n0),sk.esym,1.0e-5,"esym density 2");
  sk.parent_method=false;

  n0=0.50;
  sk.esym=sk.fesym(n0);
  sk.parent_method=true;
  t.test_rel(sk.fesym(n0),sk.esym,1.0e-5,"esym density 3");
  sk.parent_method=false;
  cout << endl;

  cout << "See if the compressibility works at all densities:" << endl;

  n0=0.08;
  sk.comp=sk.fcomp(n0);
  sk.parent_method=true;
  t.test_rel(sk.fcomp(n0),sk.comp,1.0e-5,"comp density 1");
  sk.parent_method=false;

  n0=0.32;
  sk.comp=sk.fcomp(n0);
  sk.parent_method=true;
  t.test_rel(sk.fcomp(n0),sk.comp,1.0e-5,"comp density 2");
  sk.parent_method=false;

  n0=0.50;
  sk.comp=sk.fcomp(n0);
  sk.parent_method=true;
  t.test_rel(sk.fcomp(n0),sk.comp,1.0e-5,"comp density 3");
  sk.parent_method=false;
  cout << endl;

  cout << "See if the effective mass works at all densities:" << endl;

  n0=0.08;
  sk.msom=sk.fmsom(n0);
  sk.parent_method=true;
  t.test_rel(sk.fmsom(n0),sk.msom,1.0e-5,"msom density 1");
  sk.parent_method=false;

  n0=0.16;
  sk.msom=sk.fmsom(n0);
  sk.parent_method=true;
  t.test_rel(sk.fmsom(n0),sk.msom,1.0e-5,"msom density 2");
  sk.parent_method=false;

  n0=0.50;
  sk.msom=sk.fmsom(n0);
  sk.parent_method=true;
  t.test_rel(sk.fmsom(n0),sk.msom,1.0e-5,"msom density 3");
  sk.parent_method=false;
  cout << endl;

  cout << "See if the binding energy works at all densities:" << endl;

  n0=0.08;
  sk.eoa=sk.feoa(n0);
  sk.parent_method=true;
  t.test_rel(sk.feoa(n0),sk.eoa,1.0e-5,"eoa density 1");
  sk.parent_method=false;

  n0=0.16;
  sk.eoa=sk.feoa(n0);
  sk.parent_method=true;
  t.test_rel(sk.feoa(n0),sk.eoa,1.0e-5,"eoa density 2");
  sk.parent_method=false;

  n0=0.50;
  sk.eoa=sk.feoa(n0);
  sk.parent_method=true;
  t.test_rel(sk.feoa(n0),sk.eoa,1.0e-5,"eoa density 3");
  sk.parent_method=false;
  cout << endl;

  cout << "See if the skewness works at all densities:" << endl;

  n0=0.08;
  sk.kprime=sk.fkprime(n0);
  sk.parent_method=true;
  t.test_rel(sk.fkprime(n0),sk.kprime,2.0e-3,"kprime density 1");
  sk.parent_method=false;

  n0=0.16;
  sk.kprime=sk.fkprime(n0);
  sk.parent_method=true;
  t.test_rel(sk.fkprime(n0),sk.kprime,5.0e-4,"kprime density 2");
  sk.parent_method=false;

  n0=0.50;
  sk.kprime=sk.fkprime(n0);
  sk.parent_method=true;
  t.test_rel(sk.fkprime(n0),sk.kprime,5.0e-4,"kprime density 3");
  sk.parent_method=false;
  cout << endl;

#ifdef O2SCL_NEVER_DEFINED

  cout << "Compare model PeHF with known saturation values: " << endl;

  sk.load("PeHF");
  
  n0=sk.fn0(0.0,eoa2);
  cout << n0 << " " << 2.0*pow(1.34,3.0)/3.0/pi2 << endl;
  cout << -sk.feoa(n0)*hc_mev_fm << " " << -eoa2*hc_mev_fm << " " << 16.00 
       << endl;
  cout << sk.fesym(n0)*hc_mev_fm << " " << 32.0 << endl;
  t.test_rel(-sk.feoa(n0)*hc_mev_fm,-eoa2*hc_mev_fm,1.0e-4,"eoa1 PeHF");
  t.test_rel(-eoa2*hc_mev_fm,16.0,1.0e-3,"eoa2 PeHF");
  t.test_rel(sk.fesym(n0)*hc_mev_fm,32.0,1.0e-3,"esym PeHF");
  cout << endl;
  
  cout << "Compare model SkSC10 with known saturation values: " << endl;
  
  sk.load("SkSC10");
  
  n0=sk.fn0(0.0,eoa2);
  cout << n0 << endl;
  cout << -sk.feoa(n0)*hc_mev_fm << " " << -eoa2*hc_mev_fm << endl;
  cout << sk.fmsom(n0) << " " << 1.0 << endl;
  cout << sk.fcomp(n0)*hc_mev_fm << " " << 235.8 << endl;
  cout << sk.fesym(n0)*hc_mev_fm << " " << 32.0 << endl;
  t.test_rel(-sk.feoa(n0)*hc_mev_fm,-eoa2*hc_mev_fm,1.0e-4,"eoa SkSC10");
  t.test_rel(sk.fmsom(n0),1.0,1.0e-4,"msom SkSC10");
  t.test_rel(sk.fcomp(n0)*hc_mev_fm,235.8,5.0e-4,"n0 SkSC10");
  t.test_rel(sk.fesym(n0)*hc_mev_fm,32.0,1.0e-4,"esym SkSC10");
  cout << endl;

#endif

  cout << "Test calpar():" << endl;
  load_skms(sk);
  cout << sk.t0 << " " << sk.t1 << " " << sk.t2 << " " << sk.t3 
       << " " << sk.alpha << endl;

  double t0old, t1old, t2old, t3old, alphaold;
  t0old=sk.t0;
  t1old=sk.t1;
  t2old=sk.t2;
  t3old=sk.t3;
  alphaold=sk.alpha;

  sk.n0=0.1603;
  sk.eoa=-15.77/hc_mev_fm;
  sk.msom=0.79;
  sk.comp=216.7/hc_mev_fm;
  sk.esym=30.03/hc_mev_fm;
  sk.calpar();
  cout << sk.t0 << " " << sk.t1 << " " << sk.t2 << " " << sk.t3 
       << " " << sk.alpha << endl;
  t.test_rel(sk.t0,t0old,1.0e-1,"calpar t0");
  t.test_rel(sk.t1,t1old,1.0e-1,"calpar t1");
  t.test_rel(sk.t2,t2old,1.0e-1,"calpar t2");
  t.test_rel(sk.t3,t3old,1.0e-1,"calpar t3");
  t.test_rel(sk.alpha,alphaold,1.0e-1,"calpar alpha");
  cout << endl;

  // -----------------------------------------------------------
  // Test finite-temperature functions
  // -----------------------------------------------------------

  cout << "Compare finite T to zero T code: " << endl;
  load_skms(sk);

  n.n=0.16;
  p.n=0.16;

  // Try zero-temperature version
  sk.calc_e(n,p,th);
    
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  double tmp10, tmp11, tmp12, tmp13;

  tmp1=th.ed*hc_mev_fm/0.32-939.0;
  tmp2=n.n;
  tmp3=p.n;
  tmp4=n.ms;
  tmp5=p.ms;
  tmp6=n.mu;
  tmp7=p.mu;
  tmp8=th.ed;
  tmp9=th.pr;
  tmp10=n.pr;
  tmp11=p.pr;
  tmp12=n.ed;
  tmp13=p.ed;
    
  // Compare with finite temperature version:
  sk.calc_temp_e(n,p,0.0001,th);

  t.test_rel(tmp1,th.ed*hc_mev_fm/0.32-939.0,1.0e-2,"test(1)");
  t.test_rel(tmp2,n.n,1.0e-2,"test(2)");
  t.test_rel(tmp3,p.n,1.0e-2,"test(3)");
  t.test_rel(tmp4,n.ms,1.0e-2,"test(4)");
  t.test_rel(tmp5,p.ms,1.0e-2,"test(5)");
  t.test_rel(tmp6,n.mu,1.0e-2,"test(6)");
  t.test_rel(tmp7,p.mu,1.0e-2,"test(7)");
  t.test_rel(tmp8,th.ed,1.0e-2,"test(8)");
  t.test_rel(tmp9,th.pr,1.0e-2,"test(9)");
  t.test_rel(tmp10,n.pr,1.0e-2,"test(10)");
  t.test_rel(tmp11,p.pr,1.0e-2,"test(11)");
  t.test_rel(tmp12,n.ed,1.0e-2,"test(12)");
  t.test_rel(tmp13,p.ed,1.0e-2,"test(13)");
  cout << endl;

  // Test chemical potentials 
  cout << "Test chemical potentials at finite T: " << endl;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);
  double mun=n.mu, mup=p.mu, ed0=th.ed;
  n.n+=1.0e-4;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);
  double edn=th.ed;
  n.n-=1.0e-4;
  p.n+=1.0e-4;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);
  double edp=th.ed;
  p.n-=1.0e-4;
  t.test_rel((edn-ed0)/1.0e-4,mun,1.0e-4,"mun");
  t.test_rel((edp-ed0)/1.0e-4,mup,1.0e-4,"mup");
  cout << endl;

  // Test EOS with inc_rest_mass=false at zero T
  n.n=0.08;
  p.n=0.08;
  sk.calc_e(n,p,th);
  double edi1=th.ed;
  double nuni1=n.nu;
  double nupi1=p.nu;
  n.inc_rest_mass=false;
  p.inc_rest_mass=false;
  sk.calc_e(n,p,th);

  t.test_rel(edi1-n.n*n.m-p.n*p.m,th.ed,1.0e-6,"inc_rest_mass=false E T=0");
  t.test_rel(nuni1-n.m,n.nu,1.0e-6,"inc_rest_mass=false mun T=0");
  t.test_rel(nupi1-p.m,p.nu,1.0e-6,"inc_rest_mass=false mup T=0");
  
  // Test EOS with inc_rest_mass=false at finite T
  n.inc_rest_mass=true;
  p.inc_rest_mass=true;
  n.n=0.08;
  p.n=0.08;
  // Give a good guess for the effective chemical potentials
  n.nu+=n.m;
  p.nu+=p.m;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);
  edi1=th.ed;
  nuni1=n.nu;
  nupi1=p.nu;
  n.inc_rest_mass=false;
  p.inc_rest_mass=false;
  // Give a good guess for the effective chemical potentials
  n.nu-=n.m;
  p.nu-=p.m;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);

  t.test_rel(edi1-n.n*n.m-p.n*p.m,th.ed,1.0e-6,"inc_rest_mass=false E");
  t.test_rel(nuni1-n.m,n.nu,1.0e-6,"inc_rest_mass=false mun");
  t.test_rel(nupi1-p.m,p.nu,1.0e-6,"inc_rest_mass=false mup");

  cout << sk.fesym(0.16)*hc_mev_fm << endl;
  cout << sk.fesym_T(0.16,0.3/hc_mev_fm)*hc_mev_fm << endl;
  cout << sk.fesym_T(0.16,1.0/hc_mev_fm)*hc_mev_fm << endl;
  cout << sk.fesym_T(0.16,3.0/hc_mev_fm)*hc_mev_fm << endl;
  cout << sk.fesym_T(0.16,10.0/hc_mev_fm)*hc_mev_fm << endl;

  // -----------------------------------------------------------
  // Test alt_params_saturation
  // -----------------------------------------------------------

  // Test alt_params_saturation() with the UNEDF2 couplings from
  // Kortelainen et al. (2014)
  sk.alt_params_saturation(0.15631,(-15.8)/hc_mev_fm,239.93/hc_mev_fm,
			   1.0/1.074,29.131/hc_mev_fm,40.0/hc_mev_fm,
			   1.0/1.249,-46.831/hc_mev_fm,-113.164/hc_mev_fm,
			   -64.309/hc_mev_fm,-38.650/hc_mev_fm);
  
  sk.saturation();
  cout << sk.n0 << " " << sk.eoa*hc_mev_fm << " " << sk.comp*hc_mev_fm
       << " " << sk.esym*hc_mev_fm << " " << 1.0/sk.msom << " "
       << sk.fesym_slope(sk.n0)*hc_mev_fm << endl;
  t.test_rel(sk.n0,0.15631,1.0e-4,"n0");
  t.test_rel(sk.eoa*hc_mev_fm,-15.8,1.0e-4,"eoa");
  t.test_rel(sk.comp*hc_mev_fm,239.93,1.0e-4,"comp");
  t.test_rel(sk.msom,1.0/1.074,1.0e-4,"msom");
  t.test_rel(sk.esym*hc_mev_fm,29.131,1.0e-4,"esym");
  t.test_rel(sk.fesym_slope(sk.n0)*hc_mev_fm,40.0,1.0e-4,"L");
  t.test_rel(sk.f_effm_vector(sk.n0),1.0/1.249,1.0e-4,"Mv*");

  // -----------------------------------------------------------
  // Test calc_deriv_temp_e
  // -----------------------------------------------------------

  fermion_deriv ne(939.0/hc_mev_fm,2.0), pr(938.0/hc_mev_fm,2.0);
  ne.non_interacting=false;
  pr.non_interacting=false;
  
  n.n=0.25;
  p.n=0.35;
  sk.calc_temp_e(n,p,4.0/hc_mev_fm,th);
  cout << n.mu << " " << p.mu << endl;
  ne.mu=n.mu;
  pr.mu=p.mu;
  ne.nu=n.nu;
  pr.nu=p.nu;
  
  ne.n=0.25;
  pr.n=0.35;
  thermo_np_f_deriv tnfd;

  cout << "Here." << endl;
  sk.calc_deriv_temp_e(ne,pr,4.0/hc_mev_fm,th,tnfd);
  cout << "Here2." << endl;
  double dmundnn1=ne.mu;
  double dmupdnp1=pr.mu;
  double dmundnp1=ne.mu;
  double dmupdnn1=pr.mu;
  double dmundT1=ne.mu;
  double dmupdT1=pr.mu;
  double dsdT1=th.en;

  ne.n+=1.0e-4;
  sk.calc_deriv_temp_e(ne,pr,4.0/hc_mev_fm,th,tnfd);
  double dmundnn2=ne.mu;
  double dmupdnn2=pr.mu;
  ne.n-=1.0e-4;

  pr.n+=1.0e-4;
  sk.calc_deriv_temp_e(ne,pr,4.0/hc_mev_fm,th,tnfd);
  double dmundnp2=ne.mu;
  double dmupdnp2=pr.mu;
  pr.n-=1.0e-4;

  sk.calc_deriv_temp_e(ne,pr,(4.0+1.0e-3)/hc_mev_fm,th,tnfd);
  double dmundT2=ne.mu;
  double dmupdT2=pr.mu;
  double dsdT2=th.en;

  sk.calc_deriv_temp_e(ne,pr,4.0/hc_mev_fm,th,tnfd);

  // Test the six derivatives
  cout << tnfd.dsdT << " " << (dsdT2-dsdT1)/(1.0e-3/hc_mev_fm) << endl;
  cout << tnfd.dmundT << " " << (dmundT2-dmundT1)/(1.0e-3/hc_mev_fm) << endl;
  cout << tnfd.dmupdT << " " << (dmupdT2-dmupdT1)/(1.0e-3/hc_mev_fm) << endl;
  cout << tnfd.dmundnn << " " << (dmundnn2-dmundnn1)/1.0e-4 << endl;
  cout << tnfd.dmupdnp << " " << (dmupdnp2-dmupdnp1)/1.0e-4 << endl;
  cout << tnfd.dmudn_mixed << " "
       << (dmundnp2-dmundnp1)/1.0e-4 << " " 
       << (dmupdnn2-dmupdnn1)/1.0e-4 << endl;

  t.report();

  return 0;
}

