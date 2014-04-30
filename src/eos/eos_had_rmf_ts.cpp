/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/rmf_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

// Create a namespace to avoid warnings about shadowing
// global variables

namespace rmf_eos_ts_ns {

  int load_nl3(rmf_eos &rmf) {

    rmf.ms=508.194;
    rmf.mw=782.501;
    rmf.mr=763.0;
    rmf.mnuc=939.0;
    rmf.ms/=hc_mev_fm; 
    rmf.mw/=hc_mev_fm; 
    rmf.mr/=hc_mev_fm; 
    rmf.mnuc/=hc_mev_fm;
    
    double gs, gw, gr;
    gs=10.217;
    gw=12.868;
    gr=4.474;
    rmf.b=-10.431;
    rmf.c=-28.885;
    rmf.b/=-rmf.mnuc*pow(fabs(gs),3.0);
    rmf.c/=pow(gs,4.0);
    gr*=2.0;
    rmf.cs=gs/rmf.ms;
    rmf.cw=gw/rmf.mw;
    rmf.cr=gr/rmf.mr;
    
    rmf.xi=0.0; 
    rmf.zeta=0.0;
    rmf.a1=0.0;
    rmf.a2=0.0;
    rmf.a3=0.0;
    rmf.a4=0.0;
    rmf.a5=0.0;
    rmf.a6=0.0;
    rmf.b1=0.0;
    rmf.b2=0.0;
    rmf.b3=0.0;
    
    return 0;
  }

  fermion nferm(939.0/hc_mev_fm,2.0), p(939.0/hc_mev_fm,2.0);

  fermion e(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  fermion mu(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_muon),2.0);

  fermion nue(0.0,1.0), numu(0.0,1.0);
  rmf_eos rmf;
  fermion_eff eff;
  thermo thx;
  double nb;

  int nstar_mat(size_t nv, const ubvector &x, ubvector &y) {
    nferm.mu=x[0];
    p.mu=x[1];
    e.mu=x[2];
    mu.mu=x[3];
    nue.mu=x[4];
    double sig=x[5];
    double ome=x[6];
    double rho=x[7];
    double f1, f2, f3;
    numu.mu=mu.mu-e.mu+nue.mu;
  
    rmf.calc_eq_p(nferm,p,sig,ome,rho,f1,f2,f3,thx);
    eff.calc_mu_zerot(e);
    eff.calc_mu_zerot(mu);
    eff.calc_mu_zerot(nue);
    eff.calc_mu_zerot(numu);
    
    y[0]=p.n-e.n-mu.n;
    y[1]=nferm.n+p.n-nb;
    y[2]=nferm.mu+nue.mu-p.mu-e.mu;
    y[3]=e.n/nb-0.35;
    y[4]=mu.n-numu.n;
    y[5]=f1;
    y[6]=f2;
    y[7]=f3;

    return 0;
  }
  
  mm_funct_fptr<> mff(nstar_mat);

}

using namespace rmf_eos_ts_ns;

int main(void) {
  cout.setf(ios::scientific);
  rmf_eos re;
  mroot_hybrids<mm_funct<> > mrp;
  test_mgr t;
  t.set_output_level(1);

  nferm.non_interacting=false;
  p.non_interacting=false;
  double sig, ome, rho, f1, f2, f3, barn;
  thermo th;

  // -----------------------------------------------------------------
  // Check calc_e() for both NL3 and RAPR
  // -----------------------------------------------------------------

  load_nl3(re);

  // Nuclear matter
  for(double nbx=0.04;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx/2.0;
    p.n=nbx/2.0;
    int ret=re.calc_e(nferm,p,th);
    t.test_gen(ret==0,"NL3 nuclear matter ret. val.");
    t.test_rel(nferm.n,nbx/2.0,1.0e-6,"NL3 nuclear matter n.n");
    t.test_rel(p.n,nbx/2.0,1.0e-6,"NL3 nuclear matter p.n");
  }
  
  // This calls the error handler for pure neutron matter, which
  // should be fixed.
  // Neutron matter
  for(double nbx=0.08;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx;
    p.n=0.0;
    int ret=re.calc_e(nferm,p,th);
    cout << nbx << " " << ret << endl;
    t.test_gen(ret==0,"NL3 neutron matter ret. val.");
    t.test_rel(nferm.n,nbx,1.0e-6,"NL3 neutron matter n.n");
    t.test_rel(p.n,0.0,1.0e-6,"NL3 neutron matter p.n");
  }

  // Neutron-rich matter
  for(double nbx=0.04;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx*0.75;
    p.n=nbx/4.0;
    int ret=re.calc_e(nferm,p,th);
    t.test_gen(ret==0,"NL3 neutron-rich matter ret. val.");
    t.test_rel(nferm.n,nbx*0.75,1.0e-6,"NL3 neutron-rich matter n.n");
    t.test_rel(p.n,nbx/4.0,1.0e-6,"NL3 neutron-rich matter p.n");
  }

#ifdef O2SCL_NEVER_DEFINED

  re.load("RAPR");

  // Nuclear matter
  for(double nbx=0.04;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx/2.0;
    p.n=nbx/2.0;
    int ret=re.calc_e(nferm,p,th);
    t.test_gen(ret==0,"RAPR nuclear matter ret. val.");
    t.test_rel(nferm.n,nbx/2.0,1.0e-6,"RAPR nuclear matter n.n");
    t.test_rel(p.n,nbx/2.0,1.0e-6,"RAPR nuclear matter p.n");
  }
  
  // This calls the error handler for pure neutron matter, which
  // should be fixed.
  // Neutron matter
  for(double nbx=0.04;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx;
    p.n=0.0;
    int ret=re.calc_e(nferm,p,th);
    t.test_gen(ret==0,"RAPR neutron matter ret. val.");
    t.test_rel(nferm.n,nbx,1.0e-6,"RAPR neutron matter n.n");
    t.test_rel(p.n,0.0,1.0e-6,"RAPR neutron matter p.n");
  }

  // Neutron-rich matter
  for(double nbx=0.04;nbx<=1.29;nbx*=2.0) {
    nferm.n=nbx*0.75;
    p.n=nbx/4.0;
    int ret=re.calc_e(nferm,p,th);
    t.test_gen(ret==0,"RAPR neutron-rich matter ret. val.");
    t.test_rel(nferm.n,nbx*0.75,1.0e-6,"RAPR neutron-rich matter n.n");
    t.test_rel(p.n,nbx/4.0,1.0e-6,"RAPR neutron-rich matter p.n");
  }

#endif

  // -----------------------------------------------------------------
  // Check FSU Gold and compare with NL3
  // -----------------------------------------------------------------

#ifdef O2SCL_NEVER_DEFINED

  if (false) {
    cout << "Here." << endl;
    re.load("NL3");
    cout << re.ms*hc_mev_fm << endl;
    cout << re.cs*re.ms*re.cs*re.ms << " "
	 << re.cw*re.mw*re.cw*re.mw << " "
	 << re.cr*re.mr*re.cr*re.mr << endl;
    cout << re.b*2.0*re.mnuc*hc_mev_fm << endl;
    cout << re.c*6.0 << endl;
    cout << re.zeta << endl;
    cout << re.b1/pow(re.mw*re.cw,2.0) << endl;

    re.load("FSUGold");
    cout << re.ms*hc_mev_fm << endl;
    cout << re.cs*re.ms*re.cs*re.ms << " "
	 << re.cw*re.mw*re.cw*re.mw << " "
	 << re.cr*re.mr*re.cr*re.mr << endl;
    cout << re.b*2.0*re.mnuc*hc_mev_fm << endl;
    cout << re.c*6.0 << endl;
    cout << re.zeta << endl;
    cout << re.b1/pow(re.mw*re.cw,2.0) << endl;
    nferm.m=939.0/hc_mev_fm;
    p.m=939.0/hc_mev_fm;
    nferm.mu=4.9;
    p.mu=4.9;
    re.set_fields(0.1,0.07,0.0);
    re.saturation();
    cout << re.n0 << " " << hc_mev_fm*re.eoa << " " 
	 << re.comp*hc_mev_fm << " "
	 << re.esym*hc_mev_fm << " " << re.msom << endl;
    cout << re.fesym_slope(re.n0)*hc_mev_fm << endl;
    double fgeoa;
    re.set_fields(0.1,0.07,0.0);
    double fgn0=re.fn0(0.0,fgeoa);
    cout << fgn0 << " " << hc_mev_fm*fgeoa << endl;

    exit(-1);
  }

#endif

  // -----------------------------------------------------------------

  load_nl3(re);
  re.set_n_and_p(nferm,p);
  re.set_thermo(th);
  nferm.mu=4.8;
  p.mu=4.8;
  re.saturation();
  cout << "NL3: " << endl;
  cout << "  Saturation density: " << re.n0 << endl;
  t.test_rel(re.n0,0.148,1.0e-2,"sat density");
  cout << "  Effective mass: " << re.msom << " " << nferm.ms/nferm.m << endl;
  t.test_rel(re.msom,nferm.ms/nferm.m,1.0e-6,"msom");
  cout << "  Zero pressure: " << th.pr << endl;
  t.test_rel(th.pr,0.0,1.0e-8,"zero press");
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (th.ed/re.n0-re.mnuc)*hc_mev_fm << endl;
  t.test_rel(re.eoa,th.ed/re.n0-re.mnuc,1.0e-6,"eoa");
  cout << "  Thermodynamic identity: " 
       << th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu << endl;
  t.test_rel(th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu,0.0,1.0e-9,"TI");
  cout << endl;
  
  re.set_n_and_p(nferm,p);
  re.set_thermo(th);

  // Test calc_e(). 
  re.saturation();
  nferm.n=re.n0/2.0;
  p.n=re.n0/2.0;
  nferm.mu=nferm.m*1.1;
  p.mu=p.m*1.1;
  re.calc_e(nferm,p,th);
  t.test_rel((th.ed/(nferm.n+p.n)-nferm.m)*hc_mev_fm,re.eoa*hc_mev_fm,
	     1.0e-5,"calc_e");
  
  cout << "\n1. Testing fix_saturation()\n" << endl;
  cout << "  From PRL 86, 5647 - NL3" << endl;

  re.zeta=0.0;
  re.xi=0.0;

  cout << "  m_n=939.0 MeV, m_{sigma}=508.194 MeV, m_{omega}=782.5 MeV, "
       << "m_{rho}=763.0 MeV" << endl;
  re.mnuc=939.0/hc_mev_fm;
  re.ms=508.194/hc_mev_fm;
  re.mw=782.5/hc_mev_fm;
  re.mr=763.0/hc_mev_fm;

  cout << "  M^{*}=0.59, K=271 MeV, E_{b}=-16.25 MeV, n_0=0.1484 fm^{-3}" 
       << endl;
  re.esym=35.0/hc_mev_fm;
  re.msom=0.59;
  re.comp=271.0/hc_mev_fm;
  re.eoa=-16.25/hc_mev_fm;
  re.n0=0.1484;

  cout << "  => g_{sigma}^2=104.387, g_{omega}^2=165.585, kappa=3.860, "
       << "lambda=-0.01591" << endl;
  
  re.fix_saturation();

  t.test_rel(re.cs*re.cs*re.ms*re.ms,104.387,2.0e-2,"gs2");
  t.test_rel(re.cw*re.cw*re.mw*re.mw,165.585,3.0e-2,"gw2");
  t.test_rel(re.b*hc_mev_fm*2.0*re.mnuc,3.860,5.0e-2,"b");
  t.test_rel(re.c*6.0,-0.01591,1.0e-2,"c");
  cout << "  " << re.cs*re.cs*re.ms*re.ms << " " 
       << re.cw*re.cw*re.mw*re.mw << " " 
       << re.b*hc_mev_fm*2.0*re.mnuc << " " << re.c*6.0 << endl;
  cout << "  " << 104.387 << " " << 165.585 << " " << 3.860 << " " 
       << -0.01591 << endl;
  cout << endl;

  cout << endl << "  From PRL 86, 5647 - Z271" << endl;

  re.zeta=0.06;

  cout << "  m_n=939.0 MeV, m_{sigma}=465.0 MeV, m_{omega}=783.0 MeV, "
       << "m_{rho}=763.0 MeV" << endl;
  re.mnuc=939.0/hc_mev_fm;
  re.ms=465.0/hc_mev_fm;
  re.mw=783.0/hc_mev_fm;
  re.mr=763.0/hc_mev_fm;

  cout << "  M^{*}=0.8, K=271 MeV, E_{b}=-16.25 MeV, n_0=0.1484 fm^{-3}" 
       << endl;
  re.esym=35.0/hc_mev_fm;
  re.msom=0.8;
  re.comp=271.0/hc_mev_fm;
  re.eoa=-16.25/hc_mev_fm;
  re.n0=0.1484;

  cout << "  => g_{sigma}^2=49.44, g_{omega}^2=70.669, kappa=6.17, "
       << "lambda=-0.15634" << endl;

  re.fix_saturation();
  t.test_rel(re.cs*re.cs*re.ms*re.ms,49.44,1.0e-2,"gs2");
  t.test_rel(re.cw*re.cw*re.mw*re.mw,70.669,1.0e-2,"gw2");
  t.test_rel(re.b*hc_mev_fm*2.0*re.mnuc,6.17,2.0e-2,"b");
  t.test_rel(re.c*6.0,0.15634,1.0e-2,"c");
  cout << "  " << re.cs*re.cs*re.ms*re.ms << " " 
       << re.cw*re.cw*re.mw*re.mw << " " 
       << re.b*hc_mev_fm*2.0*re.mnuc << " " << re.c*6.0 << " " 
       << re.cr*re.cr*re.mr*re.mr << endl;
  cout << "  " << 49.44 << " " << 70.669 << " " << 6.17 << " " 
       << 0.15634 << endl;

  cout << endl << "  From PRD 46, 1274" << endl;

  re.zeta=0.0;
  re.xi=0.0;

  cout << "  m_n=? MeV, m_{sigma}=? MeV, m_{omega}=? MeV, m_{rho}=? MeV" 
       << endl;
  re.mnuc=939.0/hc_mev_fm;
  re.ms=3.3446;
  re.mw=3.9679;
  re.mr=3.9021;

  cout << "  E_{sym}=32.5 MeV, M^{*}=0.855, K=225 MeV, E_{b}=-16.0 MeV, "
       << "n_0=0.16 fm^{-3}" << endl;
  re.esym=32.5/hc_mev_fm;
  re.msom=0.855;
  re.comp=225.0/hc_mev_fm;
  re.eoa=-16.0/hc_mev_fm;
  re.n0=0.16;

  cout << "  => C_{sigma}^2=7.487, C_{omega}^2=2.615, C_{rho}^2=4.774, " 
       << "b=0.0, c=0.0 " << endl;
  re.fix_saturation(sqrt(7.5),sqrt(2.6),1.0e-4,1.0e-4);
  cout << "  " << re.cs*re.cs << " " << re.cw*re.cw << " " << re.cr*re.cr 
       << " " << re.b << " " << re.c << endl;
  cout << endl;

  t.test_rel(re.comp*hc_mev_fm,225.0,1.0e-5,"fcomp");
  //  t.test_rel(re.fcomp(re.n0)*hc_mev_fm,225.0,1.0e-5,"fcomp");

  cout << "2. calc_p() - Nuclear matter with zeta=0.0, xi=0.0, "
       << "lamv=0.0, lam4=0.0" << endl << endl;
  nferm.n=re.n0/2.0;
  p.n=re.n0/2.0;
  fermion_eff eff;
  eff.kf_from_density(nferm);
  eff.kf_from_density(p);
  sig=re.mnuc*(1.0-re.msom)/re.ms/re.cs;
  ome=re.n0*re.cw/re.mw;
  rho=0.0;
  nferm.mu=sqrt(nferm.kf*nferm.kf+re.mnuc*re.mnuc*re.msom*re.msom)+
    re.mw*re.cw*ome;
  p.mu=sqrt(p.kf*p.kf+re.mnuc*re.mnuc*re.msom*re.msom)+re.mw*re.cw*ome;
  re.calc_eq_p(nferm,p,sig,ome,rho,f1,f2,f3,th);
  barn=nferm.n+p.n;
  
  t.test_rel(barn,re.n0,1.0e-6,"sat. den.");
  t.test_rel(re.msom,nferm.ms/nferm.m,1.0e-6,"msom");
  t.test_abs(th.pr,0.0,1.0e-4,"pressure 0");
  t.test_rel(re.eoa*hc_mev_fm,(th.ed/barn-re.mnuc)*hc_mev_fm,
	     5.0e-4,"eoa");
  t.test_abs(f1,0.0,1.0e-3,"sigma field");
  t.test_abs(f2,0.0,1.0e-10,"omega field");
  t.test_abs(f3,0.0,1.0e-10,"rho field");
  t.test_abs(th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu,0.0,1.0e-10,
	     "thermo ident.");
  cout << "  Saturation density: " << barn << " " << re.n0 << endl;
  cout << "  Effective mass: " << re.msom << " " << nferm.ms/nferm.m << endl;
  cout << "  Zero pressure: " << th.pr << endl;
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (th.ed/barn-re.mnuc)*hc_mev_fm << endl;
  cout << "  Field equations: " << f1 << " " << f2 << " " << f3 << endl;
  cout << "  Thermodynamic identity: " 
       << th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu << endl;
  cout << endl;

  re.saturation();

  t.test_rel(re.n0,0.16,1.0e-4,"sat n0");
  t.test_rel(re.msom,0.855,1.0e-4,"sat msom");
  t.test_rel(re.eoa*hc_mev_fm,-16.0,5.0e-4,"sat eoa");
  t.test_rel(re.comp*hc_mev_fm,225.0,1.0e-3,"sat comp");
  cout << "  Results of sat()" << endl;
  cout << "  Saturation density: " << re.n0 << endl;
  cout << "  Effective mass: " << re.msom << endl;
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << endl;
  cout << "  Compressiblity: " << re.comp*hc_mev_fm << endl;
  cout << "  Kprime: " << re.kprime*hc_mev_fm << endl;
  cout << endl;
  
  if (false) {
    load_nl3(re);
    re.saturation();
    cout << re.comp*hc_mev_fm << endl;
    cout << re.fcomp(re.n0)*hc_mev_fm << endl;
    
    cout << "--------------These don't agree yet-----------------" << endl;
    cout << re.kprime*hc_mev_fm << endl;
    cout << re.fkprime(re.n0)*hc_mev_fm << endl;
    cout << "--------------These don't agree yet-----------------" << endl;
    exit(-1);
  }
  
  cout << "3. Testing finite values of zeta, xi, a's, b's:" << endl << endl;

  cout << "  E_{sym}=32.5 MeV, M^{*}=0.855, K=225 MeV, E_{b}=-16.0 MeV, "
       << "n_0=0.16 fm^{-3}" << endl;
  re.esym=32.5/hc_mev_fm;
  re.msom=0.855;
  re.comp=225.0/hc_mev_fm;
  re.eoa=-16.0/hc_mev_fm;
  re.n0=0.16;
  
  cout << "  a1=0.1, a2=1.0, a6=10.0, b2=0.5, zeta=0.02, xi=1.0" << endl;
  re.a1=0.1;
  re.a2=1.0;
  re.a6=10.0;
  re.b2=0.5;
  re.zeta=0.02;
  re.xi=1.0;
  re.fix_saturation();
  
  t.test_rel(barn,re.n0,1.0e-6,"sat. den.");
  t.test_rel(re.msom,nferm.ms/nferm.m,5.0e-6,"msom");
  t.test_abs(th.pr,0.0,1.0e-5,"pressure 0");
  t.test_rel(re.eoa*hc_mev_fm,(th.ed/barn-re.mnuc)*hc_mev_fm,
	     5.0e-4,"eoa");
  cout << "  Saturation density: " << barn << " " << re.n0 << endl;
  cout << "  Effective mass: " << re.msom << " " << nferm.ms/nferm.m << endl;
  cout << "  Zero pressure: " << th.pr << endl;
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (th.ed/barn-re.mnuc)*hc_mev_fm << endl;
  cout << endl;

  nferm.n=re.n0/2.0;
  p.n=re.n0/2.0;
  eff.kf_from_density(nferm);
  eff.kf_from_density(p);
  sig=re.mnuc*(1.0-re.msom)/re.ms/re.cs;

  // Manually msolve the cubic for omega:
  double cq, lr, lr13, sqt, fir, sec;
  cq=8.0/(9.0*pow(re.cw*re.cw,3.0)*re.n0*re.n0*re.zeta);
  lr=3*re.n0/re.zeta;
  lr13=pow(lr,1.0/3.0);
  sqt=sqrt(1.0+cq);
  fir=pow((1.0+sqt),1.0/3.0);
  sec=pow(fabs(1.0-sqt),1.0/3.0);
  ome=lr13*(fir-sec);
  ome/=re.cw*re.mw;

  rho=0.0;

  nferm.mu=sqrt(nferm.kf*nferm.kf+re.mnuc*re.mnuc*re.msom*re.msom)+
    re.mw*re.cw*ome;
  p.mu=sqrt(p.kf*p.kf+re.mnuc*re.mnuc*re.msom*re.msom)+re.mw*re.cw*ome;

  barn=nferm.n+p.n;
  re.calc_eq_p(nferm,p,sig,ome,rho,f1,f2,f3,th);
  t.test_abs(th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu,0.0,1.0e-10,"ti");
  t.test_abs(f1,0.0,1.0e-5,"sigma field");
  t.test_abs(f2,0.0,1.0e-10,"omega field");
  t.test_abs(f3,0.0,1.0e-10,"rho field");
  cout << "  Thermodynamic identity: " 
       << th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu << endl;
  cout << "  Field equations: " << f1 << " " << f2 << " " << f3 << endl;
  cout << endl;

  re.saturation();
  t.test_rel(re.n0,0.16,1.0e-4,"sat2 n0");
  t.test_rel(re.msom,0.855,1.0e-4,"sat2 msom");
  t.test_rel(re.eoa*hc_mev_fm,-16.0,1.0e-4,"sat2 eoa");
  t.test_rel(re.comp*hc_mev_fm,225.0,1.0e-4,"sat2 comp");
  cout << "  Results of sat()" << endl;
  cout << "  Saturation density: " << re.n0 << endl;
  cout << "  Effective mass: " << re.msom << endl;
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << endl;
  cout << "  Compressiblity: " << re.comp*hc_mev_fm << endl;
  cout << "  Kprime: " << re.kprime*hc_mev_fm << endl;
  cout << endl;
  
#ifdef O2SCL_NEVER_DEFINED

  re.load("NL4");
  nferm.m=939.0/hc_mev_fm;
  p.m=939.0/hc_mev_fm;
  nferm.mu=4.8;
  p.mu=4.8;
  re.saturation();
  cout << re.n0 << " " << hc_mev_fm*re.eoa << " " << re.comp*hc_mev_fm << " "
       << re.esym*hc_mev_fm << " " << re.msom << endl;

#endif

    load_nl3(rmf);
  double nsig=0.1;
  double nome=0.07;
  double nrho=0.0;
  rmf.def_neutron.mu=4.9;
  rmf.def_proton.mu=4.9;
  rmf.set_fields(nsig,nome,nrho);
  double eoa;
  double n0=rmf.fn0(0.0,eoa);
  cout << n0 << " " << eoa*hc_mev_fm << endl;
  
  rmf.zeta=0.06;
  rmf.xi=1.5;
  rmf_eos cn;
  rmf.check_naturalness(cn);
  cout << cn.cs << " " << cn.cw << " " << cn.cr << " " << cn.b << "\n"
       << cn.c << " " << cn.zeta << " " << cn.xi << endl;

  t.report();

  return 0;
}


