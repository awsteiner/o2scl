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
#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/eos_had_apr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  {
    eos_had_apr ap;
    fermion n(939.0/hc_mev_fm,2.0), pr(939.0/hc_mev_fm,2.0);
    thermo hb;
    double barn, dtemp;

    ap.select(eos_had_apr::a18_uix_deltav);

    n.non_interacting=false;
    pr.non_interacting=false;

    {
      // Test the liquid-gas phase transition
      fermion n2(939.0/197.33,2.0), p2(939.0/197.33,2.0);
      n2.non_interacting=false;
      p2.non_interacting=false;
      thermo th2, th;
      n.n=0.09;
      pr.n=n.n;
      n2.n=1.0e-2;
      p2.n=1.0e-2;
      ap.def_mroot.ntrial=1000;
      ap.calc_liqgas_dens_temp_e(n,pr,n2,p2,2.0/197.33,th,th2);
      t.test_rel(n.mu,n2.mu,1.0e-8,"mun");
      t.test_rel(pr.mu,p2.mu,1.0e-8,"mup");
      cout << n.n << " " << pr.n << " " << n2.n << " " << p2.n << endl;

      double chi=0.1;
      double nB=(n.n+pr.n)*chi+(n2.n+p2.n)*(1.0-chi)+1.0e-4;
      double Ye=(pr.n*chi+p2.n*(1.0-chi))/nB;
      ap.calc_liqgas_temp_e(n,pr,n2,p2,nB,Ye,2.0/197.33,th,th2,chi);
      t.test_rel(n.mu,n2.mu,1.0e-8,"mun");
      t.test_rel(pr.mu,p2.mu,1.0e-8,"mup");
      cout << n.n << " " << pr.n << " " << n2.n << " " << p2.n << endl;
    }

    n.n=0.08;
    pr.n=0.08;
    ap.calc_e(n,pr,hb);
    cout << (hb.ed/0.16-n.m)*hc_mev_fm << endl;
    double mund, mupd, munde, mupde;
    thermo th3;
    ap.check_mu(n,pr,th3,mund,mupd,munde,mupde);
    t.test_abs(n.mu,mund,fabs(munde),"neutron chem pot.");
    t.test_abs(pr.mu,mupd,fabs(mupde),"proton chem pot.");
    double end, ende;
    ap.check_en(n,pr,5.0/hc_mev_fm,th3,end,ende);
    t.test_abs(th3.en,end,ende,"entropy");

    if (true) {
      // Testing low-density Neutron matter for APR
      // Show that the energy per baryon is bad
      n.n=2.6e-5;
      pr.n=0.0;
      ap.calc_e(n,pr,hb);
      cout << (hb.ed/n.n-n.m)*hc_mev_fm << endl;
    }

    cout.setf(ios::scientific);
    cout.precision(4);

    cout << "A18+UIX*+deltav" << endl << endl;

    cout << "Table VI on page 1809 of PRC 58 (1998) 1804:" << endl;
    for(barn=0.04;barn<=0.96;barn+=0.04) {
      if (barn>0.25) barn+=0.04;
      if (barn>0.65) barn+=0.08;

      n.n=barn/2.0;
      pr.n=n.n;
      ap.calc_e(n,pr,hb);
      cout << barn << " " << (hb.ed/barn-n.m)*hc_mev_fm << endl;

    }

    cout << "Table VII on page 1809 of PRC 58 (1998) 1804:" << endl;
    pr.n=0.0;
    for(barn=0.02;barn<=0.96;barn+=0.04) {
      if (barn>0.25) barn+=0.04;
      if (barn>0.65) barn+=0.08;

      n.n=barn;
      ap.calc_e(n,pr,hb);
      cout << barn << " " << (hb.ed/barn-n.m)*hc_mev_fm << endl;

      if (barn<0.03) barn=0.0;
    }

    cout.precision(6);
  
    ap.set_n_and_p(n,pr);
    ap.fn0(0.0,barn,dtemp);
    cout << "Saturation density: " << barn << endl;
    cout << "Binding energy: " << ap.feoa(barn)*hc_mev_fm << endl;
    cout << "Symmetry energy: " << ap.fesym(barn)*hc_mev_fm << endl;
    cout << "Compressibility: " << ap.fcomp(barn)*hc_mev_fm << endl;
    cout << "Effective mass: " << ap.fmsom(barn) << endl;
    cout << "Skewness: " << ap.fkprime(barn)*hc_mev_fm << endl;
    cout << endl;

    t.test_rel(barn,0.16,5.0e-3,"Saturation density.");
    t.test_rel(ap.feoa(barn)*hc_mev_fm,-16.0,1.0e-4,"Binding energy.");
    t.test_rel(ap.fesym(barn)*hc_mev_fm,32.5,1.0e-2,"Symmetry energy.");
    t.test_rel(ap.fcomp(barn)*hc_mev_fm,266.0,1.0e-3,"Compressibility.");

    // Testing chemical potentials in LDP
    double edold, ednew;
    ap.pion=1;
    for(n.n=0.03;n.n<=0.301;n.n+=0.05) {
      pr.n=n.n/4.0;
      ap.calc_e(n,pr,hb);
      edold=hb.ed;

      n.n+=1.0e-4;
      ap.calc_e(n,pr,hb);
      ednew=hb.ed;

      t.test_rel(n.mu,(ednew-edold)/1.0e-4,1.0e-3,"mun");
      n.n-=1.0e-4;
      pr.n+=1.0e-4;
      ap.calc_e(n,pr,hb);
      ednew=hb.ed;

      t.test_rel(pr.mu,(ednew-edold)/1.0e-4,1.0e-3,"mup");
    }

    // Testing chemical potentials in HDP
    ap.pion=2;
    for(n.n=0.14;n.n<=1.0;n.n+=0.1) {
      pr.n=n.n/2.0;
      ap.calc_e(n,pr,hb);
      edold=hb.ed;

      n.n+=1.0e-4;
      ap.calc_e(n,pr,hb);
      ednew=hb.ed;

      t.test_rel(n.mu,(ednew-edold)/1.0e-4,1.0e-3,"mun2");
      n.n-=1.0e-4;
      pr.n+=1.0e-4;
      ap.calc_e(n,pr,hb);
      ednew=hb.ed;

      t.test_rel(pr.mu,(ednew-edold)/1.0e-4,1.0e-3,"mup2");
    }

    // Testing I/O
    ap.n0=barn;
    ap.eoa=ap.feoa(barn);
    ap.esym=ap.fesym(barn);
    ap.comp=ap.fcomp(barn);
    ap.msom=ap.fmsom(barn);
    cout << endl;

    // Testing qij's
    n.n=0.09;
    pr.n=0.07;
    double qnn, qnp, qpp, d1, d2, d3, d4, d5, d6;
    double qqnn, qqnp, qqpp, dd1, dd2, dd3, dd4, dd5, dd6;
    ap.gradient_qij(n,pr,hb,qnn,qnp,qpp,d1,d2,d3,d4,d5,d6);
    cout << qnn << " " << qnp << " " << qpp << endl;
    cout << d1 << " " << d2 << " " << d3 << "\n"
	 << d4 << " " << d5 << " " << d6 << endl;
    ap.gradient_qij2(n.n,pr.n,qqnn,qqnp,qqpp,dd1,dd2,dd3,dd4,dd5,dd6);
    t.test_rel(qnn,qqnn,5.0e-4,"qnn");
    t.test_rel(qnp,qqnp,5.0e-4,"qnp");
    t.test_rel(qpp,qqpp,5.0e-4,"qpp");
    t.test_rel(d1,dd1,5.0e-4,"d1");
    t.test_rel(d2,dd2,5.0e-4,"d2");
    t.test_rel(d3,dd3,5.0e-4,"d3");
    t.test_rel(d4,dd4,5.0e-4,"d4");
    t.test_rel(d5,dd5,5.0e-4,"d5");
    t.test_rel(d6,dd6,5.0e-4,"d6");
    cout << endl;
  
    cout << "Check compressibility." << endl;
    ap.pion=eos_had_apr::ldp;
    for(double nb=0.08;nb<=0.32001;nb+=0.08) {
      ap.parent_method=true;
      double tx1=ap.fcomp(nb);
      ap.parent_method=false;
      double tx2=ap.fcomp(nb);
      t.test_rel(tx1,tx2,1.0e-8,"comp ldp");
    }
    ap.pion=eos_had_apr::hdp;
    for(double nb=0.17;nb<=0.60001;nb+=0.08) {
      ap.parent_method=true;
      double tx1=ap.fcomp(nb);
      ap.parent_method=false;
      double tx2=ap.fcomp(nb);
      t.test_rel(tx1,tx2,1.0e-8,"comp hdp");
    }
    cout << endl;

    cout << "Check symmetry energy." << endl;
    ap.pion=eos_had_apr::ldp;
    for(double nb=0.16;nb<=0.32001;nb+=0.08) {
      ap.parent_method=true;
      double tx1=ap.fesym_diff(nb);
      ap.parent_method=false;
      double tx2=ap.fesym_diff(nb);
      t.test_rel(tx1,tx2,1.0e-8,"esym_diff ldp");
    }
    ap.pion=eos_had_apr::hdp;
    for(double nb=0.17;nb<=0.60001;nb+=0.08) {
      ap.parent_method=true;
      double tx1=ap.fesym_diff(nb);
      ap.parent_method=false;
      double tx2=ap.fesym_diff(nb);
      t.test_rel(tx1,tx2,1.0e-8,"esym_diff hdp");
    }
    cout << endl;

  }

  // ----------------------------------------------------
  // Finite temperature testing
  // ----------------------------------------------------

  {
    fermion n(939.0/hc_mev_fm,2.0), pr(939.0/hc_mev_fm,2.0);
    thermo hb;
    double nb, littlea, eth, pth, eta, eoa0, pr0;
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    double tmp10, tmp11, tmp12;
    double thc=1.0, tsubf;
    eos_had_apr at;
  
    n.non_interacting=false;
    pr.non_interacting=false;

    at.select(1);

    cout.setf(ios::scientific);

    n.n=0.16;
    pr.n=0.16;
    at.calc_e(n,pr,hb);
    tmp1=hb.ed;
    tmp2=n.n;
    tmp3=pr.n;
    tmp4=n.ms;
    tmp5=pr.ms;
    tmp6=n.mu;
    tmp7=pr.mu;
    tmp8=n.ed;
    tmp9=pr.ed;
    tmp10=n.pr;
    tmp11=pr.pr;
    tmp12=hb.pr;
    at.calc_temp_e(n,pr,0.001,hb);
    t.test_rel(tmp1,hb.ed,1.0e-6,"test(1)");
    t.test_rel(tmp2,n.n,1.0e-6,"test(2)");
    t.test_rel(tmp3,pr.n,1.0e-6,"test(3)");
    t.test_rel(tmp4,n.ms,1.0e-6,"test(4)");
    t.test_rel(tmp5,pr.ms,1.0e-6,"test(5)");
    t.test_rel(tmp6,n.mu,1.0e-6,"test(6)");
    t.test_rel(tmp7,pr.mu,1.0e-6,"test(7)");
    t.test_rel(tmp8,n.ed,1.0e-6,"test(8)");
    t.test_rel(tmp9,pr.ed,1.0e-6,"test(9)");
    t.test_rel(tmp10,n.pr,1.0e-3,"test(10)");
    t.test_rel(tmp11,pr.pr,1.0e-3,"test(11)");
    t.test_rel(tmp12,hb.pr,1.0e-4,"test(12)");
    t.test_rel(hb.ed,-hb.pr+n.mu*n.n+pr.mu*pr.n+0.001*hb.en,1.0e-6,
	       "Thermodynamic identity");

    t.set_output_level(1);
    cout << "Symmetric nuclear matter" << endl;
  
    for(nb=0.04;nb<=0.8;nb+=0.2) {
      if (nb>0.25) nb+=0.04;
      if (nb>0.65) nb+=0.08;
    
      n.n=nb/2.0;
      pr.n=n.n;
      n.kf=pow(3.0*pi2*n.n,1.0/3.0);
      littlea=pi2*n.ms/2.0/n.kf/n.kf;
      tsubf=n.kf*n.kf/2.0/n.ms;
    
      at.calc_e(n,pr,hb);
      eoa0=hb.ed/nb-n.m;
      pr0=hb.pr;
      cout << "nb=" << nb << " e(T=0)=" << eoa0*hc_mev_fm 
	   << " p(T=0)=" << pr0*hc_mev_fm << " tsubf=" 
	   << tsubf*hc_mev_fm << endl;
    
      //--------------------------------------------
      // Degenerate
      thc=tsubf*hc_mev_fm/5.0;
    
      cout << "T=" << thc << " MeV" << " T/T_F=" 
	   << thc/hc_mev_fm/tsubf << endl;
    
      at.calc_temp_e(n,pr,thc/hc_mev_fm,hb);
    
      eth=pow(thc/hc_mev_fm,2.0)*littlea;
      pth=nb*2.0/3.0*littlea*pow(thc/hc_mev_fm,2.0)*(1.0-1.5*(n.ms/n.m-1.0));
    
      cout << "Deg. e " << (hb.ed/nb-n.m)*hc_mev_fm << " " 
	   << (eth+eoa0)*hc_mev_fm << endl;
      cout << "     p " << hb.pr*hc_mev_fm << " " 
	   << (pth+pr0)*hc_mev_fm << endl;
      cout << "     s " << hb.en << " " 
	   << 2.0/3.0*thc/hc_mev_fm*n.kf*n.ms << endl;
      //t.test_rel((hb.ed/nb-n.m)*hc_mev_fm,(eth+eoa0)*hc_mev_fm,5.0e-1,
      //"Deg. e");
      t.test_rel(hb.pr*hc_mev_fm,(pth+pr0)*hc_mev_fm,5.0e-1,"Deg. p");
      //t.test_rel(hb.en,2.0/3.0*thc/hc_mev_fm*n.kf*n.ms,5.0e-1,"Deg. s");

      //--------------------------------------------
      // Very degenerate

      thc=tsubf*hc_mev_fm/20.0;

      cout << "T=" << thc << " MeV" << " T/T_F=" 
	   << thc/hc_mev_fm/tsubf << endl;
    
      at.calc_temp_e(n,pr,thc/hc_mev_fm,hb);
    
      eth=pow(thc/hc_mev_fm,2.0)*littlea;
      pth=nb*2.0/3.0*littlea*pow(thc/hc_mev_fm,2.0)*(1.0-1.5*(n.ms/n.m-1.0));
    
      cout << "Deg. e " << (hb.ed/nb-n.m)*hc_mev_fm << " " 
	   << (eth+eoa0)*hc_mev_fm << endl;
      cout << "     p " << hb.pr*hc_mev_fm << " " 
	   << (pth+pr0)*hc_mev_fm << endl;
      cout << "     s " << hb.en << " " 
	   << 2.0/3.0*thc/hc_mev_fm*n.kf*n.ms << endl;
      //t.test_rel((hb.ed/nb-n.m)*hc_mev_fm,(eth+eoa0)*hc_mev_fm,1.0e-1,
      //"Deg. e");
      t.test_rel(hb.pr*hc_mev_fm,(pth+pr0)*hc_mev_fm,1.0e-1,"Deg. p");
      //t.test_rel(hb.en,2.0/3.0*thc/hc_mev_fm*n.kf*n.ms,1.0e-1,"Deg. s");

      //--------------------------------------------
      // Non-degenerate

      thc=tsubf*5.0*hc_mev_fm;
    
      cout << "T=" << thc << " MeV" << " T/T_F=" 
	   << thc/hc_mev_fm/tsubf << endl;
    
      at.calc_temp_e(n,pr,thc/hc_mev_fm,hb);
    
      eth=pow(thc/hc_mev_fm,2.0)*littlea;
      pth=nb*2.0/3.0*littlea*pow(thc/hc_mev_fm,2.0)*(1.0-1.5*(n.ms/n.m-1.0));
    
      pth=nb*thc/hc_mev_fm*(1.0+nb*pow(pi/4.0/n.ms/thc*hc_mev_fm,1.5))-
	0.2*nb*n.kf*n.kf/n.ms;
      eth=1.5*thc/hc_mev_fm*(1.0+nb*pow(pi/4.0/n.ms/thc*hc_mev_fm,1.5))-
	0.3*n.kf*n.kf/n.ms;
      eta=4.0/3.0/sqrt(pi)*pow(hc_mev_fm/thc*n.kf*n.kf/2.0/n.ms,1.5);
    
      cout << "Ndg. e " << (hb.ed/nb-n.m)*hc_mev_fm << " " 
	   << (eth+eoa0)*hc_mev_fm << endl;
      cout << "     p " << hb.pr*hc_mev_fm << " " 
	   << (pth+pr0)*hc_mev_fm << endl;
      cout << "     s " << hb.en/nb << " " 
	   << 2.5-log(eta)+eta/pow(2,3.5) << endl;
      //t.test_rel((hb.ed/nb-n.m)*hc_mev_fm,(eth+eoa0)*hc_mev_fm,1.0e-1,
      //"Ndg. e");
      // t.test_rel(hb.pr*hc_mev_fm,(pth+pr0)*hc_mev_fm,1.0e-0,"Ndg. p");
      //t.test_rel(hb.en/nb,2.5-log(eta)+eta/pow(2.0,3.5),1.0e-1,"Ndg. s");
    
      //--------------------------------------------
      // Very non-degenerate

      thc=tsubf*20.0*hc_mev_fm;

      cout << "T=" << thc << " MeV" << " T/T_F=" 
	   << thc/hc_mev_fm/tsubf << endl;
    
      at.calc_temp_e(n,pr,thc/hc_mev_fm,hb);
    
      eth=pow(thc/hc_mev_fm,2.0)*littlea;
      pth=nb*2.0/3.0*littlea*pow(thc/hc_mev_fm,2.0)*(1.0-1.5*(n.ms/n.m-1.0));
    
      pth=nb*thc/hc_mev_fm*(1.0+nb*pow(pi/4.0/n.ms/thc*hc_mev_fm,1.5))-
	0.2*nb*n.kf*n.kf/n.ms;
      eth=1.5*thc/hc_mev_fm*(1.0+nb*pow(pi/4.0/n.ms/thc*hc_mev_fm,1.5))-
	0.3*n.kf*n.kf/n.ms;
      eta=4.0/3.0/sqrt(pi)*pow(hc_mev_fm/thc*n.kf*n.kf/2.0/n.ms,1.5);
    
      cout << "Ndg. e " << (hb.ed/nb-n.m)*hc_mev_fm << " " 
	   << (eth+eoa0)*hc_mev_fm << endl;
      cout << "     p " << hb.pr*hc_mev_fm << " " 
	   << (pth+pr0)*hc_mev_fm << endl;
      cout << "     s " << hb.en/nb << " " 
	   << 2.5-log(eta)+eta/pow(2,3.5) << endl;
      //t.test_rel((hb.ed/nb-n.m)*hc_mev_fm,(eth+eoa0)*hc_mev_fm,1.0e-1,
      //"Ndg. e");
      //  t.test_rel(hb.pr*hc_mev_fm,(pth+pr0)*hc_mev_fm,1.0e0,"Ndg. p");
      //t.test_rel(hb.en/nb,2.5-log(eta)+eta/pow(2.0,3.5),1.0e-1,"Ndg. s");
    
    }
    
    n.n=0.08;
    pr.n=0.08;
    n.nu=n.m;
    pr.nu=pr.m;
    // Test chemical potentials 
    cout << "Test chemical potentials at finite T: " << endl;
    thermo th;
    at.calc_temp_e(n,pr,4.0/hc_mev_fm,th);
    double mun=n.mu, mup=pr.mu, ed0=th.ed;
    n.n+=1.0e-4;
    at.calc_temp_e(n,pr,4.0/hc_mev_fm,th);
    double edn=th.ed;
    n.n-=1.0e-4;
    pr.n+=1.0e-4;
    at.calc_temp_e(n,pr,4.0/hc_mev_fm,th);
    double edp=th.ed;
    pr.n-=1.0e-4;
    
    t.test_rel((edn-ed0)/1.0e-4,mun,5.0e-4,"mun T>0");
    t.test_rel((edp-ed0)/1.0e-4,mup,5.0e-4,"mup T>0");

    n.n=0.08;
    pr.n=0.08;
    // Test chemical potentials 
    cout << "Test chemical potentials at zero T: " << endl;
    at.calc_e(n,pr,th);
    mun=n.mu; mup=pr.mu; ed0=th.ed;
    n.n+=1.0e-4;
    at.calc_e(n,pr,th);
    edn=th.ed;
    n.n-=1.0e-4;
    pr.n+=1.0e-4;
    at.calc_e(n,pr,th);
    edp=th.ed;
    pr.n-=1.0e-4;
    
    t.test_rel((edn-ed0)/1.0e-4,mun,1.0e-4,"mun T=0");
    t.test_rel((edp-ed0)/1.0e-4,mup,1.0e-4,"mup T=0");

    /*
      cout << at.fesym(0.16)*hc_mev_fm << endl;
      cout << at.fesym_T(0.16,0.3/hc_mev_fm)*hc_mev_fm << endl;
      cout << at.fesym_T(0.16,1.0/hc_mev_fm)*hc_mev_fm << endl;
      cout << at.fesym_T(0.16,3.0/hc_mev_fm)*hc_mev_fm << endl;
      cout << at.fesym_T(0.16,10.0/hc_mev_fm)*hc_mev_fm << endl;
    */
  }
  

  t.report();
  return 0;
}

