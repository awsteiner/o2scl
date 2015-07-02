/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/root_cern.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

class beta_temp {

public:
  
  o2scl::eos_had_rmf &re;
  o2scl::fermion &n;
  o2scl::fermion &p;
  
  beta_temp(o2scl::eos_had_rmf &ret, o2scl::fermion &nt,
	    o2scl::fermion &pt) : re(ret), n(nt), p(pt) {
    e.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  }
  
  fermion e;
  fermion_zerot fzt;
  double barn;
  
  double solve_fun(double x) {

    //cout << "x: " << x << " " << barn << endl;
    p.n=x*barn;
    n.n=barn-p.n;
    re.calc_temp_e(n,p,8.0/hc_mev_fm,re.def_thermo);
    
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    //cout << "y: " << p.n-e.n << endl;
    return p.n-e.n;
  }

};


int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  string cmd=((string)"cp ")+o2scl_settings.get_data_dir()+
    "/rmfdata/FSUGold.o2 temp.o2";
  system(cmd.c_str());
  cmd=((string)"cp ")+o2scl_settings.get_data_dir()+
    "/rmfdata/FSUGold.o2 temp2.o2";
  system(cmd.c_str());

  eos_had_rmf re;

  re.def_sat_mroot.def_jac.set_epsrel(1.0e-6);
  re.def_sat_mroot.def_jac.set_epsmin(1.0e-15);
  re.def_sat_mroot.ntrial*=10;
  
  rmf_load(re,"FSUGold");
  cout << re.ms << " " << re.mw << " " << re.mr << endl;
  cout << re.cs << " " << re.cw << " " << re.cr << endl;
  cout << re.b << " " << re.c << " " << re.mnuc << endl;
  cout << re.zeta << " " << re.b1 << endl;

  // It turns out that FSUGold needs a better initial
  // guess than the default to get saturation right
  re.set_fields(0.2,0.1,0.01);
  re.saturation();

  cout << "FSUGold: " << re.n0 << endl;
  cout << "FSUGold: " << re.esym*hc_mev_fm << endl;
  cout << endl;
  
  hdf_file hf;
  hf.open_or_create("temp.o2");
  double gs2=112.1996;
  double gw2=204.5469;
  double gr2=138.4701;
  hf.setd("gs",sqrt(gs2));
  hf.setd("gw",sqrt(gw2));
  hf.setd("gr",sqrt(gr2)/2.0);
  hf.seti("oakstyle",1);
  hf.setd("zeta",0.06);
  hf.setd("g2",-pow(gs2,3.0/2.0)*1.4203/hc_mev_fm/2.0);
  hf.setd("g3",pow(gs2,2.0)*0.023762/6.0);
  hf.setd("b1",0.03*gw2);
  hf.close();

  rmf_load(re,"temp.o2",true);
  cout << re.ms << " " << re.mw << " " << re.mr << endl;
  cout << re.cs << " " << re.cw << " " << re.cr << endl;
  cout << re.b << " " << re.c << " " << re.mnuc << endl;
  cout << re.zeta << " " << re.b1 << endl;

  // It turns out that FSUGold needs a better initial
  // guess than the default to get saturation right
  re.set_fields(0.2,0.1,0.01);
  re.saturation();
  
  cout << "FSUGold: " << re.n0 << endl;
  cout << "FSUGold: " << re.esym*hc_mev_fm << endl;
  cout << endl;

  rmf_load(re,"IUFSU");
  cout << re.ms << " " << re.mw << " " << re.mr << endl;
  cout << re.cs << " " << re.cw << " " << re.cr << endl;
  cout << re.b << " " << re.c << " " << re.mnuc << endl;
  cout << re.zeta << " " << re.b1 << endl;
  re.set_fields(0.175,0.1,0.0001);
  re.saturation();
  cout << "IUFSU: " << re.n0 << endl;
  cout << "IUFSU: " << re.esym*hc_mev_fm << endl;
  cout << endl;

  hf.open_or_create("temp2.o2");
  gs2=99.4266;
  gw2=169.8349;
  gr2=184.6877;
  hf.setd("gs",sqrt(gs2));
  hf.setd("gw",sqrt(gw2));
  hf.setd("gr",sqrt(gr2)/2.0);
  hf.seti("oakstyle",1);
  hf.sets("reference",((string)"F.J. Fattoyev, C.J. Horowitz,")+
	  "J. Piekarewicz, and G. Shen, Phys. Rev. C 82 (2010) 055803.");
  hf.setd("zeta",0.03);
  hf.setd("g2",-pow(gs2,3.0/2.0)*3.3808/hc_mev_fm/2.0);
  hf.setd("g3",pow(gs2,2.0)*0.000296/6.0);
  hf.setd("b1",0.046*gw2);
  hf.close();

  rmf_load(re,"temp2.o2",true);
  cout << re.ms << " " << re.mw << " " << re.mr << endl;
  cout << re.cs << " " << re.cw << " " << re.cr << endl;
  cout << re.b << " " << re.c << " " << re.mnuc << endl;
  cout << re.zeta << " " << re.b1 << endl;

  re.set_fields(0.175,0.1,0.0001);
  re.saturation();
  cout << "IUFSU: " << re.n0 << endl;
  cout << "IUFSU: " << re.esym*hc_mev_fm << endl;
  cout << endl;

  re.def_mroot.def_jac.set_epsrel(1.0e-6);
  re.def_mroot.def_jac.set_epsmin(1.0e-15);
  re.def_mroot.ntrial*=10;

  // Test beta equilibrium at 8 MeV
  double sigma, omega, rho;
  double xp=0.02;
  double nb=0.02;

  re.def_neutron.m=939.0/hc_mev_fm;
  re.def_proton.m=939.0/hc_mev_fm;

  beta_temp bt(re,re.def_neutron,re.def_proton);
  bt.barn=0.02;

  funct11 bf=std::bind(std::mem_fn<double(double)>
		       (&beta_temp::solve_fun),
		       &bt,std::placeholders::_1);
  root_cern<> rt;
  rt.solve(xp,bf);
  bt.solve_fun(xp);

  re.get_fields(sigma,omega,rho);
  cout << re.def_neutron.n+re.def_proton.n << " "
       << re.def_proton.n/(re.def_neutron.n+re.def_proton.n) << endl;
  cout << re.def_neutron.mu*hc_mev_fm << " "
       << re.def_proton.mu*hc_mev_fm << endl;
  cout << re.def_neutron.nu*hc_mev_fm << " "
       << re.def_proton.nu*hc_mev_fm << endl;
  cout << sigma*hc_mev_fm << " " << omega*hc_mev_fm << " "
       << rho*hc_mev_fm << endl;
  cout << 1.0/re.def_neutron.ms*(re.def_neutron.ed-3.0*re.def_neutron.pr)
       << endl;
  cout << 1.0/re.def_proton.ms*(re.def_proton.ed-3.0*re.def_proton.pr)
       << endl;
  cout << re.def_thermo.ed/nb*hc_mev_fm-939.0 << endl;
  cout << endl;

  // Verify from EOS
  re.def_neutron.n=nb*(1.0-xp);
  re.def_proton.n=nb*xp;
  re.calc_temp_e(re.def_neutron,re.def_proton,8.0/hc_mev_fm,re.def_thermo);
  re.get_fields(sigma,omega,rho);
  cout << re.def_neutron.n+re.def_proton.n << " "
       << re.def_proton.n/(re.def_neutron.n+re.def_proton.n) << endl;
  cout << re.def_neutron.mu*hc_mev_fm << " "
       << re.def_proton.mu*hc_mev_fm << endl;
  cout << re.def_neutron.nu*hc_mev_fm << " "
       << re.def_proton.nu*hc_mev_fm << endl;
  cout << sigma*hc_mev_fm << " " << omega*hc_mev_fm << " "
       << rho*hc_mev_fm << endl;
  cout << 1.0/re.def_neutron.ms*(re.def_neutron.ed-3.0*re.def_neutron.pr)
       << endl;
  cout << 1.0/re.def_proton.ms*(re.def_proton.ed-3.0*re.def_proton.pr)
       << endl;
  cout << re.def_thermo.ed/nb*hc_mev_fm-939.0 << endl;
  cout << endl;

  t.report();

  return 0;
}


