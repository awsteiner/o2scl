/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2022, Andrew W. Steiner
  
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
#include <o2scl/eos_crust_virial.h>
#include <o2scl/fermion_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  double hc=hc_mev_fm;

  fermion n(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  n.inc_rest_mass=false;
  fermion p(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  p.inc_rest_mass=false;
  fermion e(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_electron),2.0);

  thermo th;
  double T;
  eos_crust_virial ve;
  fermion_rel rf;

  vector<double> nba, pra, soaa;
  interp<vector<double> > it(itp_cspline);

  /* 
     This reproduces the pressure at n_B=0.05 fm^{-3} in 
     Figure 10 in Horowitz et al. (2005)
  */
  T=10.0/hc;
  for(double mu=-50.0;mu<-10.0;mu/=1.03) {
    n.mu=mu/hc;
    p.mu=mu/hc;
    ve.calc_temp_p(n,p,T,th);
    nba.push_back(n.n+p.n+4.0*ve.alpha.n);
    pra.push_back(th.pr*hc);
  }
  t.test_rel(it.eval(0.05,nba.size(),nba,pra),0.162,0.01,"Fig 10.");
  nba.clear();
  pra.clear();

  T=4.0/hc;
  for(double mu=-50.0;mu<-5.0;mu/=1.03) {
    n.mu=mu/hc;
    p.mu=mu/hc;
    ve.calc_temp_p(n,p,T,th);
    e.n=p.n+2.0*ve.alpha.n;
    // Need to provide an initial guess to electron chemical potential
    e.mu=e.m;
    rf.calc_density(e,T);
    nba.push_back(n.n+p.n+4.0*ve.alpha.n);
    pra.push_back((th.pr+e.pr)*hc);
    soaa.push_back((th.en+e.en)/(n.n+p.n+4.0*ve.alpha.n));
  }
  t.test_rel(it.eval(0.0016,nba.size(),nba,pra),0.015,5.0e-2,"Fig 11 right.");
  t.test_rel(it.eval(0.05,nba.size(),nba,soaa),1.35,1.0e-1,"Fig 13 right.");
  t.test_rel(it.eval(1.0e-4,nba.size(),nba,soaa),7.5,2.0e-2,"Fig 13 left.");
  
  t.report();
  
  return 0;
}
