/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2014, Andrew W. Steiner
  
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

#include <o2scl/eos_had_gogny.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/fermion.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  // Ensure that this works without GNU units
  o2scl_settings.get_convert_units().use_gnu_units=false;
  
  fermion n(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  fermion p(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  thermo th;
    
  eos_had_gogny ge;

  cout.setf(ios::showpos);

  cout << "D1S: " << endl;
  gogny_load(ge,"d1s");
  cout << " nb            E_nuc         P_nuc         E_neut        P_neut"
       << endl;
  for(double nb=0.02;nb<0.2001;nb+=0.02) {
    n.n=nb/2.0;
    p.n=nb/2.0;
    ge.calc_e(n,p,th);
    cout << nb << " " 
	 << (th.ed-n.n*n.m-p.n*p.m)/nb*o2scl_const::hc_mev_fm << " "
	 << th.pr*o2scl_const::hc_mev_fm << " ";
    n.n=nb;
    p.n=0.0;
    ge.calc_e(n,p,th);
    cout << (th.ed-n.n*n.m-p.n*p.m)/nb*o2scl_const::hc_mev_fm << " "
	 << th.pr*o2scl_const::hc_mev_fm << endl;
  }
  cout << endl;

  for(double a=0.0;a<1.01;a+=0.1) {
    p.n=0.16*(1.0-a)/2.0;
    n.n=0.16-p.n;
    cout << n.n << " " << p.n << " ";
    ge.calc_e(n,p,th);
    cout << (th.ed-n.n*n.m-p.n*p.m)/0.16*o2scl_const::hc_mev_fm << " "
	 << th.pr*o2scl_const::hc_mev_fm << endl;
  }
  cout << endl;

  ge.saturation();
  cout << "n0,B,S,K: " << ge.n0 << " " << ge.eoa*hc_mev_fm << " " 
       << ge.esym*hc_mev_fm << " " << ge.comp*hc_mev_fm << endl;
  cout << ge.fesym_diff(ge.n0)*hc_mev_fm << endl;
  cout << endl;
  
  cout << "D1N: " << endl;
  gogny_load(ge,"d1n");
  cout << " nb            E_nuc         P_nuc         E_neut        P_neut"
       << endl;
  for(double nb=0.02;nb<0.2001;nb+=0.02) {
    n.n=nb/2.0;
    p.n=nb/2.0;
    ge.calc_e(n,p,th);
    cout << nb << " " 
	 << (th.ed-n.n*n.m-p.n*p.m)/nb*o2scl_const::hc_mev_fm << " "
	 << th.pr*o2scl_const::hc_mev_fm << " ";
    n.n=nb;
    p.n=0.0;
    ge.calc_e(n,p,th);
    cout << (th.ed-n.n*n.m-p.n*p.m)/nb*o2scl_const::hc_mev_fm << " "
	 << th.pr*o2scl_const::hc_mev_fm << endl;
  }
  cout << endl;
  
  ge.saturation();
  cout << "n0,B,S,K: " << ge.n0 << " " << ge.eoa*hc_mev_fm << " " 
       << ge.esym*hc_mev_fm << " " << ge.comp*hc_mev_fm << endl;
  cout << ge.fesym_diff(ge.n0)*hc_mev_fm << endl;
  cout << endl;

  t.report();
  return 0;
}
