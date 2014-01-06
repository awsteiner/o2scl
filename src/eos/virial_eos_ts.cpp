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

#include <o2scl/test_mgr.h>

#include "virial_eos.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  double hc=hc_mev_fm;

  // Ensure that this works without GNU units
  o2scl_settings.get_convert_units().use_gnu_units=false;
  
  fermion n(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  fermion p(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  n.inc_rest_mass=false;
  p.inc_rest_mass=false;

  thermo th;
  double T;
  virial_eos ve;

  T=10.0/hc;
  cout << "mu nn_np n_B P" << endl;
  for(double mu=-50.0;mu<-10.0;mu/=1.2) {
    n.mu=mu/hc;
    p.mu=mu/hc;
    ve.calc_temp_p(n,p,T,th);
    cout << mu << " " << n.n+p.n << " " 
	 << n.n+p.n+4.0*ve.alpha.n << " " << th.pr*hc << endl;
  }
  cout << endl;

  cout << "mu nn_np n_B P xa SoA" << endl;
  for(double mu=-50.0;mu<-5.0;mu/=1.03) {
    n.mu=mu/hc;
    p.mu=mu/hc;
    ve.calc_temp_p(n,p,T,th);
    cout << mu << " " << n.n+p.n << " " 
	 << n.n+p.n+4.0*ve.alpha.n << " " << th.pr*hc << " "
	 << 4.0*ve.alpha.n/(n.n+p.n+4.0*ve.alpha.n) << " " 
	 << th.en/(n.n+p.n+4.0*ve.alpha.n) << endl;
  }
  cout << endl;
  
  t.report();
  
  return 0;
}
