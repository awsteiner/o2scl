/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2024, Andrew W. Steiner
  
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
// sphinx-example-start
/* Example: ex_eos_had_rmf.cpp
   -------------------------------------------------------------------
   A simple example for an RMF EOS. See "License Information" 
   section of the documentation for license information.
*/
#include <o2scl/test_mgr.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/fermion_deriv_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

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
  re.err_nonconv=false;
  re.def_sat_mroot.err_nonconv=false;
  re.def_sat_mroot.def_jac.err_nonconv=false;
  re.def_mroot.err_nonconv=false;
  re.def_mroot.def_jac.err_nonconv=false;
  
  int sret=re.saturation();
  
  if (sret!=0 || re.n0<0.08) {
    re.set_fields(0.2,0.1,0.01);
    re.def_neutron.mu=5.0;
    re.def_proton.mu=5.0;
    sret=re.saturation();
    if (sret!=0) {
      O2SCL_ERR("RMF EOS failed.",o2scl::exc_efailed);
    }
  }

  // Return the convergence error flag to the default value
  re.err_nonconv=true;
  re.def_sat_mroot.err_nonconv=true;
  re.def_sat_mroot.def_jac.err_nonconv=true;
  re.def_mroot.err_nonconv=true;
  re.def_mroot.def_jac.err_nonconv=true;
  
  // Compute the saturation density and the symmetry energy
  cout << "FSUGold: " << re.n0 << endl;
  cout << "FSUGold: " << re.esym*hc_mev_fm << endl;
  t.test_rel(re.n0,0.1481,4.0e-4,"FSUGold saturation density.");

  // Commpute the beta-equilibrium EOS and solve the TOV equations
  nstar_cold nc;
  nc.set_eos(re);
  nc.calc_eos();
  nc.calc_nstar();

  // Output gravitational mass-radius curve
  cout << endl;
  shared_ptr<table_units<>> tov=nc.get_tov_results();
  for(size_t i=0;i<tov->get_nlines();i++) {
    cout << tov->get("gm",i) << " " << tov->get("r",i) << endl;
  }

  t.report();

  return 0;
}


