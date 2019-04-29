/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/* Example: ex_nstar_cold.cpp
   -------------------------------------------------------------------
   This example solves the TOV equations using class nstar_cold using a
   relativistic mean-field EOS from class eos_had_rmf.
*/

#include <o2scl/nstar_cold.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

// For hc_mev_fm
using namespace o2scl_const;

// A simple function to load the NL3 model
int load_nl3(eos_had_rmf &rmf);

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  nstar_cold nst;

  // Initialize EOS
  eos_had_rmf rmf;
  load_nl3(rmf);

  rmf.def_neutron.mu=5.0;
  rmf.def_proton.mu=5.0;
  rmf.set_fields(0.2,0.1,-0.05);
  
  rmf.err_nonconv=false;
  rmf.def_sat_mroot.err_nonconv=false;
  rmf.def_sat_mroot.def_jac.err_nonconv=false;
  rmf.def_mroot.err_nonconv=false;
  rmf.def_mroot.def_jac.err_nonconv=false;
  
  int sret=rmf.saturation();

  if (sret!=0) {
    
    rmf.set_fields(0.1,0.07,0.0);
    rmf.def_neutron.mu=5.0;
    rmf.def_proton.mu=5.0;

    rmf.err_nonconv=true;
    rmf.def_sat_mroot.err_nonconv=true;
    rmf.def_sat_mroot.def_jac.err_nonconv=true;
    rmf.def_mroot.err_nonconv=true;
    rmf.def_mroot.def_jac.err_nonconv=true;
    
    rmf.saturation();
    
  } else {
    
    rmf.err_nonconv=true;
    rmf.def_sat_mroot.err_nonconv=true;
    rmf.def_sat_mroot.def_jac.err_nonconv=true;
    rmf.def_mroot.err_nonconv=true;
    rmf.def_mroot.def_jac.err_nonconv=true;

  }
  
  cout << "Saturation density: " << rmf.n0 << endl;
  cout << "Binding energy: " << rmf.eoa*hc_mev_fm << endl;
  cout << "Effective mass: " << rmf.msom << endl;
  cout << "Symmetry energy: " << rmf.esym*hc_mev_fm << endl;
  cout << "Compressibility: " << rmf.comp*hc_mev_fm << endl;
  
  // Compute EOS, include muons
  nst.include_muons=true;
  nst.set_eos(rmf);
  nst.verbose=1;
  nst.def_root.verbose=1;
  rmf.def_mroot.verbose=1;
  cout << "Hx." << endl;
  rmf.calc_e_steps=200;
  nst.calc_eos();
  std::shared_ptr<table_units<> > te=nst.get_eos_results();

  // Compute mass vs. radius
  nst.calc_nstar();
  std::shared_ptr<table_units<> > tr=nst.get_tov_results();
  cout << "Maximum mass: " << tr->max("gm") << endl;
  cout << "Radius of maximum mass star: " 
       << tr->get("r",tr->lookup("gm",tr->max("gm"))) << endl;
  cout << "Central baryon density of maximum mass star: ";
  cout << tr->get("nb",tr->lookup("gm",tr->max("gm"))) << endl;

  // Output EOS and TOV results to files
  hdf_file hf;
  hf.open_or_create("ex_nstar_cold_eos.o2");
  hdf_output(hf,*te,"eos");
  hf.close();
  hf.open_or_create("ex_nstar_cold_tov.o2");
  hdf_output(hf,*tr,"tov");
  hf.close();

  t.report();
  return 0;
}
// End of example

int load_nl3(eos_had_rmf &rmf) {

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
