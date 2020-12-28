/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <o2scl/eos_had_potential.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/nstar_cold.h>

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

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  eos_had_potential go;
  go.B=106.35/hc_mev_fm;
  go.sigma=4.0/3.0;
  go.rho0=0.16;
  go.Cu=-103.40/hc_mev_fm;
  go.Cl=-11.70/hc_mev_fm;
  go.form=go.mdi_form;
  go.Lambda=cbrt(1.5*pi2*go.rho0);
  
  go.x=0.0;
  
  go.Au=-95.98/hc_mev_fm-2.0*go.B*go.x/(go.sigma+1.0);
  go.Al=-120.75/hc_mev_fm+2.0*go.B*go.x/(go.sigma+1.0);

  nstar_cold nst;
  nst.set_eos(go);
  nst.calc_eos();
  nst.calc_nstar();
  
  std::shared_ptr<table_units<> > te=nst.get_eos_results();
  cout << "EOS results: " << endl;
  double ed1=te->interp("nb",0.16,"ed");
  double pr1=te->interp("nb",0.16,"pr");
  te->summary(&cout);
  cout << endl;

  std::shared_ptr<table_units<> > tr=nst.get_tov_results();
  tr->summary(&cout);
  cout << endl;
  cout << "M_{max} = " << tr->max("gm") << " R_{max} = "
       << tr->get("r",tr->lookup("gm",tr->max("gm"))) << " cent. density = "
       << tr->get("nb",tr->lookup("gm",tr->max("gm"))) << endl;
  t.test_rel(tr->max("gm"),1.90373,2.0e-4,"M_max");
  t.test_rel(tr->get("r",tr->lookup("gm",tr->max("gm"))),
	     9.956775,1.0e-3,"R_max");
  cout << endl;

  // Check that EOS corresponds to result in M vs. R table
  t.test_rel(ed1,tr->interp("nb",0.16,"ed"),2.0e-3,"ed");
  t.test_rel(pr1,tr->interp("nb",0.16,"pr"),4.0e-3,"pr");

  nstar_hot nh;
  eos_had_skyrme sk;
  load_sly4(sk);
  nh.set_eos_T(sk);

  nh.calc_eos_T(10.0/hc_mev_fm);
  te=nh.get_eos_results();
  te->summary(&cout);
  
  t.report();
  return 0;
}

