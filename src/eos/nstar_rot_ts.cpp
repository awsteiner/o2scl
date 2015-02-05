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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/nstar_rot.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/nstar_cold.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

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

  nstar_rot nst;
  nst.constants_rns();
  nst.eosC();
  nst.test1(t);
  nst.test2(t);
  nst.test3(t);
  nst.test4(t);
  nst.test5(t);
  nst.test6(t);
  nst.test7(t);
  nst.test8(t);

  eos_had_skyrme sk;
  load_sly4(sk);
  nstar_cold nco;
  nco.def_tov.verbose=0;
  nco.set_eos(sk);
  nco.calc_eos();
  nco.calc_nstar();
  
  o2_shared_ptr<table_units<> >::type eos=nco.get_eos_results();
  o2_shared_ptr<table_units<> >::type mvsr=nco.get_tov_results();

  nst.set_eos_fm(eos->get_nlines(),(*eos)["ed"],(*eos)["pr"],(*eos)["nb"]);

  //cout << "Here." << endl;
  //nst.fix_cent_eden_non_rot(1.0e15);
  //cout << nst.Mass << endl;

  t.report();
  
  return 0;
}
