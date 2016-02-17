/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2016, Andrew W. Steiner
  
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
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  nstar_rot nst;
  nst.constants_rns();
  nst.test1(t);
  nst.test2(t);
  nst.test3(t);
  nst.test4(t);
  nst.test5(t);
  nst.test6(t);
  nst.test7(t);
  nst.test8(t);

  eos_had_skyrme sk;
  o2scl_hdf::skyrme_load(sk,"Gs");

  nstar_cold nco;
  nco.def_tov.verbose=0;
  nco.set_eos(sk);

  // Compute the Skyrme EOS in beta-equilibrium
  nco.calc_eos();
  std::shared_ptr<table_units<> > eos=nco.get_eos_results();

  // Evalulate the mass-radius curve
  nco.calc_nstar();
  std::shared_ptr<table_units<> > mvsr=nco.get_tov_results();

  // Lookup the central energy density of a 1.4 Msun neutron star
  // in g/cm^3
  convert_units &cu=o2scl_settings.get_convert_units();
  double ed_cent=mvsr->get("ed",mvsr->lookup("gm",1.4));
  ed_cent=cu.convert("1/fm^4","g/cm^3",ed_cent);

  // Send the EOS to the nstar_rot object
  eos_nstar_rot_interp p;
  table_units<> new_eos;
  new_eos.line_of_names("ed pr nb");
  new_eos.set_unit("ed","1/fm^4");
  new_eos.set_unit("pr","1/fm^4");
  new_eos.set_unit("nb","1/fm^3");
  for(size_t i=0;i<eos->get_nlines();i+=2) {
    if (eos->get("nb",i)>0.0801) {
      double line[3]={eos->get("ed",i),eos->get("pr",i),eos->get("nb",i)};
      new_eos.line_of_data(3,line);
    }
  }
  p.set_eos_fm(new_eos.get_nlines(),new_eos["ed"],
	       new_eos["pr"],new_eos["nb"]);
  nst.set_eos(p);
  
  // Compute the mass of the non-rotating configuration with the
  // same energy density
  nst.fix_cent_eden_non_rot(ed_cent);
  // Compare with with the answer from nstar_rot
  t.test_rel(nst.Mass/nst.MSUN,1.4,0.015,"correct mass");

  t.report();
  
  return 0;
}
