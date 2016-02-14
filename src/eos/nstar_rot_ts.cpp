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
  if (false) {
    table_units<> new_eos;
    new_eos.line_of_names("ed pr nb");
    new_eos.set_unit("ed","1/fm^4");
    new_eos.set_unit("pr","1/fm^4");
    new_eos.set_unit("nb","1/fm^3");
    {
      // Read default EOS. 
      static const size_t nlines=73;
      double ed_arr[nlines]=
	{3.89999984e-18,3.93000002e-18,3.95000001e-18,4.07499982e-18,
	 5.80000020e-18,8.20000023e-18,2.25500006e-17,1.05999999e-16,
	 5.75000005e-16,5.22000020e-15,1.31099999e-14,3.29349991e-14,
	 8.27000008e-14,2.07800004e-13,5.21999978e-13,1.31100003e-12,
	 3.29400006e-12,4.14649998e-12,8.27499961e-12,1.65100000e-11,
	 3.29450009e-11,6.57500027e-11,1.31199995e-10,1.65200006e-10,
	 2.61849986e-10,4.15049994e-10,5.22499988e-10,6.57999988e-10,
	 8.28499991e-10,1.31300004e-09,2.08199991e-09,3.30050010e-09,
	 4.15599999e-09,5.22999999e-09,6.59000010e-09,8.29500024e-09,
	 1.04500000e-08,1.31549998e-08,1.65650000e-08,2.08599999e-08,
	 2.62699995e-08,3.30850014e-08,4.16600017e-08,5.24500017e-08,
	 6.61000001e-08,8.31999998e-08,9.21999970e-08,1.04800002e-07,
	 1.31999997e-07,1.66250004e-07,2.09400000e-07,//2.14950006e-07,
	 2.23000001e-07,2.61400004e-07,3.30500001e-07,3.98200001e-07,
	 4.86399983e-07,5.97999986e-07,7.35500009e-07,//8.41528747e-07,
	 //4.21516539e-06,
	 8.43854053e-06,1.26672671e-05,1.69004320e-05,
	 2.11374665e-05,2.53779855e-05,2.96217149e-05,3.38684539e-05,
	 3.81180526e-05,4.23703981e-05,4.66254054e-05,5.08830106e-05,
	 5.51431670e-05,5.94058410e-05,6.36710102e-05,6.79386612e-05};
      double pr_arr[nlines]=
	{5.64869998e-32,5.64869986e-31,5.64869986e-30,5.64870017e-29,
	 6.76729990e-28,7.82989977e-27,9.50780029e-26,3.25500004e-24,
	 1.06260006e-22,5.44959997e-21,2.77850000e-20,1.35959994e-19,
	 6.43729996e-19,2.94519994e-18,1.29639997e-17,5.45580009e-17,
	 2.18730006e-16,2.94129994e-16,8.02570017e-16,2.14369999e-15,
	 5.62639983e-15,1.45639997e-14,3.73379984e-14,4.88699996e-14,
	 9.11069993e-14,1.69410005e-13,2.30929994e-13,2.81650002e-13,
	 3.83670011e-13,7.11400014e-13,1.31769996e-12,2.43959995e-12,
	 3.16660001e-12,4.30759985e-12,5.86130016e-12,7.96970042e-12,
	 1.08389998e-11,1.39989999e-11,1.90380003e-11,2.58830006e-11,
	 3.32719997e-11,4.52399992e-11,6.15209966e-11,8.36119993e-11,
	 1.13700001e-10,1.45250006e-10,1.61739996e-10,1.83999996e-10,
	 2.50169996e-10,3.25279997e-10,4.38359987e-10,//4.36519987e-10,
	 4.41269993e-10,4.67110017e-10,5.08829978e-10,5.49830015e-10,
	 6.05699990e-10,6.81199985e-10,7.82430010e-10,//4.70257815e-10,
	 //7.04778782e-09,
	 1.45139718e-08,2.62697827e-08,4.05674724e-08,
	 5.69532689e-08,7.52445638e-08,9.53657839e-08,1.17299621e-07,
	 1.41064470e-07,1.66702091e-07,1.94270264e-07,2.23838165e-07,
	 2.55483362e-07,2.89289800e-07,3.25346430e-07,3.63172009e-07};
      double nb_arr[nlines]=
	{4.00000001e-15,4.73000011e-15,4.75999990e-15,4.91000012e-15,
	 6.99000006e-15,9.89999996e-15,2.71999999e-14,1.27000000e-13,
	 6.93000019e-13,6.29500011e-12,1.58099991e-11,3.97200016e-11,
	 9.97599989e-11,2.50600013e-10,6.29399977e-10,1.58100000e-09,
	 3.97200006e-09,4.99999997e-09,9.97599958e-09,1.98999999e-08,
	 3.97199997e-08,7.92400030e-08,1.58099994e-07,1.98999999e-07,
	 3.15499989e-07,4.99999999e-07,6.29400006e-07,7.92399987e-07,
	 9.97599955e-07,1.58099999e-06,2.50600010e-06,3.97199983e-06,
	 4.99999987e-06,6.29399983e-06,7.92399987e-06,9.97600000e-06,
	 1.25600000e-05,1.58099992e-05,1.99000006e-05,2.50600006e-05,
	 3.15500001e-05,3.97199983e-05,4.99999987e-05,6.29400020e-05,
	 7.92400024e-05,9.97600000e-05,1.10499997e-04,1.25599996e-04,
	 1.58099996e-04,1.99000002e-04,2.50599987e-04,//2.57200008e-04,
	 2.66999996e-04,3.12599994e-04,3.95100011e-04,4.75899986e-04,
	 5.81200002e-04,7.14300026e-04,8.78599996e-04,//1.00000001e-03,
	 //5.00000001e-03,
	 1.00000000e-02,1.50000000e-02,2.00000000e-02,
	 2.50000000e-02,3.00000000e-02,3.50000000e-02,4.00000000e-02,
	 4.50000000e-02,5.00000000e-02,5.50000000e-02,6.00000000e-02,
	 6.50000000e-02,7.00000000e-02,7.50000000e-02,8.00000000e-02};
      for(size_t i=80;i<nlines;i+=2) {
	double line[3]={cu.convert("Msun/km^3","1/fm^4",ed_arr[i]),
			cu.convert("Msun/km^3","1/fm^4",pr_arr[i]),
			nb_arr[i]};
	new_eos.line_of_data(3,line);
      }
      for(size_t i=0;i<eos->get_nlines();i+=2) {
	if (eos->get("nb",i)>0.0801) {
	  double line[3]={eos->get("ed",i),eos->get("pr",i),eos->get("nb",i)};
	  new_eos.line_of_data(3,line);
	}
      }
    }
    p.set_eos_fm(new_eos.get_nlines(),new_eos["ed"],
		 new_eos["pr"],new_eos["nb"]);
  } else {
    p.set_eos_fm(eos->get_nlines(),(*eos)["ed"],(*eos)["pr"],(*eos)["nb"]);
  }
  nst.set_eos(p);
  
  // Compute the mass of the non-rotating configuration with the
  // same energy density
  nst.fix_cent_eden_non_rot(ed_cent);
  // Compare with with the answer from nstar_rot
  t.test_rel(nst.Mass/nst.MSUN,1.4,0.015,"correct mass");

  t.report();
  
  return 0;
}
