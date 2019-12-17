/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2019, Andrew W. Steiner
  
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
#include <o2scl/eos_had_rmf.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  nstar_rot nst;

#ifndef O2SCL_FAST_TEST

  if (true) {
    // Perform the RNS tests
    nst.constants_rns();
    nst.test1(t);
    nst.test2(t);
    nst.test3(t);
    nst.test4(t);
    nst.test5(t);
    nst.test6(t);
    nst.test7(t);
    nst.test8(t);
    nst.constants_o2scl();
  }

  if (true) {

    // Test running with SLy4
    eos_had_skyrme sk;
    o2scl_hdf::skyrme_load(sk,"../../data/o2scl/skdata/SLy4.o2",1);

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
    convert_units<double> &cu=o2scl_settings.get_convert_units();
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

    // Compute a configuration with a fixed ratio of radii
    nst.fix_cent_eden_axis_rat(ed_cent,0.7);
    t.test_rel(nst.r_ratio,0.7,1.0e-6,"correct ratio");

    // Create an output table
    table3d t;
    nst.output_table(t);

    // Now reinterpolate to construct a new table in cartesian
    // coordinates
    table3d t2;
    t.set_interp_type(itp_linear);
    // Equatorial radius in km
    double rad_eq=nst.R_e/1.0e5;
    uniform_grid_end<double> coord_grid(0,rad_eq*1.1,100);
    t2.set_xy("x",coord_grid,"z",coord_grid);
    t2.line_of_names("ed pr h vsq rho gamma omega alpha");
    for(size_t i=0;i<t2.get_nx();i++) {
      for(size_t j=0;j<t2.get_ny();j++) {
	double r=sqrt(pow(coord_grid[i],2.0)+pow(coord_grid[j],2.0));
	double theta=atan(-coord_grid[j]/coord_grid[i])+acos(-1)/2.0;
	if (i==0 && j==0) theta=0.0;
	t2.set(i,j,"ed",t.interp(r/(r+rad_eq),cos(theta),"ed"));
	t2.set(i,j,"pr",t.interp(r/(r+rad_eq),cos(theta),"pr"));
	t2.set(i,j,"h",t.interp(r/(r+rad_eq),cos(theta),"h"));
	t2.set(i,j,"vsq",t.interp(r/(r+rad_eq),cos(theta),"vsq"));
	t2.set(i,j,"rho",t.interp(r/(r+rad_eq),cos(theta),"rho"));
	t2.set(i,j,"gamma",t.interp(r/(r+rad_eq),cos(theta),"gamma"));
	t2.set(i,j,"omega",t.interp(r/(r+rad_eq),cos(theta),"omega"));
	t2.set(i,j,"alpha",t.interp(r/(r+rad_eq),cos(theta),"alpha"));
      }
    }

    // Write both tables to a file
    o2scl_hdf::hdf_file hf;
    hf.open_or_create("nstar_rot.o2");
    o2scl_hdf::hdf_output(hf,(const table3d &)t,"nstar_rot");
    o2scl_hdf::hdf_output(hf,(const table3d &)t2,"nstar_rot2");
    hf.close();
  }

  {
    eos_had_rmf rmf;
    o2scl_hdf::rmf_load(rmf,"../../data/o2scl/rmfdata/RAPR.o2",true);

    nstar_cold nco;
    nco.def_tov.verbose=0;
    nco.set_eos(rmf);

    // Compute the RMF EOS in beta-equilibrium
    nco.calc_eos();
    std::shared_ptr<table_units<> > eos=nco.get_eos_results();

    // Evalulate the mass-radius curve
    nco.calc_nstar();
    std::shared_ptr<table_units<> > mvsr=nco.get_tov_results();

    // Lookup the central energy density of a 1.4 Msun neutron star
    // in g/cm^3
    convert_units<double> &cu=o2scl_settings.get_convert_units();
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

    // Compute a configuration with a fixed ratio of radii
    nst.fix_cent_eden_axis_rat(ed_cent,0.7);
    t.test_rel(nst.r_ratio,0.7,1.0e-6,"correct ratio");

    // Compare the original RNS solvers with the new alt()
    // functions 
    nst.fix_cent_eden_with_kepler(ed_cent);
    t.test_rel(nst.Omega_K,nst.Omega,1.0e-4,"kepler 1");

    nst.fix_cent_eden_with_kepler_alt(ed_cent);
    t.test_rel(nst.Omega_K,nst.Omega,1.0e-8,"kepler 2");

    nst.fix_cent_eden_grav_mass(ed_cent,3.072e33/nst.MSUN);
    t.test_rel(3.072e33,nst.Mass,1.0e-4,"gmass 1");
    
    nst.fix_cent_eden_grav_mass_alt(ed_cent,3.072e33/nst.MSUN);
    t.test_rel(3.072e33,nst.Mass,1.0e-9,"gmass 2");
    
    nst.fix_cent_eden_bar_mass(ed_cent,3.373e33/nst.MSUN);
    t.test_rel(3.373e33,nst.Mass,1.0e-1,"bmass 1");
    
    nst.fix_cent_eden_bar_mass_alt(ed_cent,3.373e33/nst.MSUN);
    t.test_rel(3.373e33,nst.Mass,1.0e-1,"bmass 2");
    
    nst.fix_cent_eden_ang_vel(ed_cent,1.131e4);
    t.test_rel(1.131e4,nst.Omega,1.0e-4,"ang_vel 1");
    
    nst.fix_cent_eden_ang_vel_alt(ed_cent,1.131e4);
    t.test_rel(1.131e4,nst.Omega,1.0e-9,"ang_vel 2");
    
    double J0=nst.G*nst.MSUN*nst.MSUN/nst.C;
    nst.fix_cent_eden_ang_mom(ed_cent,1.059e49/J0);
    t.test_rel(1.059e49,nst.J,1.0e-4,"ang_mom 1");

    nst.fix_cent_eden_ang_mom_alt(ed_cent,1.059e49/J0);
    t.test_rel(1.059e49,nst.J,1.0e-8,"ang_mom 2");
  }
#endif

  t.report();

  return 0;
}
