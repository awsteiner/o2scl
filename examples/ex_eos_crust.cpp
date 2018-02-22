/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

/* Example: ex_eos_crust.cpp
   -------------------------------------------------------------------
   Compute the Baym-Pethick-Sutherland equation of state
*/

#include <fstream>

#include <o2scl/test_mgr.h>
#include <o2scl/lib_settings.h>
#include <o2scl/table_units.h>
#include <o2scl/eos_crust.h>
#include <o2scl/eos_tov.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

// A simple function to load the BPS data without
// HDF support, in case O2scl was compiled without it.
int bps_load(table_units<> &bps);

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  eos_crust be;
  // Energy, pressure, etc.
  thermo th;
  // Proton and mass numbers
  int Z, A;

  // The baryon number density
  double barn;
  // Energy density from BPS table in units of g/cm^3 
  double ed_bps;
  // Pressure from BPS table in units of dyn/cm^2
  double pr_bps;

  // Read the original BPS results 
  table_units<> bps;
  bps_load(bps);
  
  // Ensure that this works without GNU units
  convert_units &cu=o2scl_settings.get_convert_units();
  cu.use_gnu_units=false;

  // Cubic spline interpolation doesn't do very well here, so we use
  // linear interpolation to interpolate the BPS table
  bps.set_interp_type(itp_linear);

  // Tabulate the o2scl results using eos_crust::calc_pressure() and 
  // compare them to the original table. The energy density
  // and pressure are output in MeV/fm^3, the baryon number density
  // is given in fm^{-3}. The last two columns are the fractional
  // deviations from the original BPS table. 

  cout << "Default mass formula, using calc_pressure()." << endl;
  cout << "eden         pres         barn         A    Z   ";
  cout << "ed_err        pr_err" << endl;

  for(double pr=8.9e-4/hc_mev_fm;pr>=2.12e-10/hc_mev_fm;pr/=2.0) {

    th.pr=pr;
    be.calc_pressure(th,barn,Z,A);

    // Output energy and pressure in MeV/fm^3
    cout << th.ed*hc_mev_fm << " " << th.pr*hc_mev_fm
	 << " " << barn << " ";
    cout.width(3);
    cout << A << " ";
    cout.width(3);
    cout << Z << " ";
    
    // Convert the BPS results from g/cm^3 and dyn/cm^2 to MeV/fm^3
    ed_bps=cu.convert("g/cm^3","MeV/fm^3",
		      bps.interp("nb",barn*1.0e39,"rho"));
      
    pr_bps=cu.convert("dyne/cm^2","MeV/fm^3",
		      bps.interp("nb",barn*1.0e39,"P"));

    cout.setf(ios::showpos);
    cout << (ed_bps-th.ed*hc_mev_fm)/ed_bps << " " 
	 << (pr_bps-th.pr*hc_mev_fm)/pr_bps << endl;
    cout.unsetf(ios::showpos);

    t.test_rel((ed_bps-th.ed*hc_mev_fm)/ed_bps,0.0,2.0e-3,"Energy");
    t.test_rel((pr_bps-th.pr*hc_mev_fm)/pr_bps,0.0,1.5e-1,"Pressure");
  }
  cout << endl;

  // Tabulate the o2scl results using eos_crust::calc_density() and 
  // compare them to the original table

  cout << "Default mass formula, using calc_density()." << endl;
  cout << "eden         pres         barn         A    Z   ";
  cout << "ed_err        pr_err" << endl;

  for(barn=4.19e-4;barn>=6.39e-9;barn/=2.0) {

    be.calc_density(barn,th,Z,A);
    
    ed_bps=cu.convert("g/cm^3","MeV/fm^3",
		      bps.interp("nb",barn*1.0e39,"rho"));
      
    pr_bps=cu.convert("dyne/cm^2","MeV/fm^3",
		      bps.interp("nb",barn*1.0e39,"P"));
    
    cout << th.ed*hc_mev_fm << " " << th.pr*hc_mev_fm
	 << " " << barn << " ";
    cout.width(3);
    cout << A << " ";
    cout.width(3);
    cout << Z << " ";
    cout.setf(ios::showpos);
    cout << (ed_bps-th.ed*hc_mev_fm)/ed_bps << " " 
	 << (pr_bps-th.pr*hc_mev_fm)/pr_bps << endl;
    cout.unsetf(ios::showpos);
    
    t.test_rel((ed_bps-th.ed*hc_mev_fm)/ed_bps,0.0,2.0e-3,"Energy");
    t.test_rel((pr_bps-th.pr*hc_mev_fm)/pr_bps,0.0,1.5e-1,"Pressure");
  }
  cout << endl;

  // ---------------------------------------------------------
  // Create a table comparing the crust EOSs from eos_tov_interp

  // Read the APR EOS table
  table_units<> tab;
  hdf_file hf;
  string name;
  hf.open("ex_eos_had_apr_nstar.o2");
  hdf_input(hf,tab,name);
  hf.close();
  
  // Send the EOS to the eos_tov_interp object
  eos_tov_interp te;
  te.read_table(tab,"ed","pr","nb");

  // Define a grid of pressures over which to evaluate the EOS
  size_t ngrid=200;
  uniform_grid_log_end<double> pr_grid(2.0e-14,4.0e-3,ngrid-1);

  // Create a table to store the data
  table_units<> crust_comp;
  crust_comp.line_of_names("pr ed_NVBPS ed_SHO ed_PNM_L40 ed_PNM_L100");
  crust_comp.line_of_names("ed_J35_L40 ed_J35_L100 ed_SLy4 ed_APR ed_Rs");
  crust_comp.line_of_names("nb_NVBPS nb_SHO nb_PNM_L40 nb_PNM_L100");
  crust_comp.line_of_names("nb_J35_L40 nb_J35_L100 nb_SLy4 nb_APR nb_Rs");

  crust_comp.set_unit("pr","1/fm^4");

  crust_comp.set_unit("ed_NVBPS","1/fm^4");
  crust_comp.set_unit("ed_SHO","1/fm^4");
  crust_comp.set_unit("ed_PNM_L40","1/fm^4");
  crust_comp.set_unit("ed_PNM_L100","1/fm^4");
  crust_comp.set_unit("ed_J35_L40","1/fm^4");
  crust_comp.set_unit("ed_J35_L100","1/fm^4");
  crust_comp.set_unit("ed_APR","1/fm^4");
  crust_comp.set_unit("ed_SLy4","1/fm^4");
  crust_comp.set_unit("ed_Rs","1/fm^4");

  crust_comp.set_unit("nb_NVBPS","1/fm^3");
  crust_comp.set_unit("nb_SHO","1/fm^3");
  crust_comp.set_unit("nb_PNM_L40","1/fm^3");
  crust_comp.set_unit("nb_PNM_L100","1/fm^3");
  crust_comp.set_unit("nb_J35_L40","1/fm^3");
  crust_comp.set_unit("nb_J35_L100","1/fm^3");
  crust_comp.set_unit("nb_APR","1/fm^3");
  crust_comp.set_unit("nb_SLy4","1/fm^3");
  crust_comp.set_unit("nb_Rs","1/fm^3");

  // Energy density and baryon density
  double ed, nb;

  cout << "NVBPS crust." << endl;
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_NVBPS",i,ed);
    crust_comp.set("nb_NVBPS",i,nb);
    if (i==199) {
      cout << pr_grid[i] << " " << ed << " " << nb << endl;
    }
  }
  
  cout << "SHO crust." << endl;
  te.sho11_low_dens_eos();
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_SHO",i,ed);
    crust_comp.set("nb_SHO",i,nb);
    if (i==199) {
      cout << pr_grid[i] << " " << ed << " " << nb << endl;
    }
  }

  cout << "NGL13, L=40 crust." << endl;
  te.ngl13_low_dens_eos(40.0);
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_PNM_L40",i,ed);
    crust_comp.set("nb_PNM_L40",i,nb);
  }

  cout << "NGL13, L=100 crust." << endl;
  te.ngl13_low_dens_eos(100.0);
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_PNM_L100",i,ed);
    crust_comp.set("nb_PNM_L100",i,nb);
  }

  cout << "NGL13, L=40 crust." << endl;
  te.ngl13_low_dens_eos(40.0,"J35");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_J35_L40",i,ed);
    crust_comp.set("nb_J35_L40",i,nb);
  }

  cout << "NGL13, L=100 crust." << endl;
  te.ngl13_low_dens_eos(100.0,"J35");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_J35_L100",i,ed);
    crust_comp.set("nb_J35_L100",i,nb);
  }

  cout << "S12, SLy4 crust." << endl;
  te.s12_low_dens_eos("SLy4");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_SLy4",i,ed);
    crust_comp.set("nb_SLy4",i,nb);
  }

  cout << "S12, APR crust." << endl;
  te.s12_low_dens_eos("APR");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_APR",i,ed);
    crust_comp.set("nb_APR",i,nb);
  }

  cout << "S12, Rs crust." << endl;
  te.s12_low_dens_eos("Rs");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.ed_nb_from_pr(cu.convert("1/fm^4","Msun/km^3",pr_grid[i]),ed,nb);
    ed=cu.convert("Msun/km^3","1/fm^4",ed);
    crust_comp.set("ed_Rs",i,ed);
    crust_comp.set("nb_Rs",i,nb);
  }
  cout << endl;

  hf.open_or_create("ex_crust_comp.o2");
  hdf_output(hf,crust_comp,"crust_comp");
  hf.close();

  t.report();

  return 0;
}
// End of example

int bps_load(table_units<> &bps) {
  bps.line_of_names("rho P nb");
  ifstream fin;
  fin.open("ex_bps_input.txt");
  string tmp;
  fin >> tmp >> tmp >> tmp;
  for(size_t i=0;i<43;i++) {
    double line[3];
    fin >> line[0] >> line[1] >> line[2];
    bps.line_of_data(3,line);
  }
  fin.close();
  return 0;
}
