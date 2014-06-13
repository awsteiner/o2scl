/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#include <o2scl/eos_tov.h>
#include <o2scl/test_mgr.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/tov_solve.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

class simple_apr {

public:
  
  fermion n, p, e;
  thermo ht, hb, l;
  eos_had_apr ap;
  fermion_zerot fzt;
  double barn;

  simple_apr() {
    n.init(939.0/hc_mev_fm,2.0);
    p.init(939.0/hc_mev_fm,2.0);
    e.init(0.511/hc_mev_fm,2.0);
    ap.select(1);
    ap.pion=0;
  }

  int nstarfun(size_t nv, const ubvector &x, ubvector &y) {
    
    n.n=x[0];
    p.n=x[1];
    e.mu=x[2];
    
    if (x[0]<0.0 || x[1]<0.0) return 1;
    
    ap.calc_e(n,p,hb);
    fzt.calc_mu_zerot(e);
    
    l.ed=e.ed;
    l.pr=e.pr;
    ht=l+hb;
    
    y[0]=n.n+p.n-barn;
    y[1]=p.n-e.n;
    y[2]=n.mu-p.mu-e.mu;
    
    return 0;
  }
  
};

/*
void test_crust(eos_tov_interp &te, convert_units &cu, double pr_low,
		double pr_high, bool new_units, test_mgr &t) {
  
  double ed_old;
  vector<double> ed_bench;
  double ed, nb;
  size_t i;

  cout << "High-density (user units: e,p,nb): " << endl;
  ed_old=0.0;
  for(double pr=pr_low;pr<pr_high;pr*=1.07) {
    te.get_eden_high(pr,ed,nb);
    t.test_gen(ed_old<ed,"Energy density increasing");
    ed_old=ed;
    cout << ed << " " << pr << " " << nb << endl;
  }
  cout << endl;
  cout << "Low-density (user units: e,p,nb): " << endl;
  ed_old=0.0;
  for(double pr=pr_low;pr<pr_high;pr*=1.07) {
    te.get_eden_low(pr,ed,nb);
    t.test_gen(ed_old<ed,"Energy density increasing");
    ed_old=ed;
    cout << ed << " " << pr << " " << nb << endl;
  }
  cout << endl;
  cout << "Final (user units: e,p,nb): " << endl;
  ed_old=0.0;
  for(double pr=pr_low;pr<pr_high;pr*=1.07) {
    int ip=te.get_eden_full(pr,ed,nb);
    ed_bench.push_back(ed);
    t.test_gen(ed_old<ed,"Energy density increasing");
    ed_old=ed;
    cout << ed << " " << pr << " " << nb << " " << ip << endl;
  }
  cout << endl;
  cout << "Final (converted to user units: e,p,nb): " << endl;
  ed_old=0.0;
  i=0;
  for(double pr=pr_low;pr<pr_high;pr*=1.07) {
    if (new_units) {
      te.get_eden(cu.convert("1/fm^4","Msun/km^3",pr),ed,nb);
      t.test_rel(cu.convert("Msun/km^3","1/fm^4",ed),
		 ed_bench[i],1.0e-12,"Two eden functions match");
      i++;
      t.test_gen(ed_old<ed,"Energy density increasing");
      ed_old=ed;
      cout << cu.convert("Msun/km^3","1/fm^4",ed) << " " 
	   << pr << " " << nb << endl;
    } else {
      te.get_eden(pr,ed,nb);
      t.test_rel(ed,ed_bench[i],1.0e-12,"Two eden functions match");
      i++;
      t.test_gen(ed_old<ed,"Energy density increasing");
      ed_old=ed;
      cout << ed << " " << pr << " " << nb << endl;
    }
  }
  ed_bench.clear();
  cout << endl;

  return;
}
*/

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // --------------------------------------------------------------
  // Calculate APR EOS, with naive phase transition and no muons

  simple_apr sa;

  ubvector x(3);
  x[0]=0.09663;
  x[1]=0.003365;
  x[2]=0.4636;

  table_units<> eos;
  eos.line_of_names("ed pr nb");
  eos.set_unit("ed","Msun/km^3");
  eos.set_unit("pr","Msun/km^3");
  eos.set_unit("nb","1/fm^3");

  mroot_hybrids<mm_funct11> gmh;
  mm_funct11 nf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&simple_apr::nstarfun),&sa,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  convert_units &cu=o2scl_settings.get_convert_units();

  gmh.err_nonconv=false;

  for(sa.barn=0.04;sa.barn<=2.0001;sa.barn+=0.01) {
    gmh.msolve(3,x,nf);
    double line[3]={cu.convert("1/fm^4","Msun/km^3",sa.ht.ed),
		    cu.convert("1/fm^4","Msun/km^3",sa.ht.pr),sa.barn};
    eos.line_of_data(3,line);
  }
  
  // Read APR EOS 
  eos_tov_interp te;
  te.verbose=2;
  te.default_low_dens_eos();
  te.read_table(eos,"ed","pr","nb");
  cout << endl;

  // Make sure baryon_column is set correctly
  t.test_gen(te.baryon_column==true,"baryon column");

  // Double check that get_names_units() correctly reports 0
  vector<string> auxp, auxu;
  size_t np;
  te.get_names_units(np,auxp,auxu);
  t.test_gen(np==0,"no aux parms");
  
  double pr_low=2.0e-7;
  double pr_high=6.0e-7;
  double ed, nb, ed_old;
  size_t i;

  vector<double> ed_bench;

  cout << "-------------------------------------------------------------- "
       << endl;
  cout << "Show interpolation results near transition" << endl;
  cout << endl;

  //test_crust(te,cu,pr_low,pr_high,false,t);
  cout << endl;

  cout << "-------------------------------------------------------------- "
       << endl;
  cout << "Switch to different units" << endl;
  cout << endl;

  cout << "Switching to units of inverse fm." << endl;
  eos.convert_to_unit("ed","1/fm^4");
  eos.convert_to_unit("pr","1/fm^4");
  te.read_table(eos,"ed","pr","nb");
  pr_low=cu.convert("Msun/km^3","1/fm^4",pr_low);
  pr_high=cu.convert("Msun/km^3","1/fm^4",pr_high);
  cout << endl;
  
  //test_crust(te,cu,pr_low,pr_high,true,t);

  cout << "-------------------------------------------------------------- "
       << endl;
  cout << "Now try with new transition density and width" << endl;
  cout << endl;

  /*
    double prt_low, prt, prt_high;
    te.get_transition(prt_low,prt,prt_high);
    cout << "Pressures near transition: " << endl;
    cout << prt_low << " " << prt << " " << prt_high << endl;
    te.set_transition(2.0e-3,1.2);
    cout << endl;
  */

  //test_crust(te,cu,pr_low,pr_high,true,t);

  cout << "-------------------------------------------------------------- "
       << endl;
  cout << "Now try with different transition method" << endl;
  cout << endl;

  te.transition_mode=eos_tov_interp::match_line;

  //test_crust(te,cu,pr_low,pr_high,true,t);

#ifdef O2SCL_NEVER_DEFINED

  cout << "GCP10: " << endl;
  te.gcp10_low_dens_eos("BSk19");
  te.get_transition(prt_low,prt,prt_high);
  cout << "Pressures near transition: " << endl;
  cout << prt_low << " " << prt << " " << prt_high << endl;
  te.gcp10_low_dens_eos("BSk20");
  te.get_transition(prt_low,prt,prt_high);
  cout << "Pressures near transition: " << endl;
  cout << prt_low << " " << prt << " " << prt_high << endl;
  te.gcp10_low_dens_eos("BSk21");
  te.get_transition(prt_low,prt,prt_high);
  cout << "Pressures near transition: " << endl;
  cout << prt_low << " " << prt << " " << prt_high << endl;

  tov_solve tsn;
  tsn.verbose=0;
  tsn.set_eos(te);
  tsn.fixed(1.4);
  cout << tsn.rad << endl;

  // ---------------------------------------------------------
  // Create a table comparing the crust EOSs
  
  size_t ngrid=200;
  uniform_grid_log_end<double> pr_grid(2.0e-14,4.0e-3,ngrid-1);
  cout << pr_grid.get_npoints() << endl;

  table_units<> crust_comp;
  crust_comp.line_of_names("pr ed_NVBPS ed_SHO ed_PNM_L40 ed_PNM_L100");
  crust_comp.line_of_names("ed_J35_L40 ed_J35_L100 ed_SLy4 ed_APR ed_Gs");
  crust_comp.line_of_names("ed_Rs nb_NVBPS nb_SHO nb_PNM_L40 nb_PNM_L100");
  crust_comp.line_of_names("nb_J35_L40 nb_J35_L100 nb_SLy4 nb_APR nb_Gs");
  crust_comp.line_of_names("nb_Rs");

  crust_comp.set_unit("pr","1/fm^4");

  crust_comp.set_unit("ed_NVBPS","1/fm^4");
  crust_comp.set_unit("ed_SHO","1/fm^4");
  crust_comp.set_unit("ed_PNM_L40","1/fm^4");
  crust_comp.set_unit("ed_PNM_L100","1/fm^4");
  crust_comp.set_unit("ed_J35_L40","1/fm^4");
  crust_comp.set_unit("ed_J35_L100","1/fm^4");
  crust_comp.set_unit("ed_APR","1/fm^4");
  crust_comp.set_unit("ed_SLy4","1/fm^4");
  crust_comp.set_unit("ed_Gs","1/fm^4");
  crust_comp.set_unit("ed_Rs","1/fm^4");

  crust_comp.set_unit("nb_NVBPS","1/fm^3");
  crust_comp.set_unit("nb_SHO","1/fm^3");
  crust_comp.set_unit("nb_PNM_L40","1/fm^3");
  crust_comp.set_unit("nb_PNM_L100","1/fm^3");
  crust_comp.set_unit("nb_J35_L40","1/fm^3");
  crust_comp.set_unit("nb_J35_L100","1/fm^3");
  crust_comp.set_unit("nb_APR","1/fm^3");
  crust_comp.set_unit("nb_SLy4","1/fm^3");
  crust_comp.set_unit("nb_Gs","1/fm^3");
  crust_comp.set_unit("nb_Rs","1/fm^3");

  // NVBPS
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_NVBPS",i,ed);
    crust_comp.set("nb_NVBPS",i,nb);
    cout << pr_grid[i] << " " << ed << " " << nb << endl;
  }
  
  // SHO
  te.sho11_low_dens_eos();
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_SHO",i,ed);
    crust_comp.set("nb_SHO",i,nb);
  }

  te.ngl13_low_dens_eos(40.0);
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_PNM_L40",i,ed);
    crust_comp.set("nb_PNM_L40",i,nb);
  }

  te.ngl13_low_dens_eos(100.0);
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_PNM_L100",i,ed);
    crust_comp.set("nb_PNM_L100",i,nb);
  }

  te.ngl13_low_dens_eos(40.0,"J35");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_J35_L40",i,ed);
    crust_comp.set("nb_J35_L40",i,nb);
  }

  te.ngl13_low_dens_eos(100.0,"J35");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_J35_L100",i,ed);
    crust_comp.set("nb_J35_L100",i,nb);
  }

  te.s12_low_dens_eos("SLy4");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_SLy4",i,ed);
    crust_comp.set("nb_SLy4",i,nb);
  }

  te.s12_low_dens_eos("APR");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_APR",i,ed);
    crust_comp.set("nb_APR",i,nb);
  }

  te.s12_low_dens_eos("Rs");
  for(size_t i=0;i<200;i++) {
    crust_comp.set("pr",i,pr_grid[i]);
    te.get_eden_low(pr_grid[i],ed,nb);
    crust_comp.set("ed_Rs",i,ed);
    crust_comp.set("nb_Rs",i,nb);
  }

  /*
    hdf_file hf;
    hf.open("crust_comp.o2");
    hdf_output(hf,crust_comp,"crust_comp");
    hf.close();
  */

#endif
  
  t.report();

  return 0;
}


