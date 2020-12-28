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
/** \file ex_eos_had_apr.cpp
    \brief File defining \ref ex_eos_had_apr class
*/
/* Example: ex_eos_had_apr.cpp
   -------------------------------------------------------------------
*/
#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/eos_had_apr.h>
#include <o2scl/table.h>
#include <o2scl/constants.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/tov_solve.h>
#include <o2scl/deriv_cern.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace o2scl_cgs;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/** \brief Compute the APR EOS with a Gibbs construction and the mass
    versus radius curve [Example class]
    
    In succession, calculates nuclear matter, neutron matter, and then
    neutron star matter with Maxwell and Gibbs constructions.
    
    We could use the more accurate masses in
    <tt>o2scl/constants.h</tt> here, but APR appears to have been
    designed to be used with neutron and protons masses equal
    to 939 MeV. 
*/
class ex_eos_had_apr {
  
protected:

  /// \name Fermions
  //@{
  /// Compute zero-temperature thermodynamics
  fermion_zerot fzt;
  /// Neutron for low-density phase
  fermion n;
  /// Proton for low-density phase
  fermion p;
  /// Neutron for high-density phase
  fermion n2;
  /// Proton for high-density phase
  fermion p2;
  /// Electron for low-density phase
  fermion e;
  /// Muon for low-density phase
  fermion mu;
  /// Electron for high-density phase
  fermion e2;
  /// Muon for high-density phase
  fermion mu2;
  //@}
  
  /// \name 'Thermo' objects
  //@{
  /// Baryon thermodynamics for low-density phase
  thermo hb;
  /// Leptonic thermodynamics for low-density phase
  thermo l;
  /// Baryon thermodynamics for high-density phase
  thermo hb2;
  /// Total thermodynamics
  thermo tot;
  /// Leptonic thermodynamics for high-density phase
  thermo l2;
  //@}

  /// \name Numerical methods
  //@{
  /// General solver
  mroot_hybrids<> solver;
  /// Solver for transition densities (lower tolerances)
  mroot_hybrids<> solver_trans_density;
  /// Derivative object
  deriv_cern<> cd;
  //@}

  /// Baryon density
  double nb;
  /// Volume fraction of low-density phase
  double chi;
  /// Baryon chemical potential
  double mub;
  /// Charge chemical potential
  double muq;
  /// Proton fraction for Fig. 7
  double f7x;
  /// Choice of model from APR
  int choice;
  /// \name Phase specification
  //@{
  int phase;
  static const int low_phase=1;
  static const int mixed_phase=2;
  static const int high_phase=3;
  //@}

  /// Base APR EOS
  eos_had_apr ap;
  
  /// Table for output
  table_units<> at;
  /// HDF file for output
  hdf_file hf;

  /// Function for the Maxwell construction in Fig. 7
  int maxwell_fig7(size_t nv, const ubvector &x, ubvector &y) {

    n.n=x[0];
    p.n=x[1];
    n2.n=x[2];
    p2.n=x[3];

    if (x[0]<0.0 || x[1]<0.0 || x[2]<0.0 || x[3]<0.0) return 1;

    // Low-density phase
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);

    e.mu=n.mu-p.mu;
    mu.mu=e.mu;
    fzt.calc_mu_zerot(e);
    fzt.calc_mu_zerot(mu);
    l.ed=e.ed+mu.ed;
    l.pr=e.pr+mu.pr;
    l.en=e.en+mu.en;

    // High-density phase
    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);

    e2.mu=n2.mu-p2.mu;
    mu2.mu=e2.mu;
    fzt.calc_mu_zerot(e2);
    fzt.calc_mu_zerot(mu2);
    l2.ed=e2.ed+mu2.ed;
    l2.pr=e2.pr+mu2.pr;
    l2.en=e2.en+mu2.en;

    // Charge neutrality for low-density phase
    y[0]=p.n-e.n-mu.n;
    // Equal neutron chemical potentials
    y[1]=n.mu-n2.mu;
    // Equal pressures
    y[2]=hb.pr-hb2.pr;
    // Charge neutrality for high-density phase
    y[3]=p2.n-e2.n-mu2.n;

    return 0;
  }

  /// Maxwell construction of the nuclear matter mixed phase 
  int mixedmaxwell(size_t nv, const ubvector &x, ubvector &y) {
    n.n=x[0];
    p.n=x[1];
    n2.n=x[2];
    p2.n=x[3];
    chi=x[4];
  
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);

    e.mu=n.mu-p.mu;
    mu.mu=e.mu;
    fzt.calc_mu_zerot(e);
    fzt.calc_mu_zerot(mu);
    l.ed=e.ed+mu.ed;
    l.pr=e.pr+mu.pr;
    l.en=e.en+mu.en;

    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);

    e2.mu=n2.mu-p2.mu;
    mu2.mu=e2.mu;
    fzt.calc_mu_zerot(e2);
    fzt.calc_mu_zerot(mu2);
    l2.ed=e2.ed+mu2.ed;
    l2.pr=e2.pr+mu2.pr;
    l2.en=e2.en+mu2.en;

    y[0]=nb-chi*(n.n+p.n)-(1.0-chi)*(n2.n+p2.n);
    y[1]=n.mu-n2.mu;
    y[2]=hb.pr-hb2.pr;
    y[3]=p.n-e.n-mu.n;
    y[4]=p2.n-e2.n-mu2.n;
  
    return 0;
  }

  /// Function to construct Fig. 7
  int fig7fun(size_t nv, const ubvector &x, ubvector &y) {

    nb=x[0];
    p.n=f7x*nb;
    n.n=nb-p.n;

    p2.n=f7x*nb;
    n2.n=nb-p2.n;
  
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);

    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);

    y[0]=hb.ed-hb2.ed;

    return 0;
  }

  /// Solve for neutron star matter (low-density phase)
  int nstar_low(size_t nv, const ubvector &x, ubvector &y) {
  
    n.n=x[0];
    p.n=nb-n.n;
  
    if (n.n<0.0 || p.n<0.0) {
      return 1;
    }

    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);

    e.mu=n.mu-p.mu;
    mu.mu=e.mu;
    fzt.calc_mu_zerot(e);
    fzt.calc_mu_zerot(mu);
    l.ed=e.ed+mu.ed;
    l.pr=e.pr+mu.pr;
    l.en=e.en+mu.en;
  
    y[0]=p.n-e.n-mu.n;

    tot.pr=hb.pr+l.pr;
    tot.ed=hb.ed+l.ed;

    return 0;
  }

  /// Solve for neutron star matter (high-density phase)
  int nstar_high(size_t nv, const ubvector &x, ubvector &y) {

    n2.n=x[0];
    p2.n=nb-n2.n;
    
    if (n2.n<0.0 || p2.n<0.0) return 1;

    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);

    e.mu=n2.mu-p2.mu;
    mu.mu=e.mu;
    fzt.calc_mu_zerot(e);
    fzt.calc_mu_zerot(mu);
    l.ed=e.ed+mu.ed;
    l.pr=e.pr+mu.pr;
    l.en=e.en+mu.en;
  
    y[0]=p2.n-e.n-mu.n;

    tot.pr=hb2.pr+l.pr;
    tot.ed=hb2.ed+l.ed;
  
    return 0;
  }

  /// Solve for neutron star matter (mixed phase)
  int nstar_mixed(size_t nv, const ubvector &x, ubvector &y) {

    n.n=x[0];
    p.n=x[1];
    e.mu=x[2];
    n2.n=x[3];
    p2.n=x[4];
    mu.mu=e.mu;

    if (phase==low_phase) chi=1.0;
    else if (phase==high_phase) chi=0.0;
    else chi=x[5];

    if (n.n<0.0 || n2.n<0.0) return 1;
    if (p.n<0.0 || p2.n<0.0) return 1;

    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);
  
    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);

    fzt.calc_mu_zerot(e);
    fzt.calc_mu_zerot(mu);
    l.ed=e.ed+mu.ed;
    l.pr=e.pr+mu.pr;
    l.en=e.en+mu.en;
  
    y[0]=nb-chi*(n.n+p.n)-(1.0-chi)*(n2.n+p2.n);
    y[1]=chi*p.n+(1.0-chi)*p2.n-e.n-mu.n;
    y[2]=n.mu-p.mu-e.mu;
    y[3]=p2.mu-p.mu;
    y[4]=n.mu-n2.mu;

    if (phase==mixed_phase) y[5]=hb.pr-hb2.pr;

    if (phase==low_phase) {
      tot.pr=hb.pr+l.pr;
      tot.ed=hb.ed+l.ed;
    } else if (phase==mixed_phase) {
      tot.pr=hb.pr+l.pr;
      tot.ed=hb.ed*chi+hb2.ed*(1.0-chi)+l.ed;
    } else {
      tot.pr=hb2.pr+l.pr;
      tot.ed=hb2.ed+l.ed;
    }

    return 0;
  }

  /// Write a line of data to the table
  void store_data() {

    std::vector<double> line;
    line.push_back(tot.ed);
    line.push_back(tot.pr);
    line.push_back(nb);
    line.push_back((chi*n.n+(1.0-chi)*n2.n));
    line.push_back((chi*p.n+(1.0-chi)*p2.n));
    line.push_back(n.n);
    line.push_back(p.n);
    line.push_back(n2.n);
    line.push_back(p2.n);
    line.push_back(chi);
    line.push_back(e.n);
    line.push_back(mu.n);
    line.push_back(n.mu);
    line.push_back(p.mu);
    line.push_back(e.mu);
    line.push_back(mu.mu);
    line.push_back(n.ms);
    line.push_back(p.ms);
    line.push_back(n2.ms);
    line.push_back(p2.ms);
    line.push_back(n.kf);
    line.push_back(p.kf);
    line.push_back(n2.kf);
    line.push_back(p2.kf);
    line.push_back(e.kf);
    line.push_back(mu.kf);
    if (line.size()!=at.get_ncolumns()) {
      O2SCL_ERR("Table misaligned.",exc_efailed);
    }
  
    at.line_of_data(line.size(),line);
    return;
  }

  /// Solve for nuclear matter (mixed phase)
  int nucmixed(size_t nv, const ubvector &x, ubvector &y) {

    n.n=x[0];
    n2.n=x[1];

    if (phase==low_phase) chi=1.0;
    else if (phase==high_phase) chi=0.0;
    else chi=x[2];

    if (n.n<0.0 || n2.n<0.0) return 1;

    p.n=n.n;
    p2.n=n2.n;

    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);
  
    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);
  
    y[0]=n.mu-n2.mu;
    y[1]=nb-chi*(n.n+p.n)-(1.0-chi)*(n2.n+p2.n);

    if (phase==mixed_phase) y[2]=hb.pr-hb2.pr;

    if (phase==low_phase) {
      tot.pr=hb.pr;
      tot.ed=hb.ed;
    } else if (phase==mixed_phase) {
      tot.pr=hb.pr;
      tot.ed=hb.ed*chi+hb2.ed*(1.0-chi);
    } else {
      tot.pr=hb2.pr;
      tot.ed=hb2.ed;
    }
 

    return 0;
  }

  /// Solve for neutron matter (mixed phase)
  int neutmixed(size_t nv, const ubvector &x, ubvector &y) {

    n.n=x[0];
    n2.n=x[1];

    if (phase==low_phase) chi=1.0;
    else if (phase==high_phase) chi=0.0;
    else chi=x[2];

    if (n.n<0.0 || n2.n<0.0) return 1;

    p.n=0.0;
    p2.n=0.0;

    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,hb);
  
    ap.pion=eos_had_apr::hdp;
    ap.calc_e(n2,p2,hb2);
  
    y[0]=n.mu-n2.mu;
    y[1]=nb-chi*(n.n+p.n)-(1.0-chi)*(n2.n+p2.n);

    if (phase==mixed_phase) y[2]=hb.pr-hb2.pr;

    if (phase==low_phase) {
      tot.pr=hb.pr;
      tot.ed=hb.ed;
    } else if (phase==mixed_phase) {
      tot.pr=hb.pr;
      tot.ed=hb.ed*chi+hb2.ed*(1.0-chi);
    } else {
      tot.pr=hb2.pr;
      tot.ed=hb2.ed;
    }
 
    return 0;
  }

  /// Solve for phase transition to nuclei
  int nucleimat(size_t nv, const ubvector &ex, ubvector &ey) {

    double nn1,np1,mun1,mup1;
    double nn2,np2,mun2,mup2;
    double u=0.0;
    thermo th1, th2;

    n.n=ex[0];
    p.n=ex[1];
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th1);
    nn1=n.n;
    np1=p.n;
    mun1=n.mu;
    mup1=p.mu;
  
    n.n=ex[2];
    p.n=0.0;
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th2);
    nn2=n.n;
    np2=p.n;
    mun2=n.mu;
    mup2=p.mu;
  
    e.mu=mun1-mup1;
    fzt.calc_mu_zerot(e);
    l.ed=e.ed;
    l.pr=e.pr;
    l.en=e.en;
    if (nv>3) u=ex[3];
    else u=1.0;
  
    ey[0]=(mun1-mun2)/1.0; 
    ey[1]=(th1.pr-th2.pr)/1.0; 
    ey[2]=u*np1-e.n;
    if (nv>3)ey[3]= barn-u*(np1+nn1)-(1.0-u)*nn2;
  
    tot.pr=th1.pr+e.pr;
    hb.ed=u*th1.ed+(1.0-u)*th2.ed;
    tot.ed=u*th1.ed+(1.0-u)*th2.ed+e.ed;

    n.n=ex[0];
    p.n=ex[1];
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th1);
  
    return 0;
  }

  /// Solve for phase transition to nuclei with a proton drip
  int nucleimat_pdrip(size_t nv, const ubvector &ex, ubvector &ey) {

    double nn1,np1,mun1,mup1;
    double nn2,np2,mun2,mup2;
    double u=0.0;
    thermo th1, th2;

    if (ex[0]<0.0 || ex[1]<0.0 || ex[2]<0.0 || ex[3]<0.0) {
      return 1;
    }

    n.n=ex[0];
    p.n=ex[1];
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th1);
    nn1=n.n;
    np1=p.n;
    mun1=n.mu;
    mup1=p.mu;
  
    n.n=ex[2];
    p.n=ex[3];
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th2);
    nn2=n.n;
    np2=p.n;
    mun2=n.mu;
    mup2=p.mu;
  
    e.mu=mun1-mup1;
    fzt.calc_mu_zerot(e);
    l.ed=e.ed;
    l.pr=e.pr;
    l.en=e.en;
    if (nv>4) u=ex[4];
    else u=1.0;
  
    ey[0]=(mun1-mun2)/1.0; 
    ey[1]=(th1.pr-th2.pr)/1.0; 
    ey[2]=e.n-u*np1-(1.0-u)*np2;
    ey[3]=(mup1-mup2)/1.0;
    if (nv>4) ey[4]=barn-u*(np1+nn1)-(1.0-u)*(nn2+np2);
  
    tot.pr=th1.pr+e.pr;
    hb.ed=u*th1.ed+(1.0-u)*th2.ed;
    tot.ed=u*th1.ed+(1.0-u)*th2.ed+e.ed;

    n.n=ex[0];
    p.n=ex[1];
    ap.pion=eos_had_apr::ldp;
    ap.calc_e(n,p,th1);
  
    return 0;
  }

public:

  ex_eos_had_apr() {

    n.init(939.0/hc_mev_fm,2.0);
    p.init(939.0/hc_mev_fm,2.0);
    n2.init(939.0/hc_mev_fm,2.0);
    p2.init(939.0/hc_mev_fm,2.0);

    // Ensure that this works without GNU units
    o2scl_settings.get_convert_units().use_gnu_units=false;
    e.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_electron),2.0);
    mu.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_muon),2.0);
    e2.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_electron),2.0);
    mu2.init(o2scl_settings.get_convert_units().convert
	     ("kg","1/fm",o2scl_mks::mass_muon),2.0);

    n.non_interacting=false;
    p.non_interacting=false;
    n2.non_interacting=false;
    p2.non_interacting=false;

    at.inc_maxlines(2000);
  }

  /// Main driver, computing the APR EOS and the associated M vs. R curve
  void run() {

    test_mgr t;
    t.set_output_level(1);

    // Output directory
    std::string prefix="ex_eos_had_apr_";
    cout << "Set output prefix to '" << prefix << "' ." << endl;
    cout << endl;

    // Density grid
    double nbstart, nb_end, dnb;

    // Density at which to start looking for a mixed phase
    double nb_switch;

    // Temporary filename
    string filenm;

    // Error code value
    int ret;

    // Output file stream
    ofstream fout;

    // Variables and function values for solvers
    ubvector x(7), y(7), xx(3), yy(3);
  
    // Function objects
    mm_funct f_nucmixed=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nucmixed),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_neutmixed=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::neutmixed),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_nstar_mixed=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nstar_mixed),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_nstar_low=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nstar_low),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_nstar_high=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nstar_high),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_fig7fun=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::fig7fun),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_maxwell_fig7=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::maxwell_fig7),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct f_mixedmaxwell=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::mixedmaxwell),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

    // Init density grid
    nbstart=0.005;
    dnb=0.002;
    nb_end=2.0;

    // Select APR EOS
    choice=1;
    ap.select(choice);

    // Init solver tolerances
    solver.tol_abs=1.0e-10;
    solver.tol_rel=1.0e-12;
    solver_trans_density.tol_abs=1.0e-10;
    solver_trans_density.tol_rel=1.0e-12;

    // Initialize table
    at.line_of_names(((string)("ed pr nb "))+
		     "nn np nn1 np1 nn2 np2 chi ne nmu "+
		     "mun mup mue mumu "+
		     "msn msp msn2 msp2 kfn kfp kfn2 kfp2 kfe kfmu");
    at.set_unit("ed","1/fm^4");
    at.set_unit("pr","1/fm^4");
    at.set_unit("nb","1/fm^3");
    at.set_unit("nn","1/fm^3");
    at.set_unit("np","1/fm^3");
    at.set_unit("nn1","1/fm^3");
    at.set_unit("np1","1/fm^3");
    at.set_unit("nn2","1/fm^3");
    at.set_unit("np2","1/fm^3");
    at.set_unit("ne","1/fm^3");
    at.set_unit("nmu","1/fm^3");
    at.set_unit("mun","1/fm");
    at.set_unit("mup","1/fm");
    at.set_unit("mue","1/fm");
    at.set_unit("mumu","1/fm");
    at.set_unit("msn","1/fm");
    at.set_unit("msp","1/fm");
    at.set_unit("msn2","1/fm");
    at.set_unit("msp2","1/fm");
    at.set_unit("kfn","1/fm");
    at.set_unit("kfp","1/fm");
    at.set_unit("kfn2","1/fm");
    at.set_unit("kfp2","1/fm");
    at.set_unit("kfe","1/fm");
    at.set_unit("kfmu","1/fm");
    
    //--------------------------------------------------------------------
    // Saturation properties

    ap.set_thermo(hb);
    ap.set_mroot(solver);
    ap.set_n_and_p(n,p);
    ap.saturation();

    cout << "--------- Saturation properties --------------------\n" << endl;
    cout << "n_0: " << ap.n0 << " fm^{-3}\nE_B: " << ap.eoa*hc_mev_fm 
	 << " MeV\nK: " << ap.comp*hc_mev_fm << " MeV\nS: " 
	 << ap.esym*hc_mev_fm << " MeV\nM^{*}/M: " << ap.msom << endl;
    t.test_rel(ap.n0,0.1598371,1.0e-5,"n0");
    t.test_rel(ap.eoa*hc_mev_fm,-16.00066,1.0e-5,"eoa");
    t.test_rel(ap.esym*hc_mev_fm,32.5887,1.0e-5,"esym");
    cout << endl;
    
    //--------------------------------------------------------------------
    // Init lepton fields to zero to start

    e.mu=0.0;
    e.n=0.0;
    e.kf=0.0;
    mu.mu=0.0;
    mu.n=0.0;
    mu.kf=0.0;

    //--------------------------------------------------------------------
    // Nuclear matter

    cout << "--------- Nuclear matter ---------------------------\n" << endl;

    // Begin with just the low density phase

    nb_switch=0.19;

    ap.pion=eos_had_apr::ldp;
    chi=1.0;
    at.clear_data();
    for(nb=nbstart;nb<nb_switch;nb+=dnb) {
    
      n.n=nb/2.0;
      p.n=n.n;
      
      ret=ap.calc_e(n,p,hb);
      if (ret!=0) {
	O2SCL_ERR("Failed to compute nuclear matter.",exc_efailed);
      }
      tot.pr=hb.pr;
      tot.ed=hb.ed;
      
      store_data();
    }
  
    phase=low_phase;
    x[0]=nb_switch/2.0;
    x[1]=x[0];

    // Now start searching for the mixed phase
  
    for(nb=nb_switch;nb<=nb_end;nb+=dnb) {
    
      if (phase!=mixed_phase) ret=solver.msolve(2,x,f_nucmixed);
      else ret=solver.msolve(3,x,f_nucmixed);
      if (ret!=0) {
	cout << nb << endl;
	O2SCL_ERR("Solving nuclear matter failed.",
		  exc_efailed);
      }
    
      if (hb.pr<hb2.pr && phase==low_phase) {
	cout << "Mixed phase begins near nb=" << nb << " fm^{-3}." << endl;
	phase=mixed_phase;
	// Pick a value of chi close to, but less than, one
	x[2]=0.90;
	nb-=dnb;
      } else if (phase==mixed_phase && x[2]<0.0) {
	cout << "Mixed phase ends near nb=" << nb << " fm^{-3}." << endl;
	phase=high_phase;
	nb-=dnb;
      } else {
	store_data();
      }
    }
    
    string fn1=((string)prefix)+"nuc.o2";
    hf.open_or_create(fn1);
    hdf_output(hf,at,"nuc");
    hf.close();
    cout << "Generated file '" << fn1 << "'." << endl;

    cout << endl;
    
    //--------------------------------------------------------------------
    // Neutron matter

    cout << "--------- Neutron matter ---------------------------\n" << endl;

    ap.pion=eos_had_apr::ldp;
    chi=1.0;
    at.clear_data();

    nb_switch=0.16;

    // Begin with just the low density phase

    for(nb=nbstart;nb<nb_switch;nb+=dnb) {
    
      n.n=nb;
      p.n=0.0;

      ret=ap.calc_e(n,p,hb);
      if (ret!=0) {
	O2SCL_ERR("Failed to compute neutron matter.",exc_efailed);
      }

      tot.pr=hb.pr;
      tot.ed=hb.ed;
    
      store_data();
    }

    phase=low_phase;
    x[0]=nb_switch;
    x[1]=nb_switch;
  
    // Now start searching for the mixed phase
  
    for(nb=nb_switch;nb<=nb_end;nb+=dnb) {

      if (phase!=mixed_phase) {
	ret=solver.msolve(2,x,f_neutmixed);
      } else {
	ret=solver.msolve(3,x,f_neutmixed);
      }
      if (ret!=0) {
	cout << nb << endl;
	O2SCL_ERR("Solving neutron matter failed.",
		  exc_efailed);
      }

      if (hb.pr<hb2.pr && phase==low_phase) {
	cout << "Mixed phase begins near nb=" << nb << " fm^{-3}." << endl;
	phase=mixed_phase;
	// Pick a value of chi close to, but less than, one
	x[2]=0.90;
	nb-=dnb;
      } else if (phase==mixed_phase && x[2]<0.0) {
	cout << "Mixed phase ends near nb=" << nb << " fm^{-3}." << endl;
	phase=high_phase;
	nb-=dnb;
      } else {
	store_data();
      }
    }

    string fn2=((string)prefix)+"neut.o2";
    hf.open_or_create(fn2);
    hdf_output(hf,at,"neut");
    hf.close();
    cout << "Generated file '" << fn2 << "'." << endl;

    cout << endl;

    //--------------------------------------------------------------------
    // Neutron star matter - Maxwell construction

    cout << "--------- Neutron star matter ----------------------\n" << endl;

    // Removing this for now, as it caused negative densities
    if (false) {
      
      cout << "Maxwell construction." << endl;
      
      // Verify Figure 7

      filenm=prefix;
      filenm+="fig7.1.txt";
      fout.open(filenm.c_str());
      fout.setf(ios::scientific);
      x[0]=0.32;
      for(f7x=0.5;f7x>=-0.001;f7x-=0.025) {
	ret=solver.msolve(1,x,f_fig7fun);
	if (ret!=0) {
	  O2SCL_ERR("Failed to perform maxwell.",exc_efailed);
	}
	fout << f7x << " " << x[0] << endl;
      }
      fout.close();
      cout << "Generated file '" << filenm << "'." << endl;

      // Compute the beginning and ending densities
      // for the Maxwell construction

      x[0]=0.94*0.2;
      x[1]=0.06*0.2;
      x[2]=0.95*0.24;
      x[3]=0.05*0.24;
      solver.msolve(4,x,f_maxwell_fig7);

      double nb_low=x[0]+x[1], nb_high=x[2]+x[3];
      cout << "Mixed phase begins at nb=" << nb_low << " fm^{-3}." << endl;
      cout << "Mixed phase ends at nb=" << nb_high << " fm^{-3}." << endl;
      x[0]=nbstart/1.1;

      filenm=prefix; 
      filenm+="fig7.2.txt";
      fout.open(filenm.c_str());
    
      /// Compute matter at densities below the maxwell construction

      for(nb=0.02;nb<nb_low*1.00001;nb+=(nb_low-nbstart)/20.0) {
	ret=solver.msolve(1,x,f_nstar_low);
	if (ret!=0) {
	  cout << nb << endl;
	  f_nstar_low(1,x,y);
	  cout << x[0] << " " << y[0] << endl;
	  O2SCL_ERR("Solving Maxwell construction failed.",
		    exc_efailed);
	}
	fout << p.n/nb << " " << nb << endl;
      }

      // Compute matter at densities inside the Maxwell construction

      x[0]=(1.0-0.07)*nb_low;
      x[1]=0.07*nb_low;
      x[2]=(1.0-0.06)*nb_high;
      x[3]=0.06*nb_high;
      x[4]=1.0;
      dnb=(nb_high-nb_low)/20.0;
      for(nb=nb_low+dnb;nb<=nb_high*1.00001;nb+=dnb) {
	ret=solver.msolve(5,x,f_mixedmaxwell);
	if (ret!=0) {
	  cout << nb << endl;
	  O2SCL_ERR("Solving Maxwell construction (part 2) failed.",
		    exc_efailed);
	}
	fout << (chi*p.n+(1.0-chi)*p2.n)/nb << " " << nb << endl;
      }
    
      // Compute matter at densities above the Maxwell construction
    
      x[0]=0.23;
      for(nb=nb_high;nb<nb_end;nb+=(nb_end-nb_high)/40.0) {
	ret=solver.msolve(1,x,f_nstar_high);
	if (ret!=0) {
	  cout << nb << endl;
	  O2SCL_ERR("Solving Maxwell construction (part 3) failed.",
		    exc_efailed);
	}
	fout << p2.n/nb << " " << nb << endl;
      }
      fout.close();
      cout << "Generated file '" << filenm << "'." << endl;
    }
    
    //--------------------------------------------------------------------
    // Neutron star matter - Gibbs construction

    cout << "\nGibbs construction." << endl;
  
    dnb=0.002;
    ap.pion=eos_had_apr::ldp;
    chi=1.0;
    at.clear_data();

    x[0]=nbstart/1.1;
    for(nb=nbstart;nb<=0.1701;nb+=dnb) {
      ret=solver.msolve(1,x,f_nstar_low);
      if (ret!=0) {
	cout << nb << endl;
	O2SCL_ERR("Solving Gibbs construction failed.",
		  exc_efailed);
      }
      store_data();
    }
    nb-=dnb;
  
    phase=low_phase;
    x[0]=n.n;
    x[1]=p.n;
    x[2]=n.mu-p.mu;
    x[3]=n.n;
    x[4]=p.n;
  
    for(nb=0.17;nb<=nb_end;nb+=dnb) {
    
      if (phase!=mixed_phase) {
	ret=solver.msolve(5,x,f_nstar_mixed);
	nstar_mixed(5,x,y);
      } else {
	ret=solver.msolve(6,x,f_nstar_mixed);
      }
      if (ret!=0) {
	O2SCL_ERR("Gibbs construction (part 2) failed.",exc_esanity);
      }
    
      if (hb.pr<hb2.pr && phase==low_phase) {
	cout << "Mixed phase begins at nb=" << nb << " fm^{-3}." << endl;
	phase=mixed_phase;
	x[5]=0.90;
	nb-=dnb;
      } else if (phase==mixed_phase && x[5]<0.0) {
	cout << "Mixed phase ends at nb=" << nb << " fm^{-3}." << endl;
	phase=high_phase;
	nb-=dnb;
      } else {
	store_data();
      }
    }
  
    //--------------------------------------------------------------------
    // Estimate transition density 
    
    cout << "\nEstimate of transition density." << endl;
  
    ubvector newx(12), newy(12);
    mm_funct nuclei_f=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nucleimat),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
    mm_funct nucleip_f=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&ex_eos_had_apr::nucleimat_pdrip),this,
       std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

    solver_trans_density.tol_abs/=1.0e4;
    solver_trans_density.tol_rel/=1.0e4;
  
    newx[0]=0.1;
    newx[1]=0.001;
    newx[2]=0.1;
    ret=solver_trans_density.msolve(3,newx,nuclei_f);
    if (ret!=0) {
      nucleimat(3,newx,newy);
      if (fabs(newy[0])>1.0e-8 || fabs(newy[1])>1.0e-8 ||
	  fabs(newy[2])>1.0e-8) {
	cout << endl;
	cout << "Problem in transition density (1)." << endl;
	cout << err_hnd->get_str() << endl;
	cout << newx[0] << " " << newy[0] << endl;
	cout << newx[1] << " " << newy[1] << endl;
	cout << newx[2] << " " << newy[2] << endl;
	cout << endl;
      }
    }
    cout << "Without proton drip: density: " << newx[0]+newx[1] 
	 << " fm^{-3},\n  pressure: " 
	 << at.interp("nb",newx[0]+newx[1],"pr") << " fm^{-4}." << endl;

    solver_trans_density.tol_abs=5.0e-8;
    solver_trans_density.tol_rel=5.0e-8;
  
    newx[3]=0.001;
    ret=solver_trans_density.msolve(4,newx,nucleip_f);
    if (ret!=0) {
      nucleimat_pdrip(4,newx,newy);
      if (fabs(newy[0])>1.0e-8 || fabs(newy[1])>1.0e-8 ||
	  fabs(newy[2])>1.0e-8 || fabs(newy[3])>1.0e-8) {
	cout << endl;
	cout << "Problem in transition density (2)." << endl;
	cout << err_hnd->get_str() << endl;
	cout << newx[0] << " " << newy[0] << endl;
	cout << newx[1] << " " << newy[1] << endl;
	cout << newx[2] << " " << newy[2] << endl;
	cout << newx[3] << " " << newy[3] << endl;
	cout << endl;
      }
    }
    nucleimat_pdrip(4,newx,newy);
    cout << newx[0] << " " << newy[0] << endl;
    cout << newx[1] << " " << newy[1] << endl;
    cout << newx[2] << " " << newy[2] << endl;
    cout << newx[3] << " " << newy[3] << endl;
    cout << "With proton drip: density: " << newx[0]+newx[1] 
	 << " fm^{-3},\n  pressure: " 
	 << at.interp("nb",newx[0]+newx[1],"pr") << " fm^{-4}." << endl;
    cout << endl;
  
    solver_trans_density.tol_abs=1.0e-16;
    solver_trans_density.tol_rel=1.0e-16;
  
    // Output to file

    string fn3=((string)prefix)+"nstar.o2";
    hf.open_or_create(fn3);
    hdf_output(hf,at,"nstar");
    hf.close();
    cout << "Generated file '" << fn3 << "'." << endl;

    cout << endl;

    //--------------------------------------------------------------------
    // Integrate TOV equations

    cout << "--------- TOV solver, M vs. R. ---------------------\n" << endl;

    tov_solve atov;
    eos_tov_interp teos;
    atov.verbose=0;
    teos.verbose=0;
    atov.set_units("1/fm^4","1/fm^4","1/fm^3");
    teos.default_low_dens_eos();
    teos.read_table(at,"ed","pr","nb");
    atov.set_eos(teos);
    atov.calc_gpot=true;

    // Get results
    
    atov.mvsr();
    std::shared_ptr<table_units<> > tov_tmp=atov.get_results();
    cout << "Maximum mass is " << tov_tmp->max("gm") 
	 << " solar masses." << endl;
    t.test_rel(tov_tmp->max("gm"),2.191338,5.0e-4,"max mass.");
    
    // Output to file

    string fn4=((string)prefix)+"mvsr.o2";
    hf.open_or_create(fn4);
    hdf_output(hf,*tov_tmp,"mvsr");
    hf.close();
    cout << "Generated file '" << fn4 << "'." << endl;
    cout << endl;

    cout << "--------- TOV solver, 1.4 Msun ---------------------\n" << endl;

    atov.fixed(1.4);
    std::shared_ptr<table_units<> > tov_tmp2=atov.get_results();
    cout << "Aprpoximate radius of a 1.4 solar mass neutron star is " 
	 << tov_tmp2->get("r",tov_tmp2->lookup("gm",1.4)) << " km." << endl;
    t.test_rel(tov_tmp2->get("r",tov_tmp2->lookup("gm",1.4)),11.46630,
	       5.0e-3,"r14");

    string fn5=((string)prefix)+"m14.o2";
    hf.open_or_create(fn5);
    hdf_output(hf,*tov_tmp2,"m14");
    hf.close();
    cout << "Generated file '" << fn5 << "'." << endl;
    cout << endl;

    t.report();

    return;
  }

};

int main(void) {

  cout.setf(ios::scientific);

  ex_eos_had_apr ac;
  ac.run();

  return 0;
}
// End of example
