/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2020, Andrew W. Steiner
  
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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/convert_units.h>
#include <o2scl/fermion.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/anneal_gsl.h>
#include <o2scl/mm_funct.h>
#include <o2scl/tov_solve.h>
#include <o2scl/test_mgr.h>
#include <o2scl/eos_quark_bag.h>
#include <o2scl/eos_quark_njl.h>
#include <o2scl/eos_had_rmf_hyp.h>
#include <o2scl/cli.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

/** \brief Compute the equation of state of matter with a Gibbs phase
    transition between hadronic and quark phase

    Following \ref Spinella16, the surface and Coulomb energies can be
    expressed in terms of \f$ \chi \f$ which is the volume fraction of
    matter in the quark phase. The quantity \f$ x \equiv
    \mathrm{min}(\chi,1-\chi) \f$ is the volume fraction in the rare
    phase. The surface energy density is written as
    \f[
    \varepsilon_{\mathrm{surf}} = \frac{d x \sigma}{r}
    \f]
    where \f$ r \f$ is the radius of the structure in the
    rare phase, \f$ d \f$ is the dimensionality, and \f$ \sigma \f$
    is the surface tension. 
    The Coulomb energy is written as
    \f[
    \varepsilon_{\mathrm{Coul}} = 2 \pi e^2 
    \left[ q_H - q_Q \right]^2 r^2 x f_{d}(x) 
    \f]
    where 
    \f[
    f_d(x) = \frac{1}{d+2} \left[ \frac{1}{d-2} 
    \left( 2-d x^{1-2/d}\right)+x\right] \, ,
    \f]
    \f$ q_H \f$ is the electric charge in the hadronic phase, and \f$
    q_Q \f$ is the electric charge in the quark phase.

*/
class ex_eos_gibbs {

public:

  typedef boost::numeric::ublas::vector<double> ubvector;

  /// Hadronic EOS pointer
  eos_had_base *ptr_h;

  /// Skyrme model
  eos_had_skyrme sk;

  /// RMF model
  eos_had_rmf rmf;

  /// RMF model with hyperons
  eos_had_rmf_hyp rmf_hyp;

  /// Quark EOS pointer
  eos_quark *ptr_q;

  /// NJL model
  eos_quark_njl njl;

  /// Bag model
  eos_quark_bag bag;

  /// \name Particle objects
  //@{
  fermion n;
  fermion p;
  fermion e;
  fermion mu;
  quark u;
  quark d;
  quark s;
  //@}
  
  /// \name Thermodynamic quantities
  //@{
  thermo hth;
  thermo qth;
  thermo tot;
  //@}

  /// Electron EOS
  fermion_zerot fzt;

  /// Solver
  mroot_hybrids<> mh;

  /// Minimizer
  mmin_simp2<> mmin;
  
  /// Alternate minimizer (unused)
  anneal_gsl<> mmin2;

  /// Bag constant
  double B;

  /// Surface tension (default \f$ (1~\mathrm{MeV})/(\hbar c) \f$
  double sigma;

  /** \brief If non-zero, then the bag constant will be adjusted 
      to ensure this value is the beginning of the mixed phase
  */
  double mp_start_fix;
  
  /** \brief Determine the bag constant by fixing the density at which
      the mixed phase begins
  */
  int f_bag_constant(size_t nv, const ubvector &x, ubvector &y,
		     double &nB) {

    n.n=x[0];
    p.n=nB-n.n;

    if (ptr_q==&bag) {
      bag.bag_constant=x[1];
    } else if (ptr_q==&njl) {
      njl.B0=x[1];
    }

    ptr_h->calc_e(n,p,hth);
    
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);

    tot.pr=hth.pr;
    tot.ed=hth.ed+e.ed;
    
    y[0]=hth.pr-qth.pr;
    y[1]=p.n-e.n-mu.n;

    return 0;
  }

  /** \brief Solve for the hadronic phase
   */
  int f_had_phase(size_t nv, const ubvector &x, ubvector &y,
		  double &nB) {
    n.n=x[0];
    p.n=nB-n.n;
    
    ptr_h->calc_e(n,p,hth);
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);

    y[0]=e.n-p.n;

    return 0;
  }
  
  /** \brief Solve for the mixed phase (with no Coulomb or surface)
   */
  int f_mixed_phase(size_t nv, const ubvector &x, ubvector &y,
		    double &nB, double &chi) {
    n.n=x[0];
    p.n=x[1];
    
    ptr_h->calc_e(n,p,hth);
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);

    double quark_nqch, quark_nQ;
    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);
    
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;
      
    chi=(e.n-quark_nqch)/(p.n-quark_nqch);

    // Update the energy density and pressure
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+e.ed;
    tot.pr=hth.pr+e.pr;

    y[0]=hth.pr-qth.pr;
    y[1]=nB-(n.n+p.n)*chi-quark_nQ*(1.0-chi)/3.0;

    return 0;
  }

  /** \brief Solve for the quark phase
   */
  int f_quark_phase(size_t nv, const ubvector &x, ubvector &y,
		    double &nB) {
    
    double muQ=x[0];
    e.mu=x[1];
    
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);
    
    double quark_nqch, quark_nQ;
    u.mu=muQ-e.mu*2.0/3.0;
    d.mu=muQ+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);
    
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;

    y[0]=e.n-quark_nqch;
    y[1]=quark_nQ/3.0-nB;

    return 0;
  }

  /** \brief Solve for the beginning of the mixed phase
   */
  int f_beg_mixed_phase(size_t nv, const ubvector &x, ubvector &y,
			double &nB) {

    n.n=x[0];
    p.n=x[1];
    
    ptr_h->calc_e(n,p,hth);
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);
    
    // Compute total baryon density
    nB=n.n+p.n;
    
    y[0]=hth.pr-qth.pr;
    // Ensure electric neutrality in the hadronic phase
    y[1]=p.n-e.n-mu.n;

    return 0;
  }
  
  /** \brief Solve for the end of the mixed phase
   */
  int f_end_mixed_phase(size_t nv, const ubvector &x, ubvector &y,
			double &nB) {

    n.n=x[0];
    p.n=x[1];
    
    ptr_h->calc_e(n,p,hth);
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);

    double quark_nqch, quark_nQ;

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;
    
    // Update the energy density and pressure
    tot.pr=hth.pr+e.pr;

    // Compute total baryon density
    nB=quark_nQ/3.0;
    
    y[0]=hth.pr-qth.pr;
    // Ensure electric neutrality from quarks alone
    y[1]=quark_nqch-e.n-mu.n;

    return 0;
  }
  
  int f_min_densities(size_t nv, const ubvector &x, ubvector &y,
		      double &nB) {
    
    double chi=x[0];
    p.n=x[1];

    if (chi<0.0 || p.n<0.0 || n.n<0.0) return 1;

    if (false) {
      double t1, t2;
      sk.eff_mass(n,p,t1,t2);
    }
    
    if (n.ms<0.0 || p.ms<0.0) return 2;

    ptr_h->calc_e(n,p,hth);
    e.mu=n.mu-p.mu;
    fzt.calc_mu_zerot(e);
    mu.mu=e.mu;
    fzt.calc_mu_zerot(mu);
    
    double quark_nqch, quark_nQ;

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    
    ptr_q->calc_p(u,d,s,qth);
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;
    
    y[0]=e.n-p.n*chi-quark_nqch*(1.0-chi);
    y[1]=nB-(n.n+p.n)*chi-quark_nQ*(1.0-chi)/3.0;
    return 0;
  }
  
  double f_mixed_phase_min(size_t nv, const ubvector &x,
			   double &nB, double &chi) {
    n.n=x[0];

    mm_funct fp_min_densities=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double &)>
       (&ex_eos_gibbs::f_min_densities),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::ref(nB));

    ubvector x2(2), y2(2);
    x2[0]=chi;
    x2[1]=p.n;
    mh.err_nonconv=false;
    mh.def_jac.err_nonconv=false;
    int ret=mh.msolve(2,x2,fp_min_densities);
    mh.def_jac.err_nonconv=true;
    mh.err_nonconv=true;
    chi=x2[0];
    p.n=x2[1];
    f_min_densities(2,x2,y2,nB);
    
    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);

    // Update the energy density and pressure
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+e.ed;
    tot.pr=hth.pr+e.pr;

    if (ret!=0 || chi>1.0) tot.ed+=1.0e6;
    
    return tot.ed;
  }

  /** \brief The Coulomb function
   */
  double fcoul(double d, double chi) {
    double ret;
    if (d==2.0) {
      return 0.25*(chi-1.0-log(chi));
    }
    ret=(2.0/(d-2.0)*(1.0-0.5*d*pow(chi,1.0-2.0/d))+chi)/(d+2.0);
    return ret;
  }

  double f_mixed_phase_min_r(size_t nv, const ubvector &x,
			     double &nB, double &chi, double &dim,
			     double &esurf, double &ecoul) {
    n.n=x[0];
    double r_rare=x[1];

    // Set up the functor
    mm_funct fp_min_densities=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double &)>
       (&ex_eos_gibbs::f_min_densities),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,
       std::ref(nB));

    // Determine the proton density and chi by solving 
    ubvector x2(2), y2(2);
    x2[0]=chi;
    x2[1]=p.n;
    mh.err_nonconv=false;
    mh.def_jac.err_nonconv=false;
    int sret=mh.msolve(2,x2,fp_min_densities);
    mh.def_jac.err_nonconv=true;
    mh.err_nonconv=true;
    chi=x2[0];
    p.n=x2[1];
    f_min_densities(2,x2,y2,nB);

    double quark_nqch, quark_nQ;
    
    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;

    // Surface and Coulomb energy
    double xglen=chi;
    if (chi>0.5) xglen=1.0-chi;
    bool andrew_mode=true;
    if (andrew_mode) {
      esurf=dim*xglen*sigma/r_rare;
      ecoul=2.0*pi*pow(p.n-quark_nqch,2.0)*
	o2scl_const::fine_structure*r_rare*r_rare*xglen*fcoul(dim,xglen);
    } else {
      esurf=dim*chi*sigma/r_rare;
      ecoul=2.0*pi*pow(p.n-quark_nqch,2.0)*
	o2scl_const::fine_structure*r_rare*r_rare*chi*fcoul(dim,xglen);
    }

    // Update the energy density and pressure
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+e.ed+esurf+ecoul;

    // Use the thermodynamic identity for the pressure
    tot.pr=-tot.ed+n.n*chi*n.mu+p.n*chi*p.mu+
      quark_nQ*(1.0-chi)*n.mu/3.0;

    //cout << 1.0/cbrt(4.0*pi*fabs(p.n*p.n-quark_nqch*quark_nqch)*
    //o2scl_const::fine_structure/
    //		       sigma/dim*fcoul(dim,chi)) << endl;
    /*
      cout << "z: " << x[0] << " " << x[1] << " "
      << fabs(p.n*p.n-quark_nqch*quark_nqch) << " "
      << fcoul(dim,chi) << " " 
      << tot.ed << endl;
    */
    if (chi<0.0 || chi>1.0 || sret!=0 || r_rare>50.0 || r_rare<0.0) tot.ed+=1e6;
    return tot.ed;
  }
  
  ex_eos_gibbs() {

    skyrme_load(sk,"NRAPR");
    ptr_h=&sk;

    ptr_q=&bag;
    
    n.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    p.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_proton),2.0);
    e.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_electron),2.0);
    mu.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_muon),2.0);
    u.init(0.0,6.0);
    d.init(0.0,6.0);
    s.init(95.0/hc_mev_fm,6.0);
    
    n.non_interacting=false;
    p.non_interacting=false;
    n.inc_rest_mass=true;
    p.inc_rest_mass=true;
    
    // Make the minimizer a bit more accurate
    mmin.tol_rel/=1.0e2;
    mmin.tol_abs/=1.0e2;
    mmin.ntrial*=100.0;

    // Set a small stepsize for the minimizer
    double step[1]={0.01};
    mmin.set_step(1,step);

    sigma=1.0/hc_mev_fm;
    
    mp_start_fix=0.0;
  }

  /** \brief Compute the EOS and M-R curve at T=0
   */
  int mvsr(vector<string> &sv, bool itive_com) {

    double mp_end, chi, nB, dim=3.0, esurf, ecoul, mp_start;
    
    mm_funct fp_had_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double &)>
       (&ex_eos_gibbs::f_had_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::ref(nB));
    mm_funct fp_mixed_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double &, double &)>
       (&ex_eos_gibbs::f_mixed_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,
       std::ref(nB),std::ref(chi));
    mm_funct fp_quark_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &, double &)>
       (&ex_eos_gibbs::f_quark_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,
       std::ref(nB));
    mm_funct fp_beg_mixed_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double &)>
       (&ex_eos_gibbs::f_beg_mixed_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,
       std::ref(nB));
    mm_funct fp_end_mixed_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double &)>
       (&ex_eos_gibbs::f_end_mixed_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,
       std::ref(nB));
    multi_funct fp_mixed_phase_min=std::bind
      (std::mem_fn<double(size_t,const ubvector &, double &, double &)>
       (&ex_eos_gibbs::f_mixed_phase_min),this,std::placeholders::_1,
       std::placeholders::_2,std::ref(nB),std::ref(chi));
    multi_funct fp_mixed_phase_min_r=std::bind
      (std::mem_fn<double(size_t,const ubvector &, double &, double &,
			  double &, double &, double &)>
       (&ex_eos_gibbs::f_mixed_phase_min_r),this,std::placeholders::_1,
       std::placeholders::_2,std::ref(nB),std::ref(chi),std::ref(dim),
       std::ref(esurf),std::ref(ecoul));

    ubvector x(2), y(2);
    
    cout << "Masses (n,p,e): " << n.m*hc_mev_fm << " "
	 << p.m*hc_mev_fm << " " << e.m*hc_mev_fm << " MeV" << endl;
    cout << endl;

    // -----------------------------------------------------------------
    // Determine bag constant

    if (mp_start_fix>0 && (ptr_q==&bag || ptr_q==&njl)) {
      
      mp_start=mp_start_fix;
      
      mm_funct fp_bag_constant=std::bind
	(std::mem_fn<int(size_t,const ubvector &,ubvector &,
			 double &)>
	 (&ex_eos_gibbs::f_bag_constant),this,std::placeholders::_1,
	 std::placeholders::_2, std::placeholders::_3,
	 std::ref(nB));

      cout << "Determine B by fixing the "
	   << "beginning of the mixed phase to\n n_B=" << mp_start
	   << " fm^{-3}:" << endl;
      nB=mp_start;
      x[0]=mp_start*0.9;
      x[1]=1.0;
      mh.msolve(2,x,fp_bag_constant);
      B=x[1];
      cout << "B: " << B*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "Densities (n,p,e): " << n.n << " " << p.n << " " << e.n
	   << " fm^{-3}" << endl;
      cout << "Quark densities (u,d,s): " << u.n << " " << d.n << " "
	   << s.n << " fm^{-3}" << endl;
      cout << "Chem pots. (n,p,e): " << n.mu*hc_mev_fm << " "
	   << p.mu*hc_mev_fm << " " << e.mu*hc_mev_fm << " MeV" << endl;
      cout << "Quark chem. pots. (u,d,s): " << u.mu*hc_mev_fm << " "
	   << d.mu*hc_mev_fm << " " << s.mu*hc_mev_fm << " fm^{-3}" << endl;
      cout << "Electron energy density: " << e.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Total pressure: " << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "Total energy density: " << tot.ed*hc_mev_fm
	   << " MeV/fm^3" << endl;
      {
	cout << "Skyrme energy density: " << hth.ed*hc_mev_fm << " MeV/fm^3"
	     << endl;
	cout << "Skyrme pressure: " << hth.pr*hc_mev_fm << " MeV/fm^3"
	     << endl;
      }
      {
	cout << "Quark energy dens. & pressure: "
	     << qth.ed*hc_mev_fm << " MeV/fm^3, "
	     << qth.pr*hc_mev_fm << " MeV/fm^3"
	     << endl;
	//cout << "Number density of quarks: " << quark_nQ << " 1/fm^3" << endl;
	//cout << "Charge density in quark matter: " << quark_nqch
	//<< " 1/fm^3" << endl;
      }
      cout << endl;

    } else {

      // -----------------------------------------------------------------
      // Find the beginning of the mixed phase with fp_beg_mixed_phase()
      
      cout << "Beginning of mixed phase: " << endl;
      x[0]=0.3;
      x[1]=0.1;
      mh.msolve(2,x,fp_beg_mixed_phase);
      n.n=x[0];
      p.n=x[1];
      f_beg_mixed_phase(2,x,y,nB);
      mp_start=nB;
      cout << "Baryon density: " << nB << " fm^{-3}" << endl;
      cout << endl;
      
    }

    // -----------------------------------------------------------------
    // Find the end of the mixed phase with fp_end_mixed_phase()

    if (true) {
      cout << "End of mixed phase: " << endl;
      x[0]=0.9;
      x[1]=0.2;
      mh.msolve(2,x,fp_end_mixed_phase);
      n.n=x[0];
      p.n=x[1];
      f_end_mixed_phase(2,x,y,nB);
      mp_end=nB;
      cout << "Baryon density: " << nB << " fm^{-3}" << endl;
      cout << endl;
    }

    // -----------------------------------------------------------------
    // Test solver and minimizer inside mixed phase

    cout << "Mixed phase at n_B=" << mp_start+0.01
	 << " fm^{-3} from solver:" << endl;
    nB=mp_start+0.01;
    x[0]=nB*0.9;
    x[1]=0.02;
    mh.msolve(2,x,fp_mixed_phase);
    n.n=x[0];
    p.n=x[1];
    f_mixed_phase(2,x,y,nB,chi);
    cout << "Chi: " << chi << endl;
    cout << "Densities (n,p,e): " << n.n << " " << p.n << " " << e.n 
	 << " fm^{-3}" << endl;
    cout << "Chem pots. (n,p,e): " << n.mu*hc_mev_fm << " "
	 << p.mu*hc_mev_fm << " " << e.mu*hc_mev_fm << " MeV" << endl;
    cout << "Electron energy density: " << e.ed*hc_mev_fm << " MeV/fm^3"
	 << endl;
    cout << "Total pressure: " << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
    cout << "Total energy density: " << tot.ed*hc_mev_fm << " MeV/fm^3" << endl;
    {
      cout << "Skyrme energy density: " << hth.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Skyrme pressure: " << hth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
    }
    {
      cout << "Quark energy dens. & pressure: "
	   << qth.ed*hc_mev_fm << " MeV/fm^3, "
	   << qth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
      //cout << "Number density of quarks: " << quark_nQ << " 1/fm^3" << endl;
      //cout << "Charge density in quark matter: " << quark_nqch
      //<< " 1/fm^3" << endl;
    }
    {
      cout << "Check thermodynamic identities:" << endl;
      cout << "\tSkyrme: " << fabs(hth.pr+hth.ed-n.n*n.mu-p.n*p.mu)/
	(hth.pr) << endl;
      cout << "\tElectrons: " << fabs(e.pr+e.ed-e.n*e.mu)/e.pr << endl;
      /*
	cout << "\tQuarks: " << fabs(qth.pr+qth.ed-quark_nQ*n.mu/3.0+
	quark_nqch*e.mu)/(qth.pr) << endl;
	cout << "\tTotal: " << fabs(tot.pr+tot.ed-chi*n.n*n.mu-
	chi*p.n*p.mu-e.n*e.mu-
	(1.0-chi)*quark_nQ*n.mu/3.0+
	(1.0-chi)*quark_nqch*e.mu)/
	fabs(tot.pr) << endl;
      */
    }
    cout << endl;
      exit(-1);
    
    cout << "Mixed phase at n_B=" << nB << " fm^{-3} from minimizer:" << endl;
    cout << "  Notice this is slightly different from the above" << endl;
    cout << "  due to numerical errors." << endl;
    x[0]=0.31;
    mmin.mmin(1,x,y[0],fp_mixed_phase_min);
    cout << "Chi: " << chi << endl;
    cout << "Densities (n,p,e): " << n.n << " " << p.n << " " << e.n 
	 << " fm^{-3}" << endl;
    cout << "Chem pots. (n,p,e): " << n.mu*hc_mev_fm << " "
	 << p.mu*hc_mev_fm << " " << e.mu*hc_mev_fm << " MeV" << endl;
    cout << "Electron energy density: " << e.ed*hc_mev_fm << " MeV/fm^3"
	 << endl;
    cout << "Total pressure: " << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
    cout << "Total energy density: " << tot.ed*hc_mev_fm << " MeV/fm^3" << endl;
    {
      cout << "Skyrme energy density: " << hth.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Skyrme pressure: " << hth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
    }
    {
      cout << "Quark energy dens. & pressure: "
	   << qth.ed*hc_mev_fm << " MeV/fm^3, "
	   << qth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
      //cout << "Number density of quarks: " << quark_nQ << " 1/fm^3" << endl;
      //cout << "Charge density in quark matter: " << quark_nqch
      //<< " 1/fm^3" << endl;
    }
    {
      cout << "Check thermodynamic identities:" << endl;
      cout << "\tSkyrme: " << fabs(hth.pr+hth.ed-n.n*n.mu-p.n*p.mu)/
	(hth.pr) << endl;
      cout << "\tElectrons: " << fabs(e.pr+e.ed-e.n*e.mu)/e.pr << endl;
      /*
	cout << "\tQuarks: " << fabs(qth.pr+qth.ed-quark_nQ*n.mu/3.0+
	quark_nqch*e.mu)/(qth.pr) << endl;
	cout << "\tTotal: " << fabs(tot.pr+tot.ed-chi*n.n*n.mu-
	chi*p.n*p.mu-e.n*e.mu-
	(1.0-chi)*quark_nQ*n.mu/3.0+
	(1.0-chi)*quark_nqch*e.mu)/
	fabs(tot.pr) << endl;
      */
    }
    cout << endl;

    // -----------------------------------------------------------------
    // Find the end of the mixed phase with fp_mixed_phase() with
    // solver and minimizer

    cout << "Mixed phase at n_B=" << mp_end-0.01
	 << " fm^{-3} from solver:" << endl;
    nB=mp_end-0.01;
    x[0]=1.0;
    x[1]=0.12;
    mh.msolve(2,x,fp_mixed_phase);
    n.n=x[0];
    p.n=x[1];
    f_mixed_phase(2,x,y,nB,chi);
    cout << "Chi: " << chi << endl;
    cout << "Densities (n,p,e): " << n.n << " " << p.n << " " << e.n 
	 << " fm^{-3}" << endl;
    cout << "Chem pots. (n,p,e): " << n.mu*hc_mev_fm << " "
	 << p.mu*hc_mev_fm << " " << e.mu*hc_mev_fm << " MeV" << endl;
    cout << "Electron energy density: " << e.ed*hc_mev_fm << " MeV/fm^3"
	 << endl;
    cout << "Total pressure: " << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
    cout << "Total energy density: " << tot.ed*hc_mev_fm << " MeV/fm^3" << endl;
    {
      cout << "Skyrme energy density: " << hth.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Skyrme pressure: " << hth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
    }
    {
      cout << "Quark energy dens. & pressure: "
	   << qth.ed*hc_mev_fm << " MeV/fm^3, "
	   << qth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
      //cout << "Number density of quarks: " << quark_nQ << " 1/fm^3" << endl;
      //cout << "Charge density in quark matter: " << quark_nqch
      //<< " 1/fm^3" << endl;
    }
    {
      cout << "Check thermodynamic identities:" << endl;
      cout << "\tSkyrme: " << fabs(hth.pr+hth.ed-n.n*n.mu-p.n*p.mu)/
	(hth.pr) << endl;
      cout << "\tElectrons: " << fabs(e.pr+e.ed-e.n*e.mu)/e.pr << endl;
      /*
	cout << "\tQuarks: " << fabs(qth.pr+qth.ed-quark_nQ*n.mu/3.0+
	quark_nqch*e.mu)/(qth.pr) << endl;
	cout << "\tTotal: " << fabs(tot.pr+tot.ed-chi*n.n*n.mu-
	chi*p.n*p.mu-e.n*e.mu-
	(1.0-chi)*quark_nQ*n.mu/3.0+
	(1.0-chi)*quark_nqch*e.mu)/
	fabs(tot.pr) << endl;
      */
    }
    cout << endl;

    nB=mp_end-0.01;
    chi=0.00572;
    p.n=0.35748;
    cout << "Mixed phase at n_B=" << nB
	 << " fm^{-3} from minimizer:" << endl;
    x[0]=0.40223;
    mmin.mmin(1,x,y[0],fp_mixed_phase_min);
    cout << "Chi: " << chi << endl;
    cout << "Densities (n,p,e): " << n.n << " " << p.n << " " << e.n 
	 << " fm^{-3}" << endl;
    cout << "Chem pots. (n,p,e): " << n.mu*hc_mev_fm << " "
	 << p.mu*hc_mev_fm << " " << e.mu*hc_mev_fm << " MeV" << endl;
    cout << "Electron energy density: " << e.ed*hc_mev_fm << " MeV/fm^3"
	 << endl;
    cout << "Total pressure: " << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
    cout << "Total energy density: " << tot.ed*hc_mev_fm << " MeV/fm^3" << endl;
    {
      cout << "Skyrme energy density: " << hth.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Skyrme pressure: " << hth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
    }
    {
      cout << "Quark energy dens. & pressure: "
	   << qth.ed*hc_mev_fm << " MeV/fm^3, "
	   << qth.pr*hc_mev_fm << " MeV/fm^3"
	   << endl;
      //cout << "Number density of quarks: " << quark_nQ << " 1/fm^3" << endl;
      //cout << "Charge density in quark matter: " << quark_nqch
      //<< " 1/fm^3" << endl;
    }
    {
      cout << "Check thermodynamic identities:" << endl;
      cout << "\tSkyrme: " << fabs(hth.pr+hth.ed-n.n*n.mu-p.n*p.mu)/
	(hth.pr) << endl;
      cout << "\tElectrons: " << fabs(e.pr+e.ed-e.n*e.mu)/e.pr << endl;
      /*
	cout << "\tQuarks: " << fabs(qth.pr+qth.ed-quark_nQ*n.mu/3.0+
	quark_nqch*e.mu)/(qth.pr) << endl;
	cout << "\tTotal: " << fabs(tot.pr+tot.ed-chi*n.n*n.mu-
	chi*p.n*p.mu-e.n*e.mu-
	(1.0-chi)*quark_nQ*n.mu/3.0+
	(1.0-chi)*quark_nqch*e.mu)/
	fabs(tot.pr) << endl;
      */
    }
    cout << endl;

    // -----------------------------------------------------------------
    // Tabulate full hadronic phase

    table_units<> thad;
    if (true) {
      thad.line_of_names("nB nn np ne ede edh edt pre prh prt");
      cout << "Hadronic phase: " << endl;
      x[0]=0.07;
      for(nB=0.08;nB<1.5001;nB+=0.01) {
	mh.msolve(1,x,fp_had_phase);
	std::vector<double> line={nB,n.n,p.n,e.n,e.ed*hc_mev_fm,
				  hth.ed*hc_mev_fm,(hth.ed+e.ed)*hc_mev_fm,
				  e.pr*hc_mev_fm,hth.pr*hc_mev_fm,
				  (hth.pr+e.pr)*hc_mev_fm};
	thad.line_of_data(line);
      }
      hdf_file hf;
      hf.open_or_create("ex_eos_gibbs.o2");
      hdf_output(hf,thad,"had");
      hf.close();
    }

    // -----------------------------------------------------------------
    // Tabulate full quark phase

    table_units<> tq;
    {
      tq.line_of_names("nB muQ mue ede pre edq prq edt prt");
      cout << "Quark phase: " << endl;
      x[0]=1.0;
      x[1]=0.1;
      for(nB=0.3;nB<1.5001;nB+=0.01) {
	mh.msolve(2,x,fp_quark_phase);
	std::vector<double> line={nB,x[0],x[1],e.ed*hc_mev_fm,
				  e.pr*hc_mev_fm,qth.ed*hc_mev_fm,
				  qth.pr*hc_mev_fm,(e.ed+qth.ed)*hc_mev_fm,
				  (e.pr+qth.pr)*hc_mev_fm};
	tq.line_of_data(line);
      }
      hdf_file hf;
      hf.open_or_create("ex_eos_gibbs.o2");
      hdf_output(hf,tq,"quark");
      hf.close();
    }

    // -----------------------------------------------------------------
    // Tabulate full mixed phase without Coulomb or surface

    if (true) {
      cout << "Mixed phase: " << endl;
      cout.precision(4);
      table<> t;
      t.line_of_names("nB nn ede pre edt prt edh prh edq prq mun mup chi");
      
      double nB_start=(mp_start+mp_end)/2.0;

      // initial guesses
      chi=0.5;
      p.n=nB_start*0.2;
      x[0]=nB_start*0.8;

      bool high_done=false;
      for(nB=nB_start;nB<2.0 && high_done==false;nB+=0.01) {
	mmin.mmin(1,x,y[0],fp_mixed_phase_min);
	std::vector<double> line;
	line.push_back(nB);
	line.push_back(n.n);
	line.push_back(e.ed*hc_mev_fm);
	line.push_back(e.pr*hc_mev_fm);
	line.push_back(tot.ed*hc_mev_fm);
	line.push_back(tot.pr*hc_mev_fm);
	{
	  line.push_back(hth.ed*hc_mev_fm);
	  line.push_back(hth.pr*hc_mev_fm);
	}
	{
	  line.push_back(qth.ed*hc_mev_fm);
	  line.push_back(qth.pr*hc_mev_fm);
	}
	line.push_back(n.mu*hc_mev_fm);
	line.push_back(p.mu*hc_mev_fm);
	line.push_back(chi);
	t.line_of_data(line.size(),line);
	if (chi<0.05) high_done=true;
      }

      // initial guesses
      chi=0.5;
      p.n=nB_start*0.2;
      x[0]=nB_start*0.8;

      bool low_done=false;
      for(nB=nB_start-0.01;nB>0.0 && low_done==false;nB-=0.01) {
	mmin.mmin(1,x,y[0],fp_mixed_phase_min);
	std::vector<double> line;
	line.push_back(nB);
	line.push_back(n.n);
	line.push_back(e.ed*hc_mev_fm);
	line.push_back(e.pr*hc_mev_fm);
	line.push_back(tot.ed*hc_mev_fm);
	line.push_back(tot.pr*hc_mev_fm);
	{
	  line.push_back(hth.ed*hc_mev_fm);
	  line.push_back(hth.pr*hc_mev_fm);
	}
	{
	  line.push_back(qth.ed*hc_mev_fm);
	  line.push_back(qth.pr*hc_mev_fm);
	}
	line.push_back(n.mu*hc_mev_fm);
	line.push_back(p.mu*hc_mev_fm);
	line.push_back(chi);
	t.line_of_data(line.size(),line);
	if (chi>0.95) low_done=true;
      }
      
      hdf_file hf;
      t.sort_table("nB");
      hf.open_or_create("ex_eos_gibbs.o2");
      hdf_output(hf,t,"mixed");
      hf.close();

      cout << endl;
    }

    // -----------------------------------------------------------------
    // Tabulate full mixed phase with Coulomb and surface

    table_units<> t[4];

    if (true) {
      sigma=1.0/hc_mev_fm;
      for(size_t id=0;id<3;id++) {
	t[id].line_of_names(((string)"nB nn np ede pre edt prt edh ")+
			    "prh edq prq mun mup chi r_rare esurf ecoul "+
			    "nqch ne r_ws");
	dim=((double)id)+1.0;

	double nB_start=(mp_start+mp_end)/2.0;
	
	// initial guesses
	chi=0.5;
	p.n=nB_start*0.2;
	x[0]=nB_start*0.8;
	x[1]=1.3;

	bool high_done=false;
	for(nB=nB_start;nB<2.0 && high_done==false;nB+=0.001) {
	  //cout << "nB: " << id << " " << nB << " " << chi << endl;
	  
	  mmin.mmin(2,x,y[0],fp_mixed_phase_min_r);
	  double r_rare=x[1];
	  //cout << id << " " << nB << " " << n.n << " " << p.n << " "
	  //<< x[1] << " " << tot.ed << endl;
	  double xglen=chi;
	  if (chi>0.5) xglen=1.0-chi;
	  double r_ws=r_rare/pow(xglen,1.0/dim);
	  std::vector<double> line={nB,n.n,p.n,e.ed*hc_mev_fm,
				    e.pr*hc_mev_fm,tot.ed*hc_mev_fm,
				    tot.pr*hc_mev_fm,hth.ed*hc_mev_fm,
				    hth.pr*hc_mev_fm,qth.ed*hc_mev_fm,
				    qth.pr*hc_mev_fm,n.mu*hc_mev_fm,
				    p.mu*hc_mev_fm,chi,x[1],
				    esurf*hc_mev_fm,ecoul*hc_mev_fm,
				    0.0,e.n,r_ws};
	  t[id].line_of_data(line);
	  if (chi<0.05) high_done=true;
	  
	}

	// initial guesses
	chi=0.5;
	p.n=nB_start*0.2;
	x[0]=nB_start*0.8;
	x[1]=1.3;

	bool low_done=false;
	for(nB=nB_start-0.01;nB>0.0 && low_done==false;nB-=0.001) {
	  //cout << "nB: " << id << " " << nB << " " << chi << endl;
	  
	  mmin.mmin(2,x,y[0],fp_mixed_phase_min_r);
	  double r_rare=x[1];
	  //cout << id << " " << nB << " " << n.n << " " << p.n << " "
	  //<< x[1] << " " << tot.ed << endl;
	  double xglen=chi;
	  if (chi>0.5) xglen=1.0-chi;
	  double r_ws=r_rare/pow(xglen,1.0/dim);
	  std::vector<double> line={nB,n.n,p.n,e.ed*hc_mev_fm,
				    e.pr*hc_mev_fm,tot.ed*hc_mev_fm,
				    tot.pr*hc_mev_fm,hth.ed*hc_mev_fm,
				    hth.pr*hc_mev_fm,qth.ed*hc_mev_fm,
				    qth.pr*hc_mev_fm,n.mu*hc_mev_fm,
				    p.mu*hc_mev_fm,chi,x[1],
				    esurf*hc_mev_fm,ecoul*hc_mev_fm,
				    0.0,e.n,r_ws};
	  t[id].line_of_data(line);
	  if (chi>0.95) low_done=true;
	  
	}

	if (true) {
	  t[id].sort_table("nB");
	  hdf_file hf;
	  hf.open_or_create("ex_eos_gibbs.o2");
	  hdf_output(hf,t[id],"mixed"+o2scl::itos(id+1));
	  hf.close();
	}
      }
      
      if (true) {
	t[3].line_of_names(((string)"nB nn np ede pre edt prt edh ")+
			   "prh edq prq mun mup chi r_rare esurf ecoul "+
			   "nqch ne r_ws");
	size_t count=0;
	double nB_start=t[0].min("nB");
	if (t[1].min("nB")>nB_start) nB_start=t[1].min("nB");
	if (t[2].min("nB")>nB_start) nB_start=t[2].min("nB");
	double nB_end=t[0].max("nB");
	if (t[1].max("nB")<nB_end) nB_end=t[1].max("nB");
	if (t[2].max("nB")<nB_end) nB_end=t[2].max("nB");
	for(nB=nB_start;nB<nB_end;nB+=0.001) {
	  std::vector<double> line;
	  if (t[0].interp("nB",nB,"edt")<t[1].interp("nB",nB,"edt") &&
	      t[0].interp("nB",nB,"edt")<t[2].interp("nB",nB,"edt")) {
	    for(size_t j=0;j<t[3].get_ncolumns();j++) {
	      line.push_back(t[0].interp("nB",nB,t[3].get_column_name(j)));
	    }
	    if (count%20==0) std::cout << nB << " d=1" << std::endl;
	  } else if (t[1].interp("nB",nB,"edt")<t[0].interp("nB",nB,"edt") &&
		     t[1].interp("nB",nB,"edt")<t[2].interp("nB",nB,"edt")) {
	    if (count%20==0) std::cout << nB << " d=2" << std::endl;
	    for(size_t j=0;j<t[3].get_ncolumns();j++) {
	      line.push_back(t[1].interp("nB",nB,t[3].get_column_name(j)));
	    }
	  } else {
	    if (count%20==0) std::cout << nB << " d=3" << std::endl;
	    for(size_t j=0;j<t[3].get_ncolumns();j++) {
	      line.push_back(t[2].interp("nB",nB,t[3].get_column_name(j)));
	    }
	  }
	  t[3].line_of_data(line.size(),line);
	  count++;
	}
      }

      // Output the mixed phase with optimized dimensionality
      if (true) {
	hdf_file hf;
	hf.open_or_create("ex_eos_gibbs.o2");
	hdf_output(hf,t[3],"mixed_min");
	hf.close();
      }
      
    }

    // -----------------------------------------------------------------
    // Construct the full neutron star EOS

    t[3].set_interp_type(itp_linear);
    double had_end=t[3].interp("chi",1.05,"nB");
    double quark_start=t[3].interp("chi",-0.05,"nB");
    cout << "H: " << had_end << " " << quark_start << endl;

    table_units<> ns;
    ns.line_of_names("ed pr");
    ns.set_unit("ed","MeV/fm^3");
    ns.set_unit("pr","MeV/fm^3");

    for(double nB=0.08;nB<had_end;nB+=0.01) {
      double line[2]={thad.interp("nB",nB,"edt"),
		      thad.interp("nB",nB,"prt")};
      ns.line_of_data(2,line);
    }
    for(size_t i=0;i<t[3].get_nlines();i++) {
      double line[2]={t[3].get("edt",i),t[3].get("prt",i)};
      ns.line_of_data(2,line);
    }
    for(nB=quark_start;nB<1.4;nB+=0.01) {
      double line[2]={tq.interp("nB",nB,"edt"),
		      tq.interp("nB",nB,"prt")};
      ns.line_of_data(2,line);
    }

    hdf_file hf;
    hf.open_or_create("ex_eos_gibbs.o2");
    hdf_output(hf,ns,"nstar");
    hf.close();

    // tov solve
    eos_tov_interp eti;
    eti.default_low_dens_eos();
    eti.read_table(ns,"ed","pr");
    cout << "Going to tov_solve." << endl;
    tov_solve ts;
    ts.set_eos(eti);
    ts.mvsr();
    std::shared_ptr<table_units<> > tov=ts.get_results();
    cout << "M_max: " << tov->max("gm") << endl;
    cout << tov->get_unit("ed") << endl;
    cout << tov->get("ed",tov->lookup("gm",tov->max("gm"))) << endl;

    hf.open_or_create("ex_eos_gibbs.o2");
    hdf_output(hf,*tov,"tov");
    hf.close();
    
    return 0;
  }

  /** \brief Select a hadronic model
   */
  int model_hadrons(vector<string> &sv, bool itive_com) {

    if (sv.size()<2) {
      cerr << "No model specified in model_hadrons()." << endl;
      return 1;
    }
    
    if (sv[1]=="sk") {
      ptr_h=&sk;
      skyrme_load(sk,sv[2]);
    } else if (sv[1]=="rmf") {
      ptr_h=&rmf;
    } else if (sv[1]=="SLB00") {
      ptr_h=&rmf;

      n.m=939.0/hc_mev_fm;
      p.m=939.0/hc_mev_fm;
      
      rmf.set_n_and_p(n,p);
      
      rmf.mnuc=939.0/hc_mev_fm;
      rmf.n0=0.16;
      rmf.eoa=16.0/hc_mev_fm;
      rmf.comp=250.0/hc_mev_fm;
      rmf.msom=0.6;
      rmf.esym=36.0/hc_mev_fm;
      rmf.ms=500.0/hc_mev_fm;
      rmf.mw=763.0/hc_mev_fm;
      rmf.mr=770.0/hc_mev_fm;
      rmf.zeta=0.0;
      rmf.xi=0.0;
      rmf.a1=0.0;
      rmf.a2=0.0;
      rmf.a3=0.0;
      rmf.a4=0.0;
      rmf.a5=0.0;
      rmf.a6=0.0;
      rmf.b1=0.0;
      rmf.b2=0.0;
      rmf.b3=0.0;
      n.mu=n.m*1.1;
      p.mu=p.m*1.1;
      rmf.set_fields(0.2,0.15,-0.001);
      cout << "Going to fix_saturation()." << endl;
      rmf.fix_saturation();
      cout << "Done." << endl;
      rmf.saturation();
      cout << rmf.n0 << " " << rmf.msom << endl;
      cout << "Selected the RMF model from SLB00 (no hyperons)." << endl;
    } else if (sv[1]=="SLB00_hyp") {
      ptr_h=&rmf_hyp;
      rmf_hyp.n0=0.16;
      rmf_hyp.eoa=16.0/hc_mev_fm;
      rmf_hyp.comp=250.0/hc_mev_fm;
      rmf_hyp.msom=0.6;
      rmf_hyp.esym=36.0/hc_mev_fm;
      rmf_hyp.zeta=0.0;
      rmf_hyp.xi=0.0;
      rmf_hyp.a1=0.0;
      rmf_hyp.a2=0.0;
      rmf_hyp.a3=0.0;
      rmf_hyp.a4=0.0;
      rmf_hyp.a5=0.0;
      rmf_hyp.a6=0.0;
      rmf_hyp.b1=0.0;
      rmf_hyp.b2=0.0;
      rmf_hyp.b3=0.0;
      rmf_hyp.fix_saturation();
      rmf_hyp.xs=0.8;
      rmf_hyp.xw=0.895;
      rmf_hyp.xr=0.8;
    } else if (sv[1]=="rmfh") {
      ptr_h=&rmf_hyp;
    } else {
      cerr << "Hadronic model specification " << sv[1]
	   << " not understood." << endl;
      return 2;
    }
    return 0;
  }
  
  /** \brief Select a quark model
   */
  int model_quarks(vector<string> &sv, bool itive_com) {
    
    if (sv.size()<2) {
      cerr << "No model specified in model_quarks()." << endl;
      return 1;
    }
    
    if (sv[1]=="bag") {
      ptr_q=&bag;
    } else if (sv[1]=="SLB00_bag") {
      ptr_q=&bag;
      u.m=5.5/hc_mev_fm;
      d.m=5.5/hc_mev_fm;
      s.m=140.7/hc_mev_fm;
      bag.bag_constant=200.0/hc_mev_fm;
      mp_start_fix=0.0;
      cout << "Selected the bag model from SLB00." << endl;
    } else if (sv[1]=="SLB00_njl") {
      ptr_q=&njl;
      u.m=5.5/hc_mev_fm;
      d.m=5.5/hc_mev_fm;
      s.m=140.7/hc_mev_fm;
      njl.set_quarks(u,d,s);
      njl.up_default_mass=5.5/hc_mev_fm;
      njl.down_default_mass=5.5/hc_mev_fm;
      njl.strange_default_mass=140.7/hc_mev_fm;
      njl.set_parameters(602.3/hc_mev_fm,
			 1.835/njl.L/njl.L,
			 12.36/pow(njl.L,5.0));
      mp_start_fix=0.0;
    } else if (sv[1]=="njl") {
      ptr_q=&njl;
    } else {
      cerr << "Quark model specification " << sv[1]
	   << " not understood." << endl;
      return 2;
    }
    
    return 0;
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  ex_eos_gibbs ehg;
  
  cli cl;
  cl.prompt="ex_eos_gibbs>";
  int comm_option_both=2;
  
  static const int narr=3;
  comm_option_s options_arr[narr]={
    {0,"hadrons","Select hadronic model.",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::model_hadrons),
     cli::comm_option_both},
    {0,"quarks","Select quark model.",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::model_quarks),
     cli::comm_option_both},
    {0,"mvsr","Compute EOS and M-R curve",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::mvsr),
     cli::comm_option_both}
  };
  
  cl.set_comm_option_vec(narr,options_arr);
  cl.cmd_name="ex_eos_gibbs";
  
  cl.run_auto(argc,argv);
  
  return 0;
}
