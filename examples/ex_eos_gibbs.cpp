/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2019-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/convert_units.h>
#include <o2scl/fermion.h>
#include <o2scl/mroot_hybrids.h>
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

#define XCHECK if (x.size()<=ix) { cout << __LINE__ << endl; exit(-1); }
#define YCHECK if (y.size()<=ix) { cout << __LINE__ << endl; exit(-1); }

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
  eos_had_temp_base *ptr_h;

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

  /// NJL model
  eos_quark_njl_vec njlv;

  /// Bag model
  eos_quark_bag bag;

  /// \name Particle objects
  //@{
  fermion n;
  fermion p;
  fermion e;
  fermion mu;
  fermion nu_e;
  fermion nu_mu;
  fermion lam;
  fermion sigz;
  fermion sigm;
  fermion sigp;
  fermion casz;
  fermion casm;
  quark u;
  quark d;
  quark s;
  //@}
  
  /// \name Thermodynamic quantities
  //@{
  thermo hth;
  thermo lep;
  thermo qth;
  thermo tot;
  //@}

  /// Lepton thermodynamics
  fermion_rel fr;

  /// Solver
  mroot_hybrids<> mh;

  /// Bag constant
  double B;

  /// Surface tension (default \f$ (1~\mathrm{MeV})/(\hbar c) \f$
  double surften;

  /// The entropy per baryon (default 0.0)
  double sonB;

  /// The electron lepton fraction (default -1.0)
  double YLe;
  
  /// The muon lepton fraction (default -1.0)
  double YLmu;
  
  /** \brief If non-zero, then the bag constant will be adjusted 
      to ensure this value is the beginning of the mixed phase
  */
  double mp_start_fix;

  /** \brief Guess for the hadronic phase

      Baryon number density, neutron density, temperature, electron
      neutrino chemical potential, and muon neutrino chemical
      potential, in order
  */
  ubvector had_phase_guess;

  /** \brief Guess for the quark phase

      Quark number chemical potential, electron chemical potential,
      temperature, electron neutrino chemical potential, and muon
      neutrino chemical potential, in order
   */
  ubvector quark_phase_guess;

  /** \brief Guess for the beginning of the mixed phase
      
      neutron number density, proton number density, temperature,
      electron neutrino chemical potential, and muon neutrino chemical
      potential
   */
  ubvector beg_mixed_phase_guess;
  
  /** \brief Determine the bag constant by fixing the density at which
      the mixed phase begins at \f$ T=0 \f$
  */
  int f_bag_constant(size_t nv, const ubvector &x, ubvector &y,
		     double &nB) {

    //cout << "x: ";
    //vector_out(cout,nv,x,true);

    n.n=x[0];
    p.n=nB-n.n;

    if (ptr_q==&bag) {
      bag.bag_constant=x[1];
    } else if (ptr_q==&njl) {
      njl.B0=x[1];
    }

    ptr_h->calc_e(n,p,hth);
    
    e.mu=n.mu-p.mu;
    fr.calc_mu_zerot(e);
    mu.mu=e.mu;
    fr.calc_mu_zerot(mu);

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_p(u,d,s,qth);

    y[0]=hth.pr-qth.pr;
    y[1]=p.n-e.n-mu.n;

    tot=hth+qth+lep;
    
    //cout << "y: ";
    //vector_out(cout,nv,y,true);
    
    return 0;
  }

  /** \brief Given the current index in the vector \c x,
      set the neutrino chemical potentials
  */
  void leptons_in(size_t &ix, const ubvector &x) {
    
    if (YLe>-0.5) {
      XCHECK;
      nu_e.mu=x[ix++];
    } else {
      nu_e.mu=0.0;
    }
    if (YLmu>-0.5) {
      XCHECK;
      nu_mu.mu=x[ix++];
    } else {
      nu_mu.mu=0.0;
    }

    return;
  }

  /** \brief Given the current index in the vector \c y and presuming
      the electron chemical potential is set, compute all leptons and
      add two new equations
  */
  void leptons_out(size_t &ix, ubvector &y, double T, double nB) {

    fr.pair_mu(e,T);
    
    mu.mu=e.mu-nu_e.mu+nu_mu.mu;
    fr.pair_mu(mu,T);
    
    lep.ed=e.ed+mu.ed;
    lep.pr=e.pr+mu.pr;
    lep.en=e.en+mu.en;
    
    if (nu_e.mu!=0.0) {
      fr.massless_pair_mu(nu_e,T);
      lep.ed+=nu_e.ed;
      lep.pr+=nu_e.pr;
      lep.en+=nu_e.en;
    }
    if (nu_mu.mu!=0.0) {
      fr.massless_pair_mu(nu_mu,T);
      lep.ed+=nu_mu.ed;
      lep.pr+=nu_mu.pr;
      lep.en+=nu_mu.en;
    }

    if (YLe>-0.5) {
      YCHECK;
      y[ix++]=(e.n+nu_e.n)/nB-YLe;
    }
    if (YLmu>-0.5) {
      YCHECK;
      y[ix++]=(mu.n+nu_mu.n)/nB-YLmu;
    }
    
    return;
  }
  
  /** \brief Solve for the hadronic phase at \f$ T=0 \f$ or \f$
      s=\mathrm{constant} \f$
  */
  int f_had_phase(size_t nv, const ubvector &x, ubvector &y,
		  double &nB) {

    //cout << "x: ";
    //vector_out(cout,nv,x,true);
    
    size_t ix=0;
    XCHECK;
    n.n=x[ix++];
    p.n=nB-n.n;
    
    double T=0.0;
    if (sonB>0.0) {
      XCHECK;
      T=x[ix++];
    }
    leptons_in(ix,x);
    
    int hret=ptr_h->calc_temp_e(n,p,T,hth);
    
    e.mu=n.mu-p.mu+nu_e.mu;

    ix=0;
    leptons_out(ix,y,T,nB);
    
    tot=hth+lep;
    
    YCHECK;
    y[ix++]=e.n+mu.n-p.n;
    if (sonB>0.0) {
      YCHECK;
      y[ix++]=sonB-tot.en/nB;
    }
    
    //cout << "y: ";
    //vector_out(cout,nv,y,true);
    //exit(-1);

    return 0;
  }
  
  /** \brief Solve for the hadronic phase with hyperons
   */
  int f_had_phase_hyp(size_t nv, const ubvector &x, ubvector &y,
		      double &nB) {

    size_t ix=0;
    n.mu=x[ix++];
    e.mu=x[ix++];
    double sigma=x[ix++];
    double omega=x[ix++];
    double rho=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      T=x[ix++];
    }
    leptons_in(ix,x);

    p.mu=n.mu-e.mu+nu_e.mu;
    sigm.mu=n.mu+e.mu-nu_e.mu;
    sigz.mu=n.mu;
    sigp.mu=p.mu;
    casm.mu=sigm.mu;
    casz.mu=sigz.mu;

    mu.mu=e.mu-nu_e.mu+nu_mu.mu;

    double f1, f2, f3;
    rmf_hyp.calc_eq_hyp_p(n,p,lam,sigp,sigz,sigm,casz,casm,
			  sigma,omega,rho,f1,f2,f3,hth);

    ix=0;
    leptons_out(ix,y,T,nB);
    
    tot=hth+lep;
    
    y[ix++]=n.n+p.n+lam.n+sigp.n+sigz.n+sigm.n+casz.n+casm.n-nB;
    y[ix++]=p.n+sigp.n-sigm.n-casm.n-e.n-mu.n;
    y[ix++]=f1;
    y[ix++]=f2;
    y[ix++]=f3;
    
    if (sonB>0.0) {
      y[ix++]=sonB-tot.en/nB;
    }
    
    return 0;
  }
  
  /** \brief Solve for the mixed phase at fixed chi
      (with no Coulomb or surface)
   */
  int f_mixed_phase_chi(size_t nv, const ubvector &x, ubvector &y,
                        const double &chi) {

    //cout.precision(10);
    //cout << "x: ";
    //vector_out(cout,nv,x,true);
    
    size_t ix=0;
    XCHECK;
    n.n=x[ix++];
    XCHECK;
    p.n=x[ix++];
    XCHECK;
    double nB=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      XCHECK;
      T=x[ix++];
    }
    // Add neutrino chemical potential variables
    leptons_in(ix,x);

    ptr_h->calc_temp_e(n,p,T,hth);
    
    e.mu=n.mu-p.mu;

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    
    ptr_q->calc_temp_p(u,d,s,T,qth);
    
    double quark_nqch, quark_nQ;
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;

    // Use the electron and neutrino chemical potentials to compute
    // the lepton thermodynamics, store in lep, and add y values for
    // electron or muon lepton number
    ix=0;
    leptons_out(ix,y,T,nB);
    
    // Update the total thermodynamics
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+lep.ed;
    tot.pr=hth.pr+lep.pr;
    if (sonB>0.0) {
      tot.en=hth.en*chi+qth.en*(1.0-chi)+lep.en;
    } else {
      tot.en=0.0;
    }

    // pressure equality
    YCHECK;
    y[ix++]=hth.pr-qth.pr;
    // baryon number 
    YCHECK;
    y[ix++]=nB-(n.n+p.n)*chi-quark_nQ*(1.0-chi)/3.0;
    // charge
    YCHECK;
    y[ix++]=chi*p.n+(1.0-chi)*quark_nqch-e.n-mu.n;

    if (sonB>0.0) {
      y[ix++]=sonB-tot.en/nB;
    }

    //cout << "y: ";
    //vector_out(cout,nv,y,true);
    
    return 0;
  }

  /** \brief Solve for the mixed phase (with no Coulomb or surface)
   */
  int f_mixed_phase_Ye(size_t nv, const ubvector &x, ubvector &y,
                       double nB, double Ye, double &chi) {
    
    size_t ix=0;
    XCHECK;
    n.n=x[ix++];
    XCHECK;
    p.n=x[ix++];
    XCHECK;
    e.mu=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      XCHECK;
      T=x[ix++];
    }
    // Add neutrino chemical potential variables
    leptons_in(ix,x);
    
    ptr_h->calc_temp_e(n,p,T,hth);
    
    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    
    ptr_q->calc_temp_p(u,d,s,T,qth);
    
    double quark_nqch, quark_nQ;
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;
      
    chi=(e.n+mu.n-quark_nqch)/(p.n-quark_nqch);

    // Use the electron and neutrino chemical potentials to compute
    // the lepton thermodynamics, store in lep, and add y values for
    // electron or muon lepton number
    ix=0;
    leptons_out(ix,y,T,nB);
    
    // Compute the total thermodynamics
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+lep.ed;
    tot.pr=hth.pr+lep.pr;
    if (sonB>0.0) {
      tot.en=hth.en*chi+qth.en*(1.0-chi)+lep.en;
    } else {
      tot.en=0.0;
    }

    YCHECK;
    y[ix++]=hth.pr-qth.pr;
    YCHECK;
    y[ix++]=nB-(n.n+p.n)*chi-quark_nQ*(1.0-chi)/3.0;
    YCHECK;
    y[ix++]=e.n/nB-Ye;

    if (sonB>0.0) {
      YCHECK;
      y[ix++]=sonB-tot.en/nB;
    }
    
    return 0;
  }

  /** \brief Solve for the quark phase
   */
  int f_quark_phase(size_t nv, const ubvector &x, ubvector &y,
		    double &nB) {

    //cout << "x: ";
    //vector_out(cout,nv,x,true);
    
    size_t ix=0;
    XCHECK;
    double muQ=x[ix++];
    XCHECK;
    e.mu=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      XCHECK;
      T=x[ix++];
    }
    
    leptons_in(ix,x);
    
    u.mu=muQ-e.mu*2.0/3.0;
    d.mu=muQ+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_temp_p(u,d,s,T,qth);

    ix=0;
    leptons_out(ix,y,T,nB);
    
    tot=hth+e+mu+qth;

    YCHECK;
    y[ix++]=(2.0*u.n-d.n-s.n)/3.0-e.n-mu.n;
    YCHECK;
    y[ix++]=(u.n+d.n+s.n)/3.0-nB;
    
    if (sonB>0.0) {
      YCHECK;
      y[ix++]=sonB-tot.en/nB;
    }
    
    //cout << "y: ";
    //vector_out(cout,nv,y,true);
    
    return 0;
  }

  /** \brief Solve for the beginning of the mixed phase
   */
  int f_beg_mixed_phase(size_t nv, const ubvector &x, ubvector &y,
			double &nB) {

    //cout << "x: ";
    //vector_out(cout,nv,x,true);
    
    size_t ix=0;
    XCHECK;
    n.n=x[ix++];
    XCHECK;
    p.n=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      XCHECK;
      T=x[ix++];
    }
    // Add neutrino chemical potential variables
    leptons_in(ix,x);

    ptr_h->calc_temp_e(n,p,T,hth);

    // Compute total baryon density
    nB=n.n+p.n;

    // Compute electron and quark chemical potentials
    e.mu=n.mu-p.mu+nu_e.mu;
    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_temp_p(u,d,s,T,qth);

    // Use the electron and neutrino chemical potentials to compute
    // the lepton thermodynamics, store in lep, and add y values for
    // electron or muon lepton number
    ix=0;
    leptons_out(ix,y,T,nB);
    
    YCHECK;
    y[ix++]=hth.pr-qth.pr;
    // Ensure electric neutrality in the hadronic phase
    YCHECK;
    y[ix++]=p.n-e.n-mu.n;
    if (sonB>0.0) {
      tot.en=hth.en+e.en+mu.en;
      if (nu_e.mu!=0.0) tot.en+=nu_e.en;
      if (nu_mu.mu!=0.0) tot.en+=nu_mu.en;
      YCHECK;
      y[ix++]=sonB-tot.en/nB;
    }

    cout << "bmp: " << n.n << " " << p.n << " " << u.n << " "
         << d.n << " " << s.n << endl;
    //cout << "y: ";
    //vector_out(cout,nv,y,true);

    return 0;
  }
  
  /** \brief Solve for the beginning of the mixed phase at fixed Ye
   */
  int f_beg_mixed_phase_Ye(size_t nv, const ubvector &x, ubvector &y,
                           double Ye, double &nB) {
    
    size_t ix=0;
    n.n=x[ix++];
    p.n=x[ix++];
    e.mu=x[ix++];
    
    double T=0.0;
    if (sonB>0.0) {
      T=x[ix++];
    }
    leptons_in(ix,x);

    ptr_h->calc_temp_e(n,p,T,hth);

    // Compute total baryon density
    nB=n.n+p.n;

    nu_e.mu=n.mu-p.mu-e.mu;
    mu.mu=e.mu-nu_e.mu+nu_mu.mu;

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    ptr_q->calc_temp_p(u,d,s,T,qth);
    
    ix=0;
    leptons_out(ix,y,T,nB);

    y[ix++]=hth.pr-qth.pr;
    // Ensure electric neutrality in the hadronic phase
    y[ix++]=p.n-e.n-mu.n;
    y[ix++]=e.n/nB-Ye;
    if (sonB>0.0) {
      tot.en=hth.en+e.en+mu.en;
      if (nu_e.mu!=0.0) tot.en+=nu_e.en;
      if (nu_mu.mu!=0.0) tot.en+=nu_mu.en;
      y[ix++]=sonB-tot.en/nB;
    }

    //cout << "y: ";
    //vector_out(cout,nv,y,true);

    return 0;
  }
  
  /** \brief Minimize the free energy with respect to the densities
   */
  int f_min_densities(size_t nv, const ubvector &x, ubvector &y,
		      double &nB) {
    
    double chi=x[0];
    p.n=x[1];

    double T=0.0;
    if (sonB>0.0) {
      T=x[2];
    }

    if (chi<0.0 || p.n<0.0 || n.n<0.0) return 1;

    if (false) {
      double t1, t2;
      sk.eff_mass(n,p,t1,t2);
    }
    
    if (n.ms<0.0 || p.ms<0.0) return 2;

    ptr_h->calc_temp_e(n,p,T,hth);
    e.mu=n.mu-p.mu;
    fr.calc_mu(e,T);
    mu.mu=e.mu;
    fr.calc_mu(mu,T);
    
    double quark_nqch, quark_nQ;

    u.mu=n.mu/3.0-e.mu*2.0/3.0;
    d.mu=n.mu/3.0+e.mu/3.0;
    s.mu=d.mu;
    
    ptr_q->calc_temp_p(u,d,s,T,qth);
    quark_nqch=(2.0*u.n-d.n-s.n)/3.0;
    quark_nQ=u.n+d.n+s.n;
    
    y[0]=e.n+mu.n-p.n*chi-quark_nqch*(1.0-chi);
    y[1]=nB-(n.n+p.n)*chi-quark_nQ*(1.0-chi)/3.0;
    return 0;
  }

  /** \brief Construct the mixed phase by minimizing the free energy
   */
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
    mh.ntrial=1000;
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
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+lep.ed;
    tot.pr=hth.pr+lep.pr;

    if (ret!=0 || chi>1.0) tot.ed+=1.0e6;
    
    return tot.ed;
  }

  /** \brief The Coulomb energy function
   */
  double fcoul(double dim, double chi) {
    double ret;
    if (dim==2.0) {
      return 0.25*(chi-1.0-log(chi));
    }
    ret=(2.0/(dim-2.0)*(1.0-0.5*dim*pow(chi,1.0-2.0/dim))+chi)/(dim+2.0);
    return ret;
  }

  /** \brief Construct the mixed phase by minimizing the free energy
      and minimize over the size of the nucleus
  */
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
      esurf=dim*xglen*surften/r_rare;
      ecoul=2.0*pi*pow(p.n-quark_nqch,2.0)*
	o2scl_const::fine_structure_f<double>()*r_rare*r_rare*
        xglen*fcoul(dim,xglen);
    } else {
      esurf=dim*chi*surften/r_rare;
      ecoul=2.0*pi*pow(p.n-quark_nqch,2.0)*
	o2scl_const::fine_structure_f<double>()*r_rare*r_rare*
        chi*fcoul(dim,xglen);
    }

    // Update the energy density and pressure
    tot.ed=hth.ed*chi+qth.ed*(1.0-chi)+lep.ed+esurf+ecoul;

    // Use the thermodynamic identity for the pressure
    tot.pr=-tot.ed+n.n*chi*n.mu+p.n*chi*p.mu+
      quark_nQ*(1.0-chi)*n.mu/3.0;

    //cout << 1.0/cbrt(4.0*pi*fabs(p.n*p.n-quark_nqch*quark_nqch)*
    //o2scl_const::fine_structure_f<double>()/
    //		       surften/dim*fcoul(dim,chi)) << endl;
    /*
      cout << "z: " << x[0] << " " << x[1] << " "
      << fabs(p.n*p.n-quark_nqch*quark_nqch) << " "
      << fcoul(dim,chi) << " " 
      << tot.ed << endl;
    */
    if (chi<0.0 || chi>1.0 || sret!=0 ||
        r_rare>50.0 || r_rare<0.0) tot.ed+=1e6;
    return tot.ed;
  }
  
  ex_eos_gibbs() {

    skyrme_load(sk,"NRAPR");
    ptr_h=&sk;
    bag.bag_constant=200.0/hc_mev_fm;
    ptr_q=&bag;

    // Nucleon initialization
    n.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_const::mass_neutron_f<double>()),2.0);
    p.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_const::mass_proton_f<double>()),2.0);
    n.non_interacting=false;
    p.non_interacting=false;
    n.inc_rest_mass=true;
    p.inc_rest_mass=true;

    // Lepton masses
    e.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_const::mass_electron_f<double>()),2.0);
    mu.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_const::mass_muon_f<double>()),2.0);

    // For the neutrinos, we set g=1, since we will include
    // antineutrinos using the fermion_rel::pair_mu() functions which
    // thus already includes both spin states.
    nu_e.init(0.0,1.0);
    nu_mu.init(0.0,1.0);

    // Default quark masses
    u.init(mass_up_MeV_f<double>()/hc_mev_fm,6.0);
    d.init(0.0,6.0);
    s.init(95.0/hc_mev_fm,6.0);
    u.inc_rest_mass=true;
    d.inc_rest_mass=true;
    s.inc_rest_mass=true;

    // Hyperon initialization
    lam.init(mass_lambda_MeV_f<double>()/hc_mev_fm,2.0);
    sigp.init(mass_sigma_plus_MeV_f<double>()/hc_mev_fm,2.0);
    sigz.init(mass_sigma_zero_MeV_f<double>()/hc_mev_fm,2.0);
    sigm.init(mass_sigma_minus_MeV_f<double>()/hc_mev_fm,2.0);
    casz.init(mass_cascade_zero_MeV_f<double>()/hc_mev_fm,2.0);
    casm.init(mass_cascade_minus_MeV_f<double>()/hc_mev_fm,2.0);
    lam.non_interacting=false;
    sigp.non_interacting=false;
    sigz.non_interacting=false;
    sigm.non_interacting=false;
    casz.non_interacting=false;
    casm.non_interacting=false;
    lam.inc_rest_mass=true;
    sigp.inc_rest_mass=true;
    sigz.inc_rest_mass=true;
    sigm.inc_rest_mass=true;
    casz.inc_rest_mass=true;
    casm.inc_rest_mass=true;
    
    surften=1.0/hc_mev_fm;
    
    mp_start_fix=0.0;

    sonB=0.0;
    YLe=-1.0;
    YLmu=-1.0;

    had_phase_guess.resize(5);
    had_phase_guess[0]=0.24;
    had_phase_guess[1]=0.16;
    had_phase_guess[2]=0.05;
    had_phase_guess[3]=0.05;
    had_phase_guess[4]=0.05;

    beg_mixed_phase_guess.resize(5);
    beg_mixed_phase_guess[0]=0.16;
    beg_mixed_phase_guess[1]=0.08;
    beg_mixed_phase_guess[2]=0.07;
    beg_mixed_phase_guess[3]=0.95;
    beg_mixed_phase_guess[4]=-0.06;

    quark_phase_guess.resize(5);
  }

  /** \brief Desc
   */
  int guess(vector<string> &sv, bool itive_com) {
    if (sv.size()>3 && sv[1]=="beg-mixed") {
      size_t nv=sv.size()-2;
      beg_mixed_phase_guess.resize(nv);
      for(size_t i=2;i<sv.size();i++) {
        beg_mixed_phase_guess[i-2]=function_to_double(sv[i]);
      }
      cout << "Guess for beginning of mixed phase is: ";
      vector_out(cout,beg_mixed_phase_guess,true);
      return 0;
    }
    cerr << "Invalid input in guess." << endl;
    return 0;
  }
  
  /** \brief Compute the EOS and M-R curve at T=0
   */
  int mvsr(vector<string> &sv, bool itive_com) {

    if (sv.size()<2) {
      cerr << "File not specified in mvsr()." << endl;
      return 1;
    }
    
    double mp_end, chi, nB, dim=3.0, esurf, ecoul, mp_start;
    
    mm_funct fp_had_phase=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,double &)>
       (&ex_eos_gibbs::f_had_phase),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::ref(nB));
    mm_funct fp_mixed_phase_chi=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,const double &)>
       (&ex_eos_gibbs::f_mixed_phase_chi),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::cref(chi));
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

    ubvector x(8), y(8);
    
    cout << "Masses (n,p,e): " << n.m*hc_mev_fm << " "
	 << p.m*hc_mev_fm << " "
         << e.m*hc_mev_fm << " MeV" << endl;
    cout << endl;

    /*
      cout << rmf.b << endl;
      cout << rmf.c << endl;
      cout << rmf.cr << endl;
      cout << rmf.ms << endl;
      cout << rmf.zeta << endl;
      cout << rmf.cs << endl;
      cout << rmf.cw << endl;
      cout << rmf.b1 << endl;
      cout << rmf.xi << endl;
      cout << rmf.a1 << endl;
      cout << rmf.a2 << endl;
      cout << rmf.a3 << endl;
      cout << rmf.a4 << endl;
      cout << rmf.a5 << endl;
      cout << rmf.a6 << endl;
      cout << rmf.b2 << endl;
      cout << rmf.b3 << endl;
      exit(-1);
    */

    // ───────────────────────────────────────────────────────────────────
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
	   << p.mu*hc_mev_fm << " "
           << e.mu*hc_mev_fm << " MeV" << endl;
      cout << "Quark chem. pots. (u,d,s): "
           << u.mu*hc_mev_fm << " "
	   << d.mu*hc_mev_fm << " "
           << s.mu*hc_mev_fm << " fm^{-3}" << endl;
      cout << "Electron energy density: "
           << e.ed*hc_mev_fm << " MeV/fm^3"
	   << endl;
      cout << "Total pressure: "
           << tot.pr*hc_mev_fm << " MeV/fm^3" << endl;
      cout << "Total energy density: "
           << tot.ed*hc_mev_fm
	   << " MeV/fm^3" << endl;
      cout << "Hadronic energy density: "
           << hth.ed*hc_mev_fm << " MeV/fm^3"
           << endl;
      cout << "Hadronic pressure: "
           << hth.pr*hc_mev_fm << " MeV/fm^3"
           << endl;
      cout << "Quark energy dens. & pressure: "
           << qth.ed*hc_mev_fm << " MeV/fm^3, "
           << qth.pr*hc_mev_fm << " MeV/fm^3"
           << endl;
      cout << endl;

      beg_mixed_phase_guess[0]=n.n;
      beg_mixed_phase_guess[1]=p.n;
      beg_mixed_phase_guess[2]=0.05;
      beg_mixed_phase_guess[3]=0.05;
      beg_mixed_phase_guess[4]=0.05;

    } else {

      // ───────────────────────────────────────────────────────────────────
      // Find the beginning of the mixed phase with fp_beg_mixed_phase()
      
      cout << "Beginning of mixed phase: " << endl;
      size_t ix=0;
      XCHECK;
      x[ix++]=beg_mixed_phase_guess[0];
      XCHECK;
      x[ix++]=beg_mixed_phase_guess[1];
      if (sonB>0.0) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[2];
      }
      if (YLe>-0.5) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[3];
      }
      if (YLmu>-0.5) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[4];
      }
      size_t nvar=ix;

      mh.msolve(nvar,x,fp_beg_mixed_phase);
      
      f_beg_mixed_phase(nvar,x,y,nB);
      mp_start=nB;
      cout << "Mixed phase begins at nB: " << nB << " fm^{-3}" << endl;

      // Copy successful results back over to beg_mixed_phase_guess
      ix=0;
      XCHECK;
      beg_mixed_phase_guess[0]=x[ix++];
      XCHECK;
      beg_mixed_phase_guess[1]=x[ix++];
      if (sonB>0.0) {
        XCHECK;
        beg_mixed_phase_guess[2]=x[ix++];
      }
      if (YLe>-0.5) {
        XCHECK;
	beg_mixed_phase_guess[3]=x[ix++];
      }
      if (YLmu>-0.5) {
        XCHECK;
	beg_mixed_phase_guess[4]=x[ix++];
      }

      cout << endl;
      
    }

    // ───────────────────────────────────────────────────────────────────
    // Tabulate full hadronic phase

    table_units<> thad;
    if (true) {
      thad.line_of_names("nB nn np ne nmu ede edh ed pre prh pr");
      thad.line_of_units(((string)"1/fm^3 1/fm^3 1/fm^3 1/fm^3 ")+
                         "1/fm^3 MeV/fm^3 MeV/fm^3 MeV/fm^3 "+
                         "MeV/fm^3 MeV/fm^3 MeV/fm^3 ");
      size_t ix=0;
      
      x[ix++]=beg_mixed_phase_guess[0];
      if (sonB>0.0) {
	x[ix++]=beg_mixed_phase_guess[2];
      }
      if (YLe>-0.5) {
	x[ix++]=beg_mixed_phase_guess[3];
      }
      if (YLmu>-0.5) {
	x[ix++]=beg_mixed_phase_guess[4];
      }
      size_t nvar=ix;
      
      double nb_had_last=((double)((int)(mp_start*100.0)))/100.0;
                                         
      cout << "nb_had_last is: " << nb_had_last << endl;

      cout << "Hadronic phase (nB,nn,np):" << endl;
      for(nB=nb_had_last;nB>0.07999;nB-=0.01) {
	mh.msolve(nvar,x,fp_had_phase);
	cout << nB << " " << n.n << " " << p.n << endl;
	std::vector<double> line=
          {nB,n.n,p.n,e.n,mu.n,e.ed*hc_mev_fm,
           hth.ed*hc_mev_fm,
           (hth.ed+lep.ed)*hc_mev_fm,
           e.pr*hc_mev_fm,hth.pr*hc_mev_fm,
           (hth.pr+lep.pr)*hc_mev_fm};
	thad.line_of_data(line);
      }
      cout << endl;

      thad.sort_table("nB");
    }

    // ───────────────────────────────────────────────────────────────────
    // Tabulate full mixed phase without Coulomb or surface

    table_units<> tmixed;
    if (true) {
      cout << "Mixed phase (nB,nn,np,nu,nd,ns,chi): " << endl;
      tmixed.line_of_names(((string)"nB nn np nu nd ns ede pre ed ")+
			   "pr edh prh edq prq mun mup chi");
      tmixed.line_of_units(((string)"1/fm^3 1/fm^3 1/fm^3 1/fm^3 ")+
                         "1/fm^3 1/fm^3 MeV/fm^3 MeV/fm^3 MeV/fm^3 "+
                         "MeV/fm^3 MeV/fm^3 MeV/fm^3 MeV/fm^3 "+
                         "MeV/fm^3 MeV MeV .");

      size_t ix=0;
      XCHECK;
      // neutron number density
      x[ix++]=beg_mixed_phase_guess[0];
      XCHECK;
      // proton number density
      x[ix++]=beg_mixed_phase_guess[1];
      XCHECK;
      // baryon number density
      x[ix++]=beg_mixed_phase_guess[0]+beg_mixed_phase_guess[1];
      if (sonB>0.0) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[2];
      }
      if (YLe>-0.5) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[3];
      }
      if (YLmu>-0.5) {
        XCHECK;
	x[ix++]=beg_mixed_phase_guess[4];
      }
      size_t nvar=ix;
      
      double mp_end_muQ=1.0;
      double mp_end_mue=0.1;
      
      for (chi=1.0;chi>-1.0e-4;chi-=0.05) {

	mh.msolve(nvar,x,fp_mixed_phase_chi);
        nB=x[2];

        mp_end=nB;
        mp_end_muQ=(2.0*d.mu+u.mu)/3.0;
        mp_end_mue=d.mu-u.mu;

        quark_phase_guess[0]=mp_end_muQ;
        quark_phase_guess[1]=mp_end_mue;
        if (sonB>0.0) {
          quark_phase_guess[2]=x[3];
        }
        quark_phase_guess[3]=nu_e.mu;
        quark_phase_guess[4]=nu_mu.mu;

        cout.precision(4);
	cout << nB << " " << n.n << " " << p.n << " " << u.n << " "
             << d.n << " " << s.n << " " << chi << endl;
        //cout << "  " << e.n << " " << mu.n << " "
        //<< nu_e.n << " " << nu_mu.n << endl;
        cout.precision(6);
        
	std::vector<double> line;
	line.push_back(nB);
	line.push_back(n.n);
	line.push_back(p.n);
	line.push_back(u.n);
	line.push_back(d.n);
	line.push_back(s.n);
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
	tmixed.line_of_data(line.size(),line);
      }
      
      cout << endl;
    }

    // ───────────────────────────────────────────────────────────────────
    // Tabulate full quark phase

    table_units<> tq;
    if (true) {
      tq.line_of_names("nB muQ mue ede pre edq prq ed pr");
      tq.line_of_units(((string)"1/fm^3 MeV MeV ")+
                         "MeV/fm^3 MeV/fm^3 MeV/fm^3 "+
                         "MeV/fm^3 MeV/fm^3 MeV/fm^3 ");

      double nb_had_next=((double)((int)(mp_end*100.0)))/100.0+0.01;
      cout << "nb_had_next: " << nb_had_next << endl;
      
      cout << "Quark phase (nB,nu,nd,ns): " << endl;
      size_t ix=0;
      XCHECK;
      //x[ix++]=mp_end_muQ;
      x[ix++]=quark_phase_guess[0];
      XCHECK;
      //x[ix++]=mp_end_mue;
      x[ix++]=quark_phase_guess[1];
      if (sonB>0.0) {
        XCHECK;
	x[ix++]=quark_phase_guess[2];
      }
      if (YLe>-0.5) {
        XCHECK;
	x[ix++]=quark_phase_guess[3];
      }
      if (YLmu>-0.5) {
        XCHECK;
	x[ix++]=quark_phase_guess[4];
      }
      size_t nvar=ix;
      
      for(nB=nb_had_next;nB<1.5001;nB+=0.01) {
	mh.msolve(nvar,x,fp_quark_phase);
	cout << nB << " " << u.n << " " << d.n << " " << s.n << endl;
	std::vector<double> line=
          {nB,x[0],x[1],e.ed*hc_mev_fm,
           e.pr*hc_mev_fm,qth.ed*hc_mev_fm,
           qth.pr*hc_mev_fm,
           (qth.ed+lep.ed)*hc_mev_fm,
           (qth.pr+lep.pr)*hc_mev_fm};
	tq.line_of_data(line);
      }
      cout << endl;
    }

    // ───────────────────────────────────────────────────────────────────
    // Construct the neutron star EOS
    
    table_units<> ns;
    ns.line_of_names("ed pr nB");
    ns.set_unit("ed","MeV/fm^3");
    ns.set_unit("pr","MeV/fm^3");
    ns.set_unit("nB","1/fm^3");

    for(size_t i=0;i<thad.get_nlines();i++) {
      double line[3]={thad.get("ed",i),
                      thad.get("pr",i),
                      thad.get("nB",i)
      };
      ns.line_of_data(3,line);
    }
    for(size_t i=0;i<tmixed.get_nlines();i++) {
      double line[3]={tmixed.get("ed",i),
                      tmixed.get("pr",i),
                      tmixed.get("nB",i)};
      ns.line_of_data(3,line);
    }
    for(size_t i=0;i<tq.get_nlines();i++) {
      double line[3]={tq.get("ed",i),
                      tq.get("pr",i),
                      tq.get("nB",i)};
      ns.line_of_data(3,line);
    }

    // ───────────────────────────────────────────────────────────────────
    // Solve the TOV equations
    
    eos_tov_interp eti;
    eti.default_low_dens_eos();
    eti.read_table(ns,"ed","pr","nB");
    cout << "Going to tov_solve." << endl;
    tov_solve ts;

    // This additional accuracy is necessary to get the plots
    // right
    ts.princ=1.01;
    ts.set_eos(eti);
    ts.mvsr();

    // The tov_solve class uses 'nb', so we just rename
    std::shared_ptr<table_units<> > tov=ts.get_results();
    tov->rename_column("nb","nB");
    cout << "M_max: " << tov->max("gm") << endl;
    cout << "Central energy density of maximum mass star:\n  "
         << tov->get("ed",tov->lookup("gm",tov->max("gm"))) << " "
         << tov->get_unit("ed") << " ." << endl;
    cout << "Central baryon density of maximum mass star:\n  "
         << tov->get("nB",tov->lookup("gm",tov->max("gm"))) << " "
         << tov->get_unit("nB") << " ." << endl;

    // ───────────────────────────────────────────────────────────────────
    // "Insert" the TOV table into the tables for the individual
    // phases.

    thad.insert_table_units(*tov,"nB");
    tmixed.insert_table_units(*tov,"nB");
    tq.insert_table_units(*tov,"nB");
    
    // ───────────────────────────────────────────────────────────────────
    // Output results to a file

    hdf_file hf;
    hf.open_or_create(sv[1]);
    hdf_output(hf,thad,"hadrons");
    hdf_output(hf,tmixed,"mixed");
    hdf_output(hf,tq,"quarks");
    hdf_output(hf,ns,"nstar");
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
    } else if (sv[1]=="SPL00") {
      ptr_h=&rmf;

      n.m=939.0/hc_mev_fm;
      p.m=939.0/hc_mev_fm;
      
      rmf.set_n_and_p(n,p);
      
      rmf.mnuc=939.0/hc_mev_fm;
      rmf.n0=0.16;
      rmf.eoa=16.0/hc_mev_fm;
      rmf.comp=250.0/hc_mev_fm;
      rmf.msom=0.6;
      rmf.esym=35.0/hc_mev_fm;
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
      cout << "n0,m*/m: " << rmf.n0 << " " << rmf.msom << " "
           << rmf.eoa*hc_mev_fm << " " << rmf.comp*hc_mev_fm << " "
           << rmf.esym*hc_mev_fm << endl;
      
      cout << "Selected the RMF model from SPL00 "
           << "(no hyperons, zeta=xi=0)." << endl;
      
    } else if (sv[1]=="SPL00_hyp") {
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
      u.non_interacting=true;
      d.non_interacting=true;
      s.non_interacting=true;
    } else if (sv[1]=="SPL00_bag") {
      ptr_q=&bag;
      u.m=5.5/hc_mev_fm;
      d.m=5.5/hc_mev_fm;
      s.m=140.7/hc_mev_fm;
      u.non_interacting=true;
      d.non_interacting=true;
      s.non_interacting=true;
      bag.bag_constant=200.0/hc_mev_fm;
      mp_start_fix=0.0;
      cout << "Selected the bag model from SPL00." << endl;
    } else if (sv[1]=="SPL00_njl") {
      ptr_q=&njl;
      u.m=5.5/hc_mev_fm;
      d.m=5.5/hc_mev_fm;
      s.m=140.7/hc_mev_fm;
      u.non_interacting=false;
      d.non_interacting=false;
      s.non_interacting=false;
      njl.set_quarks(u,d,s);
      njl.up_default_mass=5.5/hc_mev_fm;
      njl.down_default_mass=5.5/hc_mev_fm;
      njl.strange_default_mass=140.7/hc_mev_fm;
      double L=602.3/hc_mev_fm;
      njl.set_parameters(L,1.835/L/L,12.36/pow(L,5.0));
      mp_start_fix=0.32;
    } else if (sv[1]=="njl") {
      ptr_q=&njl;
    } else if (sv[1]=="none") {
      ptr_q=0;
    } else {
      cerr << "Quark model specification " << sv[1]
	   << " not understood." << endl;
      return 2;
    }
    
    return 0;
  }

  /** \brief Desc
   */
  int spl00(vector<string> &sv, bool itive_com) {

    {
      vector<string> sv2={"hadrons","SPL00"};
      model_hadrons(sv2,itive_com);
    }
    {
      vector<string> sv2={"quarks","SPL00_bag"};
      model_quarks(sv2,itive_com);
    }
    {
      vector<string> sv2={"mvsr","ex_eos_gibbs_spl00.o2"};
      mvsr(sv2,itive_com);
    }

    sonB=1.0;
    YLe=0.4;
    YLmu=0.0;
    {
      vector<string> sv2={"mvsr","ex_eos_gibbs_spl00.o2"};
      mvsr(sv2,itive_com);
    }
    
    sonB=2.0;
    YLe=-1.0;
    YLmu=-1.0;
    {
      vector<string> sv2={"mvsr","ex_eos_gibbs_spl00.o2"};
      mvsr(sv2,itive_com);
    }
    
    return 0;
  }

  /** \brief Test
   */
  int test(vector<string> &sv, bool itive_com) {
    test_mgr t;
    t.report();
    return 0;
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  ex_eos_gibbs ehg;
  
  cli cl;
  cl.prompt="ex_eos_gibbs>";
  int comm_option_both=2;
  
  static const int narr=6;
  comm_option_s options_arr[narr]={
    {0,"hadrons","Select hadronic model.",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::model_hadrons),
     cli::comm_option_both},
    {0,"quarks","Select quark model.",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::model_quarks),
     cli::comm_option_both},
    {0,"mvsr","Compute EOS and M-R curve",1,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::mvsr),
     cli::comm_option_both},
    {0,"spl00","Desc",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::spl00),
     cli::comm_option_both},
    {0,"test","Desc",0,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::test),
     cli::comm_option_both},
    {0,"guess","Desc",-1,-1,"","",
     new comm_option_mfptr<ex_eos_gibbs>(&ehg,&ex_eos_gibbs::guess),
     cli::comm_option_both}
  };
  
  cl.set_comm_option_vec(narr,options_arr);
  cl.cmd_name="ex_eos_gibbs";

  cli::parameter_double p_sonB;
  p_sonB.d=&ehg.sonB;
  p_sonB.help="The entropy per baryon (default 0.0)";
  cl.par_list.insert(make_pair("sonB",&p_sonB));

  cli::parameter_double p_YLe;
  p_YLe.d=&ehg.YLe;
  p_YLe.help="Electron lepton fraction (default is -1 for no neutrinos)";
  cl.par_list.insert(make_pair("YLe",&p_YLe));

  cli::parameter_double p_YLmu;
  p_YLmu.d=&ehg.YLmu;
  p_YLmu.help="Muon lepton fraction (default is -1 for no neutrinos)";
  cl.par_list.insert(make_pair("YLmu",&p_YLmu));

  cl.run_auto(argc,argv);
  
  return 0;
}
