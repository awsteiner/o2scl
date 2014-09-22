/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014, Andrew W. Steiner
  
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
/** \file eos_nse_full.h
    \brief File defining \ref o2scl::dense_matter, 
    \ref o2scl::nucmass_densmat, and \ref o2scl::eos_nse_full
*/
#ifndef EOS_NSE_FULL_H
#define EOS_NSE_FULL_H 

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>

#include <o2scl/classical.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>

#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/mroot_hybrids.h>

#include <o2scl/eos_had_skyrme.h>

#include <o2scl/mmin_simp2.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A container for the properties of dense matter at a 
      specified baryon density, electron fraction and temperature

      This class is experimental.

      By default, the rest mass is not included in the neutron,
      or the proton. It is, however, included in the electron. 
  */
  class dense_matter {
  
  public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
  
    /// Temperature
    double T;
    /// Baryon number density
    double nB;
    /// Electron fraction 
    double Ye;

    /// Neutrons
    o2scl::fermion n;
    /// Protons
    o2scl::fermion p;
    /// Electrons
    o2scl::fermion e;
    /// Muons
    o2scl::fermion mu;
    /// Photons
    o2scl::boson photon;
    /// Distribution of nuclei
    std::vector<o2scl::nucleus> dist;

    /// Neutron chemical potential
    double eta_n;
    /// Proton chemical potential
    double eta_p;
  
    /// Nuclear chemical potentials
    ubvector eta_nuc;

    /// Total thermodynamic quantities
    o2scl::thermo th;

    /// Thermodynamic quantities for dripped particles
    o2scl::thermo drip_th;

    /** \brief Constructor

	This constructor automatically sets the neutron, proton, and
	electron masses.
    */
    dense_matter();

    /// Copy constructor
    dense_matter(const dense_matter &dm) {
      n=dm.n;
      p=dm.p;
      e=dm.e;
      mu=dm.mu;
      photon=dm.photon;
      
      drip_th=dm.drip_th;
      th=dm.th;
      
      for (size_t i=0;i<dm.dist.size();i++) {
	dist.push_back(dm.dist[i]);
      }
      
      eta_n=dm.eta_n;
      eta_p=dm.eta_p;
      eta_nuc=dm.eta_nuc;
      
      T=dm.T;
      nB=dm.nB;
      Ye=dm.Ye;
    }

    /// Copy constructor with operator=()
    dense_matter &operator=(const dense_matter &dm) {

      if (this!=&dm) {
	n=dm.n;
	p=dm.p;
	e=dm.e;
	mu=dm.mu;
	photon=dm.photon;
	
	drip_th=dm.drip_th;
	th=dm.th;
	
	for (size_t i=0;i<dm.dist.size();i++) {
	  dist.push_back(dm.dist[i]);
	}
	
	eta_n=dm.eta_n;
	eta_p=dm.eta_p;
	eta_nuc=dm.eta_nuc;
	
	T=dm.T;
	nB=dm.nB;
	Ye=dm.Ye;
      }

      return *this;
    }

    /** \brief Compute an average inter-ionic spacing

	This function returns 
	\f[
	\left<a\right> \equiv \left(\frac{4 \pi}{3}\sum_i n_i
	\right)^{-1/3} \, .
	\f]
    */
    double average_a();

    /** \brief Compute the number-averaged mass number

	This function returns 
	\f$ \left<A\right> \equiv \sum_i n_i A_i / \sum_i n_i \f$ .
    */
    double average_A();

    /** \brief Compute the number-averaged neutron number

	This function returns 
	\f$ \left<N\right> \equiv \sum_i n_i N_i / \sum_i n_i \f$ .
    */
    double average_N();

    /** \brief Compute the number-averaged proton number

	This function returns 
	\f$ \left<Z\right> \equiv \sum_i n_i Z_i / \sum_i n_i \f$ .
    */
    double average_Z();

    /** \brief Compute the impurity parameter 

	This function returns the impurity parameter, 
	\f[
	\left<Q\right> \equiv \left[\sum_i n_i \left(Z_i-\left<Z\right>
	\right)^2 \right]
	\left(\sum_i n_i\right)^{-1}
	\f]
    */
    double impurity();

    /** \brief Compute the baryon density in nuclei
      
	This function returns \f$ \sum_i n_i A_i \f$ .
    */
    double baryon_density_nuclei();

    /** \brief Compute the total baryon density
      
	This function returns \f$ n_n + n_p + \sum_i n_i A_i \f$ .
    */
    double baryon_density();

    /** \brief Compute the electron fraction
      
	This function returns 
	\f[
	\frac{1}{n_B} \left(n_p + \sum_i Z_i n_i \right)
	\f]
	where \f$ n_B \f$ is the value returned by \ref baryon_density() .
    */
    double electron_fraction();

    /** \brief Return true if nucleus (Z,N) is in the distribution and
	store it's index in \c index
     */
    bool nuc_in_dist(int Z, int N, size_t &index) {
      for(size_t i=0;i<dist.size();i++) {
	if (dist[i].Z==Z && dist[i].N==N) {
	  index=i;
	  return true;
	}
      }
      return false;
    }
    
  };

  /** \brief A nuclear mass formula for dense matter
      
      This class is experimental.
  */
  class nucmass_densmat {

  protected:

    /** \brief Pointer to the nuclear mass formula (points to \ref ame 
	by default)
    */
    nucmass *massp;

  public:
    
    /// Return the type, \c "nucmass_densmat".
    virtual const char *type() { return "nucmass_densmat"; }
    
    nucmass_densmat();

    /// Default nuclear masses
    nucmass_ame_exp ame;
    
    /// Set base nuclear masses
    void set_mass(nucmass &nm);

    /** \brief Test the derivatives for 
	\ref binding_energy_densmat_derivs()
    */
    virtual void test_derivatives(double eps, double &t1, double &t2, 
				  double &t3, double &t4);

    /** \brief Compute the binding energy of a nucleus in dense matter
	and derivatives

	This function computes the binding energy of a nucleus in a
	sea of protons, neutrons, and negative charges (usually
	electrons) at a fixed temperature, relative to homogeneous
	nucleonic matter with the same number densities of protons,
	neutrons, and negative charges. The proton number Z and
	neutron number N should also be counted relative homogeneous
	nucleonic matter, not relative to the vacuum.

	As in \ref o2scl::nucmass::binding_energy_d(), the binding
	energy returned in \c E has units of MeV. All densities are
	expected to be in \f$ \mathrm{fm}^{-3} \f$, and the
	temperature should be in MeV. 

	\future Extend to negative N and Z?
    */
    virtual void binding_energy_densmat_derivs
      (double Z, double N, double npout, double nnout, 
       double nneg, double T, double &E, double &dEdnp, double &dEdnn,
       double &dEdnneg, double &dEdT);

  };

  /** \brief EOS for nuclear statistical equilibrium with interactions

      This class is experimental.

      \future Add positrons, muons, and anti-muons
      \future Add fermion and boson statistics to the nuclei in the
      distribution
  */
  class eos_nse_full {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Compute particle properties assuming classical thermodynamics
    o2scl::classical cla;

    /// The integer which indicates an invalid configuration
    int invalid_config;
    
    /// Relativistic fermions
    o2scl::fermion_rel relf;

    /// Relativistic fermions with derivatives
    o2scl::fermion_deriv_rel snf;

    /// Mass formula (points to \ref nuc_dens by default)
    o2scl::nucmass_densmat *massp;

    /// The default distribution
    std::vector<o2scl::nucleus> def_dist;

    /// The full distribution of all nuclei to consider
    std::vector<o2scl::nucleus> *ad;

    /** \brief Compute the free energy from a vector of densities 
	of the nuclei
    */
    double free_energy(const ubvector &n_nuc, dense_matter &dm);

    /// Nucleonic EOS (0 by default)
    o2scl::eos_had_temp_base *ehtp;

    /// Default solver
    mroot_hybrids<> def_mroot;

#endif

  public:

    eos_nse_full();
    
    /// The minimizer
    o2scl::mmin_simp2<> def_mmin;

    /// Compute nuclei in dense matter
    o2scl::nucmass_densmat nuc_dens;

    /// Default EOS ("SLy4")
    o2scl::eos_had_skyrme sk;

    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /// If true, include electrons and photons (default true)
    bool inc_lept_phot;

    /** \brief If true, include dripped protons in the Coulomb energy 
	(default true)
    */
    bool inc_prot_coul;
    
    /** \brief Set nuclear mass formula
     */
    void set_mass(o2scl::nucmass_densmat &m) {
      massp=&m;
      return;
    }

    /** \brief Set homogeneous matter EOS
     */
    void set_eos(o2scl::eos_had_temp_base &e) {
      ehtp=&e;
      return;
    }

    /** \brief Get homogeneous matter EOS
     */
    o2scl::eos_had_temp_base &get_eos() {
      return *ehtp;
    }

    /** \brief Set distribution of nuclei
     */
    void set_dist(std::vector<o2scl::nucleus> &dist) {
      ad=&dist;
      return;
    }
  
    /** \brief Compute the properties of matter from the densities,
	not presuming equilibrium
    */
    int calc_density_noneq(dense_matter &dm, int verbose=0);

    /** \brief Compute the properties of matter from 
	neutron and proton densities, using NSE
    */
    int calc_density_fixnp(dense_matter &dm, int verbose=0);
  
    /** \brief Compute the free energy for a fixed composition

	Given a fixed baryon density (dm.nB), electron fraction
	(dm.Ye), temperature (dm.T), this minimizes the free energy
	over the densities of the nuclei currently present in the
	distribution.

	\note This function only performs a very simple minimization
	and currently works in only limited circumstances.
    */
    int calc_density_fixcomp(dense_matter &dm, int verbose=0);

    /** \brief Output properties of a \ref dense_matter object to
	std::cout

	This function was particularly designed for comparing results
	with \ref o2scl::eos_sn_base derived classes.
    */
    void output(dense_matter &dm, int verbose=0);

    /** \brief Adjust the particle densities to match specified
	density and electron fraction
    */
    int density_match(dense_matter &dm);

    /** \brief Desc
     */
    int calc_density(dense_matter &dm, int verbose);

    /** \brief Desc
     */
    int solve_fixnp(size_t n, const ubvector &x, ubvector &y,
			      dense_matter &dm);

#ifdef O2SCL_NEVER_DEFINED

    /** \brief Desc

	This was intended to be a version of calc_density_fixcomp()
	which optimized the composition, but it doesn't really work
	yet.
     */
    int calc_density(dense_matter &dm, int verbose=0);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
