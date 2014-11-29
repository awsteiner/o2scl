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
#ifndef NUCMASS_DENSMAT_H
#define NUCMASS_DENSMAT_H

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/classical.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>

#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>

/** \file nucmass_densmat.h
    \brief File defining \ref o2scl::dense_matter and
    \ref o2scl::nucmass_densmat
*/

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
  
    /// Temperature (in \f$ \mathrm{fm}^{-1} \f$)
    double T;
    /// Baryon number density (in \f$ \mathrm{fm}^{-3} \f$)
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
	
	This function performs a simple brute-force search.
    */
    bool nuc_in_dist(int Z, int N, size_t &index);

    /** \brief Desc
     */
    void prune_distribution(double factor);

    /** \brief Desc
     */
    void copy_densities_from(dense_matter &dm2);
    
  };

  /** \brief A nuclear mass formula for dense matter
      
      This class is experimental.

      The default set of nuclear masses is from the AME 2012
      mass evaluation and is automatically loaded in the 
      constructor.
  */
  class nucmass_densmat {

  protected:
    
    /** \brief Pointer to the nuclear mass formula 
     */
    nucmass *massp;

  public:
    
    /// Return the type, \c "nucmass_densmat".
    virtual const char *type() { return "nucmass_densmat"; }
    
    nucmass_densmat();

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

    /** \brief Compute the binding energy of a nucleus in dense matter
     */
    virtual void binding_energy_densmat
      (double Z, double N, double npout, double nnout, 
       double nneg, double T, double &E);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

