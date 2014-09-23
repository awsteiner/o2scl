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
    \brief File defining \ref o2scl::eos_nse_full
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
#include <o2scl/nucmass_densmat.h>
#include <o2scl/nucdist.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/mroot_hybrids.h>

#include <o2scl/eos_had_skyrme.h>

#include <o2scl/mmin_simp2.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief EOS for nuclear statistical equilibrium with interactions

      This class is experimental.

      Several of the functions have a verbose parameter. Generally,
      0 means no output, 1 means the function will output the
      composition and thermodynamics for the first 10 or so 
      nuclei in the distribution, and 2 means the function will
      output the entire distribution.

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

    /** \brief Function which is solved by \ref calc_density_saha()

	This function takes two inputs, the neutron and proton
	densities, and solves to ensure that \ref
	dense_matter::baryon_density() matches \ref
	o2scl::dense_matter::nB and \ref
	o2scl::dense_matter::electron_fraction() matches \ref
	dense_matter::Ye.

	This function calls \ref calc_density_fixnp() .
    */
    int solve_fixnp(size_t n, const ubvector &x, ubvector &y,
			      dense_matter &dm);

#endif

  public:

    eos_nse_full();
    
    /// The minimizer
    o2scl::mmin_simp2<> def_mmin;

    /// Default solver
    mroot_hybrids<> def_mroot;

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
	neutron and proton densities, using the Saha equation

	Given a fixed neutron and proton density, this computes their
	chemical potentials using the homogeneous matter EOS. Then the
	electrons are computed assuming their density is given from
	\ref o2scl::dense_matter::nB and \ref o2scl::dense_matter::Ye.
	Finally, the Saha equation is used to determine the nuclear
	chemical potentials and this gives the nuclear densities.

	This function only works when \ref inc_prot_coul is
	<tt>false</tt>.

	Note that, after this function completes, the value returned
	by \ref o2scl::dense_matter::baryon_density() will not necessarily be
	the same as that stored in \ref o2scl::dense_matter::nB (and
	similarly for the electron fraction). 
    */
    int calc_density_fixnp(dense_matter &dm, int verbose=0);
  
    /** \brief Compute the free energy for a fixed composition 
	by minimization

	Given a fixed baryon density (dm.nB), electron fraction
	(dm.Ye), temperature (dm.T), this minimizes the free energy
	over the densities of the nuclei currently present in the
	distribution. The neutron and proton drip densities are 
	determined by ensuring that the baryon density and electron
	fraction are correctly reproduced. The function which
	is minimized is \ref free_energy() .

	\note This function currently only performs a very simple
	minimization and currently works in only limited
	circumstances.
    */
    int calc_density_by_min(dense_matter &dm, int verbose=0);

    /** \brief Compute properties of matter for baryon density and
	electron fraction using the Saha equation

	This function solves the function specified by \ref
	solve_fixnp() using the current values of <tt>dm.n.n</tt> and
	<tt>dm.p.n</tt> as initial guesses.
     */
    int calc_density_saha(dense_matter &dm, int verbose);

    /** \brief Output properties of a \ref o2scl::dense_matter object to
	std::cout

	This function was particularly designed for comparing results
	with \ref o2scl::eos_sn_base derived classes.
    */
    void output(dense_matter &dm, int verbose=0);

    /** \brief Adjust the particle densities to match specified
	density and electron fraction
    */
    int density_match(dense_matter &dm);

#ifdef O2SCL_NEVER_DEFINED

    /** \brief Desc

	This was intended to be a version of calc_density_by_min()
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
