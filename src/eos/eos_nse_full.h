/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014-2015, Andrew W. Steiner
  
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

#include <o2scl/nucmass_densmat.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/root_cern.h>

#include <o2scl/eos_had_skyrme.h>

#include <o2scl/mmin_simp2.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief EOS for nuclear statistical equilibrium with interactions

      This class is experimental.

      For the verbose parameter, generally 0 means no output, 1 means
      the function will output the composition and thermodynamics for
      the first 10 or so nuclei in the distribution, and 2 means the
      function will output the entire distribution.

      This class retains the usual mechanism using \ref err_nonconv to
      handle what to do if one of the functions does not converge. In
      addition, \ref calc_density_fixnp() and \ref
      calc_density_noneq() return \ref invalid_config for invalid
      configurations, which sometimes occur during normal execution.
      Since these invalid configurations are 'normal', they do not
      cause the error handler to be called, independent of the value
      of \ref err_nonconv . Practically, this means the end-user must
      check the return value of these two functions every time they
      are called.

      This class presumes that electrons include their rest mass,
      but nucleons and nuclei do not. The error handler is called
      by some functions if this is not the case (determined 
      by the values in <tt>o2scl::part::inc_rest_mass</tt>). 

      \todo I don't think inc_lept_phot=false works because then
      all WS cells have infinite size because of no electrons.
      For the moment, this variable is protected to discourage
      the user from changing it.

      \future There is a bit of duplication between calc_density_noneq()
      and calc_density_fixnp() which could be streamlined.
      \future Add fermion and boson statistics to the nuclei in the
      distribution. 
  */
  class eos_nse_full {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

  protected:
    
    /** \brief Check the \ref o2scl::dense_matter object to 
	see if the rest masses are correctly included or not, 
	etc.
     */
    virtual void check_dm(o2scl::dense_matter &dm);
    
    /** \brief Output a \ref o2scl::dense_matter object 
	according to the setting of \ref verbose 
	for function specified in \c func_name .
     */
    virtual void verb_output(o2scl::dense_matter &dm, 
			     std::string func_name);

    /** \brief If true, include electrons and photons (default true)
     */
    bool inc_lept_phot;

    /// Compute particle properties assuming classical thermodynamics
    o2scl::classical cla;

    /// Relativistic fermions with derivatives
    o2scl::fermion_deriv_rel snf;

    /// Mass formula (points to \ref nuc_dens by default)
    o2scl::nucmass_densmat *massp;

    /** \brief The full distribution of all nuclei to consider
	
	\note Currently, the \c ad variable doesn't do much, but
	it's important to leave this in as future functions may
	want to automatically adjust the distribution
    */
    std::vector<o2scl::nucleus> *ad;

    /** \brief Solve for charge neutrality assuming the specified
	electron chemical potential and proton number density. 

	This function returns 
	\f[
	\left[n_p-n_e(\mu_e)-n_{\mu}(\mu_{\mu}=\mu_e)\right]/n_p
	\f]
	using \ref relf.
    */
    virtual double charge_neutrality(double mu_e, double np_tot, 
				     dense_matter &dm);
    
    /** \brief Compute the free energy from a vector of densities 
	of the nuclei

	This calls \ref calc_density_noneq() and then returns the free
	energy. The vector \c n_nuc and the distribution \c dm.dist
	must both have the same size. The nuclear densities are
	taken from \c n_nuc and the proton and neutron densities
	are determined automatically from subtracting the density
	contributions of nuclei from the total neutron and proton
	densities as determined in \ref o2scl::dense_matter::nB
	and \ref o2scl::dense_matter::Ye . 

	If the call to \ref calc_density_noneq() returns a non-zero
	value, e.g. because of an invalid configuration,
	then the value \f$ 10^{4} \f$ is returned. 
    */
    virtual double free_energy(const ubvector &n_nuc, dense_matter &dm);
    
    /// Nucleonic EOS (0 by default)
    o2scl::eos_had_temp_base *ehtp;

  public:

    eos_nse_full();
    
    /** \brief The integer return value which indicates an invalid 
	configuration 
     */
    static const int invalid_config=-10;

    /// \name Various settings
    //@{
    /// Verbose parameter
    int verbose;

    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /** \brief If true, include dripped protons and
	neutrons in the nuclear mass (default true)
    */
    bool inc_prot_coul;

    /// If true, include muons (default false)
    bool include_muons;
    //@}
    
    /** \brief Function which is solved by \ref calc_density_saha()

	This function takes two inputs, the neutron and proton
	densities, and solves to ensure that \ref
	dense_matter::baryon_density() matches \ref
	o2scl::dense_matter::nB and \ref
	o2scl::dense_matter::electron_fraction() matches \ref
	dense_matter::Ye.

	This function calls \ref calc_density_fixnp() .
    */
    virtual int solve_fixnp(size_t n, const ubvector &x, ubvector &y,
			    dense_matter &dm, bool from_densities=true);
    
    /** \brief Desc
     */
    virtual int bracket_mu_solve(double &mun_low, double &mun_high,
				 double &mup_low, double &mup_high,
				 dense_matter &dm);

    /** \brief Fix electron fraction by varying proton chemical
	potential

	At some fixed values of <tt>dm.Ye</tt> and <tt>dm.nB</tt>,
	given a value of \f$ \mu_p \f$, and given an initial bracket
	for \f$ \mu_n \f$ (stored in <tt>mun_low</tt> and
	<tt>mun_high</tt>), this function attempts to find the value
	of \f$ \mu_n \f$ which ensures that the baryon density in
	nuclei matches that in <tt>dm.nB</tt> using a bracketing
	solver. It then returns the difference between the value of
	the proton fraction in nuclei and the value in <tt>dm.Ye</tt>.
	
	If <tt>mun_low</tt> and <tt>mun_high</tt> do not bracket the
	correct value of \f$ \mu_n \f$, this function attempts to
	modify them to give a proper bracket for the root. The
	finaly value of \f$ \mu_n \f$ is stored in <tt>dm.n.mu</tt>. 

	Currently, the values of <tt>dm.n.n</tt> and <tt>dm.p.n</tt>
	are ignored and set to zero.
    */
    double mup_for_Ye(double mup, double &mun_low,
		      double &mun_high, dense_matter &dm);
    
    /** \brief Fix the baryon density by varying the neutron 
	chemical potential
	
	Given a value of \f$ \mu_n \f$ (the value in <tt>dm.n.mu</tt>
	is ignored), this function computes the baryon density 
	in nuclei and returns the difference between this value
	and that stored in <tt>dm.nB</tt>.
	
	Currently, the values of <tt>dm.n.n</tt> and <tt>dm.p.n</tt>
	are ignored and set to zero.
     */
    virtual double solve_mun(double mun, dense_matter &dm);

    /** \brief Compute the properties of matter from the densities,
	not presuming equilibrium

	The values of <tt>dm.nB</tt> and <tt>dm.Ye</tt> are ignored
	and unchanged by this function. The electron and muon density
	are determined by charged neutrality and assuming their
	chemical potentials are equal. Photons are always included.

	If the nuclear densities are all zero, then this just
	returns nuclear matter with leptons and photons as
	determined by charge neutrality. 

	This function is designed to return non-zero values for
	invalid configurations and can return the value
	\ref invalid_config without calling the error handler, 
	independent of the value of \ref err_nonconv .

	Possible invalid configurations are:
	- negative nucleon or nucleus densities, or
	- proton radii larger than WS cell radii, i.e.
	\f$ (0.08 - n_p) / (n_e+n_{\mu}-n_p) < 1 \f$ or 
	\f$ n_p > 0.08 \f$ . 
    */
    virtual int calc_density_noneq(dense_matter &dm);

    /** \brief Compute the properties of matter from 
	neutron and proton densities, using the Saha equation
	
	If the parameter <tt>from_densities</tt> is true, then this
	computes nucleonic matter using the neutron and proton
	densities stored in <tt>dm.n.n</tt> and <tt>dm.p.n</tt>.
	Otherwise, nucleonic matter is computed using the chemical
	potential stored in <tt>dm.n.mu</tt> and <tt>dm.p.mu</tt>.
	Either way, electrons are computed assuming their density is
	given from \ref o2scl::dense_matter::nB and \ref
	o2scl::dense_matter::Ye. Muons are added assuming their
	chemical potential is equal to the electron chemical
	potential. Finally, the Saha equation is used to determine the
	nuclear chemical potentials and this gives the nuclear
	densities.

	This function only works when \ref inc_prot_coul is
	<tt>false</tt>.

	The values in \ref o2scl::dense_matter::nB and \ref
	o2scl::dense_matter::Ye are unchanged by this function. Note
	that, after this function completes, the value returned by
	\ref o2scl::dense_matter::baryon_density() will not
	necessarily be the same as that stored in \ref
	o2scl::dense_matter::nB (and similarly for the electron
	fraction).

	This function is designed to return non-zero values for
	invalid configurations and can return the value
	\ref invalid_config without calling the error handler, 
	independent of the value of \ref err_nonconv .

	Possible invalid configurations are:
	- negative nucleon densities, or
	- proton radii larger than WS cell radii, i.e.
	\f$ (0.08 - n_p) / (n_e+n_{\mu}-n_p) < 1 \f$ or 
	\f$ n_p > 0.08 \f$ . 
    */
    virtual int calc_density_fixnp(dense_matter &dm, bool from_densities=true);
  
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
    virtual int calc_density_by_min(dense_matter &dm);

    /** \brief Compute properties of matter for baryon density and
	electron fraction using the Saha equation

	This function solves the function specified by \ref
	solve_fixnp() using the current values of <tt>dm.n.n</tt> and
	<tt>dm.p.n</tt> as initial guesses.
    */
    virtual int calc_density_saha(dense_matter &dm);

    /** \brief Output properties of a \ref o2scl::dense_matter object to
	std::cout

	This function was particularly designed for comparing results
	with \ref o2scl::eos_sn_base derived classes.

	If output level is 0, then just the basic quantities are
	output without any information about the distribution. If
	output_level is 1, then only about 10 nuclei in the
	distribution are output, and if output_level is 2,
	then all nuclei in the distribution are output. 
    */
    virtual void output(dense_matter &dm, int output_level);

    /** \brief Adjust the particle densities to match specified
	density and electron fraction

	This function attempts to match the nuclear and nucleon
	densities so that the baryon density and electron fraction are
	equal to those specified in \ref o2scl::dense_matter::nB and
	\ref o2scl::dense_matter::Ye .
    */
    virtual int density_match(dense_matter &dm);

    /** \brief Relativistic fermions

	\comment
	Must currently be public for tcan/ecn.cpp.
	\endcomment
    */
    o2scl::fermion_rel relf;

    /// \name Nuclei and nuclear masses
    //@{
    /// Compute nuclei in dense matter
    o2scl::nucmass_densmat nuc_dens;

    /** \brief Set nuclear mass formula
     */
    void set_mass(o2scl::nucmass_densmat &m) {
      massp=&m;
      return;
    }

    /** \brief Set distribution of nuclei
     */
    void set_dist(std::vector<o2scl::nucleus> &dist) {
      ad=&dist;
      return;
    }
    //@}

    /// \name Nucleonic matter EOS
    //@{
    /** \brief Set homogeneous matter EOS
     */
    void set_eos(o2scl::eos_had_temp_base &e) {
      ehtp=&e;
      return;
    }

    /** \brief Get homogeneous matter EOS

	This function calls the error handler if no EOS has been set
     */
    o2scl::eos_had_temp_base &get_eos() {
      if (ehtp==0) {
	O2SCL_ERR2("Homogeneous matter EOS not specified in ",
		   "eos_nse_full::get_eos().",exc_efailed);
      }
      return *ehtp;
    }

    /** \brief Return true if an EOS was specified
     */
    bool is_eos_set() {
      if (ehtp==0) return false;
      return true;
    }
    //@}

    /// \name Numerical methods
    //@{
    /// The default minimizer
    o2scl::mmin_simp2<> def_mmin;

    /// Default solver
    mroot_hybrids<> def_mroot;

    /// Lepton solver
    root_cern<> def_root;
    //@}

#ifdef O2SCL_NEVER_DEFINED

    /** \brief Desc

	This was intended to be a version of calc_density_by_min()
	which optimized the composition, but it doesn't really work
	yet.
    */
    int calc_density(dense_matter &dm);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
