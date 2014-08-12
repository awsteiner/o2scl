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
/** \file eos_had_sym4.h
    \brief File defining \ref o2scl::eos_had_sym4
*/
#ifndef O2SCL_SYM4_EOS_H
#define O2SCL_SYM4_EOS_H

#include <iostream>
#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_potential.h>
#include <o2scl/test_mgr.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A class to compute quartic contributions to the symmetry 
      energy [abstract base]

      The standard usage is that a child class implements the virtual
      function calc_e_sep() which is then used by calc_e_alpha()
      and calc_muhat(). These functions are employed by \ref eos_had_sym4
      to compute the EOS for an arbitrary dependence of the 
      symmetry energy on the isospin.

      \hline
      <b>References:</b>

      Created for \ref Steiner06.

      \bug Testing was disabled in HDF conversion. Fix this.
   */
  class eos_had_sym4_base {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// An electron for the computation of the \f$ \hat{\mu}\f$
    fermion e;

    /// Zero-temperature fermion thermodynamics
    fermion_zerot fzt2;

#endif

  public:

    eos_had_sym4_base();
  
    virtual ~eos_had_sym4_base() {}

    /** \brief Compute alpha at the specified density
     */
    virtual int calc_e_alpha(fermion &ne, fermion &pr, thermo &lth,
			     double &alphak, double &alphap, double &alphat,
			     double &diff_kin, double &diff_pot,
			     double &ed_kin_nuc, double &ed_pot_nuc);
    
    /** \brief Compute \f$ \hat{\mu} \f$, the out-of-whack parameter
     */
    virtual double calc_muhat(fermion &ne, fermion &pr);

    /** \brief Compute the potential and kinetic parts separately (to
	be overwritten in children)
     */
    virtual int calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			   double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			   double &mu_n_pot, double &mu_p_pot)=0;

  };

  /** \brief A version of \ref eos_had_rmf to separate potential and kinetic
      contributions

      \hline
      <b>References:</b>

      Created for \ref Steiner06.
   */
  class eos_had_sym4_rmf : public eos_had_rmf, public eos_had_sym4_base {

  public:
  
    /** \brief Compute the potential and kinetic parts separately
     */
    virtual int calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			   double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			   double &mu_n_pot, double &mu_p_pot);

  };

  /** \brief A version of \ref eos_had_apr to separate potential and kinetic
      contributions

      \hline
      <b>References:</b>

      Created for \ref Steiner06.
  */
  class eos_had_sym4_apr : public eos_had_apr, public eos_had_sym4_base {
    
  public:

    /** \brief Compute the potential and kinetic parts separately
     */
    virtual int calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			   double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			   double &mu_n_pot, double &mu_p_pot);
  };

  /** \brief A version of \ref eos_had_skyrme to separate potential and kinetic
      contributions

      \hline
      <b>References:</b>

      Created for \ref Steiner06.
  */
  class eos_had_sym4_skyrme : public eos_had_skyrme, public eos_had_sym4_base {
    
  public:

    /** \brief Compute the potential and kinetic parts separately
     */
    virtual int calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			   double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			   double &mu_n_pot, double &mu_p_pot);
  };
  
  /** \brief A version of \ref eos_had_potential to 
      separate potential and kinetic contributions

      \hline
      <b>References:</b>

      Created for \ref Steiner06.
  */
  class eos_had_sym4_mdi : public eos_had_potential, 
    public eos_had_sym4_base {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Compute the kinetic part of the energy density
    double energy_kin(double var);

    /// Compute the potential part of the energy density
    double energy_pot(double var);

#endif
    
  public:

    /** \brief Compute the potential and kinetic parts separately
     */
    virtual int calc_e_sep(fermion &ne, fermion &pr, double &ed_kin, 
			   double &ed_pot, double &mu_n_kin, double &mu_p_kin, 
			   double &mu_n_pot, double &mu_p_pot);

    /// Test the separation of the potential and kinetic energy parts
    virtual int test_separation(fermion &ne, fermion &pr, test_mgr &t);

  };

  /** \brief Construct an EOS with an arbitrary choice for the terms
      in the symmetry energy that are quartic in the isospin asymmetry

      \hline
      <b>References:</b>

      Created for \ref Steiner06.
  */
  class eos_had_sym4 : public eos_had_eden_base {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The base equation of state to use
    eos_had_sym4_base *sp;

#endif

  public:

    /// The strength of the quartic terms
    double alpha;

    /// Set the base equation of state
    int set_base_eos(eos_had_sym4_base &seb);
  
    /** \brief Test the equation of state 
    
	This compares the chemical potentials from calc_e_sep() to
	their finite-difference approximations in order to ensure that
	the separation into potential and kinetic parts is done
	properly.
     */
    virtual int test_eos(fermion &ne, fermion &pr, thermo &lth);

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lth);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
