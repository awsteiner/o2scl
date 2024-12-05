/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2025, Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef NUCMASS_LDROP_EXT_H
#define NUCMASS_LDROP_EXT_H

#include <iostream>
#include <vector>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/table.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/cli_readline.h>
#include <o2scl/convert_units.h>
#include <o2scl/min_cern.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/eos_tov.h>
#include <o2scl/tov_solve.h>
#include <o2scl/deriv_cern.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/nucmass_ldrop_shell.h>
#include <o2scl/hdf_io.h>
#include <o2scl/test_mgr.h>

namespace o2scl {
  
  /** \brief Liquid droplet model for the neutron star crust

      The excluded volume correction is performed inside the bulk energy
      part of the mass formula, by subtracting
      \f[
      \frac{4 \pi}{3} R_n^3 f(n_{n,\mathrm{drip}},n_p=0,T)
      \f]
      where \f$ f \f$ is the free energy density of nucleonic matter.

      \note The derivative \ref dcoul_dchi diverges when chi is zero,
      and so the provided derivative is not accurate in that case.
  */
  class nucmass_ldrop_ext : public nucmass_ldrop_shell {

  public:

    /// The experimental mass model
    o2scl::nucmass_ame ame;

    /// Moller et al. mass model
    o2scl::nucmass_mnmsk moller;

    /** \brief Exponent for the high-density correction to 
	the internal densities (default in acc.cpp is 5)
    */
    double hd_exp;

    /** \brief Coefficient for the high-density correction to 
	the internal densities (default in acc.cpp is 0.5)
    */
    double hd_coeff;

    /** \brief If true, include excluded volume correction 
	(default in acc.cpp is true)
    */
    bool exc_volume;

    /// Excluded volume part 
    double exc;

    /// If true, use AME masses when available (default true)
    bool use_ame;

    /// If true, use MNMSK masses when available (default true)
    bool use_moller;

    /** \brief If true, include excluded volume and Coulomb corrections
	when using the AME and MNMSK masses (default true)
     */
    bool extra_corr;

    /// Use FRDM to compute the central densities 
    o2scl::nucmass_frdm frdm;
    
    /// \name Derivatives wrt volume fraction, \f$ \chi \f$
    //@{
    /// In MeV
    double dbulk_dchi;
    /// In MeV
    double dcoul_dchi;
    /// In MeV
    double dsurf_dchi;
    /// In MeV
    double dexc_dchi;
    /// In MeV
    double dshell_dchi;
    /// In \f$ \mathrm{fm}^{-3} \f$
    double df_dchi;
    /// In fm
    double dRn_dchi;
    /// In fm
    double dRp_dchi;
    //@}

    /// \name Derivatives wrt density of quasi-free neutrons
    //@{
    /// In MeV
    double dexc_dnn;
    /// In MeV
    double dshell_dnn;
    //@}

    /// Dimensionality of the nucleus (default 3)
    double dim;

    /// \name Derivatives wrt density of quasi-free protons
    //@{
    /// In MeV
    double dexc_dnp;
    /// In MeV
    double dcoul_dnp;
    /// In MeV
    double dshell_dnp;
    //@}

    nucmass_ldrop_ext();

    /// Test the derivatives for several nuclei
    int test_derivatives();
    
    /** \brief Test the derivatives for a specified nucleus and cell

	Note that the accuracy of, e.g. the bulk derivatives, is 
	depedent on the accuracy of the chemical potentials
	and thermodynamics from the EOS
    */
    int run_test(double Z, double N, double npout, double nnout, 
		 double chi, double T, o2scl::test_mgr &t);

    /// Compute the shell energy for nucleus Z and N
    virtual double shell_energy_new(double Z, double N, double pfact, 
                                    double nfact, double T, double dpfact, 
                                    double dnfact, double dndc, double dpdc,
                                    double &dsdp, double &dsdn, double &dsdc);

    /// Return the binding energy of the nucleus
    double drip_binding_energy_full_d(double Z, double N, double npout, 
				      double nnout, double chi, double T);
    
    /** \brief Compute the binding energy of a nucleus
	embedded in an electron gas

	This function uses drip_binding_energy_full_d() to solve for
	\c chi and the Wigner-Seitz radius \c Rws given the ambient
	electron density \c ne.
      
	The external proton and neutron densities, \c np and \c nn, and
	the electron density \c ne must be given in \f$ \mathrm{fm}^{-3}
	\f$. The temperature should be given in MeV. The radius of the
	Wigner-Seitz cell is returned in \c Rws in fm and the volume
	fraction is returned in \c chi.
      
	(This function works with a proton drip too, except that the
	chemical potential derivatives may not be correct in that case.)
      
	If the radius of the nucleus is larger than the radius of the
	cell, this function returns \f$ 10^{100} \f$.
      
	\todo Change to return an error code or something 
	instead of 1.0e100 if it's unphysical. This might be
	fixed instead also by allowing nuclei to be "inside-out".
    */
    double nucleus_be(int Z, int N, double npout, double nnout, double T,
		      double ne, double &Rws, double &chi);
    
    /** \brief Compute the nuclear binding energy of a pasta shape
	by minimizing with respect to the dimensionality
    */
    double nucleus_be_pasta(int Z, int N, double npout, double nnout, 
			    double T, double ne, double &Rws, double &chi);
  
  protected:

    /** \brief Solve for chi
	
	Used in nucmass_ldrop_ext::nucleus_be_solve().
    */
    class ldrop_solve_chi {
  
    public:

      /** \brief Solve for \f$ \chi \f$ at the specified electron density
       */
      ldrop_solve_chi(nucmass_ldrop_ext *tp, int Z, int N, double np, double nn, 
		      double T, double ne);
  
      virtual ~ldrop_solve_chi() {};
  
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual double operator()(double chi) const;

    protected:
  
      /// Store the pointer to the class instance
      nucmass_ldrop_ext *tptr;
  
      /// \name The parameters
      //@{
      double np_, nn_, T_, ne_;
      int Z_, N_;
      //@}

    };

    /** \brief Class for computing derivatives of different
	parts of the mass formula

	This is used in the function run_test() to test the
	derivatives.
    */
    class ldrop_be_deriv {

    public:

      /// Construct an object for the specified nucleus
      ldrop_be_deriv(nucmass_ldrop_ext &ld, double Z, double N, double npout, 
		     double nnout, double chi, double T, size_t ix, 
		     size_t jx);

      virtual ~ldrop_be_deriv() {}

      /// The operator() for computing derivatives
      virtual double operator()(double x);

    protected:

      /// The mass formula
      nucmass_ldrop_ext *ldp;

      /// \name The parameters of the mass model
      //@{
      double Z_,N_,np_,nn_,chi_,T_;
      //@}

      /// The index of the part of the mass formula
      size_t ix_;

      /// Parameter index
      size_t jx_;

    };

  public:

    /// Solver for nucleus_be()
    o2scl::root_brent_gsl<ldrop_solve_chi> grb;
    
  };

}

#endif
