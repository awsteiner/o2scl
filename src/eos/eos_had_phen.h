/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2022, Xingfu Du and Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software: you can redistribute it and/or modify
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
#ifndef O2SCL_EOS_HAD_PHEN_H
#define O2SCL_EOS_HAD_PHEN_H

#include <gsl/gsl_sf_hyperg.h>

#include <o2scl/test_mgr.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/fit_nonlin.h>
#include <o2scl/eos_crust_virial.h>
#include <o2scl/rng.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/deriv_cern.h>

#include <o2scl/eos_had_virial.h>

namespace o2scl {

  /** \brief An updated version of \ref o2scl::eos_crust_virial
      with a better fit for the virial coefficients
  */
  class eos_crust_virial_v2 : public o2scl::eos_crust_virial {
  
  public:
  
    /** \brief The neutron-neutron virial coefficient given the
	function parameters specified in \c par
    */
    double bn_func(size_t np, const std::vector<double> &par, double T);
  
    /** \brief The neutron-proton virial coefficient given the
	function parameters specified in \c par
    */
    double bpn_func(size_t np, const std::vector<double> &par, double T);

    /** \brief The neutron-neutron virial coefficient
     */
    double bn_f(double T);
  
    /** \brief The neutron-proton virial coefficient
     */
    double bpn_f(double T);
  
    /** \brief The temperature derivative of the
	neutron-neutron virial coefficient
    */
    double dbndT_f(double T);
  
    /** \brief The temperature derivative of the
	neutron-proton virial coefficient
    */
    double dbpndT_f(double T);
  
    /** \brief The current neutron-neutron virial coefficient parameters
     */
    std::vector<double> bn_params;
  
    /** \brief The current neutron-proton virial coefficient parameters
     */
    std::vector<double> bpn_params;

    /** \brief The number of neutron-neutron virial coefficient parameters
     */
    static const size_t bn_np=10;

    /** \brief The number of neutron-proton virial coefficient parameters
     */
    static const size_t bpn_np=6;
  
    /** \brief Perform the fit to the scattering data
     */
    virtual void fit(bool show_fit=false);

    /// \name Constructor and destructor
    //@{
    eos_crust_virial_v2();

    virtual ~eos_crust_virial_v2() {}
    //@}
  
  }; 

  /** \brief Phenomenological EOS for homogeneous nucleonic matter
   */
  class eos_had_phen : public o2scl::eos_had_temp_eden_base {
  
  protected:

    /// \name Main EOS parameters [protected]
    //@{
    /// The first exponent for density in the QMC EOS (unitless)
    double qmc_alpha;

    /// The first coefficient for the QMC EOS (in MeV)
    double qmc_a;
  
    /** \brief The speed of sound in neutron star matter at 
	\f$ n_B=2~\mathrm{fm}^{-3} \f$
    */
    double phi;

    /// The symmetry energy
    double eos_S;

    /// The slope of the symmetry energy
    double eos_L;

    /// The index of the neutron star model
    int i_ns;

    /// The index of the Skyrme model
    int i_skyrme;
    //@}

    /// \name Basic EOS functions [protected]
    //@{
    /** \brief Return the total free energy density of matter
	(without the rest mass contribution for the nucleons)
    */
    double free_energy_density
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th);

    /** \brief Compute the free energy returning several 
	details as output parameters
      
	f1 is g*f_virial
	f2 is (1-g)*f_skyrme
	f3 is (1-g)*delta^2*esym
	f4 is (1-g)*delta f_hot

	so that the total homogeneous free energy is f1+f2+f3+f4
    */
    double free_energy_density_detail
    (o2scl::fermion &n, o2scl::fermion &p, double T, o2scl::thermo &th,
     double &zn, double &zp,
     double &f1, double &f2, double &f3, double &f4, 
     double &g_virial, double &dgvirialdT);

    /** \brief Compute the free energy density using the virial 
	expansion including derivative information
    */
    virtual double free_energy_density_virial
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th, double &dmundnn, double &dmundnp,
     double &dmupdnn, double &dmupdnp, double &dmundT,
     double &dmupdT);    

    /** \brief Compute the free energy density using the virial 
	expansion
    */
    virtual double free_energy_density_virial
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th) {
      double x1, x2, x3, x4, x5, x6;
      return free_energy_density_virial(n,p,T,th,x1,x2,x3,x4,x5,x6);
    }
  
    /** \brief Alternate form of \ref free_energy_density() for
	computing derivatives

	This function does not include electrons or photons.
    */
    double free_energy_density_alt(o2scl::fermion &n, o2scl::fermion &p,
				   double nn, double np, double T,
				   o2scl::thermo &th);

    /** \brief Alternate form of \ref free_energy_density() 
	which includes electrons, positrons, and photons.
    */
    double free_energy_density_ep(double nn, double np, double T);
  
    /** \brief Compute the entropy density including photons,
	electrons, and positrons
    */
    double entropy(o2scl::fermion &n, o2scl::fermion &p,
		   double nn, double np, double T, o2scl::thermo &th);

    /** \brief Compute energy density including photons and electons
	(without the rest mass energy density for the nucleons)
    */
    double ed(o2scl::fermion &n, o2scl::fermion &p,
	      double nn, double np, double T, o2scl::thermo &th);

    /** \brief Compute the squared speed of sound 
     */
    double cs2_func(o2scl::fermion &n, o2scl::fermion &p, double T,
		    o2scl::thermo &th);
    //@}

    /// \name Internal variables [protected]
    //@{
    /// The table which stores the neutron star EOS results
    o2scl::table_units<> nstar_tab;

    /// The table which stores the Skyrme fits
    o2scl::table_units<> UNEDF_tab;
  
    /** \brief If true, a model has been selected (default false)
     */
    bool model_selected;
  
    /// Random number generator
    o2scl::rng<> rg;
    //@}
  
    /// \name EOS outputs
    //@{
    /// The free energy of degenerate matter
    double f_deg;
  
    /// The virial free energy
    double f_virial;
  
    /// The virial entropy
    double s_virial;

    /** \brief The value of \f$ \bar{\Lambda} \f$ 
	for a 1.4 solar mass neutron
	star
    */
    double Lambda_bar_14;
    //@}
  
    /// \name The fit to the neutron star EOS [protected]
    //@{
    /** \brief Compute the energy density (in \f$ \mathrm{fm}^{-4} \f$)
	of neutron matter at high density from the neutron star data
	using the most recent fit (without the rest mass contribution)

	\note Currently this just returns the value of
	\ref ed_fit() .
    */
    double energy_density_ns(double nn);

    /// Parameters for the function which fits the neutron star EOS
    std::vector<double> ns_fit_parms;

    /** \brief The fit function for the energy per particle
	in units of MeV as a function of the baryon density
	(in \f$ \mathrm{fm}^{-3} \f$ )

	Note this function does not include the rest mass 
	energy density for the nucleons. 
    */
    double fit_fun(size_t np, const std::vector<double> &parms,
		   double nb);

    /** \brief The energy density (in \f$ \mathrm{fm}^{-4} \f$ )
	as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )

	Note this function does not include the rest mass 
	energy density for the nucleons. 
    */
    double ed_fit(double nb);
  
    /** \brief The inverse susceptibility (in \f$ \mathrm{fm}^{2} \f$ )
	as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )
    */
    double dmudn_fit(double nb);
  
    /** \brief The speed of sound
	as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )
    */
    double cs2_fit(double nb);
  
    /** \brief Compute the minimum and maximum speed of sound
	between 0.08 and \ref ns_nb_max
    */
    void min_max_cs2(double &cs2_min, double &cs2_max);

    /** \brief Fit neutron star data from Bamr to an analytical 
	expression 
    */
    void ns_fit(int row);

    /// The chi-squared for the neutron star fit
    double chi2_ns;

    /** \brief The maximum baryon density at which the neutron star
	EOS is causal

	This quantity is determined by \ref ns_fit()
    */
    double ns_nb_max;
  
    /** \brief The baryon number chemical potential (in \f$
	\mathrm{fm}^{-1} \f$ ) as a function of number density (in \f$
	\mathrm{fm}^{-3} \f$ )

	Note this function does not include the rest mass 
	for the nucleons. 
    */
    double mu_fit(double nb);
    //@}

    /// \name Other EOS functions [protected]
    //@{
    /** \brief Compute the energy density (in \f$ \mathrm{fm}^{-4} \f$)
	of neutron matter from quantum Monte Carlo (without the rest
	mass contribution)
    */
    double energy_density_qmc(double nn, double pn);

    /** \brief Construct a new neutron star EOS which ensures
	causality at high densities
    */
    int new_ns_eos(double nb, o2scl::fermion &n, double &e_ns,
		   double &densdnn);

    /** \brief Compute dfdnn including photons and electons
     */
    double dfdnn_total(o2scl::fermion &n, o2scl::fermion &p,
		       double nn, double pn, double T, o2scl::thermo &th);
  
    /** \brief Compute dfdnp including photons and electons
     */
    double dfdnp_total(o2scl::fermion &n, o2scl::fermion &p,
		       double nn, double pn, double T, o2scl::thermo &th);
  
    /** \brief Solve for Ye to ensure a specified value of muL at fixed T
     */
    int solve_Ye(size_t nv, const ubvector &x, ubvector &y,
		 double nb, double T, double muL);
  
    /** \brief solve for a1 and a2 when cs_ns(2.0)>cs_ns(1.28)
     */
    int solve_coeff_big(size_t nv, const ubvector &x, ubvector &y, 
			double nb_last, double cs_ns_2, double cs_ns_last);

    /** \brief solve for a1 and a2 when cs_ns(2.0)<cs_ns(1.28)
     */
    int solve_coeff_small(size_t nv, const ubvector &x, ubvector &y, 
			  double nb_last, double cs_ns_2, double cs_ns_last);
  
    /** \brief Internal select function
     */
    int select_internal(int i_ns_loc, int i_skyrme_loc,
			double qmc_alpha_loc, double qmc_a_loc,
			double eos_L_loc, double eos_S_loc, double phi_loc);
    //@}

    /// \name Particle objects [protected]
    //@{
    /** \brief Electron/positron
     */
    o2scl::fermion electron;

    /** \brief Muon/anti-muon
     */
    o2scl::fermion muon;

    /** \brief Photon
     */
    o2scl::boson photon;

    /// Neutron
    o2scl::fermion neutron;

    /// Proton
    o2scl::fermion proton;

    /// Neutron for chiral part
    o2scl::fermion n_chiral;

    /// Proton for chiral part
    o2scl::fermion p_chiral;

    /// Neutrino
    o2scl::fermion neutrino;  
    //@}

    /// \name Base physics objects [protected]
    //@{
    /// The virial equation solver
    eos_had_virial_deriv vsd;

    /// Old virial solver
    eos_had_virial vs;

    /** \brief Object for computing electron/positron thermodynamic integrals
     */
    o2scl::fermion_rel relf;

    /// Thermodynamic quantities
    o2scl::thermo th2;
  
    /// Thermodynamic quantities for chiral part
    o2scl::thermo th_chiral;
  
    /// Base EOS model
    o2scl::eos_had_skyrme sk;

    /// Skyrme model for finite-temperature correction
    o2scl::eos_had_skyrme sk_Tcorr;

    /// Pointer to EOS for finite-temperature corrections
    o2scl::eos_had_temp_eden_base *eos_Tcorr;

    /// The virial EOS
    eos_crust_virial_v2 ecv;

    /// Alternative skryme model
    o2scl::eos_had_skyrme sk_alt;

    /// Pointer to alternative model
    o2scl::eos_had_temp_base *eosp_alt;
    //@}

    /// \name The parameters for the QMC energy density [protected]
    //@{
    /// The second exponent for density in the QMC EOS (unitless)
    double qmc_beta;
    /// The second coefficient for the QMC EOS (in MeV)
    double qmc_b;
    /** \brief The saturation density of the QMC EOS, equal to 
	\f$ 0.16~\mathrm{fm}^{-3} \f$
    */
    double qmc_n0;
    //@}
  
    /// \name Output saturation properties [protected]
    //@{
    /// The binding energy per particle
    double eos_EoA;
  
    /// The incompressibility
    double eos_K;
  
    /// The saturation density
    double eos_n0;
    //@}

  public:

    /// \name Constructor and destructor
    //@{
    eos_had_phen();

    virtual ~eos_had_phen() {
    }
    //@}

    /** \brief Load the required data files
     */
    void load_files();
    
    /// \name Settings [public]
    //@{
    /** \brief If true, use the EOS from the Du et al. (2019) paper
	instead of the Du et al. (2020) update (default false)
    */
    bool old_version;
  
    /** \brief Use a Skyrme model rather than the Du et al. 
	combined EOS
    */
    bool use_skalt;
  
    /** \brief If true, test the neutron star speed of sound 
	(default true)
    */
    bool test_ns_cs2;
  
    /** \brief If true, save the results of the neutron star fit to
	a file, and immediately exit (default false)
    */
    bool ns_record;

    /** \brief If true, use the old neutron star fit (default true)

	This defaults to true because the old fit performs a bit
	better than the new one. The new fit was never used
	in a publication. 
    */
    bool old_ns_fit;

    /// Verbose parameter
    int verbose;

    /// If true, create output files for individual EOSs
    bool output_files;

    /// Coefficient for modulation of virial EOS
    double a_virial;

    /// Coefficient for modulation of virial EOS
    double b_virial;
  
    /** \brief If true, include muons (default false)
     */
    bool include_muons;

    /** \brief If true, test cs2 in the \ref select_internal() function
	(default true)
    */
    bool select_cs2_test;
    //@}

    /** \brief Equation of state as a function of densities at 
	finite temperature
    */
    virtual int calc_temp_e(fermion &n, fermion &p, double T, 
			    thermo &th) {
      double fr=free_energy_density(n,p,T,th);
      return 0;
    }
    
    /** \brief Equation of state as a function of densities at 
	zero temperature
    */
    virtual int calc_e(fermion &n, fermion &p, 
                       thermo &th) {
      return calc_temp_e(n,p,0.0,th);
    }
    
    /// \name Command-line interface functions [public]
    //@{
    /** \brief Construct a table at fixed electron fraction
     */
    int table_Ye(std::vector<std::string> &sv,
		 bool itive_com);

    /** \brief Construct a table at fixed baryon density
     */
    int table_nB(std::vector<std::string> &sv,
		 bool itive_com);

    /** \brief Construct the EOS for a proto-neutron star
     */
    int pns_eos(std::vector<std::string> &sv, bool itive_com);
  
    /** \brief Construct a full table 
     */
    int table_full(std::vector<std::string> &sv, bool itive_com);

    /** \brief Test the code
     */
    int test_deriv(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select a model by specifying the parameters
     */
    int select_model(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compare the full free energy with the free energy
	from the virial expansion
    */
    int vir_comp(std::vector<std::string> &sv, bool itive_com);

    /** \brief Evaluate the EOS at one point
     */
    int point(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select a random model
     */
    int random();

    /** \brief Compute the data for the comparison figures
     */
    int comp_figs(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the data for the Monte Carlo figures
     */
    int mcarlo_data(std::vector<std::string> &sv, bool itive_com);

    /** \brief Perform the virial fit
     */
    int vir_fit(std::vector<std::string> &sv, bool itive_com);

    /** \brief Test the electron and photon contribution
     */
    int test_eg(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the EOS from previously generated EOS
	tables at several points
    */
    int eos_sn(std::vector<std::string> &sv, bool itive_com);
    //@}

    /// \name Miscellaneous functions [public]
    //@{
    /** \brief Solve for fixed entropy per baryon and fixed
	lepton fraction
    */
    int solve_fixed_sonb_YL(size_t nv, const ubvector &x, ubvector &y,
			    double nB, double sonb, double YL);

    /** \brief Solve for T to ensure a specified value of sonb at fixed Ye
     */
    int solve_T(size_t nv, const ubvector &x, ubvector &y,
		double nb, double Ye, double sonb);
  
    //@}
  
  };

}

#endif
