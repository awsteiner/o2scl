/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
/** \file nucleus_rmf.h
    \brief File defining \ref o2scl::nucleus_rmf
*/
#ifndef RMF_NUCLEUS_H
#define RMF_NUCLEUS_H

#include <iostream>
#include <string>
#include <vector>
#include <o2scl/interp.h>
#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/table_units.h>
#include <o2scl/ode_rkck_gsl.h>
#include <o2scl/ode_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Spherical closed-shell nuclei with a relativistic
      mean-field model in the Hartree approximation

      This code is very experimental.

      This class is based on a code developed by C.J. Horowitz and
      B.D. Serot, and used in \ref Horowitz81 which was then adapted
      by P.J. Ellis and used in \ref Heide94 and \ref Prakash94. Ellis
      and A.W. Steiner adapted it for the parameterization in in \ref
      eos_had_rmf for \ref Steiner05b, and then converted to C++ by
      Steiner afterwards.

      The standard usage is something like:
      \code
      nucleus_rmf rn;
      o2scl_hdf::rmf_load(rn.rmf,"NL4");
      rn.run_nucleus(82,208,0,0);
      cout << rn.rnrp << endl;
      \endcode
      which computes the structure of \f$ ^{208}\mathrm{Pb} \f$ and
      outputs the neutron skin thickness using the model \c 'NL4'.

      Potential exceptions are
      - Failed to converge
      - Failed to solve meson field equations
      - Energy not finite (usually a problem in the equation of
      state)
      - Energy not finite in final calculation
      - Function \ref iterate() called before \ref init_run()
      - Not a closed-shell nucleus

      The initial level pattern is
      \verbatim
      1 S 1/2
      // 2 nucleons
      1 P 3/2
      1 P 1/2
      // 8 nucleus
      1 D 5/2
      1 D 3/2
      2 S 1/2
      // 20 nucleons
      1 F 7/2
      // 28 nucleons
      1 F 5/2
      2 P 3/2
      2 P 1/2
      // 40 nucleons
      1 G 9/2
      // 50 nucleus
      1 G 7/2
      2 D 5/2
      1 H 11/2
      2 D 3/2
      3 S 1/2
      // 82 nucleons
      1 H 9/2
      2 F 7/2
      1 I 13/2
      2 F 5/2
      3 P 3/2
      3 P 1/2
      // 126 nucleons
      2 G 9/2
      1 I 11/2
      1 J 15/2
      3 D 5/2
      4 S 1/2
      2 G 7/2
      3 D 3/2
      // 184 nucleons
      \endverbatim
      
      Below, \f$ \alpha \f$ is a generic index for the isospin, the
      radial quantum number \f$ n \f$ and the angular quantum numbers
      \f$ \kappa \f$ and \f$ m \f$. The meson fields are \f$
      \sigma(r), \omega(r) \f$ and \f$ \rho(r) \f$. The baryon density
      is \f$ n(r) \f$, the neutron and proton densities are \f$ n_n(r)
      \f$ and \f$ n_p(r) \f$, and the baryon scalar density is \f$
      n_s(r) \f$.
      The nucleon field equations are
      \f{eqnarray*}
      F^{\prime}_{\alpha}(r)- \frac{\kappa}{r} F_{\alpha}(r)
      + \left[ \varepsilon_{\alpha} - g_{\omega} \omega(r) 
      - t_{\alpha} g_{\rho} \rho(r) - \left( t_{\alpha}+\frac{1}{2} \right) 
      e A(r) - M + g_{\sigma} \sigma(r) \right] G_{\alpha}(r) &=& 0 \\
      G^{\prime}_{\alpha}(r)+ \frac{\kappa}{r} G_{\alpha}(r)
      - \left[ \varepsilon_{\alpha} - g_{\omega} \omega(r) 
      - t_{\alpha} g_{\rho} \rho(r) - \left( t_{\alpha}+\frac{1}{2}
      \right) e A(r) + M - g_{\sigma} \sigma(r) \right] F_{\alpha}(r) &=& 0
      \f}
      where the isospin number, \f$ t_{\alpha} \f$ is \f$ 1/2 \f$ for
      protons and \f$ -1/2 \f$ for neutrons.
      The meson field equations are
      \f{eqnarray*}
      \sigma^{\prime \prime}(r) + \frac{2}{r} \sigma^{\prime}(r) 
      - m_{\sigma}^2 \sigma &=& - g_{\sigma} n_s(r) + b M g_{\sigma}^3
      \sigma^2 + c g_{\sigma}^4 \sigma^3 - g_{\rho}^2 \rho^2 
      \frac{\partial f}{\partial \sigma} \\
      \omega^{\prime \prime}(r) + \frac{2}{r} \omega^{\prime}(r)
      - m_{\omega}^2 \omega &=& - g_{\omega} n(r) + 
      \frac{\zeta}{6} g_{\omega}^4 \omega^3 + g_{\rho}^2 \rho^2 
      \frac{\partial f}{\partial \omega} \\
      \rho^{\prime \prime}(r) + \frac{2}{r} \rho^{\prime}(r)
      - m_{\rho}^2 \rho &=& - \frac{g_{\rho}}{2} 
      \left[n_n(r)-n_p(r)\right] + 2 g_{\rho}^2 \rho f + \frac{\xi}{6}
      g_{\rho}^4 \rho^3
      \f}
      and the Coulomb field equation is 
      \f[
      A^{\prime \prime}(r) + \frac{2}{r} A^{\prime}(r) = - e n_p(r) 
      \f]
      The meson field equations plus a pair of Dirac-like nucleon
      field equations for each index \f$ \alpha \f$ must be solved to
      compute the structure of a given nucleus.

      The densities (scalar, baryon, isospin, and charge) are
      \f{eqnarray*}
      n_s(r) &=& \sum_{\alpha} \left\{ \int d^3 r \left[ G(r)^2-F(r)^2 
      \right] \right\} \\
      n(r) &=& \sum_{\alpha} \left\{ \int d^3 r \left[ G(r)^2+F(r)^2 
      \right] \right\} \\
      n_i(r) &=& \sum_{\alpha} \left\{ t_{\alpha} \int d^3 r 
      \left[ G(r)^2-F(r)^2 \right] \right\} \\
      n_c(r) &=& \sum_{\alpha} \left\{ \left[t_{\alpha}+\frac{1}{2}\right] 
      \int d^3 r \left[ G(r)^2-F(r)^2 \right] \right\}
      \f}

      The total energy is
      \f[
      E = \sum_{\alpha} \varepsilon_{\alpha} (2 j_{\alpha}+1)
      - \frac{1}{2} \int d^{3} r 
      \left[- g_{\sigma} \sigma(r) \rho_s(r) + g_{\omega} \omega(r)
      \rho(r) + \frac{1}{2} g_{\rho} \rho(r) + e A(r) n_p(r) \right]
      \f]

      The charge density is the proton density folded with the 
      charge density of the proton, i.e.
      \f[
      \rho_{\mathrm{ch}}(r) = \int d^{3} r^{\prime} 
      \rho_{\mathrm{prot}}(r-r^{\prime}) \rho_p(r)
      \f]
      where the proton charge density is assumed to be of the form
      \f[
      \rho_{\mathrm{prot}}(r) = \frac{\mu^3}{8 \pi} \exp \left(
      - \mu |r|\right)
      \f]
      and the parameter \f$ \mu = (0.71)^{1/2}~\mathrm{GeV} \f$ (see
      Eq. 20b in \ref Horowitz81). The default value of \ref a_proton
      is the value of \f$ \mu \f$ converted into \f$ \mathrm{fm}^{-1}
      \f$.

      Generally, the first array index associated with a function
      is not the value at \f$ r=0 \f$, but at \f$ r=\Delta r \f$
      (stored in \ref step_size) and so the \f$ r=0 \f$ part of the
      algorithm is handled separately. 

      \todo Better documentation
      \todo Convert energies() to use EOS and possibly
      replace sigma_rhs() and related functions by the associated
      field equation method of eos_had_rmf.

      \todo Document hw=3.923+23.265/cbrt(atot);

      \comment
      the hw=blah blah term for the CM correction is discussed
      a bit in Negele, PRC 1 (1970) 1260, esp. see Eq. 2.30 but
      the numerical coefficients are different here.
      \endcomment

      \todo I believe currently the center of mass correction
      for the binding energy per nucleon is not done and has
      to be added in after the fact

      \future Sort energy levels at the end by energy
      \future Improve the numerical methods
      \future Make the neutron and proton orbitals more configurable
      \future Generalize to \f$ m_n \neq m_p \f$ .
      \future Allow more freedom in the integrations
      \future Consider converting everything to inverse fermis.
      \future Convert to zero-indexed arrays (mostly done)
      \future Warn when the level ordering is wrong, and unoccupied levels
      are lower energy than occupied levels
      \future Connect with \ref o2scl::nucmass ?
  */
  class nucleus_rmf {

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

#ifndef DOXYGEN_INTERNAL

  protected:

    /// A convenient struct for the solution of the Dirac equations
    typedef struct {
      /// Eigenvalue
      double eigen;
      /// Quantum number \f$ \kappa \f$ .
      double kappa;
      /// The meson fields
      ubmatrix *fields;
      /// Desc
      ubvector *varr;
    } odparms;

    /// The total number of shells stored internally
    static const int n_internal_levels=29;

#endif

  public:
    
    nucleus_rmf();

    ~nucleus_rmf();

    /** \name Basic operation
     */
    //@{
    /// A shell of nucleons for \ref nucleus_rmf
    typedef struct {
      /// Degeneracy \f$ 2 j+1 \f$ .
      int twojp1;
      /// \f$ \kappa \f$ 
      int kappa;
      /// Energy eigenvalue
      double energy;
      /// Isospin ( \f$ +1/2 \f$ or \f$ -1/2 \f$ .
      double isospin;
      /// Angular momentum-spin state \f$ ^{2s+1} \ell_{j} \f$ 
      std::string state;
      /// Matching radius (in fm)
      double match_point;
      /// Desc.
      double eigen;
      /// Desc.
      double eigenc;
      /// Number of nodes in the wave function
      int nodes;
    } shell;

    /** \brief Computes the structure of a nucleus with the specified
	number of levels

	Note that \ref rmf must be set before run_nucleus() is called.

	This calls init_run(), and then iterate() until \c iconverged is
	1, and then post_converge().
    */
    int run_nucleus(int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N);

    /// Set output level
    void set_verbose(int v) { verbose=v; };
    //@}
    
    /** \name Lower-level interface 
     */
    //@{
    /** \brief Initialize a run

        Note that \ref rmf must be set before run_nucleus() is called.
    */
    void init_run(int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N);
    
    /// Perform an iteration
    void iterate(int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N,
		 int &iconverged);
    
    /// After convergence, make CM corrections, etc.
    int post_converge(int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N);
	
    //@}
    
    /// \name Results
    //@{
    /** \brief Get the radial profiles
	
	The profiles are calculated each iteration by iterate().
    */
    std::shared_ptr<table_units<> > get_profiles() { return profiles; };
    
    /** \brief The final charge densities
     */
    std::shared_ptr<table_units<> > get_chden() { return chden_table; };
    
    /// The number of levels
    int nlevels;

    /** \brief The levels (protons first, then neutrons)
    
	An array of size \ref nlevels
    */
    std::vector<shell> levels;
    
    /// The number of unoccupied levels (equal to \c unocc_Z + \c unocc_N)
    int nuolevels;

    /** \brief The unoccupied levels (protons first, then neutrons)
	
	An array of size \ref nuolevels
    */
    std::vector<shell> unocc_levels;

    /** \brief Surface tension (in \f$ \mathrm{fm}^{-3} \f$ )

	Computed in post_converge() or automatically in run_nucleus()
    */
    double stens;

    /** \brief Skin thickness (in fm)
    
	Computed every iteration in iterate()
	or automatically in run_nucleus()
    */
    double rnrp;

    /** \brief Neutron RMS radius (in fm)
    
	Computed every iteration in iterate() or automatically in
	run_nucleus()
    */
    double rnrms;

    /** \brief Proton RMS radius (in fm)
    
	Computed every iteration in iterate() or automatically in
	run_nucleus()
    */
    double rprms;

    /** \brief Total energy (in MeV)
    
	Computed every iteration in iterate() or automatically in
	run_nucleus()
    */
    double etot;

    /** \brief Charge radius (in fm)
    
	Computed in post_converge() or automatically in run_nucleus()
    */
    double r_charge;
    
    /** \brief Charge radius corrected by the center of mass (in fm)
    
	Computed in post_converge() or automatically in run_nucleus()
    */
    double r_charge_cm;
    //@}

    /** \name Equation of state
     */
    //@{
    /** \brief The default equation of state (default NL3)

	This is set in the constructor to be the default
	model, NL3, using the function \ref load_nl3().
    */
    eos_had_rmf def_rmf;
    
    /** \brief Set the base EOS to be used
	
	The equation of state must be set before run_nucleus() or
	init_fun() are called, including the value of eos_had_rmf::mnuc.
    */
    int set_eos(eos_had_rmf &r) {
      rmf=&r;
      return 0;
    }

    /** \brief \ref thermo object for the EOS

	This is just used as temporary storage.
    */
    thermo hb;

    /** \brief The neutron 

	The mass of the neutron is ignored and set by init_run() 
	to be eos_had_rmf::mnuc from \ref rmf.
    */
    fermion n;

    /** \brief The proton

	The mass of the proton is ignored and set by init_run() 
	to be eos_had_rmf::mnuc from \ref rmf.
    */
    fermion p;
    //@}

    /** \name Numeric configuration
     */
    //@{
    /** \brief If true, call the error handler if the routine does not 
	converge or reach the desired tolerance (default true)
    */
    bool err_nonconv;
  
    /// Set the stepper for the Dirac differential equation
    void set_step(ode_step<ubvector,ubvector,ubvector,
      ode_funct> &step) {
      ostep=&step;
      return;
    }

    /// Maximum number of total iterations (default 70)
    int itmax;

    /** \brief Maximum number of iterations for solving the meson
	field equations (default 10000)
    */
    int meson_itmax;

    /** \brief Maximum number of iterations for solving the Dirac
	equations (default 100)
    */
    int dirac_itmax;

    /// Tolerance for Dirac equations (default \f$ 5 \times 10^{-3} \f$ ).
    double dirac_tol;
    
    /** \brief Second tolerance for Dirac equations 
	(default \f$ 5 \times 10^{-4} \f$ ).
    */
    double dirac_tol2;

    /// Tolerance for meson field equations (default \f$ 10^{-6} \f$ ).
    double meson_tol;
    //@}

    /** \brief Initial guess structure

	The initial guess for the meson field profiles is 
	a set of Fermi-Dirac functions, i.e.
	\f[
	\sigma(r)=\mathrm{sigma0}/
	[1+\exp((r-\mathrm{fermi\_radius})/\mathrm{fermi\_width})]
	\f]
    */
    typedef struct {
      /// Scalar field at r=0
      double sigma0;
      /// Vector field at r=0
      double omega0;
      /// Isubvector field at r=0
      double rho0;
      /// Coulomb field at r=0
      double A0;
      /// The radius for which the fields are half their central value
      double fermi_radius;
      /// The "width" of the Fermi-Dirac function
      double fermi_width;
    } initial_guess;

    /** \brief Parameters for initial guess

	Default is {310,240,-6,25.9,6.85,0.6}
    */
    initial_guess ig;

    /** \brief If true, use the generic ODE solver instead of the 
	internal 4th order Runge-Kutta
     */
    bool generic_ode;

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The parameter for the charge density of the proton
	(default is about 4.27073)
    */
    double a_proton;

    /// Load the default model NL3 into the given \ref eos_had_rmf object
    int load_nl3(eos_had_rmf &r);
    
    /// The base EOS
    eos_had_rmf *rmf;

    /** \brief The radial profiles
     */
    std::shared_ptr<table_units<> > profiles;
    
    /** \brief The final charge densities
     */
    std::shared_ptr<table_units<> > chden_table;
    
    /** \brief A pointer to the current vector of levels 
	(either \ref levels or \ref unocc_levels)
    */
    std::vector<shell> *levp;

    /** \brief Control output (default 1)
     */
    int verbose;

    /// The starting neutron levels
    shell neutron_shells[n_internal_levels];
    
    /// The starting proton levels
    shell proton_shells[n_internal_levels];

    /// The grid size
    static const int grid_size=300;
    
    /// The grid step size (default 0.04)
    double step_size;

    /// The nucleon mass (automatically set in init_fun())
    double mnuc;

    /** \name The meson and photon fields and field equations (protected)
     */
    //@{  
    /// Values of the fields from the last iteration
    ubmatrix field0;

    /// The values of the fields
    ubmatrix fields;

    /// The Green's functions inside
    ubmatrix gin;

    /// The Green's functions outside
    ubmatrix gout;

    /// Scalar density RHS
    double sigma_rhs(double sig, double ome, double rho);
    
    /// Vector density RHS
    double omega_rhs(double sig, double ome, double rho);

    /// Isubvector density RHS
    double rho_rhs(double sig, double ome, double rho);

    /** \brief Calculate meson and photon 
	Green's functions \ref gin and \ref gout
    */
    void meson_init();

    /** \brief Calculate meson and photon fields
	
	The meson and photon field equations are of the 
	Sturm-Liouville form, e.g.
	\f[
	\left[r^2 \sigma^{\prime}(r) \right]^{\prime} - 
	r^2 m_{\sigma}^2 \sigma(r) = f(r)
	\f]
	where \f$ \sigma(0) = \sigma_0 \f$ and \f$ \sigma(+\infty) = 0
	\f$. A solution of the homogeneous equation with \f$ f(r) =0
	\f$ is \f$ \sigma(r) = \sigma_0 \sinh( m_{\sigma} r)/
	(m_{\sigma} r) \f$. The associated Green's function is
	\f[
	D(r,r^{\prime},m_{\sigma}) = \frac{-1}{m_{\sigma} r r^{\prime}} 
	\sinh (m_{\sigma} r_{<}) \exp (-m_{\sigma} r_{>})
	\f] 
	where \f$ r_{<} \f$ is the smaller of \f$ r \f$ and 
	\f$ r^{\prime} \f$ and \f$ r_{>} \f$ is the larger.
	One can then write the solution of the meson field
	equation given the density
	\f[
	\sigma(r) = \int_0^{\infty} r^{\prime 2}~d r^{\prime}
	\left[ - g_{\sigma} n_{s}(r) \right] 
	D\left(r,r^{\prime},m_{\sigma}\right)
	\f]

	Since radii are positive, \f$ \sinh (m r) \approx
	e^{m r}/2 \f$ and 
	\f[
	D(r,r^{\prime},m_{\sigma}) \approx \left[
	\frac{-1}{2 m_{\sigma} r_{<}} 
	\exp (m_{\sigma} r_{<}) 
	\right] \left[
	\frac{1}{r_{>}} 
	\exp (-m_{\sigma} r_{>})
	\right]
	\f] 
	The two terms in the Green's function are separated into
	\f[
	\mathrm{gin}(r) = \frac{e^{m_{\sigma} r}}{2 m_{\sigma} r}
	\f]
	and
	\f[
	\mathrm{gout}(r) = \frac{e^{-m_{\sigma} r}}{r} \, .
	\f]
	These functions are computed in \ref meson_init() . Then the
	field is given by
	\f[
	\sigma(r)=  \left(\frac{e^{-m_{\sigma} r}}{r}\right)
	\int_0^{r} r^{\prime 2} g_{\sigma} n_{s} 
	\left(\frac{e^{m_{\sigma} r^{\prime}}}{2 m_{\sigma} 
	r^{\prime}} \right)~d r^{\prime} +
	\left(\frac{e^{m_{\sigma} r}}{2 m_{\sigma} r} \right)
	\int_r^{\infty} r^{\prime 2} g_{\sigma} n_{s} 
	\left(\frac{e^{-m_{\sigma} r^{\prime}}}
	{r^{\prime}}\right)~d r^{\prime}
	\f]
	where the first integral is stored in <tt>xi2</tt> and 
	the second is in <tt>xi1</tt> in the function \ref meson_iter() .
	The part of <tt>xi2</tt> at \f$ r^{\prime}=0 \f$ is 
	stored in <tt>xi20</tt>.
    */
    void meson_iter(double ic);
    
    /** \brief Solve for the meson profiles
     */
    void meson_solve();

    /** \brief The grid index corresponding to the nuclear surface 
	(computed by init_run())
    */
    int surf_index;
    //@}

    /** \name Density information (protected)
     */
    //@{
    /// The densities times radius squared
    ubmatrix xrho;
    
    /// The proton scalar density times radius squared
    ubvector xrhosp; 
    
    /// The scalar field RHS
    ubvector xrhos; 

    /// The vector field RHS
    ubvector xrhov; 

    /// The isubvector field RHS
    ubvector xrhor; 

    /// Charge density
    ubvector chden1; 

    /// Charge density
    ubvector chdenc;

    /// Baryon density
    ubvector arho;
    //@}
    
    /// Energy integrand
    ubvector energy;

    /// Initialize the meson and photon fields, the densities, etc.
    void init_meson_density();

    /// Calculate the energy profile
    void energies(double xpro, double xnu, double e);

    /// Compute the center of mass correction
    void center_mass_corr(double atot);

    /** \name Calculating the form factor, etc. (protected)
     */
    //@{
    /// Fold in proton form factor
    void pfold(double x, double &xrhof);

    /// Function representing proton form factor
    double xpform(double x, double xp, double a);

    /// Perform integrations for form factor
    void gauss(double xmin, double xmax, double x, double &xi);

    /// Desc.
    double xrhop(double x1);

    /// Interpolation object
    interp_vec<ubvector> *gi;
    //@}

    /** \name Solving the Dirac equations (protected)
     */
    //@{
    /** \brief Solve the Dirac equations
	
	Solves the Dirac equation in from 12 fm to the match point and
	then out from .04 fm and adjusts eigenvalue with
	\f[
	\Delta \varepsilon = -g(r=\mathrm{match\_point}) 
	\times (f^{+}-f^{-})
	\f]
    */
    void dirac(int ilevel);
  
    /// Take a step in the Dirac equations
    void dirac_step(double &x, double h, double eigen,
		    double kappa, ubvector &varr);
    
    /// The form of the Dirac equations for the ODE stepper
    int odefun(double x, size_t nv, const ubvector &y,
	       ubvector &dydx, odparms &op);

    /// Compute the fields for the Dirac equations
    void field(double x, double &s, double &v, ubvector &varr);
    
    /// The default stepper
    ode_rkck_gsl<ubvector,ubvector,ubvector,
      ode_funct> def_step;
    
    /// The ODE stepper
    ode_step<ubvector,ubvector,ubvector,
      ode_funct> *ostep;

    /** \brief Integrate the Dirac equations using a simple 
	inline 4th order Runge-Kutta
    */
    double dirac_rk4(double x, double g1, double f1, double &funt, 
		     double eigen, double kappa, ubvector &varr);
    //@}
    
    /// True if init() has been called
    bool init_called;
    
    /// ODE functions
    ubvector ode_y;

    /// ODE derivatives
    ubvector ode_dydx;

    /// ODE errors
    ubvector ode_yerr;

    /// \name Gauss-Legendre integration points and weights
    //@{
    double x12[6], w12[6];
    double x100[50], w100[50];
    //@}
    
#endif
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
