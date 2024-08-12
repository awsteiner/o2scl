/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
/** \file tov_solve.h
    \brief File defining \ref o2scl::tov_solve
*/
#ifndef O2SCL_TOV_SOLVE_H
#define O2SCL_TOV_SOLVE_H

#include <o2scl/eos_tov.h>
#include <o2scl/interp.h>
#include <o2scl/table_units.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/min_brent_gsl.h>

namespace o2scl {

  /** \brief Class to solve the Tolman-Oppenheimer-Volkov equations

      \verbatim embed:rst
      See the :ref:`Solution of the Tolman-Oppenheimer-Volkov Equations`
      section of the User's Guide for
      the mathematical background. 
      \endverbatim

      This class uses adaptive integration to compute the
      gravitational mass, the radius, the baryonic mass (if the EOS
      supplies the baryon density), and the gravitational potential
      (if requested). The equation of state may be changed at any
      time, by specifying the appropriate \ref eos_tov object

      <b>Basic Usage</b>

      After specifying the EOS through \ref tov_solve::set_eos(), one
      calls either \ref tov_solve::mvsr(), \ref tov_solve::fixed(),
      \ref tov_solve::max() or \ref tov_solve::fixed_pr() to compute
      the mass versus radius curve, the profile of a star of a given
      fixed mass, the profile of the maximum mass star, or the profile
      of a star with a fixed central pressure. These functions all
      generate results in the form of a \ref table_units object, which
      can be obtained from \ref tov_solve::get_results().
      
      Screen output:
      - \c verbose=0 - Nothing
      - \c verbose=1 - Basic information
      - \c verbose=2 - For each profile computation, report solution
      information at every kilometer.
      - \c verbose=3 - Report profile information at every 20 grid
      points and output the final interpolation to zero pressure. A
      keypress is required after each profile.

      <b>Mass versus radius curve</b>

      The neutron star mass versus radius curve is constructed by
      computing profiles of all neutron stars with a range of central
      pressures. This range is from \ref prbegin to \ref prend, and
      the ratio between successive central pressures is specified in
      \ref princ (making a logarithmic central pressure grid). 
      
      \note If \ref prend is too close to (or smaller than) the
      central pressure of the maximum mass star, then the mass-radius
      curve generated by mvsr() might not include the maximum mass
      star.

      \note The table generated by mvsr() may include
      unstable configurations (including those with central
      pressures larger than that of the maximum mass star).
      
      <b>Profiles of fixed mass</b>

      Profiles for a fixed gravitational mass are computed by solving
      for the correct central pressure. In order to ensure that the
      solver does not accidentally select a central pressure beyond
      the maximum mass neutron star, the profile with the maximum mass
      is computed beforehand automatically (using \ref max()). The
      intial guess to the solver is always the value of \ref
      fixed_pr_guess, which defaults to \f$ 5.2 \times
      10^{-5}~\mathrm{Msun}/\mathrm{km}^3 \f$ . Alternatively, the
      user can specify the central pressure of the maximum mass star
      so that it does not have to be computed.

      In order to handle multiply-branched mass-radius relations, the
      value of \ref fixed_pr_guess must be specified in order to
      ensure that the correct profile is generated. Even if \ref
      fixed_pr_guess is specified, there is no guarantee that the
      Newton-Raphson will select the correct profile (though it is
      very likely if the guess for the central pressure is
      sufficiently accurate).

      The \ref fixed() function does not support computing profiles
      with central pressures with gravitational masses larger than the
      maximum value. For this, use \ref fixed_pr() .
      
      <b>Profile for the maximum mass star</b>

      This is provided by the function \ref max() . Internally, this
      uses the minimizer specified by \ref set_min() or the default
      minimizer, \ref def_min, to minimize the function \ref max_fun()
      . In order to generate a good initial guess, the \ref max()
      function begins by looping over central pressures from \ref
      max_begin to \ref max_end and choosing the best guess from that
      set of configurations.

      <b>Profile for a fixed central pressure</b>

      This is provided by the function \ref fixed_pr() and
      is relatively fast because it does not require any solving
      or minimizing. This function allows central pressures
      larger than that of the maximum mass star. 

      <b>Output tables</b>

      The functions \ref tov_solve::fixed(), \ref tov_solve::fixed_pr(), 
      and \ref tov_solve::max()
      produce output tables which represent the profile of the neutron
      star of the requested mass. The columns are
      - \c gm, the enclosed gravitational mass in \f$ \mathrm{M}_{\odot} \f$
      - \c r, the radial coordinate in \f$ \mathrm{km} \f$
      - \c gp, the gravitational potential (unitless) when
      \ref calc_gpot is true
      - \c rjw, the value of \f$ g \f$ (see definition below; present
      if \ref ang_vel is true)
      - \c omega_rat, the value of \f$ f \f$ (see definition below; 
      present if \ref ang_vel is true)
      - \c bm, the baryonic mass in \f$ \mathrm{M}_{\odot} \f$ (when 
      \ref eos_tov::baryon_column is true). 
      - \c pr, the pressure in user-specified units
      - \c ed, the energy density in user-specified units
      - \c nb, the baryon density in user-specified units 
      (if \ref eos_tov::baryon_column is true)
      - \c sg, the local surface gravity 
      (in \f$ \mathrm{g}/\mathrm{cm}^{2} \f$ )
      - \c rs, the local redshift (unitless),
      - \c dmdr, the derivative of the enclosed gravitational mass
      in \f$ \mathrm{M}_{\odot}/\mathrm{km} \f$
      - \c dlogpdr, the derivative of the natural logarithm of the 
      pressure
      - \c dgpdr, the derivative of the gravitational potential
      in \f$ 1/\mathrm{km} \f$ (if \ref calc_gpot is true)
      - \c dbmdr, the derivative of the enclosed baryonic mass
      (if \ref eos_tov::baryon_column is true).

      The function \ref tov_solve::mvsr() produces a different kind of
      output table corresponding to the mass versus radius curve. Some
      points on the curve may correspond to unstable configurations.

      - \c gm, the total gravitational mass in \f$ \mathrm{M}_{\odot} \f$
      - \c r, the radius in \f$ \mathrm{km} \f$
      - \c gp, the gravitational potential in the center (unitless) when
      \ref calc_gpot is true
      - \c rjw, the value of \f$ g \f$ at the surface
      (see definition below; present if \ref ang_vel is true)
      - \c omega_rat, the value of \f$ f \f$ at the surface
      (see definition below; present if \ref ang_vel is true)
      - \c bm, total the baryonic mass in \f$ \mathrm{M}_{\odot} \f$ (when 
      \ref eos_tov::baryon_column is true). 
      - \c pr, the central pressure in user-specified units 
      - \c ed, the central energy density in user-specified units 
      - \c nb, the central baryon density in user-specified units 
      (if \ref eos_tov::baryon_column is true)
      - \c sg, the surface gravity 
      (in \f$ \mathrm{g}/\mathrm{cm}^{2} \f$ )
      - \c rs, the redshift at the surface,
      - \c dmdr, the derivative of the gravitational mass
      - \c dlogpdr, the derivative of the natural logarithm of the 
      pressure
      - \c dgpdr, the derivative of the gravitational potential
      in \f$ 1/\mathrm{km} \f$ (if \ref calc_gpot is true)
      - \c dbmdr, the derivative of the enclosed baryonic mass
      (if \ref eos_tov::baryon_column is true).

      <b>Unit systems</b>

      By default, this class operates with energy density and
      pressure in units of \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$
      and baryon density in units of \f$ \mathrm{fm}^{-3} \f$. 

      The function \ref set_units(std::string,std::string,std::string)
      allows one to use different unit systems for energy density,
      pressure, and baryon density. The following list of units of
      energy density and pressure are hard-coded into the library and
      always work:
      - <tt>"g/cm^3"</tt>,
      - <tt>"erg/cm^3"</tt>,
      - <tt>"dyne/cm^2"</tt>,
      - <tt>"MeV/fm^3"</tt>,
      - <tt>"1/fm^4"</tt>, and
      - <tt>"Msun/km^3"</tt> (i.e. \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ )

      The list of hard-coded units for baryon density are:
      - <tt>"1/m^3"</tt>,
      - <tt>"1/cm^3"</tt>, and
      - <tt>"1/fm^3"</tt>.
      
      Other units choices will work if the conversion is either
      already added to the global unit conversion object (from
      <tt>o2scl_settings.get_convert_units()</tt> ) or the global unit
      conversion object is able to compute them by a <tt>system()</tt>
      call to GNU <tt>units</tt> (see documentation in \ref
      convert_units). Note that the choice of what units the tables
      are produced in is independent of the unit system specified in
      the associated \ref eos_tov object, i.e. the input EOS and
      output EOS units need not be the same.

      Alternatively, using \ref set_units(double,double,double) 
      allows one to specify the conversion factors directly without
      using the global unit conversion object.
      
      <b>Accuracy</b>

      The present code, as demonstrated in the tests, gives the
      correct central pressure and energy density of the analytical
      solution by Buchdahl to within less than 1 part in \f$ 10^8 \f$.

      <b>Rotation</b>

      Rotation is considered if \ref tov_solve::ang_vel is set to
      <tt>true</tt>. This adds two more differential equations to
      solve for each central pressure. 

      \verbatim embed:rst
      See also the :ref:`Moment of Inertia in the Slowly-Rotating
      Approximation` section of the User's Guide.
      \endverbatim

      The differential equation for \f$ \bar{\omega} \f$ is independent
      of the relative scale for \f$ \bar{\omega} \f$ and \f$ j \f$ .
      (Note that \f$ j \f$ is a metric function not simply related to
      the angular momentum, \f$ J \f$ .) First, one rescales \f$
      \bar{\omega} \f$ and rewrites everything in terms of \f$ f\equiv
      \bar{\omega}/\Omega \f$ and \f$ g \equiv r^4 j~df/dr \f$ . The
      quantity \f$ f \f$ is unitless and \f$ g \f$ has units of \f$
      \mathrm{km}^3 \f$ . Second, pick a central pressure, \f$
      m(r=0)=g(r=0)=0 \f$, arbitrary values for \f$ \Phi(r=0) \f$ and
      \f$ f(r=0) \f$, and integrate
      \f{eqnarray*}
      \frac{d P}{d r} &=& - \frac{G \varepsilon m}{r^2} 
      \left( 1+\frac{P}{\varepsilon}\right)
      \left( 1+\frac{4 \pi P r^3}{m} \right)
      \left( 1-\frac{2 G m}{r}\right)^{-1}
      \nonumber \\
      \frac{d m}{d r} &=& 4 \pi r^2 \varepsilon
      \nonumber \\
      \frac{d \Phi}{d r} &=& - \frac{1}{\varepsilon}
      \frac{ d P}{d r} \left(1+\frac{P}{\varepsilon}\right)^{-1}
      \nonumber \\
      \frac{d g}{dr} &=& -4 r^3 \frac{d j}{dr} f
      \nonumber \\
      \frac{d f}{dr} &=& \frac{g}{r^4 j}
      \f}
      Afterwards, shift \f$ \Phi \f$ by a constant to ensure
      the correct value at \f$ r=R \f$, and multiply \f$ g \f$
      by a constant to ensure that \f$ g=r^4 j (df/dr) \f$ holds
      for the new potential \f$ \Phi \f$. Then, multiply
      \f$ f \f$ by a constant to ensure that 
      \f[
      f(r=R) + \frac{R}{3} \left(\frac{df}{dr}\right)_{r=R} = 1
      \f]
      
      The functions \f$ f \f$ and \f$ g \f$ are stored in columns
      called <tt>"omega_rat"</tt> and <tt>"rjw"</tt>, respectively.
      One can compute the baryonic mass by integration or by adding
      one additional differential equation, bringing the total to six.

      The moment of inertia is,
      \f[
      I = \frac{R^4}{6 G} \left.\frac{df}{dr}\right|_{r=R}
      \f]
      where the last fraction is stored in \ref domega_rat .
      For an object named \c ts of type \ref tov_solve 
      after a call to \ref fixed(), \ref max(), or
      \ref fixed_pr(), the moment of inertia can be computed
      with, e.g. 
      \code
      tov_solve ts;
      ts.max();
      double schwarz_km=o2scl_cgs::schwarzchild_radius/1.0e5;
      double I=ts.domega_rat*pow(ts.rad,4.0)/3.0/schwarz_km;
      \endcode

      After a call to \ref mvsr(), the values of \f$ f(r=R) \f$
      and \f$ g(r=R) \f$ are stored in columns labeled
      <tt>"omega_rat"</tt> and <tt>"rjw"</tt>. The moment
      of inertia of a 1.4 solar mass neutron star 
      can be computed with, e.g. 
      \code
      tov_solve ts;
      ts.mvsr();
      std::shared_ptr<table_units<> > tab=ts.get_results();
      double schwarz_km=o2scl_cgs::schwarzchild_radius/1.0e5;
      double I_14=tab->interp("gm",1.4,"rjw")/3.0/schwarz_km;
      \endcode
      
      <b>Convergence details</b>

      By default, if the TOV solver fails to converge, the error
      handler is called and an exception is thrown. If \ref
      o2scl::tov_solve::err_nonconv is false, then \ref
      o2scl::tov_solve::mvsr(), \ref o2scl::tov_solve::fixed(), and
      \ref o2scl::tov_solve::max(), return an integer which gives some
      information about why the solver failed to converge.

      If \ref err_nonconv is set to false, then the \ref fixed()
      function temporarily sets the value of both \ref
      o2scl::mroot::err_nonconv for the current solver and \ref
      o2scl::jacobian::err_nonconv for the jacobian object, \ref
      o2scl::mroot_hybrids::def_jac, of \ref def_solver equal to false.

      <b>Other details</b>

      The ODE solution is stored in a buffer which can be directly
      accessed using \ref o2scl::tov_solve::get_rkx(), \ref
      o2scl::tov_solve::get_rky(), and \ref
      o2scl::tov_solve::get_rkdydx(). In the case that the ODE steps
      are small enough that there isn't enough space in the buffers,
      the ODE solution data is thinned by a factor of two to allow for
      the remaining ODE steps to be stored. The size of the buffers
      can be controlled in \ref buffer_size .

      If \ref o2scl::tov_solve::reformat_results is true (the
      default), then the results are placed in a shared pointer to the
      \ref table_units object which can be accessed using \ref
      o2scl::tov_solve::get_results(). The maximum size of the output
      table can be controlled with \ref max_table_size. The output
      table may be smaller than this, as it cannot be larger than the
      number of steps stored in the buffer.

      \note The function \ref o2scl::tov_solve::integ_star() returns
      <tt>gsl_efailed</tt> without calling the error handler in the
      case that the solver can recover gracefully from, for example, a
      negative pressure.

      \future Convert to \ref o2scl::ode_iv_solve?
  */
  class tov_solve {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
    typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

#ifndef DOXYGEN_INTERNAL

  protected:

    /// ODE function object
    ode_funct ofm;

    /// Interpolation object for listed radii in \ref mvsr()
    interp_vec<ubvector> iop;

    /** \brief Set up column names and units

	When this function is used by \ref mvsr(), \c mvsr_mode is set
	to true. 
     */
    void column_setup(bool mvsr_mode=false);

    /** \brief If true, \ref integ_star() computes all the profile info,
	otherwise it only computes the gravitational mass
     */
    bool integ_star_final;

    /// The last row index in \ref rky
    size_t ix_last;

    /// The schwarzchild radius in km
    double schwarz_km;
    
    /** \brief Target mass for integ_star()
	
	Has a value of zero, unless set to a non-zero value 
	by \ref fixed(). 
    */
    double tmass;
    
    /** \brief Ensure \c col does not match strings in \c cnames
	
	Underscores are added to \c col until it matches none of
	the strings in \c cnames.
    */
    void make_unique_name(std::string &col, 
			  std::vector<std::string> &cnames);

    /// \name User EOS
    //@{
    /// The EOS
    eos_tov *te;

    /// True if the EOS has been set
    bool eos_set;
    //@}

    /// \name Units for output table
    //@{
    /// Units for energy density
    std::string eunits;
    /// Units for pressure
    std::string punits;
    /// Units for baryon density
    std::string nunits;
    /// unit conversion factor for energy density (default 1.0)
    double efactor;
    /// unit conversion factor for pressure (default 1.0)
    double pfactor;
    /// unit conversion factor for baryon density (default 1.0)
    double nfactor;
    //@}

    /** \brief Smallest allowed pressure for integration (default: -100)
   
	This quantity can't be much smaller than -100 since we need to
	compute numbers near \f$ e^{-\mathrm{min\_log\_pres}} \f$

	\future Replace this with the proper value from std::limits?
    */
    double min_log_pres;

    /// \name Integration storage
    //@{
    /// Radial coordinate (in kilometers)
    ubvector rkx;
    /** \brief ODE functions

	If \c rky is viewed as a matrix, then the first column of each
	row is the gravitational mass in solar masses, and the second
	column is the natural logarithm of the pressure in \f$
	\mathrm{M}_{\odot}/km^3 \f$ . When \ref calc_gpot is true, the
	next column is the gravitational potential (which is
	unitless), and when \ref eos_tov::baryon_column is true, the
	next column is the baryonic mass in \f$ \mathrm{M}_{\odot}
	\f$.
    */
    std::vector<ubvector> rky;
    /// The derivatives of the ODE functions
    std::vector<ubvector> rkdydx;
    //@}

    /// The output table
    std::shared_ptr<table_units<> > out_table;

    /// \name Numerical methods
    //@{
    /// The solver for fixed gravitational masses
    mroot<mm_funct,ubvector,jac_funct> *mroot_ptr;

    /// The minimizer for maximum masses
    min_base<> *min_ptr;
    
    /// The adaptive stepper
    astep_base<ubvector,ubvector,ubvector,ode_funct> *as_ptr;
    //@}

    /// The ODE step function
    virtual int derivs(double x, size_t nv, const ubvector &y,
		       ubvector &dydx);
    
    /// The minimizer function to compute the maximum mass
    virtual double max_fun(double maxx);

    /** \brief The solver function to compute the stellar profile
     */
    virtual int integ_star(size_t ndvar, const ubvector &ndx, 
			ubvector &ndy);
    
#endif

  public:

    tov_solve();

    virtual ~tov_solve();
    
    /// Size of the ODE solution buffer (default \f$ 10^5 \f$)
    size_t buffer_size;

    /// Maximum number of lines in the output table (default 400)
    size_t max_table_size;

    /// \name Basic properties
    //@{
    /// Gravitational mass (in \f$ \mathrm{M}_{\odot} \f$)
    double mass;
    /// Radius (in km)
    double rad;
    /// Baryonic mass (when computed; in \f$ \mathrm{M}_{\odot} \f$)
    double bmass;
    /// Gravitational potential (when computed; unitless)
    double gpot;
    /// The value of \f$ r^4 j df / dr \f$ at the surface
    double last_rjw;
    /// The value of \f$ \bar{\omega} / \Omega \f$ at the surface
    double last_f;
    /** \brief The value of \f$ d (\bar{\omega}/\Omega)/dr \f$ 
	at the surface (when \ref ang_vel is true)
    */
    double domega_rat;

    /** \brief Maximum value for central pressure in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (default \f$ 10^{20} \f$ )
	
	This variable is set by the <tt>mvsr()</tt>, <tt>max()</tt>,
	<tt>fixed()</tt> and <tt>fixed_pr()</tt>
	functions and used in \ref integ_star() .
    */
    double pcent_max;
    //@}

    /** \name Solution parameters
	
	These parameters can be changed at any time.
    */
    //@{
    /** \brief If true, then fixed() and max() reformat results into
	a \ref o2scl::table_units object
	
	Note that \ref mvsr() always places its results in the
	output table independent of the value of this variable.
    */
    bool reformat_results;
    /** \brief The mass of one baryon
	
	The mass of one baryon in kg for the total baryon mass
	calculation (defaults to the proton mass).
    */
    double baryon_mass;
    /// If true, solve for the angular velocity (default false)
    bool ang_vel;
    /// Apply general relativistic corrections (default true)
    bool gen_rel;
    /** \brief calculate the gravitational potential (default false)
    */
    bool calc_gpot;
    /// smallest allowed radial stepsize in km (default 1.0e-4)
    double step_min;
    /// largest allowed radial stepsize in km (default 0.05)
    double step_max;
    /// initial radial stepsize in km (default 4.0e-3)
    double step_start;
    /// control for output (default 1)
    int verbose;
    /// Maximum number of integration steps (default 100000)
    size_t max_integ_steps;
    /** \brief If true, call the error handler if the solution does 
	not converge (default true)
    */
    bool err_nonconv;
    //@}

    /** \brief Default value of maximum pressure for maximum mass star
	in \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ 
	(default \f$ 10^{20} \f$)
     */
    double pmax_default;

    /// \name Mass versus radius parameters
    //@{
    /** \brief Beginning pressure in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (default 7.0e-7)
    */
    double prbegin;
    /** \brief Ending pressure in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (default 8.0e-3)
    */
    double prend;
    /// Increment factor for pressure (default 1.1)
    double princ;
    /** \brief List of pressures at which more information should be
	recorded
	
	If pressures (in the user-specified units) are added to this
	vector, then in mvsr(), the radial location, enclosed
	gravitational mass, and (if \ref o2scl::eos_tov::baryon_column
	is true) enclosed baryon mass are stored in the table for each
	central pressure. The associated columns are named 
	<tt>r0, gm0, bm0, r1, gm1, bm1,</tt> etc.
    */
    std::vector<double> pr_list;
    //@}

    /// \name Fixed mass parameter
    //@{
    /** \brief Guess for central pressure in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$
	(default \f$ 5.2 \times 10^{-5} \f$)

	This guess is used in the functions fixed().
    */
    double fixed_pr_guess;
    //@}

    /// \name Maximum mass profile parameters
    //@{
    /** \brief Beginning pressure for maximum mass guess in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (default 7.0e-5)
    */
    double max_begin;
    /** \brief Ending pressure for maximum mass guess in 
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (default 5.0e-3)
    */
    double max_end;
    /** \brief Increment for pressure for maximum mass guess 
	(unitless; default 1.3)
    */
    double max_inc;
    //@}

    /// \name Basic operation
    //@{
    /// Set the EOS to use
    void set_eos(eos_tov &ter) {
      te=&ter;
      eos_set=true;
      return;
    }

    /** \brief Set output units for the table
     */
    void set_units(double s_efactor=1.0, double s_pfactor=1.0, 
		   double s_nbfactor=1.0);
    
    /** \brief Set output units for the table
    */
    void set_units(std::string eunits="", std::string punits="", 
		   std::string nunits="");

    /** \brief Get output units for the table
    */
    void get_units(std::string &eunits, std::string &punits, 
		   std::string &nunits);

    /** \brief Calculate the mass vs. radius curve
     */
    virtual int mvsr();

    /** \brief Calculate the profile of a star with fixed mass

	This function computes the profile for a star with a fixed
	mass. If the target mass is negative, it is interpreted as
	subtracting from the maximum mass configuration. For a 2.15
	solar mass neutron star, <tt>target_mass=-0.2</tt> corresponds
	to 1.95 solar mass configuration.

	The variable \c pmax is the maximum allowable central pressure
	in \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ (This often, but
	not always, equal to the central pressure of the maximum mass
	star.) This ensures that the function does not unintentionally
	select an unstable configuration. If \c pmax is greater than
	or equal to the default value (\ref pmax_default), then the
	maximum mass star will be computed with \ref max() first
	in order to determine the maximum allowable central pressure.

	Note that this function will likely fail when the mass-radius
	curve has two central pressures with the same gravitational
	mass.
    */
    virtual int fixed(double target_mass, double pmax=1.0e20);

    /** \brief Calculate the profile of a star with a specified
	central pressure

	This function computes the profile of a star with a fixed
	central pressure. The central pressure, \c pcent, should be in
	the unit system specified by the user which defaults to solar
	masses per cubic kilometer "Msun/km^3" but can be changed with
	a call to one of the <tt>set_units()</tt> functions.

	The variable \c pmax is the maximum allowable central pressure
	in \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$, and must 
	be larger than the value of \c pcent converted to to
	\f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ . 
    */
    virtual int fixed_pr(double pcent, double pmax=1.0e20);
    
    /** \brief Calculate the profile of the maximum mass star
	
	Note that this function maximizes the gravitational mass. If
	the M-R curve has two stable branches, then this function does
	not necessarily give the configuration with the largest
	central pressure.

	This function may also depend on the accuracy of the initial
	interval determined by \ref max_begin and \ref max_end.
     */
    virtual int max();

    /** \brief Construct a table from the results
	
	This function constructs a \ref table_units object from the
	information in \ref rkx, \ref rky, and \ref rkdydx . Note that
	the table is always constructed by default so this function
	need not be called unless \ref tov_solve::reformat_results is
	<tt>false</tt>.>
     */
    virtual void make_table();
    
    /** \brief Return the results data table
     */
    std::shared_ptr<table_units<> > get_results() {
      return out_table;
    }

    /** \brief Return the results data table

	This function immediately adds four constants to the table,
	<tt>schwarz, Msun, pi</tt> and <tt>mproton</tt>.
     */
    void set_table(std::shared_ptr<table_units<> > t) {
      out_table=t;
      // Add native constants
      out_table->add_constant("schwarz",schwarz_km);
      out_table->add_constant("Msun",o2scl_const::solar_mass_f<double>());
      out_table->add_constant("pi",o2scl_const::pi);
      out_table->add_constant("mproton",o2scl_const::mass_proton_f<double>());
      return;
    }
    //@}
    
    /// \name Convergence information flags
    //@{
    static const int fixed_solver_failed=128;
    static const int fixed_integ_star_failed=256;
    static const int fixed_wrong_mass=512;
    static const int max_minimizer_failed=1024;
    static const int max_integ_star_failed=2048;
    static const int mvsr_integ_star_failed=4096;
    static const int ang_vel_failed=8192;
    static const int cent_press_large=16384;
    static const int cent_press_neg=32768;
    static const int over_max_steps=65536;
    static const int last_step_large=131072;
    //@}

    /// \name Get internally formatted results (in internal unit system)
    //@{
    /// Get the vector for the radial grid
    const ubvector &get_rkx() const {
      return rkx;
    }
    /// Get a list of vectors for the ODEs
    const std::vector<ubvector> &get_rky() const {
      return rky;
    }
    /// Get a list of vectors for the corresponding derivatives
    const std::vector<ubvector> &get_rkdydx() const {
      return rkdydx;
    }
    //@}

    /// \name Control numerical methods
    //@{
    /** \brief Set solver
    */
    void set_mroot(mroot<mm_funct,ubvector,jac_funct> &s_mrp) {
      mroot_ptr=&s_mrp; 
      return;
    };
    
    /** \brief Set minimizer
    */
    void set_min(min_base<> &s_mp) {
      min_ptr=&s_mp; 
      return;
    };

    /** \brief Set the adaptive stepper
     */
    void set_stepper(astep_base<ubvector,ubvector,ubvector,ode_funct> &sap) {
      as_ptr=&sap; 
      return;
    };
    //@}

    /// \name Default numerical objects
    //@{
    /// The default minimizer
    min_brent_gsl<> def_min;
    
    /// The default solver
    mroot_hybrids<mm_funct,ubvector,ubmatrix,jac_funct> def_solver;

    /// The default adaptive stepper
    astep_gsl<ubvector,ubvector,ubvector,ode_funct> def_stepper;
    //@}

  };

}

#endif


