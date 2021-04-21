/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2021, Andrew W. Steiner
  
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
#ifndef TOV_LOVE_H
#define TOV_LOVE_H

#include <gsl/gsl_math.h>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/astep_gsl.h>
#include <o2scl/table_units.h>
#include <o2scl/astep_nonadapt.h>
#include <o2scl/ode_rk8pd_gsl.h>
#include <o2scl/ode_iv_solve.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Determination of the neutron star Love number

      We use \f$ c=1 \f$ but keep factors of \f$ G \f$, which has
      units \f$ \mathrm{km}/\mathrm{M_{\odot}} \f$.

      \warning The Love number can be sensitive to inaccuracy in
      the input EOS table, thus it's important to make sure the
      EOS table has a sufficiently fine grid to give accurate
      results.

      \verbatim embed:rst
      Following the notation in [Postnikov10]_, define
      the function :math:`H(r)`, which is the solution of
      \endverbatim
      \f[
      H^{\prime \prime} (r) + H^{\prime}(r) \left\{
      \frac{2}{r} + e^{\lambda(r)} \left[ \frac{2 G m(r)}{r^2} +
      4 \pi G r P(r) - 4 \pi G r \varepsilon(r) \right]
      \right\} + H(r) Q(r) = 0 
      \f]
      where (now supressing the dependence on \f$ r \f$), 
      \f[
      \nu^{\prime} \equiv 2 G e^{\lambda} 
      \left(\frac{m+4 \pi P r^3}{r^2}\right) \,
      \f]
      which has units of \f$ 1/\mathrm{km} \f$ ,
      \f[
      e^{\lambda} \equiv \left(1-\frac{2 G m}{r}\right)^{-1} \, ,
      \f]
      and 
      \f[
      Q \equiv 4 \pi G e^{\lambda} \left( 5 \varepsilon + 9 P + 
      \frac{\varepsilon+P}{c_s^2}\right) - 6 \frac{e^{\lambda}}{r^2}
      - \nu^{\prime 2}
      \f]
      which has units of \f$ 1/\mathrm{km}^2 \f$ .
      The boundary conditions on \f$ H(r) \f$ are that \f$ H(r) = a_0
      r^2 \f$ and \f$ H^{\prime} = 2 a_0 r \f$ for an arbitrary
      constant \f$ a_0 \f$ (\f$ a_0 \f$ is chosen to be equal to 1).
      Internally, \f$ P \f$ and \f$ \varepsilon \f$ are stored in
      units of \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ .

      From this we can define another (unitless) function 
      \f$ y(r) \equiv r H^{\prime}(r)/H(r) \f$, which obeys
      \f[
      r y^{\prime} + y^2 + y e^{\lambda} \left[
      1+4 \pi G r^2 \left( P-\varepsilon \right)
      \right] + r^2 Q = 0 
      \f]
      with boundary condition is \f$ y(0) = 2 \f$ .
      Solving for \f$ y^{\prime} \f$,
      \f[
      y^{\prime} = \frac{1}{r} 
      \left\{-r^2 Q - y e^{\lambda} \left[ 1+ 4 \pi G r^2
      \left( P - \varepsilon \right) \right] -y^2 \right\}
      \f]
      Define \f$ y_R = y(r=R) \f$. This form for \f$ y^{\prime}(r) \f$
      is specified in \ref y_derivs() . 

      The unitless quantity \f$ k_2[\beta,y_R] \f$ (the Love number) 
      is defined by (this is the expression from Postnikov et al. (2010) )
      \f{eqnarray*}
      k_2[\beta,y(r=R)] &\equiv& \frac{8}{5} \beta^5 
      \left(1-2 \beta\right)^2
      \left[ 2 - y_R + 2 \beta \left( y_R - 1 \right) \right] 
      \nonumber \\
      && \times \left\{ 2 \beta \left( 6 - 3 y_R + 3 \beta ( 5 y_R - 8)
      + 2 \beta^2 \left[ 13 - 11 y_R + \beta (3 y_R-2) 
      \right.\right.\right. \nonumber \\
      && + \left.\left.\left. 2 \beta^2 (1+y_R) \right] \right) + 3 
      (1-2 \beta)^2 \left[ 2 - y_R + 2 \beta (y_R - 1) \right] 
      \log (1-2 \beta) \right\}^{-1}
      \f}

      \verbatim embed:rst
      [Hinderer10]_ writes the differential equation for 
      :math:`H(r)`
      in a slightly different (but equivalent) form,
      \endverbatim
      \f{eqnarray*}
      H^{\prime \prime}(r) &=& 2 \left( 1 - \frac{2 G m}{r}\right)^{-1} 
      H(r) \left\{ - 2 \pi G \left[ 5 \varepsilon + 9 P + \frac{\left(
      \varepsilon+P\right)}{c_s^2} \right] + \frac{3}{r^2}
      \right. \nonumber \\ && \left. + 
      2 \left( 1 - \frac{2 G m}{r}\right)^{-1} 
      \left(\frac{G m}{r^2}+4 \pi G r P\right)^2 \right\}
      +\frac{2 H^{\prime}(r)}{r} \left( 1 - \frac{2 G m}{r}\right)^{-1}
      \nonumber \\ && \times
      \left[ -1+\frac{G m}{r} + 2 \pi G r^2 \left(\varepsilon-P\right)
      \right] \, .
      \f}
      This is the form given in \ref H_derivs() .
      
      The tidal deformability is then 
      \f[
      \lambda \equiv \frac{2}{3} k_2 R^5
      \f]
      and has units of \f$ \mathrm{km}^5 \f$ or can be converted to
      \f$ \mathrm{g}~\mathrm{cm}^2~\mathrm{s}^2 \f$ .

      It is assumed that \ref tab stores a stellar profile (such as
      one computed with \ref tov_solve::fixed(), \ref
      tov_solve::fixed_pr(), or \ref tov_solve::max() ) has been
      specified before-hand and contains (at least) the following
      columns
      - <tt>ed</tt> energy density in units of 
      \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$
      - <tt>pr</tt> pressure in units of 
      \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$
      - <tt>cs2</tt> sound speed squared (unitless)
      - <tt>gm</tt> gravitational mass in \f$ \mathrm{M}_{\odot} \f$
      - <tt>r</tt> radius in \f$ \mathrm{km} \f$
      (Note that the \ref o2scl::tov_solve class doesn't automatically compute
      the column <tt>cs2</tt>.)

      This class handles the inner boundary by starting from the small
      non-zero radius stored in \ref eps instead of at \f$ r=0 \f$. The
      value of \ref eps defaults to 0.2 km.
      
      \verbatim embed:rst
      If there is a discontinuity in the EOS (i.e. a jump in 
      the energy density at some radius :math:`r_d`), then 
      the function :math:`y(r)` must satisfy (see [Damour09]_,
      [Postnikov10]_, and [Hinderer10]_).
      \endverbatim
      \f[
      y(r_d+\delta) - y(r_d-\delta) =
      \frac{ 
      \varepsilon(r_d+\delta)-\varepsilon(r_d-\delta)}{m(r_d)/(4 \pi r_d^3) + 
      p}
      \f]
      
      \note The function \ref calc_H() cannot yet handle 
      discontinuities (if there are any then the error handler
      is called). 

      \future Improve calc_H() to handle discontinuities and to 
      tabulate the EOS.
      \future Improve the handling at small r using an expansion,
      similar to that used in e.g. Detweiler and Lindblom (1985)?
      \future Allow specification of, e.g., an eos_tov like object
      rather than an EOS table.
  */
  class tov_love {

  public:

    /// The ODE function type
    typedef std::function<int(double,size_t,
			      const std::vector<double> &,
			      std::vector<double> &)> ode_funct2;
    
#ifndef DOXYGEN_INTERNAL
  
  protected:

    /// A pointer to the ODE integrator
    o2scl::ode_iv_solve<ode_funct2,std::vector<double> > *oisp;

    /** \brief The derivative \f$ y^{\prime}(r) \f$
     */
    int y_derivs(double r, size_t nv, const std::vector<double> &vals,
		 std::vector<double> &ders);

    /** \brief The derivatives \f$ H^{\prime \prime}(r) \f$ and
	\f$ H^{\prime}(r) \f$
    */
    int H_derivs(double r, size_t nv, const std::vector<double> &vals,
		 std::vector<double> &ders);

    /// Schwarzchild radius in km (set in constructor)
    double schwarz_km;
  
    /** \brief Compute \f$ k_2(\beta,y_R) \f$ using the analytic 
	expression

	Used in both \ref tov_love::calc_y() and \ref
	tov_love::calc_H().
    */
    double eval_k2(double beta, double yR);

    /// List of discontinuities
    std::vector<double> disc;
    
#endif

  public:
  
    tov_love();

    /** \brief If greater than zero, show the ODE output (default 0)
     */
    int show_ode;
    
    /** \brief Additional testing if the ODE solver fails
     */
    bool addl_testing;
    
    /** \brief If true, call the error handler if the solution does 
	not converge (default true)
    */
    bool err_nonconv;
    
    /** \brief A table containing the solution to the differential 
        equation 

        This table is filled when \ref calc_y() is called with 
        <tt>tabulate=true</tt>.
    */
    o2scl::table_units<> results;

    /** \brief The radial step for resolving discontinuities in km 
	(default \f$ 10^{-4} \f$)
    */
    double delta;
    
    /// The first radial point in \f$ \mathrm{km} \f$ (default 0.02)
    double eps;

    /// The default ODE integrator
    o2scl::ode_iv_solve<ode_funct2,std::vector<double> > def_ois;

    /// Pointer to the input profile
    std::shared_ptr<o2scl::table_units<> > tab;

    /// Set ODE integrator
    void set_ODE(o2scl::ode_iv_solve<ode_funct2,std::vector<double> >
		 &ois_new) {
      oisp=&ois_new;
    }
    
    /** \brief Compute the Love number using y
        
        The initial values of \c yR, \c beta, \c k2, \c lambda_km5,
        and \c lambda_cgs are ignored. 

        If \c tabulate is true, then the results of the ODE
        integration are stored in \ref results . The \ref
        o2scl::table::clear() function is called beforehand, so any
        data stored in in \ref results is destroyed.

        \note The final results for the tidal deformability may differ
        when \c tabulate is true versus when \c tabulate is false.
        This difference is likely within the uncertainty of the ODE
        integration.
     */
    int calc_y(double &yR, double &beta, double &k2, double &lambda_km5,
		double &lambda_cgs, bool tabulate=false);

    /** \brief Add a discontinuity at radius \c rd (in km)
     */
    void add_disc(double rd) {
      disc.push_back(rd);
    }

    /** \brief Remove all discontinuities
     */
    void clear_discs() {
      disc.clear();
    }
    
    /** \brief Compute the love number using H

        The initial values of \c yR, \c beta, \c k2, \c lambda_km5,
        and \c lambda_cgs are ignored. 
     */
    int calc_H(double &yR, double &beta, double &k2, double &lambda_km5,
		double &lambda_cgs);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
