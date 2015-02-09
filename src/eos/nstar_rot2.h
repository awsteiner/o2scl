/*
  -------------------------------------------------------------------
  
  This file is part of O2scl. It has been adapted from RNS v1.1d
  written by N. Stergioulas and S. Morsink. The modifications made in
  this version from the original are copyright (C) 2015, Andrew W.
  Steiner.
  
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
/*
  -------------------------------------------------------------------
  Relativistic models of rapidly rotating compact stars,
  using tabulated or polytropic equations of state.
  
  Author:  Nikolaos Stergioulas
  
  Current Address:
  
  Department of Physics
  University of Wisconsin-Milwaukee
  PO Box 413, Milwaukee, WI 53201, USA
  
  E-mail: niksterg@csd.uwm.edu, or
  niksterg@pauli.phys.uwm.edu
  
  Version: 1.1
  
  Date:    June, 1995
  
  Changes made to code by Sharon Morsink
   
  03-03-97: Corrected the units for polytropic stars
  10-28-98: Added the star's quadrupole moment to the output.
  
  References:
  KEH : H. Komatsu, Y. Eriguchi and I. Hachisu, Mon. Not. R. astr. Soc. 
  (1989) 237, 355-379.
  CST : G. Cook, S. Shapiro and S. Teukolsky, Ap. J (1992) 398, 203-223.
  
  -------------------------------------------------------------------
*/
/** \file nstar_rot2.h
    \brief File defining \ref o2scl::nstar_rot2
*/
#ifndef NSTAR_ROT2_H
#define NSTAR_ROT2_H

#include <cmath>
#include <iostream>

#include <o2scl/err_hnd.h>
#include <o2scl/search_vec.h>
#include <o2scl/test_mgr.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/lib_settings.h>
#include <o2scl/interp.h>
#include <o2scl/eos_tov.h>

namespace o2scl {
  
  /** \brief An EOS for \ref nstar_rot2
   */
  class eos_nstar_rot : public eos_tov {
    
  public:

    /** \brief From the pressure, return the enthalpy
     */
    virtual double enth_from_pr(double pr)=0;

    /** \brief From the enthalpy, return the pressure
     */
    virtual double pr_from_enth(double enth)=0;

    /** \brief From the baryon density, return the enthalpy
     */
    virtual double enth_from_nb(double nb)=0;

    void ed_nb_from_pr(double pr, double &ed, double &nb) {
      return;
    }
  };
  
  /** \brief
   */
  class eos_nstar_rot_interp : public eos_nstar_rot {
    
  protected:

    /// Array search object
    o2scl::search_vec<double *> sv;
    
    /// Search in array \c x of length \c n for value \c val
    int new_search(int n, double *x, double val);

    
    /** \brief number of tabulated EOS points */
    int n_tab;                           
    /** \brief rho points in tabulated EOS */
    double log_e_tab[201];               
    /** \brief p points in tabulated EOS */
    double log_p_tab[201];               
    /** \brief h points in EOS file */
    double log_h_tab[201];               
    /** \brief number density in EOS file */  
    double log_n0_tab[201];              

    /// \name Constants
    //@{
    /** \brief Speed of light in vacuum (in CGS units) */ 
    double C;
    /** \brief Gravitational constant (in CGS units) */ 
    double G;
    /** \brief Square of length scale in CGS units, 
	\f$ \kappa \equiv 10^{-15} c^2/G \f$
    */
    double KAPPA;
    /** \brief The value \f$ \kappa G c^{-4} \f$ */
    double KSCALE;
    //@}

    /** \brief Cache for interpolation
     */
    int n_nearest;
    
    /** \brief Driver for the interpolation routine. 
	
	First we find the tab. point nearest to xb, then we
	interpolate using four points around xb.
    */  
    double interp(double xp[], double yp[], int np, double xb);
    
  public:
    
    eos_nstar_rot_interp();
    
    /** \brief Set the EOS from four vectors in the native unit system
     */
    template<class vec1_t, class vec2_t, class vec3_t, class vec4_t>
      void set_eos_native(vec1_t &eden, vec2_t &pres, vec3_t &enth,
			  vec4_t &nb) {
      
      double C=o2scl_cgs::speed_of_light;
      double G=o2scl_cgs::gravitational_constant;
      double KAPPA=1.0e-15*C*C/G;
      double KSCALE=KAPPA*G/(C*C*C*C);
  
      n_tab=eden.size();

      for(int i=1;i<=n_tab;i++) {  
	log_e_tab[i]=log10(eden[i-1]*C*C*KSCALE);
	log_p_tab[i]=log10(pres[i-1]*KSCALE);
	log_h_tab[i]=log10(enth[i-1]/(C*C));
	log_n0_tab[i]=log10(nb[i-1]);
      }
      
      return;
    }

    /** \brief Set the EOS from energy density, pressure, and
	baryon density stored in powers of \f$ \mathrm{fm} \f$ .
    */
    template<class vec1_t, class vec2_t, class vec3_t>
      void set_eos_fm(size_t n, vec1_t &eden, vec2_t &pres, vec3_t &nb) {
      
      if (n>200) {
	O2SCL_ERR2("Too many EOS points in ",
		   "nstar_rot::set_eos().",o2scl::exc_einval);
      }

      // Conversion factor for energy density
      double conv1=o2scl_settings.get_convert_units().convert
	("1/fm^4","g/cm^3",1.0);
      // Conversion factor for pressure
      double conv2=o2scl_settings.get_convert_units().convert
	("1/fm^4","dyne/cm^2",1.0);
      
      n_tab=n;

      double mu0=(eden[0]+pres[0])/nb[0];
      double mu1=(eden[1]+pres[1])/nb[1];
      double mu_start=2.0*mu0-mu1;
      
      for(size_t i=0;i<n_tab;i++) {
	log_e_tab[i+1]=log10(eden[i]*conv1*C*C*KSCALE);
	log_p_tab[i+1]=log10(pres[i]*conv2*KSCALE);
	log_n0_tab[i+1]=log10(nb[i]*1.0e39);
	log_h_tab[i+1]=log10(log((eden[i]+pres[i])/nb[i])-log(mu_start));
      }

      return;
    }

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the enthalpy
     */
    virtual double enth_from_pr(double pr);

    /** \brief From the baryon density, return the enthalpy
     */
    virtual double enth_from_nb(double nb);

    /** \brief From the enthalpy, return the pressure
     */
    virtual double pr_from_enth(double enth);

  };
  
  /** \brief EOS C
   */
  class eos_nstar_rot_C : public eos_nstar_rot_interp {
  public:
    eos_nstar_rot_C();
  };
  
  /** \brief EOS L
   */
  class eos_nstar_rot_L : public eos_nstar_rot_interp {
  public:
    eos_nstar_rot_L();
  };
  
  /** \brief Rotating neutron star class based on RNS v1.1d from
      N. Stergioulas et al.
      
      \note This class is still experimental.

      Several changes have been made to the original code. The code
      using Numerical Recipes has been removed and replaced with an
      equivalent based on GSL and \o2. The overall interface has
      been changed and some code has been updated with C++
      equivalents.

      <b>Usage</b>

      <b>Initial guess</b>

      The original RNS code suggests that the initial guess is
      typically a star with a smaller angular momentum.

      <b>References</b> 

      The original RNS v1.1d can be obtained from
      http://www.gravity.phys.uwm.edu/rns/ , and you may find Nick
      Stergioulas's web page http://www.astro.auth.gr/~niksterg/ , or
      Sharon Morsink's page http://fermi.phys.ualberta.ca/~morsink/
      useful. See \ref Bonazzola73, \ref Bonazzola94, \ref Cook92,
      \ref Cook94, \ref Friedman88, \ref Gourgoulhon94, \ref
      Komatsu89, \ref Laarakkers99, \ref Nozawa98, \ref Stergioulas95,
      and \ref Stergioulas03 .

      \todo Better documentation is needed everywhere.

      \future Rework EOS interface and implement better 
      integration with the other \o2e EOSs. 
      \future Fix unit-indexed arrays.
      \future Try moving some of the storage to the heap?
      \future Some of the arrays seem larger than necessary.
      \future The function \ref o2scl::nstar_rot::new_search() is
      inefficient because it has to handle the boundary conditions
      separately. This could be improved.
      \future Give the user more control over the initial guess.
      \future Remove the CL_LOW stuff?

      \comment
      <b>Quadrupole moments</b>

      Quadrupole moments computed using the method in \ref Laarakkers99. 

      <b>Axisymmetric Instability</b>

      \ref Friedman88 shows that a secular axisymmetric instability
      sets in when the mass becomes maximum along a sequence of
      constant angular momentum. Equivalently, \ref Cook92 shows that
      the instability occurs when the angular momentum becomes minimum
      along a sequence of constant rest mass.

      A GR virial theorem for a stationary and axisymmetric system was
      found in \ref Bonazzola73. A more general two-dimensional virial
      identity was found in \ref Bonazzola94. The three-dimensional
      virial identity found in \ref Gourgoulhon94 is a generalization
      of the Newtonial virial theorem.

      Using the stationary and axisymmetric metric ( \f$ G = c = 1 \f$
      )
      \f[
      ds^2 = - e^{\gamma+\rho} dt^2 + e^{2 \alpha} \left( dr^2 + 
      r^2 d\theta^2 \right) + e^{\gamma-\rho} r^2 \sin^2 \theta
      ( d \phi - \omega dt) ^2
      \f]
      one solves for the four metric functions \f$ \rho(r,\theta) \f$,
      \f$ \gamma(r,\theta) \f$, \f$ \alpha(r,\theta) \f$ and \f$
      \omega(r,\theta) \f$ .

      It is assumed that matter is a perfect fluid, and the 
      stress-energy tensor is
      \f[
      T^{\mu \nu} = \left( \rho_0 + \rho_i + P \right) u^{\mu} u^{\nu}
      + P g^{\mu \nu}
      \f]
      
      Einstein's field equations imply four field equations for
      a specified rotation law,
      \f[
      u^{t} u_{\phi} = F(\Omega) 
      \f]
      for some function \f$ F(\omega) \f$ .

      Using Eq. (27) in \ref Cook92, one can write
      \f[
      \rho(s,\mu) = - e^{-\gamma/2} \sum_{n=0}^{\infty}
      P_{2n}(\mu) \int_0^{1}~ds^{\prime} \int_0^1~d \mu 
      f_{\rho}(n,s,s^{\prime}) P_{2n}{\mu^{\prime}} 
      \tilde{S}(s^{\prime},\mu^{\prime})
      \f]
      where the function \f$ f_{\rho} \f$ is defined by
      \f[
      f_{\rho} \equiv \Theta(s^{\prime}-s)
      \left(\frac{1-s}{s}\right)^{2 n+1} \left[\frac{s^{\prime
      2n}}{(1-s^{\prime})^{2n+2}}\right] + \Theta(s^{\prime}-s)
      \left(\frac{1-s}{s}\right)^{2 n+1} \left[\frac{s^{\prime
      2n}}{(1-s^{\prime})^{2n+2}}\right]
      \f]
      This function is stored in \ref f_rho . Similar 
      definitions are made for \ref f_gamma and \ref f_omega .

      The Keplerial orbit at the equator is 
      \f[
      \Omega_K = \frac{\omega^{\prime}}{2 \psi^{\prime}} ...
      \f]
      (eq. 31 in \ref Stergioulas03 )
      \comment
      
  */
  class nstar_rot2 {
  
 public:    
  
  /// The number of grid points in the \f$ \mu \f$ direction
  static const int MDIV=65;
  /// The number of grid points in the \f$ s \f$ direction
  static const int SDIV=129;
  /// The number of Legendre polynomials
  static const int LMAX=10;

  protected:

  /** \brief Subclass of \ref nstar_rot which specifies the function
      to invert a polytropic EOS
  */
  class polytrope_solve {

  protected:

    /** \brief The polytropic index
     */
    double _Gamma_P;

    /** \brief The energy density
     */
    double _ee;

  public:

    /** \brief Create a function object with specified 
	polytropic index and ?
    */
    polytrope_solve(double Gamma_P, double ee) {
      _Gamma_P=Gamma_P;
      _ee=ee;
    }
      
    /** \brief The function
     */
    double operator()(double rho0) {
      return pow(rho0,_Gamma_P)/(_Gamma_P-1.0)+rho0-_ee;
    }
      
  };

  /// The polytrope solver
  o2scl::root_bkt_cern<polytrope_solve> rbc;

  /// Array search object
  o2scl::search_vec<double *> sv;

  /** \brief grid point in RK integration */ 
  static const int RDIV=900;                     
    
  /** \brief Maximum value of s-coordinate (default 0.9999) */  
  double SMAX;
  /** \brief Spacing in \f$ s \f$ direction, 
      \f$ \mathrm{SMAX}/(\mathrm{SDIV}-1) \f$ 
  */
  double DS;
  /** \brief Spacing in \f$ \mu \f$ direction, \f$ 1/(\mathrm{MDIV}-1) \f$ 
   */ 
  double DM;

  /// Minimum radius for spherical stars (default \f$ 10^{-15} \f$)
  double RMIN;

  /** \brief Indicates if iteration diverged (default 0) */ 
  int a_check;                       

  /// \name Grid quantities set in make_grid()
  //@{
  /// \f$ s \f$
  double s_gp[SDIV+1];                 
  /// \f$ s (1-s) \f$
  double s_1_s[SDIV+1];
  /// \f$ 1-s \f$
  double one_s[SDIV+1];
  /// \f$ \mu \f$
  double mu[MDIV+1];                   
  /// \f$ 1-\mu^2 \f$
  double one_m2[MDIV+1];
  /// \f$ \theta \f$ defined by \f$ \mathrm{acos}~\mu \f$
  double theta[MDIV+1];
  /// \f$ \sin \theta \f$
  double sin_theta[MDIV+1];
  //@}

  /// \name Grid values computed in integrate() for spherical_star()
  //@{
  /// Isotropic radius
  double r_gp[RDIV+1];
  /// Radial coordinate
  double r_is_gp[RDIV+1];
  /// Metric function \f$ \lambda \f$
  double lambda_gp[RDIV+1];
  /// Metric function \f$ \nu \f$
  double nu_gp[RDIV+1];
  /// Enclosed gravitational mass
  double m_gp[RDIV+1];
  /// Energy density
  double e_d_gp[RDIV+1];   
  //@}

  /// \name Metric functions
  //@{
  /** \brief potential \f$ \rho \f$ */ 
  double rho[SDIV+1][MDIV+1];          
  /** \brief potential \f$ \gamma \f$ */ 
  double gamma[SDIV+1][MDIV+1];         
  /** \brief potential \f$ \omega \f$ */ 
  double omega[SDIV+1][MDIV+1];        
  /** \brief potential \f$ \alpha \f$ */ 
  double alpha[SDIV+1][MDIV+1];        
  //@}

  /// \name Initial guess computed by \ref comp()
  //@{
  /// Guess for the equatorial radius
  double r_e_guess;
  /** \brief Guess for \f$ \rho \f$ */ 
  double rho_guess[SDIV+1][MDIV+1];    
  /** \brief Guess for \f$ \gamma \f$ */
  double gamma_guess[SDIV+1][MDIV+1];   
  /** \brief Guess for \f$ \alpha \f$ */
  double omega_guess[SDIV+1][MDIV+1];  
  /** \brief Guess for \f$ \omega \f$ */
  double alpha_guess[SDIV+1][MDIV+1];  
  //@}

  /// \name EOS quantities
  //@{
  /** \brief Energy density \f$ \epsilon \f$ */
  double energy[SDIV+1][MDIV+1];       
  /** \brief Pressure */ 
  double pressure[SDIV+1][MDIV+1];     
  /** \brief Enthalpy */
  double enthalpy[SDIV+1][MDIV+1];     
  //@}

  /// \name Other quantities defined over the full two-dimensional grid
  //@{
  /** \brief Proper velocity squared */
  double velocity_sq[SDIV+1][MDIV+1];  
  /** \brief Derivative of \f$ \alpha \f$ with respect to \f$ \mu \f$ */
  double da_dm[SDIV+1][MDIV+1];
  //@}

  /// \name Quantities defined for fixed values of mu
  //@{
  /** \brief \f$ \gamma(s) \f$ at \f$ \mu=1 \f$ */
  double gamma_mu_1[SDIV+1];            
  /** \brief \f$ \gamma(s) \f$ at \f$ \mu=0 \f$ */
  double gamma_mu_0[SDIV+1];            
  /** \brief \f$ \rho(s) \f$ at \f$ \mu=1 \f$ */
  double rho_mu_1[SDIV+1];             
  /** \brief \f$ \rho(s) \f$ at \f$ \mu=0 \f$ */
  double rho_mu_0[SDIV+1];             
  /** \brief \f$ \omega(s) \f$ at \f$ \mu=0 \f$ */
  double omega_mu_0[SDIV+1];           
  //@}

  /** \brief The value of \f$ \hat{\gamma} \f$ at the pole */  
  double gamma_pole_h;                  
  /** \brief The value of \f$ \hat{\gamma} \f$ at the center */
  double gamma_center_h;                
  /** \brief The value of \f$ \hat{\gamma} \f$ at the equator */
  double gamma_equator_h;               
  /** \brief The value of \f$ \hat{\rho} \f$ at the pole */ 
  double rho_pole_h;                   
  /** \brief The value of \f$ \hat{\rho} \f$ at the center */
  double rho_center_h;                 
  /** \brief The value of \f$ \hat{\rho} \f$ at the equator */ 
  double rho_equator_h;                
  /** \brief The value of \f$ \hat{\omega} \f$ at the equator */
  double omega_equator_h;              
  /** \brief Angular velocity, \f$ \hat{\omega} \f$ */
  double Omega_h;                      
  /** \brief Central pressure */ 
  double p_center;                     
  /** \brief Central enthalpy */
  double h_center;                     

  /// \name Desc
  //@{
  /** \brief \f$ f_{\rho}(s,n,s') \f$ */
  double f_rho[SDIV+1][LMAX+1][SDIV+1];
  /** \brief \f$ f_{\gamma}(s,n,s') \f$ */
  double f_gamma[SDIV+1][LMAX+1][SDIV+1];
  /** \brief \f$ f_{\omega}(s,n,s') \f$ */
  double f_omega[SDIV+1][LMAX+1][SDIV+1];
  //@}
  
  /// \name Legendre polynomials
  //@{
  /** \brief Legendre polynomial \f$ P_{2n}(\mu) \f$ 
   */  
  double P_2n[MDIV+1][LMAX+1];         
  /** \brief Associated Legendre polynomial \f$ P^1_{2n-1}(\mu) \f$ 
   */ 
  double P1_2n_1[MDIV+1][LMAX+1];      
  //@}

  /** \brief Relative accuracy for the equatorial radius,
      \f$ r_e \f$ (default \f$ 10^{-5} \f$) 

      Used in \ref iterate() .
  */
  double eq_radius_tol_rel;                    

  /** \brief Integrated term over m in eqn for \f$ \rho \f$ */
  double D1_rho[LMAX+1][SDIV+1];  
  /** \brief Integrated term over m in eqn for \f$ \gamma \f$ */
  double D1_gamma[LMAX+1][SDIV+1]; 
  /** \brief Integ. term over m in eqn for \f$ \omega \f$ */
  double D1_omega[LMAX+1][SDIV+1];
  /** \brief Integrated term over s in eqn for \f$ \rho \f$ */
  double D2_rho[SDIV+1][LMAX+1];  
  /** \brief Integrated term over s in eqn for \f$ \gamma \f$ */
  double D2_gamma[SDIV+1][LMAX+1]; 
  /** \brief Integ. term over s in eqn for \f$ \omega \f$ */
  double D2_omega[SDIV+1][LMAX+1];

  /** \brief source term in eqn for \f$ \gamma \f$ */
  double S_gamma[SDIV+1][MDIV+1];  
  /** \brief source term in eqn for \f$ \rho \f$ */
  double S_rho[SDIV+1][MDIV+1];   
  /** \brief source term in eqn for \f$ \omega \f$ */
  double S_omega[SDIV+1][MDIV+1]; 

  /** \brief The tolerance for the functions with the prefix "fix" 
      (default \f$ 10^{-4} \f$ )
  */
  double tol_abs;

  /// \name Thermodyanmic quantities near the surface
  //@{
  /// Pressure at the surface
  double p_surface;
  /// Energy density at the surface
  double e_surface;
  /** \brief Minimum specific enthalpy
   */
  double enthalpy_min;                 
  //@}

  /// \name Polytrope parameters
  //@{
  /// Polytropic index
  double n_P;
  /// Polytropic exponent
  double Gamma_P;
  //@}

  /// \name For CL_LOW is true
  //@{
  /// Desc
  double e_match;
  /// Desc
  double p_match;
  /// Desc
  double h_match;
  /// Desc
  double n0_match;
  /** \brief Desc (default false)
   */
  bool CL_LOW;
  /// Desc
  double de_pt;
  /// Desc
  double e_cl;
  //@}

  /// \name Interpolation functions
  //@{
  /** \brief Cache for interpolation
   */
  int n_nearest;
  
  /// Search in array \c x of length \c n for value \c val
  int new_search(int n, double *x, double val);
    
  /** \brief Driver for the interpolation routine. 
	
      First we find the tab. point nearest to xb, then we
      interpolate using four points around xb.
	
      Used by \ref int_z(), \ref e_at_p(), \ref p_at_e(), 
      \ref p_at_h(), \ref h_at_p(), \ref n0_at_e(), 
      \ref comp_omega(), \ref comp_M_J(), \ref comp(), 
      \ref spherical_star(), \ref iterate().
  */  
  double interp(double xp[], double yp[], int np ,double xb);

  /** \brief Driver for the interpolation routine.

      Four point interpolation at a 
      given offset the index of the first point k. 

      Used in \ref comp() .
  */
  double interp_4_k(double xp[], double yp[], int np, double xb, int k);
  //@}

  /** \brief Integrate f[mu] from m-1 to m. 

      This implements a 8-point closed Newton-Cotes formula.
	
      Used in \ref comp() .
  */
  double int_z(double f[MDIV+1], int m);

  /// \name EOS functions
  //@{
  /** \brief Compute \f$ \varepsilon(P) \f$  
	
      Used in \ref dm_dr_is(), \ref dp_dr_is(), \ref integrate()
      and \ref iterate(). 
  */
  double e_at_p(double pp);

  /** \brief Compute \f$ P(\varepsilon) \f$  
	
      Used in \ref make_center() and \ref integrate().
  */
  double p_at_e(double ee);

  /** \brief Pressure at fixed enthalpy

      Used in \ref iterate().
  */
  double p_at_h(double hh);

  /** \brief Enthalpy at fixed pressure 

      Used in \ref make_center() and \ref integrate().
  */
  double h_at_p(double pp);
    
  /** \brief Baryon density at fixed energy density 

      Used in \ref comp_M_J() and \ref comp() .
  */
  double n0_at_e(double ee);
  //@}

  /// \name Derivatives on the grid
  //@{
  /** \brief Returns the derivative w.r.t. s of an array f[SDIV+1]. 
   */ 
  double s_deriv(double f[SDIV+1], int s);

  /** \brief Returns the derivative w.r.t. mu of an array f[MDIV+1]. 
   */ 
  double m_deriv(double f[MDIV+1], int m);

  /** \brief Returns the derivative w.r.t. s  
   */ 
  double deriv_s(double f[SDIV+1][MDIV+1], int s, int m);

  /** \brief Returns the derivative w.r.t. mu 
   */ 
  double deriv_m(double f[SDIV+1][MDIV+1], int s, int m);

  /** \brief Returns the derivative w.r.t. s and mu 
   */ 
  double deriv_sm(double f[SDIV+1][MDIV+1], int s, int m);
  //@}

  /// \name Initialization functions
  //@{
  /** \brief Returns the Legendre polynomial of degree n, evaluated at x. 

      This uses the recurrence relation and is used in \ref comp_f_P()
      which is called by the constructor.
  */
  double legendre(int n, double x);

  /** \brief Compute two-point functions
	
      This function computes the 2-point functions \f$
      f^m_{2n}(r,r') \f$ used to integrate the potentials \f$ \rho,
      \gamma \f$ and \f$ \omega \f$ (See \ref Komatsu89 for
      details). Since the grid points are fixed, we can compute the
      functions \ref f_rho, \ref f_gamma, \ref f_omega, \ref P_2n,
      and \ref P1_2n_1 once at the beginning.

      See Eqs. 27-29 of \ref Cook92 and Eqs. 33-35 of \ref
      Komatsu89. This function is called by the constructor.
  */
  void comp_f_P();

    /** \brief Create computational mesh. 

      Create the computational mesh for \f$ s=r/(r+r_e) \f$
      (where \f$ r_e \f$ is the coordinate equatorial radius) 
      and \f$ \mu = \cos \theta \f$
      using 
      \f[
      s[i]=\mathrm{SMAX}\left(\frac{i-1}{\mathrm{SDIV}-1}\right)
      \f]
      \f[
      \mu[j]=\left(\frac{i-1}{\mathrm{MDIV}-1}\right)
      \f]
      When \f$ r=0 \f$, \f$ s=0 \f$, when \f$ r=r_e \f$, 
      \f$ s=1/2 \f$, and when \f$ r = \infty \f$, \f$ s=1 \f$ .
      \comment
      (Note that some versions of the manual have a typo,
      giving \f$ 1-i \f$ rather than \f$ i-1 \f$ above.)
      \endcomment
	
      Points in the mu-direction are stored in the array
      <tt>mu[i]</tt>. Points in the s-direction are stored in the
      array <tt>s_gp[j]</tt>.

      This function sets \ref s_gp, \ref s_1_s, \ref one_s,
      \ref mu, \ref one_m2, \ref theta and \ref sin_theta .
      All of these arrays are unit-indexed. It is called by
      the constructor.
  */
  void make_grid();
  //@}

  /** \brief Compute central pressure and enthalpy from central
      energy density

      For polytropic EOSs, this also computes <tt>rho0_center</tt> .
  */
  void make_center(double e_center);

  /// \name Post-processing functions
  //@{
  /** \brief Compute Omega and Omega_K. 
   */
  void comp_omega();
  
  /** \brief Compute rest mass and angular momentum. 
   */
  void comp_M_J();

  /** \brief Compute various quantities.

      The main post-processing funciton
  */
  void comp();
  //@}

  /// \name For computing spherical stars
  //@{
  /** \brief Computes a spherically symmetric star 
	
      The metric is 
      \f[
      ds^2 = -e^{2\nu}dt^2 + e^{2\lambda} dr^2 + r^2 d\theta^2 + 
      r^2 sin^2\theta d\phi^2
      \f]
      where \f$ r \f$ is an isotropic radial coordinate 
      (corresponding to <tt>r_is</tt> in the code).
      
      This function computes \ref r_e_guess, \ref R_e, 
      \ref Mass, and \ref Z_p .
  */
  void spherical_star();

  /** \brief Derivative of gravitational mass with respect to
      isotropic radius */
  double dm_dr_is(double r_is, double r, double m, double p);
 
  /** \brief Derivative of pressure with respect to isotropic radius */
  double dp_dr_is(double r_is, double r, double m, double p);

  /** \brief Derivative of radius with respect to isotropic radius */
  double dr_dr_is(double r_is, double r, double m);
  
  /** \brief Integrate one of the differential equations for 
      spherical stars*/
  void integrate(int i_check, double &r_final, double &m_final,
		 double &r_is_final);
  //@}

  /** \brief Main iteration function
   */
  int iterate(double r_ratio);

  /// \name EOS member variables
  //@{ 
  /** \brief If true, then an EOS has been set
   */
  bool eos_set;
  
  /** \brief If true, then use a polytrope and rescale
   */
  bool scaled_polytrope;

  /** \brief Pointer to the user-specified EOS
   */
  eos_nstar_rot *eosp;
  //@}
  
  public:

  nstar_rot2();

  /** \brief Verbosity parameter
   */
  int verbose;

  /// \name Output
  //@{
  /** \brief Central energy density (in \f$ \mathrm{g}/\mathrm{cm}^3 \f$) 
   */
  double e_center;                     
  /** \brief Ratio of polar to equatorial radius
   */ 
  double r_ratio;                      
  /** \brief Coordinate equatorial radius
   */ 
  double r_e;                          
  //@}

  /// \name Quantities computed by nstar_rot2::comp() (in order)
  //@{
  /** \brief Radius at pole */      
  double r_p;                          
  /** \brief The value of the s-coordinate at the pole */
  double s_p;                          
  /** \brief The value of the s-coordinate at the equator */
  double s_e; 
  /// The velocity at the equator
  double velocity_equator;              
  /** \brief Circumferential radius (i.e. the radius defined such
      that \f$ 2 \pi R_e \f$ is the proper circumference) */
  double R_e;                          
  /// Proper mass
  double Mass_p;
  /// Gravitational mass (in g)
  double Mass;
  /// Baryonic mass (in g)
  double Mass_0;
  /// Angular momentum
  double J;
  /// Angular velocity
  double Omega;
  /// Total rotational kinetic energy
  double T;
  /// Moment of inertia
  double I;
  /// Gravitational binding energy
  double W;
  /// Polar redshift
  double Z_p;
  /// Forward equatorial redshift
  double Z_f;
  /// Backward equatorial redshift
  double Z_b;
  /** \brief Kepler rotation frequency (in 1/s) */  
  double Omega_K;                      
  /// The eccentricity
  double eccentricity;
  /// Desc
  double v_plus[SDIV+1];
  /// Desc
  double v_minus[SDIV+1];
  /// Desc
  double vel_plus;
  /// Desc
  double vel_minus;
  /** \brief Height from surface of last stable co-rotating circular 
      orbit in equatorial plane

      If this is zero then all orbits are stable.
  */
  double h_plus;
  /** \brief Height from surface of last stable counter-rotating circular 
      orbit in equatorial plane
	
      If this is zero then all orbits are stable.
  */
  double h_minus;
  /// Desc
  double Omega_plus;
  /// Desc
  double u_phi;
  /// Angular velocity of a particle in a circular orbit at the equator
  double Omega_p;
  /// Desc
  double grv2;
  /// Desc
  double grv2_new;
  /// Desc
  double grv3;
  /** \brief Ratio of potential \f$ \omega \f$ to angular 
      velocity \f$ \Omega \f$
  */
  double om_over_Om;
  /** \brief Mass quadrupole moment
   */
  double mass_quadrupole;
  //@}

  /// \name Settings
  //@{
  /// The convergence factor (default 1.0)
  double cf;
  //@}

  /// \name Internal constants
  //@{
  /** \brief Use the values of the constants from the original RNS
      code
  */
  void constants_rns();
  /** \brief Use the \o2 values
   */
  void constants_o2scl();
  /** \brief Speed of light in vacuum (in CGS units) */ 
  double C;
  /** \brief Gravitational constant (in CGS units) */ 
  double G;
  /** \brief Mass of sun (in g) */
  double MSUN;
  /** \brief Square of length scale in CGS units, 
      \f$ \kappa \equiv 10^{-15} c^2/G \f$
  */
  double KAPPA;
  /** \brief The mass of one baryon (in g)
   */
  double MB;
  /** \brief The value \f$ \kappa G c^{-4} \f$ */
  double KSCALE;
  /// The constant \f$ \pi \f$
  double PI;
  //@}

  /// \name Basic Usage
  //@{
  /** \brief Set the EOS
   */
  void set_eos(eos_nstar_rot &eos) {
    eosp=&eos;
    eos_set=true;
    scaled_polytrope=false;
    return;
  }

  /** \brief Use a polytropic EOS with a specified index
   */
  void polytrope_eos(double index) {
    n_P=index;
    scaled_polytrope=true;
    eos_set=true;
    return;
  }
    
  /** \brief Construct a configuration with a fixed central 
      energy density and a fixed axis ratio
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$ .
  */
  int fix_cent_eden_axis_rat(double cent_eden, double axis_rat);
    
  /** \brief Construct a configuration with a fixed central 
      energy density and a fixed gravitational mass
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$ and the gravitational 
      mass should be in solar masses. 
  */
  int fix_cent_eden_grav_mass(double cent_eden, double grav_mass);

  /** \brief Construct a configuration with a fixed central 
      energy density and a fixed baryonic mass
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$ and the baryonic 
      mass should be in solar masses. 
  */
  int fix_cent_eden_bar_mass(double cent_eden, double bar_mass);

  /** \brief Construct a configuration with a fixed central 
      energy density and the Keplerian rotation rate
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$ .
  */
  int fix_cent_eden_with_kepler(double cent_eden);
    
  /** \brief Construct a non-rotating configuration with a fixed central 
      energy density
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$ .
  */
  int fix_cent_eden_non_rot(double cent_eden);

  /** \brief Construct a configuration with a fixed central 
      energy density and a fixed angular velocity.
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$.
  */
  int fix_cent_eden_ang_vel(double cent_eden, double ang_vel);

  /** \brief Construct a configuration with a fixed central 
      energy density and a fixed angular momentum.
	
      The central energy density should be in \f$
      \mathrm{g}/\mathrm{cm}^3 \f$.
  */
  int fix_cent_eden_ang_mom(double cent_eden, double ang_mom);
  //@}
    
  /** \name Testing functions

      All these compare with hard-coded results obtained with
      the RNS code. 
  */
  //@{
  /** \brief Test determining configuration with fixed central
      energy density and fixed radius ratio with EOS C
  */    
  void test1(o2scl::test_mgr &t);
    
  /** \brief Test configuration rotating and Keplerian frequency
      with a fixed central energy density and EOS C
  */    
  void test2(o2scl::test_mgr &t);
    
  /** \brief Test fixed central energy density and fixed 
      gravitational mass with EOS C
  */    
  void test3(o2scl::test_mgr &t);
    
  /** \brief Test fixed central energy density and fixed baryonic 
      mass with EOS C
  */    
  void test4(o2scl::test_mgr &t);
    
  /** \brief Test fixed central energy density and fixed angular
      velocity with EOS C
  */    
  void test5(o2scl::test_mgr &t);
    
  /** \brief Test fixed central energy density and fixed angular 
      momentum with EOS C
  */    
  void test6(o2scl::test_mgr &t);

  /** \brief Test a series of non-rotating stars on a energy density
      grid with EOS C
  */    
  void test7(o2scl::test_mgr &t);
    
  /** \brief Test Keplerian frequency for a polytrope
   */    
  void test8(o2scl::test_mgr &t);
  //@}


  };

}

#endif
