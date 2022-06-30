/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/** \file eos_tov.h
    \brief File defining \ref o2scl::eos_tov
*/
#ifndef O2SCL_TOV_EOS_H
#define O2SCL_TOV_EOS_H

#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/lib_settings.h>
#include <o2scl/interp.h>
#include <o2scl/table_units.h>
#include <o2scl/vector_derint.h>
#include <o2scl/root_brent_gsl.h>

namespace o2scl {

  /** \brief A EOS base class for the TOV solver
   */
  class eos_tov {

  protected:
    
    /** \brief Set to true if the baryon density is provided in the
        EOS (default false)
    */
    bool baryon_column;

  public:
    
    eos_tov();

    virtual ~eos_tov() {}
    
    /// Control for output (default 1)
    int verbose;

    /// Return true if a baryon density is available
    bool has_baryons() {
      return baryon_column;
    }
    
    /** \brief Check that the baryon density is consistent 
        with the \f$ P(\varepsilon) \f$
    */
    void check_nb(double &avg_abs_dev, double &max_abs_dev);

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr)=0;

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed)=0;

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed)=0;

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr)=0;

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb)=0;

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb)=0;

    /** \brief Given the pressure, produce the energy and number densities
        
        The arguments \c pr and \c ed should always be in \f$
        M_{\odot}/\mathrm{km}^3 \f$ . The argument for \c nb should be
        in \f$ \mathrm{fm}^{-3} \f$ .
        
        If \ref baryon_column is false, then \c nb is unmodified.
    */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb)=0;

  };

  /** \brief The Buchdahl EOS for the TOV solver

      The Buchdahl EOS is 
      \f[
      \varepsilon = 12 \sqrt{P_{*} P}- 5 P
      \f]
      which can be inverted to give
      \f[
      P = - \frac{\varepsilon}{5} + 
      \frac{72 P^{*}}{25} \left[1+\sqrt{1-\frac{5 \varepsilon}{36 P^{*}}}
      \right]
      \f]

      Physical solutions are obtained only for \f$ P< 25 P_{*}/144 \f$
      (ensuring that the argument to the square root is positive) and
      \f$ \beta=G M/R<1/6 \f$ (ensuring that the EOS is not acausal).

      The (baryon) chemical potential is
      \f[
      \mu = \mu_1 
      \left(\sqrt{P_1}-3\sqrt{P^{*}}\right)^{1/2}
      \left(\sqrt{P}-3\sqrt{P^{*}}\right)^{-1/2}
      \f]
      where \f$ \mu_1 \f$ is the chemical potential at \f$ P=P_1 \f$.
      The baryon density, as a function of \f$ P \f$ can be obtained
      from \f$ n_B=(\varepsilon+P)/\mu \f$. In \o2, the parameters \f$
      \mu_1 \f$ and \f$ P_1 \f$ are expressed in terms of the baryon
      density (\f$ n_{B,1} \f$) and energy density (\f$ \varepsilon_1
      \f$) instead, and can be set in \ref set_baryon_density().

      Without loss of generality, one can choose \f$ P_1 = 0 \f$,
      leaving \f$ \mu_1 \f$, the value of the chemical potential when
      \f$ P=0 \f$, as the remaining parameter. In analogy to the case
      of nucleonic matter, one can choose \f$ mu_1 = m \f$ where \f$ m
      \f$ is the nucleon mass. Then, the baryon density can be
      simplified to
      \f[
      n m_n = 12 \sqrt{P p_{*}} \left( 1-\frac{1}{3} \sqrt{P/p_{*}} 
      \right)^{3/2}
      \f]
      c.f. Eq. 10 in Lattimer et al. (2001).

      The mass-radius curve is the solution of the equation
      \f[
      M = \left[ \frac{\pi}{288 p_{*} 
      \left(1-2 \beta\right)}\right]^{1/2} 
      \beta\left(1-\beta\right)
      \f]
      where \f$ \beta = GM/R \f$ and \f$ A \f$ (which
      has units of inverse km) is defined by
      \f[
      A^2 = \frac{288 \pi G P^{*}}{1-2 \beta}
      \f]
      For a fixed gravitational mass, this equation is solved (given
      \ref Pstar) to obtain the radius by the function \ref
      rad_from_gm() .

      The central pressure and energy density are
      \f[
      P_c = 36 p_{*} \beta^2 
      \f]
      \f[
      {\varepsilon}_c = 72 p_{*} \beta (1-5 \beta/2) 
      \f]

      To obtain energy density and pressure profiles can be obtained
      by solving 
      \f[
      r=r^{\prime} \left(\frac{1-\beta+u}{1-2 \beta}\right)
      \f]
      for the new coordinate \f$ r^{\prime} \f$
      where \f$ u \f$ is defined by
      \f[
      u = \beta \frac{\sin(A r^{\prime})}{A r^{\prime}}
      \f]
      Using these, the profiles are
      \f[
      P(r) = A^2 (1- 2 \beta) u^2 
      \left[ 8 \pi \left(1 - \beta+u\right)^2\right]^{-1}
      \f] 
      and
      \f[
      \varepsilon(r) = 2 A^2 (1- 2 \beta) u 
      \left( 1 - \beta + 3 u/2\right)
      \left[ 8 \pi \left(1 - \beta+u\right)^2\right]^{-1}
      \f]
      
      \verbatim embed:rst
      Based on [Lattimer01]_.
      \endverbatim

      \note The default tov_solve beginning and ending pressures
      work well for physical EOSs, but not for the Buchdahl EOS,
      so if you use this EOS in a tov_solve object and try to 
      compute the M-R curve you will need to change prbegin
      and prend.

  */
  class eos_tov_buchdahl : public eos_tov {
    
  protected:

    /** \brief The baryon density at \c ed1
     */
    double nb1;

    /** \brief The energy density for which the baryon density is known
     */
    double ed1;

    /** \brief The pressure at \c ed1
     */
    double pr1;

    /** \brief Solver
     */
    root_brent_gsl<> rbg;
    
  public:

    eos_tov_buchdahl();

    virtual ~eos_tov_buchdahl() {}

    /** \brief The parameter with units of pressure in units of solar
        masses per km cubed (default value \f$ 3.2 \times 10^{-5} \f$
        )
    */
    double Pstar;

    /** \brief The gravitational constant in kilometers per solar mass

        \note Make this static const
        
        \comment
        This is useful because the user needs to specify beta=G*M/R
        for some of the functions below.
        \endcomment
     */
    double G_km_Msun;
    
    /** \brief Set the baryon density
     */
    void set_baryon_density(double nb, double ed);

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);

    /** \brief Given the pressure, produce the energy and number densities
        
        If the baryon density is not specified, it should be set to
        zero or \ref baryon_column should be set to false
    */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb);

    /** \brief Given the gravitational mass, compute the radius

        This function solves Eq. 10.67 in Schutz' GR book 
        (see also Eq. 7 in Lattimer01)
        \f[
        M = \left[ \frac{\pi}{288 p_{*} 
        \left(1-2 \beta\right)}\right]^{1/2} 
        \beta\left(1-\beta\right)
        \f]

        This function can be used to determine the full mass-radius 
        given \c Pstar.
     */
    virtual double rad_from_gm(double gm);
    
    /** \brief Compute the energy density at radius
        \c r given the compactness (unitless)
    */
    virtual double ed_from_r_gm(double r, double beta);
    
    /** \brief Compute the pressure at radius
        \c r given the compactness (unitless)
    */
    virtual double pr_from_r_gm(double r, double beta);
    
    /** \brief Compute \f$ \exp(2 \Phi) \f$ at radius
        \c r given the compactness (unitless)
    */
    virtual double exp2phi_from_r_gm(double r, double beta);
    
    /** \brief Compute \f$ \exp(2 \Lambda) \f$ at radius
        \c r given the compactness (unitless)
    */
    virtual double exp2lam_from_r_gm(double r, double beta);
    
  protected:

    /** \brief Solve for the radius at fixed gravitational mass
     */
    virtual double solve_rad(double rad, double gm);

    /** \brief Solve for \f$ r^{\prime} \f$ as a function of
        \f$ r \f$ at fixed gravitational mass
    */
    virtual double solve_rp(double rp, double r, double beta);
    
  };

  /** \brief Standard polytropic EOS \f$ P = K \varepsilon^{1+1/n} \f$

      The quantity \f$ K \f$ must be in units of 
      \f$ \left(\mathrm{M}_{\odot}/\mathrm{km}^3\right)^{-1/n} \f$ .

      \comment
      The documentation below was taken from bamr.
      \endcomment
      For a polytrope \f$ P = K \varepsilon^{1+1/n} \f$
      beginning at a pressure of \f$ P_1 \f$, an energy
      density of \f$ \varepsilon_1 \f$ and a baryon density 
      of \f$ n_{B,1} \f$, the baryon density along the polytrope
      is 
      \f[
      n_B = n_{B,1} \left(\frac{\varepsilon}{\varepsilon_1}\right)^{1+n} 
      \left(\frac{\varepsilon_1+P_1}{\varepsilon+P}\right)^{n} \, .
      \f]
      Similarly, the chemical potential is
      \f[
      \mu_B = \mu_{B,1} \left(1 + \frac{P_1}{\varepsilon_1}\right)^{1+n}
      \left(1 + \frac{P}{\varepsilon}\right)^{-(1+n)} \, .
      \f]
      The expression for the 
      baryon density can be inverted to determine \f$ \varepsilon(n_B) \f$
      \f[
      \varepsilon(n_B) = \left[ \left(\frac{n_{B,1}}
      {n_B \varepsilon_1} \right)^{1/n}
      \left(1+\frac{P_1}{\varepsilon_1}\right)-K\right]^{-n} \, .
      \f]
      Sometimes the baryon susceptibility is also useful 
      \f[
      \frac{d \mu_B}{d n_B} = \left(1+1/n\right)
      \left( \frac{P}{\varepsilon}\right)
      \left( \frac{\mu_B}{n_B}\right) \, .
      \f]

      \future The simple formulation fo the expressions here more than
      likely break down when their arguments are sufficiently extreme.
  */
  class eos_tov_polytrope : public eos_tov {
    
  protected:

    /** \brief The baryon density at \c ed1
     */
    double nb1;

    /** \brief The energy density for which the baryon density is known
     */
    double ed1;

    /** \brief The pressure at \c ed1
     */
    double pr1;

    /** \brief Coefficient (default 1.0)
    */
    double K;

    /// Index (default 3.0)
    double n;

  public:

    eos_tov_polytrope();

    virtual ~eos_tov_polytrope() {}

    /** \brief Set the coefficient and polytropic index
     */
    void set_coeff_index(double coeff, double index);
    
    /** \brief Set the baryon density
     */
    void set_baryon_density(double nb, double ed);

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);

    /** \brief Given the pressure, produce the energy and number densities
     */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb);

  };

  /** \brief Linear EOS \f$ P = c_s^2 (\varepsilon-\varepsilon_0) \f$

      This implements a linear EOS with a fixed speed of sound and a
      fixed energy density at zero pressure. This will also compute
      the baryon density, if one calls \ref set_baryon_density() to
      set the baryon density at one fiducial energy density.

      Given a fiducial baryon density \f$ n_{B,1} \f$ at some energy
      density \f$ \varepsilon_1 \f$ and pressure \f$ P_1 =
      c_s^2(\varepsilon_1-\varepsilon) \f$, the baryon density is
      \f[
      n_B = n_{B,1} \left[ \frac{\varepsilon(1+c_s^2) - 
      c_s^2 \varepsilon_0 } {\varepsilon_1 (1 + c_s^2) - 
      c_s^2 \varepsilon_0}\right]^{1/(1+c_s^2)} = 
      n_{B,1} \left[ \frac{ \varepsilon + P }
      {\varepsilon_1 + P_1} \right]^{1/(1+c_s^2)}
      \f]

      \note AWS, 6/29/22: As can be seen from this expression, if \f$
      \epsilon_0 =0 \f$, then I don't think one can take \f$
      \varepsilon_1=P_1=0 \f$ because the baryon density does not have
      a finite value at that point, independent of the value of \f$
      n_{B,1} \f$.

      One can write the energy density as a function of the
      baryon density
      \f[
      \varepsilon = \varepsilon_0 \left(\frac{c_s^2}{1+c_s^2}\right) + 
      \left(\frac{\varepsilon_1+P_1}{1+c_s^2}\right)
      \left(\frac{n_B}{n_{B,1}}\right)^{1+c_s^2} 
      \f]
      which has the form of a constant term plus a power-law
      in the density with exponent \f$ 1+c_s^2 \f$. The pressure
      is 
      \f[
      P = c_s^2 \left[\varepsilon_1 - \left(\frac{c_s^2 \varepsilon_0}
      {1+c_s^2}\right)\right] \left(\frac{n_B}{n_{B,1}}\right)^{1+c_s^2} 
      - \varepsilon_0 \left(\frac{c_s^2}{1+c_s^2}\right) 
      \f]

      \note Experimental
   */
  class eos_tov_linear : public eos_tov {

  protected:

    /** \brief The baryon density at \c ed1
     */
    double nb1;

    /** \brief The energy density for which the baryon density is known
     */
    double ed1;

    /** \brief The pressure for which the baryon density is known
     */
    double pr1;

    /** \brief Coefficient (default 1.0)
    */
    double cs2;

    /// The energy density at zero pressure (default 0.0)
    double eps0;

  public:

    eos_tov_linear();

    virtual ~eos_tov_linear() {}

    /** \brief Set the sound speed and energy density at zero pressure
     */
    void set_cs2_eps0(double cs2_, double eps0_);

    /** \brief Set the baryon density
     */
    void set_baryon_density(double nb, double ed);

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);

    /** \brief Given the pressure, produce the energy and number densities
     */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb);

  };

  /** \brief Provide an EOS for TOV solvers based on 
      interpolation of user-supplied vectors
   */
  template<class vec_t> class eos_tov_vectors : public eos_tov {

    /** \brief Internal function to reset the interpolation
     */
    void reset_interp(size_t n) {
      pe_int.set(n,pr_vec,ed_vec,itp_linear);
      ep_int.set(n,ed_vec,pr_vec,itp_linear);
      return;
    }
    
    /** \brief Internal function to reset the interpolation
        with baryon density
     */
    void reset_interp_nb(size_t n) {
      reset_interp(n);
      pn_int.set(n,pr_vec,nb_vec,itp_linear);
      np_int.set(n,nb_vec,pr_vec,itp_linear);
      en_int.set(n,ed_vec,nb_vec,itp_linear);
      ne_int.set(n,nb_vec,ed_vec,itp_linear);
      return;
    }
    
  public:

    /** \brief Read the EOS from a set of equal length
        vectors for energy density, pressure, and baryon density

        In this version, the user-specified vectors are swapped
        with internal storage.
    */
    void read_vectors_swap(size_t user_n, vec_t &user_ed, vec_t &user_pr,
                           vec_t &user_nb) {
      std::swap(user_ed,ed_vec);
      std::swap(user_pr,pr_vec);
      std::swap(user_nb,nb_vec);
      this->baryon_column=true;
      reset_interp_nb(user_n);
      return;
    }
    
    /** \brief Read the EOS from a pair of equal length
        vectors for energy density and pressure

        In this version, the user-specified vectors are swapped
        with internal storage.
    */
    void read_vectors_swap(size_t user_n, vec_t &user_ed, vec_t &user_pr) {
      std::swap(user_ed,ed_vec);
      std::swap(user_pr,pr_vec);
      this->baryon_column=false;
      reset_interp(user_n);
      return;
    }

    /** \brief Read the EOS from a set of equal length
        vectors for energy density, pressure, and baryon density

        In this version, the user-specified vectors are copied
        to internal storage.
    */
    void read_vectors_copy(size_t user_n, vec_t &user_ed, vec_t &user_pr,
                           vec_t &user_nb) {
      if (ed_vec.size()!=user_n) ed_vec.resize(user_n);
      if (pr_vec.size()!=user_n) pr_vec.resize(user_n);
      if (nb_vec.size()!=user_n) nb_vec.resize(user_n);
      vector_copy(user_ed,ed_vec);
      vector_copy(user_pr,pr_vec);
      vector_copy(user_nb,nb_vec);
      this->baryon_column=true;
      reset_interp_nb(user_n);
      return;
    }
    
    /** \brief Read the EOS from a pair of equal length
        vectors for energy density and pressure

        In this version, the user-specified vectors are copied
        to internal storage.
    */
    void read_vectors_copy(size_t user_n, vec_t &user_ed, vec_t &user_pr) {
      if (ed_vec.size()!=user_n) ed_vec.resize(user_n);
      if (pr_vec.size()!=user_n) pr_vec.resize(user_n);
      vector_copy(user_ed,ed_vec);
      vector_copy(user_pr,pr_vec);
      this->baryon_column=false;
      reset_interp(user_n);
      return;
    }

    /// \name Basic EOS functions
    //@{
    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr) {
      return pe_int.eval(pr);
    }

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed) {
      return ep_int.eval(ed);
    }
    
    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed) {
      return en_int.eval(ed);
    }
    
    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr) {
      return pn_int.eval(pr);
    }
    
    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb) {
      return ne_int.eval(nb);
    }
    
    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb) {
      return np_int.eval(nb);
    }

    /** \brief Given the pressure, produce the energy and number densities

        The arguments \c pr and \c ed should always be in \f$
        M_{\odot}/\mathrm{km}^3 \f$ . The argument for \c nb should be
        in \f$ \mathrm{fm}^{-3} \f$ .
        
        If \ref baryon_column is false, then \c nb is unmodified.
    */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb) {      
      ed=ed_from_pr(pr);
      if (this->baryon_column) {
        nb=nb_from_pr(pr);
      }
      return;
    }
    //@}
    
  protected:
    
    /// \name EOS storage
    //@{
    /// Energy densities from full EOS
    vec_t ed_vec;
    /// Pressures from full EOS
    vec_t pr_vec;
    /// Baryon densities from full EOS
    vec_t nb_vec;
    //@}

    /// \name Interpolators
    //@{
    interp_vec<vec_t> pe_int;
    interp_vec<vec_t> pn_int;
    interp_vec<vec_t> ep_int;
    interp_vec<vec_t> en_int;
    interp_vec<vec_t> np_int;
    interp_vec<vec_t> ne_int;
    //@}

  };
  
  /** \brief An EOS for the TOV solver using simple linear
      interpolation and an optional low-density EOS

      The simplest usage of this class is simply to use \ref
      read_table() to read a tabulated EOS stored in a \ref
      table_units object and optionally specify a separate crust 
      (i.e. low-density) EOS.

      There are two methods to handle the crust-core interface. The
      method labeled <tt>smooth_trans</tt> uses the crust below
      pressure \f$ P_{\mathrm{lo}} \f$ (equal to the value of \ref
      trans_pres divided by \ref trans_width) and the core above
      pressure \f$ P_{\mathrm{hi}} \f$ (the value of \ref trans_pres
      times \ref trans_width) and then in between uses
      \f[
      \varepsilon(P) = [1-\chi(P)] \varepsilon_{\mathrm{crust}}(P) + 
      \chi(P) \varepsilon_{\mathrm{core}}(P)
      \f]
      where the value \f$ \chi(P) \f$ is determined by
      \f[
      \chi(P) = (P-P_{\mathrm{lo}})/
      (P_{\mathrm{hi}}-P_{\mathrm{lo}}) \, .
      \f]
      This method is a bit more faithful to the original EOS tables,
      but the matching can result in pressures which decrease with
      increasing energy density. Alternatively the <tt>match_line</tt>
      method uses \f$
      \varepsilon_{\mathrm{lo}}=\varepsilon_{\mathrm{crust}}
      (P_{\mathrm{lo}}) \f$ and \f$
      \varepsilon_{\mathrm{hi}}=\varepsilon_{\mathrm{core}}
      (P_{\mathrm{hi}}) \f$ and
      \f[
      \varepsilon(P) = (\varepsilon_{\mathrm{hi}} - 
      \varepsilon_{\mathrm{lo}}) \chi 
      + \varepsilon_{\mathrm{lo}} \, .
      \f]
      (using the same expression for \f$ \chi \f$ ). This method less
      frequently results in decreasing pressures, but can deviate
      further from the original tables. For the baryon density,
      expression similar to those used for \f$ \varepsilon(P) \f$
      are used for \f$ n_B(P) \f$ .

      By default, no crust EOS is used. If a crust EOS is specified
      through one of the crust EOS functions, then by default \ref
      trans_width is 1.0 and <tt>transition_mode</tt> is set equal
      <tt>smooth_trans</tt>. This creates a discontinuous energy
      density between the core and crust EOS at the transition
      pressure. A smoother transition can be chosen by increasing \ref
      trans_width to a value larger than 1. The crust EOS can be
      changed after the core EOS is specified.

      The value of \ref trans_pres is set either by \ref
      set_transition(), or any of the crust EOS functions.

      Internally, energy and pressure are stored in units of \f$
      \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ and baryon density is
      stored in units of \f$ \mathrm{fm}^{-3} \f$ . The user-specified
      EOS table is left as is, and unit conversion is performed as
      needed in ed_nb_from_pr() and other functions from the units
      specified in the input \ref table_units object.

      A faster method is provided in \ref o2scl::eos_tov_vectors where
      the user can simply specify objects of type
      <tt>std::vector<double></tt> which store the energy density,
      pressure, and baryon density (which should include the crust if
      necessary).
  */
  class eos_tov_interp : public eos_tov {
    
  public:
    
    eos_tov_interp();

    virtual ~eos_tov_interp();

    /// \name Mode of transitioning between crust and core EOS
    //@{
    int transition_mode;
    static const int smooth_trans=0;
    static const int match_line=1;
    //@}

    /** \brief If true, call the error handler if the EOS
        reports a non-finite value (default true)
    */
    bool err_nonconv;
    
    /// \name Basic EOS functions (all in the internal unit system)
    //@{
    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);
    //@}

    /// \name Basic usage
    //@{
    /** \brief Specify the EOS through a table

        If units are specified for any of the columns, then this
        function attempts to automatically determine the correct
        conversion factors using the \ref o2scl::convert_units object
        returned by \ref o2scl::o2scl_settings . If the units for any
        of the columns are blank, then they are assumed to be the
        native units for \ref o2scl::tov_solve .
        
        \note The input table must have at least 2 rows and 
        the pressure column must be strictly increasing.

        This function copies the needed information from the
        table so if it is modified then this function
        needs to be called again to read a new table.
    */
    void read_table(table_units<> &eosat, std::string s_cole, 
                    std::string s_colp, std::string s_colnb="");
    //@}
    
    /// \name Crust EOS functions
    //@{
    /** \brief Negele-Vautherin-Baym-Pethick-Sutherland crust

        \verbatim embed:rst
        This crust is taken from [Negele73ns]_ and [Baym71tg]_.
        \endverbatim
    */
    void default_low_dens_eos();

    /** \brief The EOS labeled "C" in Stergioulas' rotating 
        neutron star code

        \verbatim embed:rst
        This EOS is adapted from [Bethe74]_.
        \endverbatim
    */
    void rns_C_low_dens_eos();
    
    /** \brief Crust EOS from Shen et al. (2011)

        \verbatim embed:rst
        From [Shen11b]_.
        \endverbatim
     */
    void sho11_low_dens_eos();

    /** \brief Crust EOS from Steiner (2012)

        \verbatim embed:rst
        This function uses the neutron star crust models from \ref
        [Steiner12]_. 
        \endverbatim

        The current acceptable values for \c model are
        <tt>APR</tt>, <tt>Gs</tt>, <tt>Rs</tt> and <tt>SLy4</tt>.
    */
    void s12_low_dens_eos(std::string model="SLy4",
                              bool external=false);

    /** \brief Crust EOS from Goriely, Chamel, and Pearson
        
        \verbatim embed:rst
        From [Goriely10]_, [Pearson11]_, and [Pearson12]_.
        \endverbatim
     */
    void gcp10_low_dens_eos(std::string model="BSk20",
                          bool external=false);

    /** \brief Crust EOS from Newton et al. 2013 given L in MeV

        \verbatim embed:rst
        From [Newton13]_
        \endverbatim

        Current acceptable values for \c model are <tt>PNM</tt>
        and <tt>J35</tt>. 
    */
    void ngl13_low_dens_eos(double L, std::string model="PNM",
                             bool external=false);
    
    /** \brief Crust EOS from Newton et al. (2013) given S and L in MeV
        and a transition density

        Note that this function works only if \f$ 28 < S < 38 \f$ MeV,
        \f$ 25 < L < 115 \f$ MeV, \f$ 0.01 < n_t < 0.15 \f$, 
        and \f$ L > 5 S-65~\mathrm{MeV} \f$
        . If \c fname is a string of length 0 (the default),
        then this function looks for a file named \c newton_SL.o2
        in the \o2 data directory specified by
        \ref o2scl::lib_settings_class::get_data_dir() .

        \verbatim embed:rst
        From [Newton13]_
        \endverbatim
    */
    void ngl13_low_dens_eos2(double S, double L, double nt,
                             std::string fname="");
    
    /// Compute with no crust EOS (this is the default)
    void no_low_dens_eos() {
      use_crust=false;
      return;
    }
    //@}

    /// \name Functions used by the tov_solve class
    //@{
    /** \brief Given the pressure, produce the energy and number densities

        The arguments \c pr and \c ed should always be in \f$
        M_{\odot}/\mathrm{km}^3 \f$ . The argument for \c nb should be
        in \f$ \mathrm{fm}^{-3} \f$ .
        
        If the baryon density is not specified, it should be set to
        zero or \ref baryon_column should be set to false
    */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb);
    //@}

    /// \name Other functions
    //@{
    /** \brief Get the energy density from the pressure in the 
        user-specified unit system
    */
    virtual void get_eden_user(double pres, double &ed, double &nb);

    /** \brief Get the transition pressure (in the user-specified
        unit system) and width
    */
    void get_transition(double &ptrans, double &pwidth);
    
    /** \brief Set the transition pressure and "width"

        Sets the transition pressure and the width (specified as a
        number greater than unity in \c pw) of the transition between
        the two EOSs. The transition should be in the same units of
        the user-specified EOS. The transition is done smoothly using
        linear interpolation between \f$ P=\mathrm{ptrans}/pmathrm{pw}
        \f$ and \f$ P=\mathrm{ptrans} \times pmathrm{pw} \f$.
     */
    void set_transition(double ptrans, double pw);
    //@}

    /// \name Full EOS arrays (in internal units)
    //@{
    /// Energy density
    std::vector<double> full_vece;
    /// Pressure
    std::vector<double> full_vecp;
    /// Baryon density
    std::vector<double> full_vecnb;
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Internal function to reinterpolate if if either the
        core or crust tables are changed
     */
    void internal_read();

    /// \name Crust EOS
    //@{
    /// Set to true if we are using a crust EOS (default false)
    bool use_crust;

    /// Energy densities
    std::vector<double> crust_vece;
    /// Pressures
    std::vector<double> crust_vecp;
    /// Baryon densities
    std::vector<double> crust_vecnb;
    //@}
    
    /// \name Core EOS
    //@{
    /// Energy densities
    std::vector<double> core_vece;
    /// Pressures
    std::vector<double> core_vecp;
    /// Baryon densities
    std::vector<double> core_vecnb;
    //@}
    
    /// \name Interpolation objects
    //@{
    interp_vec<std::vector<double> > pe_int;
    interp_vec<std::vector<double> > pnb_int;
    interp<std::vector<double> > gen_int;
    //@}

    /// \name Unit conversion factors for core EOS
    //@{
    /// Unit conversion factor for energy density (default 1.0)
    double efactor;
    /// Unit conversion factor for pressure (default 1.0)
    double pfactor;
    /// Unit conversion factor for baryon density (default 1.0)
    double nfactor;
    //@}

    /// \name Properties of transition
    //@{
    /** \brief Transition pressure (in \f$ M_{\odot}/\mathrm{km}^3 \f$)
     */
    double trans_pres;
    /// Transition width (unitless)
    double trans_width;
    //@}

#endif

  };

}

#endif


