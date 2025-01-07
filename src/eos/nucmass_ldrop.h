/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
/** \file nucmass_ldrop.h
    \brief File defining \ref o2scl::nucmass_ldrop
*/
#ifndef LDROP_MASS_H
#define LDROP_MASS_H

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/constants.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/inte_qagiu_gsl.h>

namespace o2scl {

  /** \brief Simple liquid drop mass formula
      
      This class includes bulk, surface and Coulomb contributions to
      the energy. There is no pairing, no neutron skin, and no isospin
      contribution to the surface energy.

      \verbatim embed:rst
      The NRAPR equation of state from [Steiner05b]_ is used
      for the bulk energy by default.
      \endverbatim
      
      \note This class sets part::inc_rest_mass to true for the
      particle objects specified in set_n_and_p().

      \note The input parameter T (for the temperature) should be
      given in units of inverse femtometers. This is a bit confusing,
      since the binding energy is returned in MeV.
      
      <b>Definition of </b> \f$ \chi \f$ <b> and </b> \f$ n_L \f$

      The variable \f$ \chi \f$ is defined as the fractional volume
      occupied by and \f$ n_L \f$ is the density of nucleons inside
      the nucleus. If \f$ V \f$ is the total volume of the nucleus
      plus the surrounding Wigner-Seitz cell, then we have
      \f[
      A = V n_L \chi
      \f]
      
      <b>Central densities</b>

      Given a saturation density, \f$ n_0 \f$ and a transition
      density, \f$ n_t \f$, we set \f$ I = 1 - 2 Z/A \f$, and then
      assume \f$ \delta = I \f$. We further assume that the
      isospin-asymmetric saturation density \f$ n_L \f$ is
      \f[
      n_L = n_0 + n_1 \delta^2
      \f]
      and then we can compute \f$ n_{p} = n_L (1 - \delta)/2 \f$ and
      \f$ n_{n} = n_L (1 + \delta)/2 \f$ . 

      Note that \f$ \delta = I \f$ implies no neutron skin. A neutron
      skin occurs when \f$ \delta < I \f$, and \f$ \delta = 0 \f$
      implies a "maximum skin size" which is occurs when no extra
      neutrons are in center and all extra neutrons are located in the
      skin, i.e. \f$ N_{\mathrm{skin}} = N-Z \f$. See also \ref
      nucmass_ldrop_skin .

      <b>Nuclear radius</b>

      The nuclear radius is determined by
      \f{eqnarray*}
      R &=& \left( \frac{3 A}{4 \pi n_L} \right)^{1/3}
      \f}

      <b>Bulk energy contribution</b>
      
      The bulk binding energy contribution ( \f$ \sim -16 \f$
      MeV per nucleon) and the symmetry energy are computing using the
      hadronic EOS (either \ref def_had_eos or the EOS specified in
      the most recent call to set_eos_had_temp_base() ). The bulk
      energy per baryon is
      \f[
      E_{\mathrm{bulk}}/A = \frac{1}{n_{L} }
      \left[\varepsilon(n_n,n_p) - n_n m_n - n_p m_p \right]
      \f]

      <b>Surface energy contribution</b>

      The surface energy density is
      \f[
      \varepsilon_{\mathrm{surf}} = \frac{3 \sigma \chi}{R}
      \f]
      where \f$ \sigma \f$ is the surface tension.
      To compute the surface energy per baryon, using
      \f$ A = V n_L \chi \f$, we find
      \f[
      E_{\mathrm{surf}}/A = \frac{3 \sigma}{n_n + n_p} 
      \left[ \frac{3 A}{ 4 (n_n+n_p) \pi}
      \right]^{-1/3}
      \f]
      or
      \f[
      E_{\mathrm{surf}}/A = \frac{\sigma}{n_L}
      \left(\frac{36 \pi n_L}{A} \right)^{1/3}
      \f]
      where the surface tension \f$ \sigma \f$ (in \f$
      \mathrm{MeV}/\mathrm{fm}^2 \f$) is given in \ref surften.

      Taking a typical value, \f$ \sigma
      =1.1~\mathrm{MeV}/\mathrm{fm}^2 \f$ and \f$ n_L =
      0.16~\mathrm{fm}^{-3} \f$, gives the standard result, \f$
      E_{\mathrm{surf}}/A = 18~\mathrm{MeV}~A^{-1/3} \f$ (see \ref
      nucmass_semi_empirical::Ss in \ref nucmass_semi_empirical).

      \verbatim embed:rst
      See also [Ravenhall83]_.
      \endverbatim

      <b>Coulomb energy contribution</b>

      The Coulomb energy density is
      \f[
      \varepsilon_{\mathrm{Coul}} = {\cal{C}}
      \frac{4 \pi}{5} \chi n_p^2 e^2 R^2
      \f]
      where the extra coefficient \f$ {\cal{C}} \f$ is stored
      in \ref coul_coeff and defaults to a value of 1.
      The energy per baryon is
      \f[
      E_{\mathrm{Coul}}/A = {\cal{C}} \frac{4 \pi}{5 n_L} n_p^2 e^2 R^2 
      \f]
      This is the expression used in the code.

      Taking \f$ {\cal{C}} =1 \f$ and using \f$ Z = \frac{4 \pi}{3}
      R^3 n_p \f$ and \f$ R = \left[ 3 A / (4 \pi n_L) \right]^{1/3}
      \f$ gives
      \f[
      E_{\mathrm{Coul}}/A = \frac{6^{2/3}}{5} 
       \pi^{1/3} e^2 n_L^{1/3} \frac{Z^2}{A^{4/3}} 
      \f]
      and taking \f$ n_L = 0.16~\mathrm{fm}^{-3} \f$ and 
      \f$ e^2 = \hbar c/137 \f$ gives the standard result
      \f[
      E_{\mathrm{Coul}}/A = 0.76~\mathrm{MeV}~Z^2 A^{-4/3}
      \f]

      \note This class ignores the values specified in the variables
      \c npout, \c nnout, \c dim, and \c chi.
      
      <b>References</b>

      \verbatim embed:rst
      This class is based on  [Steiner08]_, which was originally
      based on [Lattimer85]_ and [Lattimer91]_.
      \endverbatim
  */
  class nucmass_ldrop : public nucmass_fit_base {

  public:
    
    /// \name Input parameters
    //@{ 
    /// Density asymmetry (default 0)
    double n1;
    
    /** \brief Saturation density ( The default is \f$
        0.16~\mathrm{fm}^{-3} \f$)
    */
    double n0;
    
    /** \brief Surface tension in \f$ \mathrm{MeV}/\mathrm{fm}^2 \f$
        (default 1.1 \f$ \mathrm{MeV}/\mathrm{fm}^2 \f$)
    */
    double surften;
    
    /// Coulomb coefficient (default 1.0)
    double coul_coeff;
    //@}    

    /// \name Output quantities
    //@{ 
    /// Internal average neutron density in \f$ \mathrm{fm}^{-3} \f$
    double nn;
    
    /// Internal average proton density in \f$ \mathrm{fm}^{-3} \f$
    double np;

    /// Neutron radius (in \f$ \mathrm{fm} \f$ )
    double Rn;

    /// Proton radius (in \f$ \mathrm{fm} \f$ )
    double Rp;

    /// Bulk energy per baryon in \f$ \mathrm{MeV} \f$
    double bulk;

    /// Surface energy per baryon in \f$ \mathrm{MeV} \f$
    double surf;

    /// Coulomb energy per baryon in \f$ \mathrm{MeV} \f$
    double coul;
    //@}

    /// \name Other settings
    //@{
    /** \brief If true, then return large mass excesses when
        unphysical parameters are selected (default false)
    */
    bool large_vals_unphys;
    //@}

    /// \name Main functions
    //@{
    /// Default constructor
    nucmass_ldrop();

    /** \brief Given \c Z and \c N, return the mass excess in MeV
 
        \comment
        We don't need to implement mass_excess() for integers because
        that's done in the parent nucmass_cont. 
        \endcomment
    */
    virtual double mass_excess_d(double Z, double N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      return mass_excess_d(Z,N);
    }

    /** \brief Given \c Z and \c N, return the binding energy of
        the nucleus in MeV

        In this class, this function is currently independent of \c
        npout, \c nnout, \c dim, and \c chi.
    */
    virtual double drip_binding_energy_d(double Z, double N,
                                         double npout, double nnout, 
                                         double ne, double dim,
                                         double T);

    /// Return the type, \c "nucmass_ldrop".
    virtual const char *type() { return "nucmass_ldrop"; }
    //@}

    /// \name EOS and particle parameters
    //@{
    /// Change the base hadronic EOS
    int set_eos_had_temp_base(eos_had_temp_base &uhe) {
      heos=&uhe;
      return 0;
    }

    /// The default hadronic EOS
    eos_had_skyrme def_had_eos;

    /// Change neutron and proton objects
    void set_n_and_p(fermion &un, fermion &up) {
      n=&un;
      p=&up;
      return;
    }

    /// Default neutron
    fermion def_neutron;

    /// Default proton
    fermion def_proton;

    /// Energy and pressure
    thermo th;
    //@}

    /// \name Fitting functions
    //@{
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);
    //@}
    
  protected:

    /// \name Base objects [protected]
    //@{
    /// Pointer to neutron 
    fermion *n;
    /// Pointer to proton
    fermion *p;
    /// The base EOS for bulk matter
    eos_had_temp_base *heos;
    //@}

  };

  /** \brief More advanced liquid drop model

      In addition to the physics in \ref nucmass_ldrop, this includes
      corrections for
      - finite temperature
      - neutron skin
      - an isospin-dependent surface energy
      - decrease in the Coulomb energy from external protons

      \note The input parameter T should be given in units of inverse
      Fermis. This is a bit confusing, since the binding energy is
      returned in MeV.
      
      <b>Nuclear radii and volume fractions</b>
      
      The nuclear, neutron and proton radii are determined by
      \f{eqnarray*}
      R &=& \left( \frac{3 A}{4 \pi n_L} \right)^{1/3} \nonumber \\
      R_n &=& \left( \frac{3 N}{4 \pi n_n} \right)^{1/3} \nonumber \\
      R_p &=& \left( \frac{3 Z}{4 \pi n_p} \right)^{1/3} 
      \f}
      where the densities \f$ n_L, n_n \f$ and \f$ n_p \f$ are
      determined in the same way as in \ref nucmass_ldrop,
      except that now \f$ \delta \equiv I \zeta \f$, where
      \f$ \zeta \f$ is stored in \ref doi .

      The volume fraction occupied by
      protons, \f$ \chi_p \f$ is
      \f[
      \chi_p = \left(\frac{R_p}{R_{\mathrm{WS}}}\right)^3
      \f]
      and similarly for neutrons. We also define \f$ \chi \f$ as
      \f[
      \chi = \left(\frac{R}{R_{\mathrm{WS}}}\right)^3 \, .
      \f]
      We need to use charge neutrality
      \f[
      \frac{4 \pi}{3} R_p^3 \left(n_p - n_{p,\mathrm{out}}\right)
      + \frac{4 \pi}{3} R_{\mathrm{WS}}^3 n_{p,\mathrm{out}} = 
      \frac{4 \pi}{3} R_{\mathrm{WS}}^3 n_{e,\mathrm{out}} 
      \f]
      thus
      \f[
      \chi_p \left(n_p - n_{p,\mathrm{out}}\right) + 
      n_{p,\mathrm{out}} = n_{e,\mathrm{out}} \, .
      \f]
     
      <b>Bulk energy</b>

      If \ref new_skin_mode is false, then the bulk energy is 
      also computed as in \ref nucmass_ldrop. Otherwise, the
      number of nucleons in the core is computed with
      \f{eqnarray*}
      A_{\mathrm{core}} = Z (n_n+n_p)/n_p~\mathrm{for}~N\geq Z \\
      A_{\mathrm{core}} = N (n_n+n_p)/n_p~\mathrm{for}~Z>N \\
      \f}
      and \f$ A_{\mathrm{skin}} = A - A_{\mathrm{core}} \f$.
      The core contribution to the bulk energy is 
      \f[
      E_{\mathrm{core}}/A = \left(\frac{A_{\mathrm{core}}}{A}\right)
      \frac{1}{n_{L} }
      \left[\varepsilon(n_n,n_p) - n_n m_n - n_p m_p \right]
      \f]
      then the skin contribution is 
      \f[
      E_{\mathrm{skin}}/A = \left(\frac{A_{\mathrm{skin}}}{A}\right)
      \frac{1}{n_{L} }
      \left[\varepsilon(n_n,0) - n_n m_n \right]
      \quad\mathrm{for}\quad N>Z
      \f]
      and
      \f[
      E_{\mathrm{skin}}/A = \left(\frac{A_{\mathrm{skin}}}{A}\right)
      \frac{1}{n_{L} }
      \left[\varepsilon(0,n_p) - n_p m_p \right]
      \quad\mathrm{for}\quad Z>N
      \f]

      <b>Surface energy</b>

      If \ref full_surface is false, then the surface energy per
      baryon is just that from \ref nucmass_ldrop, with extra
      factors for dimension and the surface symmetry energy
      \f[
      E_{\mathrm{surf}}/A = \frac{\sigma d}{3 n_L}
      \left(\frac{36 \pi n_L}{A} \right)^{1/3} 
      \left( 1- \sigma_{\delta} \delta^2 \right)
      \f]
      where \f$ \sigma_{\delta} \f$ is unitless and stored in \ref ss.

      If \ref full_surface is true, then the following
      temperature- and isospin-dependent surface energy is used.
      Taking
      \f$ x \equiv n_p /(n_n+n_p) \f$, the new surface tension is
      \f[
      \sigma(x,T) = 
      \left[ \frac{16+b}{x^{-3}+b+(1-x)^{-3}} \right]
      \left[\frac{1-T^2/T_c(x)^2}{1+a(x) T^2/T_c(x)^2}\right]^{p}
      \f]
      where
      \f[
      a(x) = a_0 + a_2 y^2 + a_4 y^4 \, ,
      \f]
      \f[
      T_c(x) = T_c(x=1/2) \left( 1-c y^2 - d y^4\right)^{1/2}
      \f]
      and \f$ y=1/2-x \f$\. The value \f$ p \f$ is stored in \ref pp and
      the value \f$ T_c(x=1/2) \f$ is stored in \ref Tchalf.
      Currently, the values of \f$ c=3.313 \f$ and \f$ d=7.362 \f$ are
      fixed and cannot be changed. The value of \f$ b \f$ is
      determined from
      \f[
      b=-16+\frac{96}{\sigma_{\delta}}
      \f]
      which is chosen to ensure that the surface energy
      is identical to the expression used when \ref full_surface
      is false for small values of \f$ \delta \f$.

      \verbatim embed:rst
      This surface energy comes from [Steiner08]_, which was
      originally based on the expression in [Lattimer85]_.
      \endverbatim
      
      <b>Coulomb energy</b>

      The Coulomb energy density is
      \f[
      \varepsilon_{\mathrm{Coul}} = 2 \pi \chi_p
      e^2 R_p^2 (n_p-n_{p,\mathrm{out}})^2 f_d(\chi_p)
      \f]
      where the function \f$ f_d(\chi_p) \f$ is 
      \f[
      f_d(\chi_p) = \frac{1}{(d+2)}
      \left[ \frac{2}{(d-2)} \left( 1 - \frac{d}{2} 
      \chi_p^{(1-2/d)} \right) + \chi_p \right] \, .
      \f]

      AWS, 12/23/24: The fit to nuclear data seems to be better
      if we use 
      \f[
      \varepsilon_{\mathrm{Coul}} = 2 \pi \chi
      e^2 R_p^2 (n_p-n_{p,\mathrm{out}})^2 f_d(\chi_p)
      \f]
      instead. Thus
      \f[
      E_{\mathrm{Coul}}/A = \left(\frac{2 \pi e^2 R_p^2}{nL}\right)
      (n_p-n_{p,\mathrm{out}})^2 f_d(\chi_p)
      \f]
      When \f$ d=3 \f$, \f$ f_3(\chi_p) \f$ reduces to
      \f[
      \frac{1}{5} \left[ 2 - 3 \chi_p^{1/3} + \chi_p \right] \, .
      \f]
      Then, using the approximation \f$
      \chi_p = \chi \f$ and the limit \f$ \chi_p \rightarrow 0 \f$
      gives the expression used in \ref nucmass_ldrop. The second term
      in square brackets above gives the Wigner-Seitz approximation to
      the so-called "lattice" contribution
      \f[
      \varepsilon_{\mathrm{L}} =
      -\frac{6 \pi e^2}{5} \chi_p^{4/3}
      e^2 R_p^2 (n_p-n_{p,\mathrm{out}})^2 
      \f]
      Taking \f$ n_{p,\mathrm{out}}=0 \f$ and using
      \f[
      R_p \Rightarrow \left( \frac{3 Z}{4 \pi n_p}\right)^{1/3}
      \f]
      we get
      \f[
      \varepsilon_{\mathrm{L}} =
      -\left(\frac{6 \cdot 3^{2/3} \pi^{1/3}}
      {5 \cdot 4^{2/3}}\right) e^2 Z^{2/3}
      \left(\chi_p n_p\right)^{4/3}
      \f]
      and noting that charge equality implies \f$ \chi_p n_p = n_e \f$,
      gives
      \f[
      \varepsilon_{\mathrm{L}} =
      -1.4508 Z^{2/3} e^2 n_e^{4/3}
      \f]
      which is approximately equal to the expression used in \ref eos_crust .

      \verbatim embed:rst
      See [Ravenhall83]_ and [Steiner12]_.
      \endverbatim
      
      <b>Todos and Future</b>

      \verbatim embed:rst
      .. todo:: 

         In class nucmass_ldrop_skin: 

         - (future) Add translational energy?
         - (future) Remove excluded volume correction and compute nuclear
           mass relative to the gas rather than relative to the vacuum.
         - (future) In principle, Tc should be self-consistently determined
           from the EOS.
         - (future) Does this class work if the nucleus is "inside-out"?
         
      \endverbatim
      
      <b>References</b>

      \verbatim embed:rst
      Designed in [Steiner08]_ and [Souza09co]_ based in part
      on [Lattimer85]_ and [Lattimer91]_.
      \endverbatim
  */
  class nucmass_ldrop_skin : public nucmass_ldrop {
    
  public:

    /// \name Basic usage
    //@{
    /// Default constructor
    nucmass_ldrop_skin();

    /** \brief Return the free binding energy of a nucleus in a many-body 
        environment
    */
    virtual double drip_binding_energy_d(double Z, double N,
                                         double npout, double nnout,
                                         double ne, double dim, double T);
    
    /// Return the type, \c "nucmass_ldrop_skin".
    virtual const char *type() { return "nucmass_ldrop_skin"; }
    //@}

    /// \name Fitting functions
    //@{
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);
    //@}

    /// \name Settings
    //@{
    /** \brief If true, define the nuclear mass relative to the vacuum
        (default true)
    */
    bool rel_vacuum;

    /** \brief If true, properly fix the surface for the pure neutron
        matter limit (default true)
    */
    bool full_surface;

    /** \brief If true, separately compute the skin for the bulk energy
        (default false)
    */
    bool new_skin_mode;
    //@}

    /// \name Surface parameters
    //@{
    /// Ratio of \f$ \delta/I \f$ (default 0.8).
    double doi;

    /// Surface symmetry energy coefficient (default 0.5)
    double ss;
    //@}

    /// \name Input parameters for temperature dependence
    //@{ 
    /// Exponent (default 1.25)
    double pp;

    /// Coefficient (default 0.935)
    double a0;

    /// Coefficient (default -5.1)
    double a2;

    /// Coefficient (default -1.1)
    double a4;
    
    /** \brief The critical temperature of isospin-symmetric matter in
        \f$ \mathrm{fm}^{-1} \f$ (default \f$ 20.085/(\hbar c)\f$.)
    */
    double Tchalf;
    //@}
    
  };

  /** \brief Liquid drop model with pairing
      
      This class adds a pairing correction
      \f[
      E_{\mathrm{pair}}/A = - \frac{\zeta}{2 A^{3/2}}
      \left[
      \cos(Z \pi) + \cos(N \pi)
      \right]
      \f]
      where \f$ \zeta \f$ is stored in \ref Epair. The trigonometric
      functions give the expected result for integer values of 
      N and Z.
  */
  class nucmass_ldrop_pair : public nucmass_ldrop_skin {
    
  public:

    /// Return the type, \c "nucmass_ldrop_pair".
    virtual const char *type() { return "nucmass_ldrop_pair"; }

    nucmass_ldrop_pair() {
      nfit=7;
      Epair=13.0;
    }

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

    /// Pairing energy coefficient (default 13 MeV)
    double Epair;

    /// Most recently computed pairing energy per baryon
    double pair;

    /** \brief Return the free binding energy of a nucleus in a many-body 
        environment
    */
    virtual double drip_binding_energy_d
      (double Z, double N, double npout, double nnout,
       double ne, double dim, double T);

  };

}

#endif
