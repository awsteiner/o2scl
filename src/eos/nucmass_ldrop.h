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
#include <o2scl/eos_had_rmf.h>
#include <o2scl/inte_qagiu_gsl.h>

namespace o2scl {

  /** \brief Simple liquid drop mass formula
      
      Includes bulk part plus surface and Coulomb (no pairing)
      without neutron skin and without any isospin contribution
      to the surface energy. 

      The NL4 EOS is loaded by default.
      
      \note This class sets part::inc_rest_mass to true 
      for the particle objects specified in set_n_and_p().

      \hline

      <b>Central densities</b>

      Given a saturation density, \f$ n_0 \f$ and a transition
      density, \f$ n_t \f$, we set \f$ I = 1 - 2 Z/A \f$, and then
      assume \f$ \delta = I \f$. We further assume that the
      isospin-asymmetric saturation density \f$ n_L \f$ is
      \f[
      n_L = n_0 + n_1 \delta^2
      \f]
      and then we can compute \f$ n_{p} = (1 - \delta)/2 n_L \f$ and
      \f$ n_{n} = (1 + \delta)/2 n_L \f$ . 

      Note that \f$ \delta = I \f$ implies no neutron skin. A neutron
      skin occurs when \f$ \delta < I \f$, and \f$ \delta = 0 \f$
      implies a "maximum skin size" which is occurs when no extra
      neutrons are in center and all extra neutrons are located in the
      skin, i.e. \f$ N_{\mathrm{skin}} = N-Z \f$.

      <b>Nuclear radii</b>

      The neutron and proton radii are determined from the 
      central densities with
      \f{eqnarray*}
      R_n &=& \left( \frac{3 N}{4 \pi n_n} \right)^{1/3} \nonumber \\
      R_n &=& \left( \frac{3 Z}{4 \pi n_p} \right)^{1/3} 
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
      \varepsilon_{\mathrm{surf}} = \frac{3 \sigma}{R}
      \f]
      where \f$ \sigma \f$ is the surface tension.
      To compute the surface energy per baryon, we divide by
      the baryon density, \f$ n_n + n_p \f$. We can rewrite this
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
      where the surface tension \f$ \sigma \f$ (in MeV) is given in 
      \ref surften. 

      Taking a typical value, \f$ \sigma =1.1~\mathrm{MeV} \f$ and 
      \f$ n_L = 0.16~\mathrm{fm}^{-3} \f$, gives the standard result,
      \f$ E_{\mathrm{surf}}/A = 18~\mathrm{MeV}~A^{-1/3} \f$. 

      \verbatim embed:rst
      See also [Ravenhall83]_.
      \endverbatim

      <b>Coulomb energy contribution</b>

      The Coulomb energy density is
      \f[
      \varepsilon_{\mathrm{Coul}} = \frac{4 \pi}{5} n_p^2 e^2 R_p^2
      \f]
      The energy per baryon is
      \f[
      E_{\mathrm{Coul}}/A = \frac{4 \pi}{5 n_L} n_p^2 e^2 R_p^2 
      \f]
      This is the expression used in the code, except for a prefactor
      \ref coul_coeff which is a fit parameter and should be close to
      unity.

      Assuming \f$ R_p = R \f$
      and \f$ Z = \frac{4 \pi}{3} R^3 n_p \f$
      and \f$ R = \left[ 3 A / (4 \pi n_L) \right]^{1/3} \f$
      gives
      \f[
      E_{\mathrm{Coul}}/A = \frac{6^{2/3}}{5} 
       \pi^{1/3} e^2 n_L^{1/3} \frac{Z^2}{A^{4/3}} 
      \f]
      and taking \f$ n_L = 0.16~\mathrm{fm}^{-3} \f$ and 
      \f$ e^2 = \hbar c/137 \f$ gives the standard result
      \f[
      E_{\mathrm{Coul}}/A = 0.76~\mathrm{MeV}~Z^2 A^{-4/3}
      \f]

      \hline

      <b>References</b>

      \verbatim embed:rst
      This class is based on  [Steiner08]_, which was originally
      based on [Lattimer85]_ and [Lattimer91]_.
      \endverbatim

      \hline
  */
  class nucmass_ldrop : public nucmass_fit_base {

  public:
    
    nucmass_ldrop();

    /// \name Input parameters
    //@{ 
    /// Density asymmetry (default 0)
    double n1;
    
    /** \brief Saturation density ( The default is \f$ 0.16
        \mathrm{fm}^{-3} \f$)
    */
    double n0;
    
    /// Surface tension in MeV (default 1.1 MeV)
    double surften;
    
    /// Coulomb coefficient (default 1.0)
    double coul_coeff;
    //@}    

    /// \name Output quantities
    //@{ 
    /// Internal average neutron density
    double nn;
    
    /// Internal average proton density
    double np;

    /// Neutron radius
    double Rn;

    /// Proton radius
    double Rp;

    /// Bulk part of energy
    double bulk;

    /// Surface part of energy
    double surf;

    /// Coulomb part of energy
    double coul;
    //@}

    /** \brief If true, then return large mass excesses when
        unphysical parameters are selected (default false)
    */
    bool large_vals_unphys;
    
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

    /** \brief Given \c Z and \c N, return the binding energy in MeV

        This function is currently independent of \c npout, \c nnout,
        \c dim, and \c chi.
    */
    virtual double drip_binding_energy_d(double Z, double N,
                                         double npout, double nnout, 
                                         double chi, double dim, double T);

    /// \name EOS and particle parameters
    //@{
    /// Change the base hadronic EOS
    int set_eos_had_temp_base(eos_had_temp_base &uhe) {
      heos=&uhe;
      return 0;
    }

    /// The default hadronic EOS
    eos_had_rmf def_had_eos;

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
    //@}

    /// \name Fitting functions
    //@{
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);
    //@}
    
    /// Return the type, \c "nucmass_ldrop".
    virtual const char *type() { return "nucmass_ldrop"; }
      
    /// Energy and pressure
    thermo th;

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Pointer to neutron 
    fermion *n;
    /// Pointer to proton
    fermion *p;
    /// The base EOS for bulk matter
    eos_had_temp_base *heos;

#endif
    
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

      <b>Bulk energy</b>

      The central densities and radii, \f$ n_n, n_p, R_n, R_p \f$
      are all determined in the same way as \ref nucmass_ldrop, 
      except that now \f$ \delta \equiv I \zeta \f$, where
      \f$ \zeta \f$ is stored in \ref doi . Note that this
      means \f$ N > Z~\mathrm{iff}~R_n>R_p \f$. 

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

      If \ref full_surface is false, then the surface energy is 
      just that from \ref nucmass_ldrop , with an extra factor
      for the surface symmetry energy
      \f[
      E_{\mathrm{surf}} = \frac{\sigma d \chi}{3 n_L}
      \left(\frac{36 \pi n_L}{A} \right)^{1/3} 
      \left( 1- \sigma_{\delta} \delta^2 \right)
      \f]
      where \f$ \sigma_{\delta} \f$ is unitless and stored in \ref ss.

      If \ref full_surface is true, then the following
      temperature- and isospin-dependent surface energy is used.
      Taking
      \f$ x \equiv n_p /(n_n+n_p) \f$, 
      \f[
      \sigma(x,T) = \left( \frac{\sigma d \chi}{3} \right)
      \frac{16+b}{x^{-3}+b+(1-x)^{-3}}
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
      b=-16+96 \frac{\sigma}{\sigma_{\delta}}
      \f]
      which is chosen to ensure that the surface energy
      is identical to the expression used when \ref full_surface
      is false for small values of \f$ \delta \f$.

      \verbatim embed:rst
      This surface energy comes from [Steiner08]_, which was
      originally based on the expression in [Lattimer85]_.
      \endverbatim
      
      <b>Coulomb energy</b>

      First, we define the volume fraction occupied by
      protons, \f$ \chi_p \f$ which is
      \f[
      \chi_p = \chi \left(\frac{R_p}{R_n}\right)^3
      \f]
      when \f$ R_n>R_p \f$ and \f$ \chi_p = \chi \f$ otherwise.
      The Coulomb energy density is
      \f[
      \varepsilon = 2 \pi \chi_p
      e^2 R_p^2 (n_p-n_{p,\mathrm{out}})^2 f_d(\chi_p)
      \f]
      where the function \f$ f_d(\chi_p) \f$ is 
      \f[
      f_d(\chi_p) = \frac{1}{(d+2)}
      \left[ \frac{2}{(d-2)} \left( 1 - \frac{d}{2} 
      \chi_p^{(1-2/d)} \right) + \chi_p \right]
      \f]
      When \f$ d=3 \f$, \f$ f_3(\chi_p) \f$ reduces to
      \f[
      \frac{1}{5} \left[ 2 - 3 \chi_p^{1/3} + \chi_p \right]
      \f]
      and the limit \f$ \chi_p \rightarrow 0 \f$ gives the expression
      used in \ref nucmass_ldrop.

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
         - (future) Does this work if the nucleus is "inside-out"?
         
      \endverbatim
      
      \hline

      Excluded volume and \ref rel_vacuum:

      Typically in a single-nucleus EOS with a neutron drip 
      (ignoring translational degrees of freedom for the nucleus) 
      \f[
      f = n_N m_N + (1-\chi_n) f_{n,\mathrm{drip}}
      \f]
      where
      \f[
      m_N = \frac{A}{n_n+n_p}(f - n_n m_n - n_p m_p)
      \f]
      Since \f$ n_N = 3/(4 \pi R_{\mathrm{ws}}^3) \f$, and 
      \f$ \chi_n = (R_n/R_{\mathrm{ws}})^3 \f$, this is 
      \f[
      f = \frac{3}{4 \pi R_{\mathrm{ws}}^3} 
      \left[ m_N - f_{n,\mathrm{drip}} \frac{4 \pi}{3} R_n^3 \right]
      + f_{n,\mathrm{drip}}
      \f]

      \endcomment

      \hline

      <b>References</b>

      \verbatim embed:rst
      Designed in [Steiner08]_ and [Souza09co]_ based in part
      on [Lattimer85]_ and [Lattimer91]_.
      \endverbatim

      \hline
  */
  class nucmass_ldrop_skin : public nucmass_ldrop {
    
  public:

    nucmass_ldrop_skin();

    /// Return the type, \c "nucmass_ldrop_skin".
    virtual const char *type() { return "nucmass_ldrop_skin"; }

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

    /** \brief If true, properly fix the surface for the pure neutron
        matter limit (default true)
    */
    bool full_surface;

    /** \brief If true, separately compute the skin for the bulk energy
        (default false)
    */
    bool new_skin_mode;

    /// Ratio of \f$ \delta/I \f$ (default 0.8).
    double doi;

    /// Surface symmetry energy (default 0.5)
    double ss;

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
    //@}
    
    /** \brief If true, define the nuclear mass relative to the vacuum
        (default true)
    */
    bool rel_vacuum;

    /** \brief The critical temperature of isospin-symmetric matter in 
        \f$ fm^{-1} \f$ (default \f$ 20.085/(\hbar c)\f$.)
    */
    double Tchalf;
    
    /** \brief Return the free binding energy of a nucleus in a many-body 
        environment
    */
    virtual double drip_binding_energy_d(double Z, double N,
                                         double npout, double nnout,
                                         double chi, double dim, double T);
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
       double chi, double dim, double T);

  };

}

#endif
