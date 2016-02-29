/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file eos_had_potential.h
    \brief File defining \ref o2scl::eos_had_potential
*/
#ifndef O2SCL_GEN_POTENTIAL_EOS_H
#define O2SCL_GEN_POTENTIAL_EOS_H

#include <iostream>
#include <string>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/part.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/fermion_nonrel.h>
#include <cstdlib>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Generalized potential model equation of state
      
      The single particle energy is defined by the functional derivative
      of the energy density with respect to the distribution function
      \f[
      e_{\tau} = \frac{\delta \varepsilon}{\delta f_{\tau}}
      \f]

      The effective mass is defined by
      \f[
      \frac{m^{*}}{m} = \left( \frac{m}{k}
      \frac{d e_{\tau}}{d k}
      \right)^{-1}_{k=k_F}
      \f]

      In all of the models, the kinetic energy density is 
      \f$\tau_n+\tau_p\f$ where 
      \f[
      \tau_i = \frac{2}{(2 \pi)^3} \int d^3 k~ 
      \left(\frac{k^2}{2 m}\right)f_i(k,T)
      \f]
      and the number density is
      \f[
      \rho_i = \frac{2}{(2 \pi)^3} \int d^3 k~f_i(k,T)
      \f]

      When \ref form is equal to \ref mdi_form or 
      \ref gbd_form, the potential energy
      density is given by \ref Das03 :
      \f[
      V(\rho,\delta) = \frac{Au}{\rho_0} \rho_n \rho_p +
      \frac{A_l}{2 \rho_0} \left(\rho_n^2+\rho_p^2\right)+
      \frac{B}{\sigma+1} \frac{\rho^{\sigma+1}}{\rho_0^{\sigma}}
      \left(1-x \delta^2\right)+V_{mom}(\rho,\delta)
      \f]
      where \f$\delta=1-2 \rho_p/(\rho_n+\rho_p)\f$.
      If \ref form is equal to \ref mdi_form, then
      \f[
      V_{mom}(\rho,\delta)=
      \frac{1}{\rho_0} 
      \sum_{\tau,\tau^{\prime}} C_{\tau,\tau^{\prime}} 
      \int \int
      d^3 k d^3 k^{\prime} 
      \frac{f_{\tau}(\vec{k}) f_{{\tau}^{\prime}} (\vec{k}^{\prime})}
      {1-(\vec{k}-\vec{k}^{\prime})^2/\Lambda^2}
      \f]
      where \f$C_{1/2,1/2}=C_{-1/2,-1/2}=C_{\ell}\f$ and
      \f$C_{1/2,-1/2}=C_{-1/2,1/2}=C_{u}\f$. Later 
      parameterizations in this form are given in \ref Chen05.

      Otherwise if \ref form is equal to \ref gbd_form, then
      \f[
      V_{mom}(\rho,\delta)=
      \frac{1}{\rho_0}\left[
      C_{\ell} \left( \rho_n g_n + \rho_p g_p \right)+
      C_u \left( \rho_n g_p + \rho_p g_n \right)
      \right]
      \f]
      where 
      \f[
      g_i=\frac{\Lambda^2}{\pi^2}\left[ k_{F,i}-\Lambda
      \mathrm{tan}^{-1} \left(k_{F,i}/\Lambda\right) \right]
      \f]

      Otherwise, if \ref form is equal to \ref bgbd_form, \ref bpal_form
      or \ref sl_form, then the potential energy density is
      given by \ref Bombaci01 :
      \f[
      V(\rho,\delta) = V_A+V_B+V_C
      \f]
      \f[
      V_A = \frac{2 A}{3 \rho_0}
      \left[ \left(1+\frac{x_0}{2}\right)\rho^2-
      \left(\frac{1}{2}+x_0\right)\left(\rho_n^2+\rho_p^2\right)\right]
      \f]
      \f[
      V_B=\frac{4 B}{3 \rho_0^{\sigma}} \frac{T}{1+4 B^{\prime} T / 
      \left(3 \rho_0^{\sigma-1} \rho^2\right)}
      \f]
      where 
      \f[
      T = \rho^{\sigma-1} \left[ \left( 1+\frac{x_3}{2} \right) \rho^2 - 
      \left(\frac{1}{2}+x_3\right)\left(\rho_n^2+\rho_p^2\right)\right]
      \f]
      The term \f$V_C\f$ is:
      \f[
      V_C=\sum_{i=1}^{i_{\mathrm{max}}} 
      \frac{4}{5} \left(C_{i}+2 z_i\right) \rho 
      (g_{n,i}+g_{p,i})+\frac{2}{5}\left(C_i -8 z_i\right)
      (\rho_n g_{n,i} + \rho_p g_{p,i})
      \f]
      where 
      \f[
      g_{\tau,i} = \frac{2}{(2 \pi)^3} \int d^3 k f_{\tau}(k,T)
      g_i(k)
      \f]

      For \ref form is equal to \ref bgbd_form or \ref form 
      is equal to \ref bpal_form, the form factor is given by
      \f[
      g_i(k) = \left(1+\frac{k^2}{\Lambda_i^2}\right)^{-1}
      \f]
      while for \ref form is equal to \ref sl_form, the form factor
      is given by
      \f[
      g_i(k) = 1-\frac{k^2}{\Lambda_i^2}
      \f]
      where \f$\Lambda_1\f$ is specified in the parameter 
      \c Lambda when necessary.

      \bug The BGBD and SL EOSs do not work.
      
      \future Calculate the chemical potentials analytically.
  */
  class eos_had_potential : public eos_had_eden_base {

  public:
    
    /// \name The parameters for the various interactions
    //@{
    double x,Au,Al,rho0,B,sigma,Cl,Cu,Lambda;
    double A,x0,x3,Bp,C1,z1,Lambda2,C2,z2,bpal_esym;
    int sym_index;
    //@}
    
    eos_had_potential();
    
    /// Equation of state as a function of density.
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lt);

    /// Form of potential
    int form;

    /** \brief The "momentum-dependent-interaction" form from 
	\ref Das03
     */
    static const int mdi_form=1;

    /// The modifed GBD form
    static const int bgbd_form=2;

    /// The form from \ref Prakash88 as formulated in \ref Bombaci01
    static const int bpal_form=3;

    /// The "SL" form. See \ref Bombaci01.
    static const int sl_form=4;

    /// The Gale, Bertsch, Das Gupta from \ref Gale87.
    static const int gbd_form=5;
    
    /// The form from \ref Prakash88.
    static const int pal_form=6;

    /** \brief Set the derivative object to calculate the 
	chemical potentials
    */
    int set_mu_deriv(deriv_base<> &de) {
      mu_deriv_set=true;
      mu_deriv_ptr=&de;
      return 0;
    }

    /// The default derivative object for calculating chemical potentials
    deriv_gsl<> def_mu_deriv;

    /// Return string denoting type ("eos_had_potential")
    virtual const char *type() { return "eos_had_potential"; }

  protected:
    
#ifndef DOXYGEN_INTERNAL

    /// Non-relativistic fermion thermodyanmics
    fermion_nonrel nrf;

    /// True of the derivative object has been set
    bool mu_deriv_set;

    /// The derivative object
    deriv_base<> *mu_deriv_ptr;
    
    /// Compute the momentum integral for \ref mdi_form
    double mom_integral(double pft, double pftp);
    
    /** \name The mode for the energy() function [protected] */
    //@{
    int mode;
    static const int nmode=1;
    static const int pmode=2;
    static const int normal=0;
    //@}
    
    /// Compute the energy
    double energy(double x);

#endif
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
