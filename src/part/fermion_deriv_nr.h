/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_DERIV_NR_H
#define O2SCL_FERMION_DERIV_NR_H

/** \file fermion_deriv_nr.h
    \brief File defining \ref o2scl::fermion_deriv_nr_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/root_cern.h>

#include <o2scl/part_deriv.h>
#include <o2scl/classical_deriv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a nonrelativistic fermion

      This does not include the rest mass energy in the chemical 
      potential or the rest mass energy density in the energy density
      to alleviate numerical precision problems at low densities

      This implements an equation of state for a nonrelativistic fermion
      using direct integration. After subtracting the rest mass from
      the chemical potentials, the distribution function is
      \f[
      \left\{1+\exp\left[\left(\frac{k^2}
      {2 m^{*}}-\nu\right)/T\right]\right\}^{-1}
      \f]
      where \f$ \nu \f$ is the effective chemical potential, \f$ m \f$ is
      the rest mass, and \f$ m^{*} \f$ is the effective mass.
      For later use, we define \f$ E^{*} = k^2/2/m^{*} \f$ .

      Uncertainties are given in \ref unc.

      \b Evaluation \b of \b the \b derivatives

      The relevant derivatives of the distribution function are
      \f[
      \frac{\partial f}{\partial T}=
      f(1-f)\frac{E^{*}-\nu}{T^2}
      \f]
      \f[
      \frac{\partial f}{\partial \nu}=
      f(1-f)\frac{1}{T}
      \f]
      \f[
      \frac{\partial f}{\partial k}=
      -f(1-f)\frac{k}{m^{*} T}
      \f]
      \f[
      \frac{\partial f}{\partial m^{*}}=
      f(1-f)\frac{k^2}{2 m^{*2} T}
      \f]

      We also need the derivative of the entropy integrand w.r.t. the 
      distribution function, which is quite simple
      \f[
      {\cal S}\equiv f \ln f +(1-f) \ln (1-f) \qquad
      \frac{\partial {\cal S}}{\partial f} = \ln 
      \left(\frac{f}{1-f}\right) = 
      \left(\frac{\nu-E^{*}}{T}\right)
      \f]
      where the entropy density is
      \f[
      s = - \frac{g}{2 \pi^2} \int_0^{\infty} {\cal S} k^2 d k
      \f]

      The derivatives can be integrated directly
      or they may be converted to integrals
      over the distribution function through an integration by parts
      \f[
      \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) g(k)\right|_{k=a}^{k=b}
      - \int_a^b g(k) \frac{d f(k)}{dk} dk 
      \f]
      using the distribution function for \f$ f(k) \f$ and 0 and \f$
      \infty \f$ as the limits, we have
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{d g(k)}{dk} f dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} g(k) f (1-f) \frac{k}{E^{*} T} dk 
      \f]
      as long as \f$ g(k) \f$ vanishes at \f$ k=0 \f$ .
      Rewriting,
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} h(k) f (1-f) dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} f \frac{T m^{*}}{k} 
      \left[ h^{\prime} - \frac{h}{k}\right] d k
      \f]
      as long as \f$ h(k)/k \f$ vanishes at \f$ k=0 \f$ .

      \b Explicit \b forms

      1) The derivative of the density wrt the chemical potential
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2}{T} f (1-f) dk
      \f]
      Using \f$ h(k)=k^2/T \f$ we get
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      m^{*} f dk
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-\nu)}{T^2} 
      f (1-f) dk
      \f]
      Using \f$ h(k)=k^2(E^{*}-\nu)/T^2 \f$ we get
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
      \left[m^{*} \left(E^{*}-\nu\right) -k^2\right] d k
      \f]

      3) The derivative of the entropy wrt the chemical potential
      \f[
      \left(\frac{d s}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-\nu)}{T^2} dk
      \f]
      This verifies the Maxwell relation
      \f[
      \left(\frac{d s}{d \mu}\right)_T =
      \left(\frac{d n}{d T}\right)_{\mu}
      \f]

      4) The derivative of the entropy wrt the temperature
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-\nu)^2}{T^3} dk
      \f]
      Using \f$ h(k)=k^2 (E^{*}-\nu)^2/T^3 \f$ 
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      f \frac{m^{*}}{T^2} \left[\left( E^{*}-\nu \right)^2 
      +\frac{2 k^2}{m^{*}} \left(E^{*}-\nu\right)\right] d k
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      \frac{k^2}{2 m^{* 2} T} f (1-f) k^2 dk
      \f]
      Using \f$ h(k)=k^4/(2 m^{* 2} T) \f$ we get
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} f
      \frac{3 k^2}{2 m^{*}} d k 
      \f]

      <b>Conversion to unitless variables:</b>

      After integrating by parts 
      \f$ u = k^2/2/m^{*}/T \f$ and \f$ y=\mu/T \f$, so
      \f[
      k d k = m^{*} T d u 
      \f]
      or
      \f[
      d k = \frac{m^{*} T}{\sqrt{2 m^{*} T u}} d u =
      \sqrt{\frac{m^{*} T}{2 u}} d u 
      \f]
      
      1) The derivative of the density wrt the chemical potential
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g m^{* 3/2} \sqrt{T}}{2^{3/2} \pi^2} \int_0^{\infty} 
      u^{-1/2} f d u
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g m^{* 3/2} \sqrt{T}}
      {2^{3/2} \pi^2} \int_0^{\infty} f d u
      \left[ 3 u^{1/2} - y u^{-1/2}\right]
      \f]

      4) The derivative of the entropy wrt the temperature
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g m^{* 3/2} T^{1/2}}{2^{3/2} \pi^2} \int_0^{\infty} 
      f \left[ 5 u^{3/2} - 6 y u^{1/2} + y^2 u^{-1/2}\right] d u
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      \frac{3 g m{* 1/2} T^{3/2}}{2^{3/2} \pi^2} 
      \int_0^{\infty} u^{1/2} f d u
      \f]
  */
  template<class fermion_deriv_t=fermion_deriv_tl<double>,
	   class fp_t=double>
    class fermion_deriv_nr_tl : public fermion_deriv_thermo_tl<fp_t> {

  public:

  /// Create a fermion with mass \c m and degeneracy \c g
  fermion_deriv_nr_tl() {
  
    flimit=20.0;

    density_root=&def_density_root;
  }

  virtual ~fermion_deriv_nr_tl() {
  }
  
  /** \brief The limit for the Fermi functions (default 20.0)
    
      fermion_deriv_nr will ignore corrections smaller than about
      \f$ \exp(-\mathrm{f{l}imit}) \f$ . 
  */
  fp_t flimit;
    
  /// Storage for the most recently calculated uncertainties 
  fermion_deriv unc;

  /** \brief Calculate properties as function of density
      at \f$ T=0 \f$
  */
  virtual void calc_density_zerot(fermion_deriv_t &f) {
    if (f.non_interacting) { f.ms=f.m; }
    f.kf=cbrt(6.0*this->pi2/f.g*f.n);
    f.nu=f.kf*f.kf/2.0/f.ms;
    f.ed=f.g*pow(f.kf,5.0)/20.0/this->pi2/f.ms;
    if (f.inc_rest_mass) {
      f.ed+=f.n*f.m;
      f.nu+=f.m;
    }
    f.pr=-f.ed+f.n*f.nu;
    f.en=0.0;
  
    if (f.non_interacting) { f.mu=f.nu; }

    f.dndT=0.0;
    if (f.inc_rest_mass) {
      f.dndmu=3.0*f.n/2.0/(f.nu-f.m);
    } else {
      f.dndmu=3.0*f.n/2.0/f.nu;
    }
    f.dsdT=0.0;
    return;
  }

    
  /** \brief Calculate properties as function of chemical potential
      at \f$ T=0 \f$
  */
  virtual void calc_mu_zerot(fermion_deriv_t &f) {
    if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
    if (f.inc_rest_mass) {
      f.kf=sqrt(2.0*f.ms*(f.nu-f.m));
    } else {
      f.kf=sqrt(2.0*f.ms*f.nu);
    }
    f.n=f.kf*f.kf*f.kf*f.g/6.0/this->pi2;
    f.ed=f.g*pow(f.kf,5.0)/20.0/this->pi2/f.ms;
    if (f.inc_rest_mass) f.ed+=f.n*f.m;
    f.pr=-f.ed+f.n*f.nu;
    f.en=0.0;

    f.dndT=0.0;
    if (f.inc_rest_mass) {
      f.dndmu=3.0*f.n/2.0/(f.nu-f.m);
    } else {
      f.dndmu=3.0*f.n/2.0/f.nu;
    }
    f.dsdT=0.0;
    return;
  }

    
  /** \brief Calculate properties as function of chemical potential
   */
  virtual int calc_mu(fermion_deriv_t &f, fp_t temper) {

    if (temper<0.0) {
      O2SCL_ERR("T<0 in fermion_deriv_nr_tl<fp_t>::calc_mu().",exc_einval);
    }
    if (temper==0.0) {
      calc_mu_zerot(f);
      return 0;
    }
    if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
    fp_t pfac2=f.g*pow(f.ms*temper/2.0/this->pi,1.5)/temper, y;
  
    if (f.inc_rest_mass) {
      y=(f.nu-f.m)/temper;
    } else {
      y=f.nu/temper;
    }
    
    fp_t half=gsl_sf_fermi_dirac_half(y);
    fp_t mhalf=gsl_sf_fermi_dirac_mhalf(y);
    fp_t thalf=gsl_sf_fermi_dirac_3half(y);
  
    // Number density:
    f.n=pfac2*half*temper;
  
    // Energy density:
    f.ed=pfac2*temper*temper*1.5*thalf;

    if (f.inc_rest_mass) {
    
      // Finish energy density
      f.ed+=f.n*f.m;
    
      // entropy density
      f.en=((f.ed-f.n*f.m)/0.6-(f.nu-f.m)*f.n)/temper;
    
      // pressure
      f.pr=(f.ed-f.n*f.m)/1.5;
    
    } else {
    
      // entropy density
      f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
    
      // pressure
      f.pr=2.0*f.ed/3.0;
    
    }
    
    f.dndT=pfac2*(1.5*half-y*mhalf);
    f.dndmu=pfac2*mhalf;
    f.dsdT=pfac2*(3.75*thalf-3.0*y*half+y*y*mhalf);
  
    return 0;
  }


  /** \brief Calculate properties as function of density
   */
  virtual int calc_density(fermion_deriv_t &f, fp_t temper) {

    if (f.m<0.0 || (f.non_interacting==false && f.ms<0.0)) {
      O2SCL_ERR2("Mass negative in ",
		 "fermion_deriv_nr_tl<fp_t>::calc_density().",exc_einval);
    }
    if (temper<0.0) {
      O2SCL_ERR2("Temperature less than zero in ",
		 "fermion_deriv_nr_tl<fp_t>::calc_density().",exc_einval);
    }
    if (f.n<0.0) {
      O2SCL_ERR2("Density less than zero in ",
		 "fermion_deriv_nr_tl<fp_t>::calc_density().",exc_einval);
    }

    if (temper==0.0) {
      calc_density_zerot(f);
      return 0;
    }
    if (f.n==0.0) {
      if (f.inc_rest_mass) {
	f.nu=f.m;
	f.mu=f.m;
      } else {
	f.nu=0.0;
	f.mu=0.0;
      }
      f.ed=0.0;
      f.en=0.0;
      f.pr=0.0;
      f.nu=0.0;
      f.dndT=0.0;
      f.dndmu=0.0;
      f.dsdT=0.0;
    }
  
    if (f.non_interacting==true) { f.ms=f.m; f.nu=f.mu; }

    nu_from_n(f,temper);
  
    if (f.non_interacting) { f.mu=f.nu; }
  
    calc_mu(f,temper);
  
    return 0;
  }

#ifdef O2SCL_NEVER_DEFINED
    // AWS, 12/29/2020, commenting this out since they are NR and pair
    // functions are not written for fermion_nonrel either
    
  /** \brief Calculate properties with antiparticles as function of
      chemical potential
  */
  virtual int pair_mu(fermion_deriv_t &f, fp_t temper) {
  
    if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

    fermion_deriv antip(f.ms,f.g);
    f.anti(antip);

    calc_mu(f,temper);

    calc_mu(antip,temper);

    f.n-=antip.n;
    f.pr+=antip.pr;
    f.ed+=antip.ed;
    f.en+=antip.en;
    f.dsdT+=antip.dsdT;
    f.dndT+=antip.dndT;
    f.dndmu+=antip.dndmu;
  
    return 0;
  }


  /** \brief Calculate properties with antiparticles as function of
      density
  */
  virtual int pair_density(fermion_deriv_t &f, fp_t temper) {
    fp_t nex;
  
    if (temper<=0.0) {
      O2SCL_ERR("T=0 not implemented in fermion_deriv_nr().",exc_eunimpl);
    }
    if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
    nex=f.nu/temper;
    funct mf=std::bind(std::mem_fn<fp_t(fp_t,fermion_deriv_t &,fp_t)>
		       (&fermion_deriv_nr_tl<fermion_deriv_t,fp_t>::pair_fun),
		       this,std::placeholders::_1,std::ref(f),temper);

    density_root->solve(nex,mf);
    f.nu=nex*temper;

    if (f.non_interacting==true) { f.mu=f.nu; }
  
    pair_mu(f,temper);

    return 0;
  }

#endif
    
  /// Calculate effective chemical potential from density
  virtual int nu_from_n(fermion_deriv_t &f, fp_t temper) {

    // Use initial value of nu for initial guess
    fp_t nex;
    if (f.inc_rest_mass) {
      nex=-(f.nu-f.m)/temper;
    } else {
      nex=f.nu/temper;
    }

    // Make a correction if nex is too small and negative
    // (Note GSL_LOG_DBL_MIN is about -708)
    if (nex>-GSL_LOG_DBL_MIN*0.9) nex=-GSL_LOG_DBL_MIN/2.0;
  
    funct mf=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
		       (&fermion_deriv_nr_tl<fermion_deriv_t,fp_t>::solve_fun),
		       this,std::placeholders::_1,f.n/f.g,f.ms*temper);
    
    // Turn off convergence errors temporarily, since we'll
    // try again if it fails
    bool enc=density_root->err_nonconv;
    density_root->err_nonconv=false;
    int ret=density_root->solve(nex,mf);
    density_root->err_nonconv=enc;

    if (ret!=0) {

      // If it failed, try to get a guess from classical_thermo particle
    
      classical_thermo cl;
      cl.calc_density(f,temper);
      if (f.inc_rest_mass) {
	nex=-(f.nu-f.m)/temper;
      } else {
	nex=-f.nu/temper;
      } 
      ret=density_root->solve(nex,mf);
    
      // If it failed again, add error information
      if (ret!=0) {
	O2SCL_ERR("Solver failed in fermion_nonrel::nu_from_n().",ret);
      }
    }

    if (f.inc_rest_mass) {
      f.nu=-nex*temper+f.m;
    } else {
      f.nu=-nex*temper;
    }

    return 0;
  }


  /** \brief Set the solver for use in calculating the chemical
      potential from the density */
  void set_density_root(root<> &rp) {
    density_root=&rp;
    return;
  }

  /// The default solver for npen_density() and pair_density()
  root_cern<> def_density_root;

  /// Return string denoting type ("fermion_deriv_nr")
  virtual const char *type() { return "fermion_deriv_nr"; };

  protected:

#ifndef DOXYGEN_INTERNAL

  /// Solver to compute chemical potential from density
  root<> *density_root;
    
  /// Function to compute chemical potential from density
  fp_t solve_fun(fp_t x, fp_t nog, fp_t msT) {

    fp_t nden;
  
    // If the argument to gsl_sf_fermi_dirac_half() is less
    // than GSL_LOG_DBL_MIN (which is about -708), then 
    // an underflow occurs. We just set nden to zero in this 
    // case, as this helps the solver find the right root.
  
    if (((-x)<GSL_LOG_DBL_MIN) || !std::isfinite(x)) nden=0.0;
    else nden=gsl_sf_fermi_dirac_half(-x)*sqrt(this->pi)/2.0;
  
    nden*=pow(2.0*msT,1.5)/4.0/this->pi2;
    return nden/nog-1.0;
  }


  /** \brief Function to compute chemical potential from density
      when antiparticles are included
  */
  fp_t pair_fun(fp_t x, fermion_deriv_t &f, fp_t T) {
    
    fp_t nden, y, yy;

    f.nu=T*x;

    // 6/6/03 - Should this be included? I think yes!
    if (f.non_interacting) f.mu=f.nu;

    if (f.inc_rest_mass) {
      y=(f.nu-f.m)/T;
    } else {
      y=f.nu/T;
    }

    nden=gsl_sf_fermi_dirac_half(y)*sqrt(this->pi)/2.0;
    nden*=f.g*pow(2.0*f.ms*T,1.5)/4.0/this->pi2;
  
    yy=nden;

    if (f.inc_rest_mass) {
      y=-(f.nu-f.m)/T;
    } else {
      y=-f.nu/T;
    }
  
    nden=gsl_sf_fermi_dirac_half(y)*sqrt(this->pi)/2.0;
    nden*=f.g*pow(2.0*f.ms*T,1.5)/4.0/this->pi2;
  
    yy-=nden;
  
    yy=yy/f.n-1.0;

    return yy;
  }


#endif

  };

  typedef fermion_deriv_nr_tl<> fermion_deriv_nr;

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
