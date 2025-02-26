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
/** \file eos_had_skyrme.h
    \brief File defining \ref o2scl::eos_had_skyrme
*/
#ifndef O2SCL_SKYRME_EOS_H
#define O2SCL_SKYRME_EOS_H

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

namespace o2scl {
  
  /** \brief Skyrme hadronic equation of state 

      \verbatim embed:rst
      Equation of state of nucleonic matter based on 
      the Skryme interaction from [Skyrme59te]_.
      \endverbatim

      \b Hamiltonian

      \verbatim embed:rst
      The Hamiltonian is defined (using the notation of 
      [Steiner05b]_ as
      \endverbatim
      \f[
      {\cal H} = 
      {\cal H}_{k1} +
      {\cal H}_{k2} +
      {\cal H}_{k3} +
      {\cal H}_{p1} +
      {\cal H}_{p2} +
      {\cal H}_{p3} +
      {\cal H}_{g1} +
      {\cal H}_{g2} \, .
      \f]
      
      The kinetic terms are:
      \f[
      {\cal H}_{k1} = \frac{\tau_n}{2 m_n} +
      \frac{\tau_p}{2 m_p} 
      \f]
      \f[
      {\cal H}_{k2} =
      n \left(\tau_n + \tau_p \right) \left[ \frac{t_1}{4} 
      \left( 1 + \frac{x_1}{2} \right)
      + \frac{t_2}{4} \left( 1 + \frac{x_2}{2} \right) \right]
      \f]
      \f[
      {\cal H}_{k3} =
      \left( \tau_n n_n + \tau_p n_p \right) \left[ \frac{t_2}{4} 
      \left( \frac{1}{2} + x_2 \right)
      - \frac{t_1}{4} \left( \frac{1}{2} + x_1 \right) \right]
      \f]
      where \f$ \tau_i \f$ is the Fermi gas energy density
      of particle \f$ i \f$ .
      
      The potential terms are:
      \f[
      {\cal H}_{p1} = 
      \frac{t_0}{2} 
      \left[ \left( 1 + \frac{x_0}{2} \right) n^2 - 
      \left( {\textstyle \frac{1}{2}} + x_0 \right) 
      \left( n_n^2 + n_p^2 \right) \right] 
      \f]
      \f[
      {\cal H}_{p2} = 
      \frac{a t_3}{6} \left[ \left( 1 + \frac{x_3}{2} \right) n^{\alpha} 
      n_n n_p + 2^{\alpha-2} \left(1 - x_3\right)
      \left(n_n^{\alpha+2} + n_p^{\alpha+2}\right) \right] 
      \f]
      \f[
      {\cal H}_{p3} = 
      \frac{b t_3}{12} \left[ \left(1 + \frac{x_3}{2} \right) n^{\alpha+2} -
      \left(\frac{1}{2} + x_3 \right) n^{\alpha} 
      \left( n_n^2+n_p^2 \right) \right]
      \f]
      
      \b Gradient \b and \b Spin-Orbit \b Terms

      The gradient terms are displayed here for completeness even though
      they are not computed in the code:
      \f[
      {\cal H}_{g1} = 
      \frac{3}{32} \left[ t_1 \left(1 - x_1 \right) - 
      t_2 \left(1 + x_2 \right) \right] \left[ \left( \nabla n_n\right)^2 + 
      \left( \nabla n_p \right)^2 \right] 
      \f]
      \f[
      {\cal H}_{g2} = 
      \frac{1}{8} \left[ 3 t_1 \left( 1 + 
      \frac{x_1}{2} \right) - t_2 \left(1 + \frac{x_2}{2} \right) \right] 
      \nabla n_n \nabla n_p
      \f]

      \verbatim embed:rst
      The values :math:`a=0, b=1` give the standard definition of the
      Skyrme Hamiltonian [Skyrme59te]_, while :math:`a=1, b=0`
      contains the modifications suggested by [Onsi94]_.
      \endverbatim

      \verbatim embed:rst
      The spin-orbit term is (following [Steiner05b]_)
      \endverbatim
      \f[
      {\cal H}_{J} = -\frac{W_0}{2} \left( n_n \vec{\nabla} \cdot 
      \vec{J}_n + n_p \vec{\nabla} \cdot \vec{J}_p + n \vec{\nabla} 
      \cdot \vec{J} \right) + \frac{t_1}{16} \left(\vec{J}_n^2 + 
      \vec{J}_p^2 - x_1 \vec{J}^2\right) - \frac{t_2}{16} 
      \left(\vec{J}_n^2 + \vec{J}_p^2 + x_2 \vec{J}^2\right)
      \f]
      where sometimes the \f$ \vec{J}^2 \f$ terms are not included. 
      Alternatively, one can separate the isoscalar and isovector
      parts in the first term
      \f[
      {\cal H}_{J} = - b_4 n \vec{\nabla} \cdot \vec{J} -
      b_4^{\prime} n_n \vec{\nabla} \cdot \vec{J}_n - 
      b_4^{\prime} n_p \vec{\nabla} \cdot \vec{J}_p
      \f]
      then the earlier Skyrme interactions have \f$ b_4 = 
      b_4^{\prime} = W_0/2 \f$. For example, for SLy4, 
      \f$ b_4 = b_4^{\prime} = W_0/2 = 61.5~\mathrm{MeV} \f$.

      \verbatim embed:rst
      Three quantities are defined in [Steiner05b]_ for
      use in computing the properties of matter at saturation
      \endverbatim
      \f[
      t_3^{\prime} = \left(a + b\right) t_3 \, ,
      \f]
      \f[
      C = \frac{3 }{10 m} \left( \frac{3 \pi^2 }
      {2} \right)^{2/3}  \, ,
      \f]
      and
      \f[
      \beta = \frac{M}{2} \left[ \frac{1}{4} 
      \left( 3 t_1 + 5 t_2 \right)
      + t_2 x_2 \right] \, . \\
      \f]
      
      \b Units

      Quantities which have units containing powers of energy are
      divided by \f$\hbar c\f$ to ensure all quantities are in units
      of \f$ \mathrm{fm} \f$. The \f$x_i\f$ and \f$\alpha\f$ are
      unitless, while the original units of the \f$t_i\f$ are:
      - \f$t_0\f$ - \f$\mathrm{MeV}\f$ \f$\mathrm{fm}^3\f$
      - \f$t_1\f$ - \f$\mathrm{MeV}\f$ \f$\mathrm{fm}^5\f$
      - \f$t_2\f$ - \f$\mathrm{MeV}\f$ \f$\mathrm{fm}^5\f$
      - \f$t_3\f$ - \f$\mathrm{MeV}\f$ \f$\mathrm{fm}^{3(1+\alpha)}\f$
      
      These are stored internally with units of:
      - \f$t_0\f$ - \f$\mathrm{fm}^2\f$
      - \f$t_1\f$ - \f$\mathrm{fm}^4\f$
      - \f$t_2\f$ - \f$\mathrm{fm}^4\f$
      - \f$t_3\f$ - \f$\mathrm{fm}^{2+3 \alpha}\f$
      
      \b Other \b Notes

      \verbatim embed:rst
      The functions for the usual saturation properties are based 
      partly on [Brack85]_.

      Models are taken from the references: [Bartel79]_,
      [Beiner75]_, [Chabanat95]_, [Chabanat97]_, [Danielewicz09]_,
      [Dobaczewski94]_, [Dutta86]_, [Friedrich86]_, [Onsi94]_,
      [Reinhard95]_, [Reinhard99]_, [Tondeur84]_, and [VanGiai81]_ .
      \endverbatim
      
      The variables \f$ \nu_n \f$ and \f$ \nu_p \f$ contain the
      expressions \f$ (-\mu_n+V_n)/T \f$ and \f$ (-\mu_p+V_p)/T \f$
      respectively, where \f$ V \f$ is the potential part of the
      single particle energy for particle i (i.e. the derivative of
      the Hamiltonian w.r.t. density while energy density held
      constant). Equivalently, \f$ \nu_n\f$ is just \f$ -k_{F_n}^2/ 2
      m_n^{*} \f$.

      \note The finite temperature code does not attempt to include
      antiparticles and uses \ref
      o2scl::fermion_nonrel_tl::calc_density(). At finite temperature,
      pure neutron matter implies a zero proton number density which
      means the proton chemical potential is \f$ - \infty \f$ and thus
      set to the additive inverse of the return value of the function
      <tt>numeric_limits<double>::infinity()</tt> The case of pure
      proton matter is handled similarly. Negative densities result in
      calling the error handler.
      
      Skyrme models can be loaded using \ref o2scl_hdf::skyrme_load() .
      The full list is given in the \o2 repository in
      <tt>o2scl/data/o2scl/skdata/model_list</tt>.

      \b Todos
      
      \verbatim embed:rst
      .. todo:: 

         In class eos_had_skyrme:

         - Convert W0 to b4 and b4p everywhere
         - Remove use of mnuc in calparfun()?
         - Update reference list.

      \endverbatim

      \future
      - This EOS typically converges very well. One exception seems
      to be using <tt>calc_temp_p()</tt> at very low densities. I have
      had problems, for example, with <tt>mun=5.0, mup=6.5</tt>
      at <tt>T=1.0/197.33</tt>. 

  */
  class eos_had_skyrme : public eos_had_temp_eden_base {

  protected:

    /// \name EOS helper functions [protected]
    //@{
    /** \brief Compute the coefficients of the potential energy which
        have unique dependence on the densities
     */
    void hamiltonian_coeffs(double &ham1, double &ham2,
                            double &ham3, double &ham4,
                            double &ham5, double &ham6);
    
    /** \brief Compute the base thermodynamic quantities

        This function computes the energy density, pressure,
        entropy, and chemical potentials.
     */
    template<class fermion_t>
      void base_thermo
      (fermion_t &ne, fermion_t &pr, double ltemper, thermo &locth,
       double term, double term2, double ham1, double ham2,
       double ham3, double ham4, double ham5, double ham6) {
      
      double nb=ne.n+pr.n;
      double na=pow(nb,alpha);
      double npa=pow(pr.n,alpha);
      double nna=pow(ne.n,alpha);
      
      double ham=ne.ed+pr.ed+ham1*nb*nb+ham2*(ne.n*ne.n+pr.n*pr.n)+
        ham3*na*ne.n*pr.n+ham4*(nna*ne.n*ne.n+npa*pr.n*pr.n)+
        ham5*nb*nb*na+ham6*(ne.n*ne.n+pr.n*pr.n)*na;
      
      double gn, gp;
      if (ne.inc_rest_mass) {
        gn=2.0*ne.ms*(ne.ed-ne.n*ne.m);
      } else {
        gn=2.0*ne.ms*ne.ed;
      }
      if (pr.inc_rest_mass) {
        gp=2.0*pr.ms*(pr.ed-pr.n*pr.m);
      } else {
        gp=2.0*pr.ms*pr.ed;
      }
      
      // Variables dhdn{n,p} are the partial derivatives of the
      // Hamiltonian wrt the neutron and proton densities
      double common=2.0*ham1*nb+ham5*(2.0+alpha)*nb*na;
      double dhdnn=common+2.0*ham2*ne.n+ham3*na*pr.n*(alpha*ne.n/nb+1.0)+
        ham4*(nna*ne.n*(2.0+alpha))+
        ham6*(2.0*ne.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/nb);
      double dhdnp=common+2.0*ham2*pr.n+ham3*na*ne.n*(alpha*pr.n/nb+1.0)+
        ham4*(npa*pr.n*(2.0+alpha))+
        ham6*(2.0*pr.n*na+(ne.n*ne.n+pr.n*pr.n)*alpha*na/nb);

      // Compute the chemical potentials
      ne.mu=ne.nu+dhdnn+(gn+gp)*term+gn*term2;
      pr.mu=pr.nu+dhdnp+(gn+gp)*term+gp*term2;
      
      // Thermodynamics
      locth.ed=ham;
      locth.en=ne.en+pr.en;
      locth.pr=ltemper*locth.en+ne.mu*ne.n+pr.mu*pr.n-locth.ed;
      
      return;
    }
    
    /** \brief Compute second derivatives of the free energy
     */
    template<class fermion_t>
      void second_deriv
      (fermion_t &ne, fermion_t &pr, double ltemper, thermo &locth,
       thermo_np_deriv_helm &locthd, double term, double term2,
       double ham1, double ham2, double ham3,
       double ham4, double ham5, double ham6) {

      double nb=ne.n+pr.n;
      double na=pow(nb,alpha);
      double npa=pow(pr.n,alpha);
      double nna=pow(ne.n,alpha);

      // Second derivatives of the potential energy
      
      double opatpa=(1.0+alpha)*(2.0+alpha);
      double common2=2.0*ham1+2.0*ham2;
      double dhdnn2=common2+ham4*nna*opatpa+ham5*na*opatpa+
        ham3*(2.0*pr.n/nb*na*alpha+
              ne.n*pr.n*na/nb/nb*(-1.0+alpha)*alpha)+
        ham6*(2.0*na+4.0*ne.n/nb*na*alpha+na/nb/nb*
              (ne.n*ne.n+pr.n*pr.n)*(-1.0+alpha)*alpha);
      double dhdnp2=common2+ham4*npa*opatpa+ham5*na*opatpa+
        ham3*(2.0*ne.n/nb*na*alpha+
              ne.n*pr.n*na/nb/nb*(-1.0+alpha)*alpha)+
        ham6*(2.0*na+4.0*pr.n/nb*na*alpha+na/nb/nb*
              (ne.n*ne.n+pr.n*pr.n)*(-1.0+alpha)*alpha);
      double dhdnndnp=2.0*ham1+na/nb/nb*
        (ham5*nb*nb*opatpa+ham6*alpha*    
         (4.0*ne.n*pr.n+ne.n*ne.n*(1.0+alpha)+pr.n*pr.n*(1.0+alpha))+
         ham3*(ne.n*ne.n*(1.0+alpha)+pr.n*pr.n*(1.0+alpha)+
               ne.n*pr.n*(2.0+alpha+alpha*alpha)));
      
      // For the kinetic part, convert from the (mu,T)
      // to (n,T) representation
      
      double n_dsdT_f=0.0, p_dsdT_f=0.0;
      double n_dmudT_f=0.0, p_dmudT_f=0.0;
      double n_dmudn_f=0.0, p_dmudn_f=0.0;
      ne.deriv_f(n_dmudn_f,n_dmudT_f,n_dsdT_f);
      pr.deriv_f(p_dmudn_f,p_dmudT_f,p_dsdT_f);

      // The product 10 mstar^2 epsilon
      
      double hn, hp;
      if (ne.inc_rest_mass) {
        hn=10.0*ne.ms*(ne.ed-ne.n*ne.m)*ne.ms;
      } else {
        hn=10.0*ne.ms*ne.ed*ne.ms;
      }
      if (pr.inc_rest_mass) {
        hp=10.0*pr.ms*(pr.ed-pr.n*pr.m)*pr.ms;
      } else {
        hp=10.0*pr.ms*pr.ed*pr.ms;
      }

      // Now combine to compute the six derivatives
      
      locthd.dsdT=n_dsdT_f+p_dsdT_f;
      
      locthd.dmundT=2.0*ltemper*ne.ms*(term+term2)*n_dsdT_f+
        2.0*ltemper*pr.ms*term*p_dsdT_f+n_dmudT_f;
      
      locthd.dmupdT=2.0*ltemper*pr.ms*(term+term2)*p_dsdT_f+
        2.0*ltemper*ne.ms*term*n_dsdT_f+p_dmudT_f;

      locthd.dmundnn=-(term+term2)*(term+term2)*hn-term*term*hp+
        pow(1.0+3.0*(term+term2)*ne.n*ne.ms,2.0)*n_dmudn_f+
        9.0*term*term*pr.n*pr.ms*pr.n*pr.ms*p_dmudn_f+dhdnn2;
      
      locthd.dmupdnp=-(term+term2)*(term+term2)*hp-term*term*hn+
        pow(1.0+3.0*(term+term2)*pr.n*pr.ms,2.0)*p_dmudn_f+
        9.0*term*term*ne.n*ne.ms*ne.n*ne.ms*n_dmudn_f+dhdnp2;
      
      locthd.dmudn_mixed=-term*(term+term2)*(hn+hp)+
        3.0*term*ne.ms*ne.n*(1.0+3.0*(term+term2)*ne.n*ne.ms)*n_dmudn_f+
        3.0*term*pr.ms*pr.n*(1.0+3.0*(term+term2)*pr.n*pr.ms)*p_dmudn_f+
        dhdnndnp;

      return;
    }

    /** Check the EOS input

        This function checks that
        - the densities and temperature are finite and positive
        - the spin denegeracies are correct
        - the masses are sensible
        - the values of 'non_interacting' are false
        - the alpha parameter is positive
        - the temperature is not negative
    */
    template<class fermion_t>
      void check_input(fermion_t &ne, fermion_t &pr, double T) {
      
      if (!std::isfinite(ne.n) || !std::isfinite(pr.n) ||
          !std::isfinite(T)) {
        O2SCL_ERR2("Nucleon densities or temperature not finite in ",
                   "eos_had_skyrme::check_input().",exc_einval);
      }
      if (ne.n<0.0 || pr.n<0.0) {
        std::string str=((std::string)"Nucleon densities negative, n_n=")+
          std::to_string(ne.n)+", n_p="+std::to_string(pr.n)+", in "+
          "eos_had_skyrme::check_input().";
        O2SCL_ERR(str.c_str(),exc_einval);
      }
      if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
        O2SCL_ERR((((std::string)"Neutron (")+std::to_string(ne.g)+
                   ") or proton ("+std::to_string(pr.g)+") spin deg"+
                   "eneracies wrong in "+
                   "eos_had_skyrme::check_input().").c_str(),
                  exc_einval);
      }
      if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
        O2SCL_ERR((((std::string)"Neutron (")+std::to_string(ne.m)+
                   ") or proton ("+std::to_string(pr.m)+") masses wrong "+
                   "in eos_had_skyrme::check_input().").c_str(),
                  exc_einval);
      }
      if (ne.non_interacting==true || pr.non_interacting==true) {
        O2SCL_ERR2("Neutron or protons non-interacting in ",
                   "eos_had_skyrme::check_input().",exc_einval);
      }
      if (alpha<=0.0) {
        O2SCL_ERR2("Parameter alpha negative in ",
                   "eos_had_skyrme::calc_e().",exc_einval);
      }
      if (T<0.0) {
        O2SCL_ERR2("Temperature negative in ",
                   "eos_had_skyrme::calc_e().",exc_einval);
      }
      
      return;
    }
    //@}
    
  public:

    /// \name Basic usage
    //@{
    /// Create a blank Skyrme EOS
    eos_had_skyrme();

    /// Destructor
    virtual ~eos_had_skyrme() {};

    /// Copy constructor
    eos_had_skyrme(const eos_had_skyrme &f) {
      this->def_thermo=f.def_thermo;
      this->eoa=f.eoa;
      this->n0=f.n0;
      this->comp=f.comp;
      this->esym=f.esym;
      this->msom=f.msom;
      this->kprime=f.kprime;
      this->err_nonconv=f.err_nonconv;
      this->t0=f.t0;
      this->t1=f.t1;
      this->t2=f.t2;
      this->t3=f.t3;
      this->x0=f.x0;
      this->x1=f.x1;
      this->x2=f.x2;
      this->x3=f.x3;
      this->alpha=f.alpha;
      this->a=f.a;
      this->b=f.b;
      this->W0=f.W0;
      this->b4=f.b4;
      this->b4p=f.b4p;
      this->reference=f.reference;
      this->parent_method=f.parent_method;
      //this->nrf=f.nrf;
      //this->nrfd=f.nrfd;
      this->def_fet=f.def_fet;
      this->def_deriv=f.def_deriv;
      this->def_deriv2=f.def_deriv2;
      this->def_neutron=f.def_neutron;
      this->def_proton=f.def_proton;
    }
    
    /// Copy construction with operator=()
    eos_had_skyrme &operator=(const eos_had_skyrme &f) {
      if (this!=&f) {
        this->def_thermo=f.def_thermo;
        this->eoa=f.eoa;
        this->n0=f.n0;
        this->comp=f.comp;
        this->esym=f.esym;
        this->msom=f.msom;
        this->kprime=f.kprime;
        this->err_nonconv=f.err_nonconv;
        this->t0=f.t0;
        this->t1=f.t1;
        this->t2=f.t2;
        this->t3=f.t3;
        this->x0=f.x0;
        this->x1=f.x1;
        this->x2=f.x2;
        this->x3=f.x3;
        this->alpha=f.alpha;
        this->a=f.a;
        this->b=f.b;
        this->W0=f.W0;
        this->b4=f.b4;
        this->b4p=f.b4p;
        this->reference=f.reference;
        this->parent_method=f.parent_method;
        //this->nrf=f.nrf;
        //this->nrfd=f.nrfd;
        this->def_fet=f.def_fet;
        this->def_deriv=f.def_deriv;
        this->def_deriv2=f.def_deriv2;
        this->def_neutron=f.def_neutron;
        this->def_proton=f.def_proton;
       }
      return *this;
    }
    
    /** \brief Equation of state as a function of densities
        at finite temperature
    */
    virtual int calc_temp_e(fermion &ne, fermion &pr, double temper, 
                            thermo &th);

    /** \brief Equation of state as a function of the densities at
        finite temperature and including second derivatives
    */
    virtual int calc_deriv_temp_e(fermion_deriv &ne, fermion_deriv &pr,
                                  double temper, thermo &th,
                                  thermo_np_deriv_helm &thd);
    
    /** \brief Equation of state as a function of densities at 
        zero temperature
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lt);

    /** \brief Equation of state as a function of densities at 
        zero temperature including second derivatives
    */
    virtual int calc_deriv_e(fermion_deriv &ne, fermion_deriv &pr,
                             thermo &th, thermo_np_deriv_helm &thd);
    
    /// Return string denoting type ("eos_had_skyrme")
    virtual const char *type() { return "eos_had_skyrme"; }
    //@}

    /// \name Basic Skyrme model parameters
    //@{
    double t0, t1, t2, t3, x0, x1, x2, x3, alpha, a, b;
    //@}

    /// \name Other parameters
    //@{
    /** \brief Spin-orbit splitting (in \f$ \mathrm{fm}^{-1} \f$)

        This is unused, but included for possible future use and
        present in the internally stored models.
    */
    double W0;

    /// Isoscalar spin-orbit term (in \f$ \mathrm{fm}^{-1} \f$)
    double b4;

    /// Isovector spin-orbit term (in \f$ \mathrm{fm}^{-1} \f$)
    double b4p;

    /// Bibliographic reference
    std::string reference;
    
    /** \brief Use eos_had_base methods for saturation properties
        
        This can be set to true to check the difference between
        the exact expressions and the numerical values from
        class eos_had_base.
    */
    bool parent_method;
    //@}

    /** \name Saturation properties

        These calculate the various saturation properties exactly from
        the parameters at any density. These routines often assume that 
        the neutron and proton masses are equal.
    */
    //@{
    /** \brief Calculate binding energy in symmetric matter
      
        \f[
        \frac{E}{A} = C n_B^{2/3} \left( 1 + \beta n_B \right) + 
        \frac{3 t_0}{8} n_B + \frac{t_3^{\prime}}{16} n_B^{\alpha+1} 
        \f]
    */
    virtual double feoa_symm(double nb);
  
    /** \brief Calculate effective mass in symmetric matter
      
        \f[
        M^{*}/M = \left(1+ \beta n_B \right)^{-1} \\
        \f]
    */
    virtual double fmsom_symm(double nb);

    /** \brief Calculate compressibility in nuclear (isospin-symmetric
        matter)

        \f[
        K = 10 C n_B^{2/3} + \frac{27}{4} t_0 n_B + 40 C \beta n_B^{5/3} + 
        \frac{9 t_3^{\prime}}{16} 
        \alpha \left( \alpha+1 \right) n_B^{1 + \alpha} +
        \frac{9 t_3^{\prime}}{8} \left( \alpha+1 \right) n_B^{1 + \alpha}
        \f]
    */
    virtual double fcomp_nuc(double nb);

    /** \brief Calculate symmetry energy

        If pf=0.5, then the exact expression below is used.
        Otherwise, the method from class eos_had_base is used.

        \f[
        E_{sym} = \frac{5}{9} C n^{2/3} + \frac{10 C m}{3}
        \left[ \frac{t_2}{6} \left(1 + \frac{5}{4} x_2 \right) - 
        \frac{1}{8} t_1 x_1 \right] n^{5/3} 
        - \frac{t_3^{\prime}}{24} 
        \left({\textstyle \frac{1}{2}} + x_3 \right) n^{1+\alpha} - 
        \frac{t_0}{4} \left( {\textstyle \frac{1}{2}} + x_0 \right) n 
        \f]
    */
    virtual double fesym(double nb, double alpha=0.0);

    /** \brief Skewness in nuclear (isospin-symetric) matter

        \f[
        2 C n_B^{2/3} \left(9-5/M^{*}/M\right)+
        \frac{27 t_3^{\prime}}{16} n^{1+\alpha} \alpha 
        \left(\alpha^2-1\right)
        \f]
    */
    virtual double fkprime_nuc(double nb);
    //@}

    /// \name Compute and test Landau parameters
    //@{
    /** \brief Check the Landau parameters for instabilities

        This returns zero if there are no instabilities.
     */
    int check_landau(double nb, double m);

    /** \brief Calculate the Landau parameters for nuclear matter

        \verbatim embed:rst
        Given ``n0`` and ``m``, this calculates the Landau parameters in
        nuclear matter as given in [Margueron02]_.

        .. todo:: 

           - This function, eos_had_skyrme::landau_nuclear() needs to 
             be checked.

        \endverbatim
        
        (Checked once on 11/05/03)
    */
    void landau_nuclear(double n0, double m,
                        double &f0, double &g0, double &f0p,
                        double &g0p, double &f1, double &g1,
                        double &f1p, double &g1p);

    /** \brief Calculate the Landau parameters for neutron matter
            
        \verbatim embed:rst
        Given ``n0`` and ``m``, this calculates the Landau parameters in
        nuclear matter as given in [Margueron02]_.
        
        .. todo:: 

           - This function, eos_had_skyrme::landau_neutron() needs to 
             be checked.
        \endverbatim
        
        (Checked once on 11/05/03)
    */
    void landau_neutron(double n0, double m, double &f0, double &g0, 
                        double &f1, double &g1);

    //@}
    
    /// \name Other functions
    //@{
    /** \brief Evaluate the effective masses for neutrons and
        protons
    */
    template<class fermion_t>
      void eff_mass(fermion_t &ne, fermion_t &pr, double &term,
                    double &term2) {
      
      // Landau effective masses
      double nb=ne.n+pr.n;
      term=0.25*(t1*(1.0+x1/2.0)+t2*(1.0+x2/2.0));
      term2=0.25*(t2*(0.5+x2)-t1*(0.5+x1));
      ne.ms=ne.m/(1.0+2.0*(nb*term+ne.n*term2)*ne.m);
      pr.ms=pr.m/(1.0+2.0*(nb*term+pr.n*term2)*pr.m);
      return;
    }

    /** \brief Calculate \f$ t_0,t_1,t_2,t_3 \f$ and \f$ \alpha \f$ from 
        the saturation properties.
      
        In nuclear matter: 
      
        \f$ E_b=E_b(n_0,M^{*},t_0,t_3,\alpha) \f$ \n
        \f$ P=P(n_0,M^{*},t_0,t_3,\alpha) \f$ \n
        \f$ K=K(n_0,M^{*},t_3,\alpha) \f$ 
        (the \f$ t_0 \f$ dependence vanishes) \n
        \f$ M^{*}=M^{*}(n_0,t_1,t_2,x_2) \f$ 
        (the \f$ x_1 \f$ dependence cancels), \n
        \f$ E_{sym}=E_{sym}(x_0,x_1,x_2,x_3,t_0,t_1,t_2,t_3,\alpha) \f$
      
        To fix the couplings from the saturation properties, we take
        \f$ n_0, M^{*}, E_b, K \f$ as inputs, and we can fix \f$
        t_0,t_3,\alpha \f$ from the first three relations, then use
        \f$ M^{*}, E_b \f$ to fix \f$ t_2 \f$ and \f$ t_1 \f$.  The
        separation into two solution steps should make for better
        convergence. All of the x's are free parameters and should be
        set before the function call.
      
        The arguments \c gt0, \c gt3, \c galpha, \c gt1, and \c gt2
        are used as initial guesses for skyme_eos::t0, eos_had_skyrme::t3,
        eos_had_skyrme::alpha, eos_had_skyrme::t1, and eos_had_skyrme::t2
        respectively.
      
        \verbatim embed:rst
        .. todo:: 

           In function eos_had_skyrme::calpar(): 

           - Does this work for both 'a' and 'b' non-zero?
           - Compare to similar formulas in [Margueron02]_.
        \endverbatim
    */
    int calpar(double gt0=-10.0, double gt3=70.0, double galpha=0.2,
               double gt1=2.0, double gt2=-1.0);

    // Unfinished.
    /* \brief 
        From \ref Margueron02
    */
    //  int calpar_new(double m);

    /** \brief Set using alternate parameterization

        \verbatim embed:rst
        From [Bender03]_ as in, e.g. [Kortelainen14]_
        \endverbatim
        \f{eqnarray*}
        C^{\rho \rho}_{00} &=& 3 t_0/8 \nonumber \\
        C^{\rho \rho}_{10} &=& -t_0/4 \left(\frac{1}{2}+x_0 \right) 
        \nonumber \\
        C^{\rho \rho}_{0D} &=& t_3/16 \nonumber \\
        C^{\rho \rho}_{1D} &=& -t_3/24 \left(\frac{1}{2}+x_3\right)
        \nonumber \\
        C^{\rho \tau}_{0} &=& 3 t_1/16+t_2/4\left(\frac{5}{4}+x_2\right)
        \nonumber \\
        C^{\rho \tau}_{1} &=& -t_1/8 \left(\frac{1}{2}+x_1\right) +
        t_2/8 \left(\frac{1}{2}+x_2\right) \nonumber \\
        C^{\rho \Delta \rho}_{0} &=& -9/64 t_1+t_2/16 
        \left(\frac{5}{4}+x_2\right) \nonumber \\
        C^{\rho \Delta \rho}_{1} &=& 3/32 t_1 \left(\frac{1}{2}+x_1\right) +
        t_2/32 \left(\frac{1}{2}+x_2\right) \nonumber \\
        C^{\rho \nabla J}_{0} &=& -b_4 -b_4^{\prime}/2 \nonumber \\
        C^{\rho \nabla J}_{1} &=& -b_4^{\prime}/2
        \f}

        The parameters should have the following units
        - <tt>Crr00</tt>: \f$ \mathrm{fm}^2 \f$
        - <tt>Crr10</tt>: \f$ \mathrm{fm}^2 \f$
        - <tt>Crr0D</tt>: \f$ \mathrm{fm}^{3 \alpha+2} \f$
        - <tt>Crr1D</tt>: \f$ \mathrm{fm}^{3 \alpha+2} \f$
        - <tt>Crt0</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>Crt1</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrDr0</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrDr1</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrnJ0</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>CrnJ1</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>alpha2</tt>: unitless

        \verbatim embed:rst
        .. todo:: 

           - In eos_had_skyrme::alt_params_set(): These expressions
             are not exactly the same as those in [Bender03]_, so I
             need to find out why and make this more clear.

        \endverbatim
    */
    void alt_params_set
      (double Crr00, double Crr10, double Crr0D, double Crr1D, 
       double Crt0, double Crt1, double CrDr0, double CrDr1, 
       double CrnJ0, double CrnJ1, double alpha2);

    /** \brief Get alternate parameterization
        
        The parameters will have the following units
        - <tt>Crr00</tt>: \f$ \mathrm{fm}^2 \f$
        - <tt>Crr10</tt>: \f$ \mathrm{fm}^2 \f$
        - <tt>Crr0D</tt>: \f$ \mathrm{fm}^{3 \alpha+2} \f$
        - <tt>Crr1D</tt>: \f$ \mathrm{fm}^{3 \alpha+2} \f$
        - <tt>Crt0</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>Crt1</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrDr0</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrDr1</tt>: \f$ \mathrm{fm}^4 \f$
        - <tt>CrnJ0</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>CrnJ1</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>alpha2</tt>: unitless

        See \ref eos_had_skyrme::alt_params_set().
    */
    void alt_params_get
      (double &Crr00, double &Crr10, double &Crr0D, double &Crr1D, 
       double &Crt0, double &Crt1, double &CrDr0, double &CrDr1, 
       double &CrnJ0, double &CrnJ1, double &alpha2);

    /** \brief Use the specified saturation properties and couplings
        and the function \ref alt_params_set() to set the 
        Skyrme coefficients

        \verbatim embed:rst
        See [Kortelainen10]_.
        \endverbatim

        This function uses the relations in Kortelainen et al. (2010),
        The parameters should have the following units
        - <tt>n0</tt>: \f$ \mathrm{fm}^{-3} \f$
        - <tt>EoA</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>K</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>Ms_star</tt>: unitless
        - <tt>a</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>L</tt>: \f$ \mathrm{fm}^{-1} \f$
        - <tt>Mv_star</tt>: unitless
        - <tt>CrDr0</tt>: \f$ \mathrm{fm}^{-3} \f$
        - <tt>CrDr1</tt>: \f$ \mathrm{fm}^{-3} \f$
        - <tt>CrnJ0</tt>: \f$ \mathrm{fm}^{-3} \f$
        - <tt>CrnJ1</tt>: \f$ \mathrm{fm}^{-3} \f$

        Kortelainen et al. (2010) assumed equal neutron and proton
        masses, so this function uses \f$ \hbar^2/(2m) =
        \hbar^2/(m_n+m_p) \f$ and the neutron and proton masses in
        \ref eos_had_base::def_neutron and \ref
        eos_had_base::def_proton, respectively. To obtain the results
        in the original paper, set neutron and proton masses to ensure
        that \f$ \hbar^2/(2m) = 20.73553~\mathrm{MeV}~\mathrm{fm}^2
        \f$ .
    */
    void alt_params_saturation
      (double n0, double EoA, double K, double Ms_star, double a, double L,
       double Mv_star, double CrDr0, double CrDr1, double CrnJ0, double CrnJ1);
    //@}
    
    /// \name Particle classes
    //@{
    /// Thermodynamics of non-relativistic fermions
    fermion_nonrel nrf;
    
    /// Thermodynamics of non-relativistic fermions with derivatives
    fermion_deriv_nr nrfd;
    //@}
    
  protected:

    /// \name Functions and parameters for calpar() [protected]
    //@{
    int calparfun(size_t nv, const ubvector &x, ubvector &y);
    int calparfun2(size_t nv, const ubvector &x, ubvector &y);
    double fixn0, fixeoa, fixesym, fixcomp, fixmsom;
    //@}

  };

}

#endif
