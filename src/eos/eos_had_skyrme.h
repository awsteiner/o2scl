/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Skyrme hadronic equation of state 

      Equation of state of nucleonic matter based on 
      the Skryme interaction from \ref Skyrme59 .

      \hline
      \b Background:

      The Hamiltonian is defined (using the notation of 
      \ref Steiner05b ) as
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

      The values \f$ a=0, b=1 \f$ give the standard definition of the
      Skyrme Hamiltonian \ref Skyrme59, while \f$a=1, b=0\f$ contains
      the modifications suggested by \ref Onsi94.

      The spin-orbit term is (following \ref Steiner05)
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

      Also, couple useful definitions
      \f[
      t_3^{\prime} = \left(a + b\right) t_3 \, ,
      \f]
      \f[
      C = \frac{3 }{10 m} \left( \frac{3 \pi^2 }{2} \right)^{2/3}  \, ,
      \f]
      and
      \f[
      \beta = \frac{M}{2} \left[ \frac{1}{4} \left( 3 t_1 + 5 t_2 \right) \, .
      + t_2 x_2 \right] \\
      \f]
      
      \hline
      \b Units:

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
      
      \hline
      \b Misc:

      The functions for the usual saturation properties are based 
      partly on \ref Brack85.
      
      Models are taken from the references: \ref Bartel79, \ref
      Beiner75, \ref Chabanat95, \ref Chabanat97, \ref Danielewicz09,
      \ref Dobaczewski94, \ref Dutta86, \ref Friedrich86, \ref Onsi94,
      \ref Reinhard95, and \ref Tondeur84, and \ref VanGiai81 .
      
      The variables \f$ \nu_n\f$ and \f$ \nu_p\f$ contain the
      expressions \f$ (-\mu_n+V_n)/T \f$ and \f$ (-\mu_p+V_p)/T \f$
      respectively, where \f$ V \f$ is the potential part of the
      single particle energy for particle i (i.e. the derivative of
      the Hamiltonian w.r.t. density while energy density held
      constant). Equivalently, \f$ \nu_n\f$ is just \f$ -k_{F_n}^2/ 2
      m^{*} \f$.

      \note The finite temperature code does not include attempt to
      include antiparticles and uses \ref
      o2scl::fermion_nonrel::calc_density(). At finite temperature,
      pure neutron matter implies a zero proton number density which
      would imply that the proton chemical potential is \f$ - \infty
      \f$ . This class handles this situation by just setting \f$
      \nu_p \f$ to zero. The case of pure proton matter is handled
      similarly.
      
      \note Since this EOS uses the effective masses and chemical
      potentials in the fermion class, the values of
      <tt>part::non_interacting</tt> for neutrons and protons are set
      to false in many of the functions.

      \hline
      
      \todo
      - Convert W0 to b4 and b4p everywhere
      - Remove use of mnuc in calparfun()?
      - Document \ref o2scl_hdf::skyrme_load() file format.
      - Update reference list.

      \future
      - There is some code duplication between calc_e() and
      calc_temp_e() which could be simplified.
      - This EOS typically converges very well. One exception seems
      to be using <tt>calc_temp_p()</tt> at very low densities. I have
      had problems, for example, with <tt>mun=5.0, mup=6.5</tt>
      at <tt>T=1.0/197.33</tt>. 

      \hline
      
  */
  class eos_had_skyrme : public eos_had_temp_eden_base {

  public:

    /// \name Basic usage
    //@{
    /// Create a blank Skyrme EOS
    eos_had_skyrme();

    /// Destructor
    virtual ~eos_had_skyrme() {};

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

  protected:
    
    void hamiltonian_coeffs(double &ham1, double &ham2,
			    double &ham3, double &ham4,
			    double &ham5, double &ham6);
    
    /** \brief Handle the zero density limit
     */
    template<class fermion_t>
      void zero_density(fermion_t &ne, fermion_t &pr,
			thermo &th) {

      ne.ms=ne.m;
      pr.ms=pr.m;
      if (ne.inc_rest_mass) {
	ne.mu=ne.m;
	ne.nu=ne.m;
      } else {
	ne.mu=0.0;
	ne.nu=0.0;
      }
      if (pr.inc_rest_mass) {
	pr.mu=pr.m;
	pr.nu=pr.m;
      } else {
	pr.mu=0.0;
	pr.nu=0.0;
      }
      ne.pr=0.0;
      pr.pr=0.0;
      ne.ed=0.0;
      pr.ed=0.0;
      ne.en=0.0;
      pr.en=0.0;
      th.pr=0.0;
      th.ed=0.0;
      th.en=0.0;
      
      return;
    }

    /** \brief Compute the base thermodynamic quantities
     */
    template<class fermion_t>
      void base_thermo
      (fermion_t &ne, fermion_t &pr, double ltemper, thermo &locth,
       double term, double term2,
       double ham1, double ham2, double ham3, double ham4, double ham5,
       double ham6) {
      
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
      
      double opatpa=(1.0+alpha)*(2.0+alpha);
      double common2=2.0*ham1+2.0*ham2;
      double dhdnn2=common2+4.0*nna*opatpa+
	na/nb/nb*(ham5*nb*nb*opatpa+ham3*pr.n*alpha*
		  (ne.n+2.0*pr.n+ne.n*alpha)+
		  ham6*(4.0*ne.n*pr.n*(1.0+alpha)+ne.n*ne.n*opatpa+
			pr.n*pr.n*(2.0+alpha*(alpha-1.0))));
      double dhdnp2=common2+4.0*npa*opatpa+
	na/nb/nb*(ham5*nb*nb*opatpa+ham3*ne.n*alpha*
		  (pr.n+2.0*ne.n+pr.n*alpha)+
		  ham6*(4.0*ne.n*pr.n*(1.0+alpha)+pr.n*pr.n*opatpa+
			ne.n*ne.n*(2.0+alpha*(alpha-1.0))));
      double dhdnndnp=2.0*ham1+na/nb/nb*
	(ham5*nb*nb*opatpa+ham6*alpha*    
	 (4.0*ne.n*pr.n+ne.n*ne.n*(1.0+alpha)+pr.n*pr.n*(1.0+alpha))+
	 ham3*(ne.n*ne.n*(1.0+alpha)+pr.n*pr.n*(1.0+alpha)+
	       ne.n*pr.n*(2.0+alpha+alpha*alpha)));
      
      // For the kinetic part, convert from (mu,T) to (n,T)
      double n_dsdT_f=0.0, p_dsdT_f=0.0;
      double n_dmudT_f=0.0, p_dmudT_f=0.0;
      double n_dmudn_f=0.0, p_dmudn_f=0.0;
      ne.deriv_f(n_dmudn_f,n_dmudT_f,n_dsdT_f);
      pr.deriv_f(p_dmudn_f,p_dmudT_f,p_dsdT_f);
      
      double X_n, X_p;
      if (ltemper>0.0) {
	X_n=2.5*ne.ed-4.5*ne.ms*ne.n*ne.n/ltemper/ne.dndmu;
	X_p=2.5*pr.ed-4.5*pr.ms*pr.n*pr.n/ltemper/pr.dndmu;
      } else {
	X_n=2.5*ne.ed-3.75*ne.n/ne.dndmu*ne.ed/ne.nu;
	X_p=2.5*pr.ed-3.75*pr.n/pr.dndmu*pr.ed/pr.nu;
      }
      
      // Now combine to compute the six derivatives
      locthd.dsdT=n_dsdT_f+p_dsdT_f;
      locthd.dmundT=2.0*ltemper*ne.ms*(term+term2)*n_dsdT_f+
	2.0*ltemper*pr.ms*term*p_dsdT_f;
      locthd.dmupdT=2.0*ltemper*pr.ms*(term+term2)*p_dsdT_f+
	2.0*ltemper*ne.ms*term*n_dsdT_f;
      locthd.dmundnn=-4.0*ne.ms*ne.ms*pow(term+term2,2.0)*X_n-
	4.0*term*term*pr.ms*pr.ms*X_p+n_dmudn_f+dhdnn2;
      locthd.dmupdnp=-4.0*pr.ms*pr.ms*pow(term+term2,2.0)*X_p-
	4.0*term*term*ne.ms*ne.ms*X_p+p_dmudn_f+dhdnp2;
      locthd.dmudn_mixed=-4.0*(term+term2)*term*
	(ne.ms*ne.ms*X_n+pr.ms*pr.ms*X_p)+dhdnndnp;
      
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
		   "eos_had_skyrme::calc_deriv_temp_e().",exc_einval);
      }
      if (ne.n<0.0 || pr.n<0.0) {
	std::string str=((std::string)"Nucleon densities negative, n_n=")+
	  std::to_string(ne.n)+", n_p="+std::to_string(pr.n)+", in "+
	  "eos_had_skyrme::calc_deriv_temp_e().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
      if (fabs(ne.g-2.0)>1.0e-10 || fabs(pr.g-2.0)>1.0e-10) {
	O2SCL_ERR((((std::string)"Neutron (")+std::to_string(ne.g)+
		   ") or proton ("+std::to_string(pr.g)+") spin deg"+
		   "eneracies wrong in "+
		   "eos_had_skyrme::calc_deriv_temp_e().").c_str(),
		  exc_einval);
      }
      if (fabs(ne.m-4.5)>1.0 || fabs(pr.m-4.5)>1.0) {
	O2SCL_ERR((((std::string)"Neutron (")+std::to_string(ne.m)+
		   ") or proton ("+std::to_string(pr.m)+") masses wrong "+
		   "in eos_had_skyrme::calc_deriv_temp_e().").c_str(),
		  exc_einval);
      }
      if (ne.non_interacting==true || pr.non_interacting==true) {
	O2SCL_ERR2("Neutron or protons non-interacting in ",
		   "eos_had_skyrme::calc_deriv_temp_e().",exc_einval);
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
    
  public:
    
    /** \brief Equation of state as a function of densities

	\note Runs the zero temperature code if \c temper is less
	than or equal to zero.
     */
    virtual int calc_temp_e(fermion &ne, fermion &pr, double temper, 
			    thermo &th);

    /** \brief Equation of state including second derivatives
	as a function of the densities
    */
    virtual int calc_deriv_temp_e(fermion_deriv &ne, fermion_deriv &pr,
				  double temper, thermo &th,
				  thermo_np_deriv_helm &thd);
    
    /// Equation of state as a function of density.
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lt);
    //@}

    /// \name Basic Skyrme model parameters
    //@{
    double t0,t1,t2,t3,x0,x1,x2,x3,alpha,a,b;
    //@}

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

    /** \name Saturation properties

	These calculate the various saturation properties exactly from
	the parameters at any density. These routines often assume that 
	the neutron and proton masses are equal.
    */
    //@{
    /** \brief Calculate binding energy
      
	\f[
	\frac{E}{A} = C n_B^{2/3} \left( 1 + \beta n_B \right) + 
	\frac{3 t_0}{8} n_B + \frac{t_3^{\prime}}{16} n_B^{\alpha+1} 
	\f]
    */
    virtual double feoa(double nb);
  
    /** \brief Calculate effective mass
      
	\f[
	M^{*}/M = \left(1+ \beta n_B \right)^{-1} \\
	\f]
    */
    virtual double fmsom(double nb);

    /** \brief Calculate compressibility

	\f[
	K = 10 C n_B^{2/3} + \frac{27}{4} t_0 n_B + 40 C \beta n_B^{5/3} + 
	\frac{9 t_3^{\prime}}{16} 
	\alpha \left( \alpha+1 \right) n_B^{1 + \alpha} +
	\frac{9 t_3^{\prime}}{8} \left( \alpha+1 \right) n_B^{1 + \alpha}
	\f]
    */
    virtual double fcomp(double nb);

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

    /** \brief skewness

	\f[
	2 C n_B^{2/3} \left(9-5/M^{*}/M\right)+
	\frac{27 t_3^{\prime}}{16} n^{1+\alpha} \alpha 
	\left(\alpha^2-1\right)
	\f]
    */
    virtual double fkprime(double nb);
    //@}

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
      
	\todo Does this work for both 'a' and 'b' non-zero?
	
	\todo Compare to similar formulas in \ref Margueron02
    */
    int calpar(double gt0=-10.0, double gt3=70.0, double galpha=0.2,
	       double gt1=2.0, double gt2=-1.0);

    // Unfinished.
    /* \brief 
	From \ref Margueron02
    */
    //  int calpar_new(double m);

    /** \brief Use eos_had_base methods for saturation properties
      
	This can be set to true to check the difference between
	the exact expressions and the numerical values from
	class eos_had_base.
    */
    bool parent_method;
  
    /** \brief Check the Landau parameters for instabilities

	This returns zero if there are no instabilities.
     */
    int check_landau(double nb, double m);

    /** \brief Calculate the Landau parameters for nuclear matter

	Given \c n0 and \c m, this calculates the Landau parameters in
	nuclear matter as given in \ref Margueron02
     
	\todo This needs to be checked.
	
	(Checked once on 11/05/03)
    */
    void landau_nuclear(double n0, double m,
		       double &f0, double &g0, double &f0p,
		       double &g0p, double &f1, double &g1,
		       double &f1p, double &g1p);

    /** \brief Calculate the Landau parameters for neutron matter
	    
	Given 'n0' and 'm', this calculates the Landau parameters in
	neutron matter as given in \ref Margueron02
	
	\todo This needs to be checked
	
	(Checked once on 11/05/03)
    */
    void landau_neutron(double n0, double m, double &f0, double &g0, 
			double &f1, double &g1);

    /// Return string denoting type ("eos_had_skyrme")
    virtual const char *type() { return "eos_had_skyrme"; }

    /** \brief Set using alternate parameterization

	From \ref Bender03 as in, e.g. \ref Kortelainen14
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

	\todo These expressions are not exactly the same
	as those in \ref Bender03, so I need to find out why
	and make this more clear.
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

	This function uses the relations in \ref Kortelainen10 .
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

	\ref Kortelainen10 assumed equal neutron and proton masses, so
	this function uses \f$ \hbar^2/(2m) = \hbar^2/(m_n+m_p) \f$
	and the neutron and proton masses in \ref
	eos_had_base::def_neutron and \ref eos_had_base::def_proton,
	respectively. To obtain the results in the original paper, set
	neutron and proton masses to ensure that \f$ \hbar^2/(2m) =
	20.73553~\mathrm{MeV}~\mathrm{fm}^2 \f$ .
    */
    void alt_params_saturation
      (double n0, double EoA, double K, double Ms_star, double a, double L,
       double Mv_star, double CrDr0, double CrDr1, double CrnJ0, double CrnJ1);
 
    /// Thermodynamics of non-relativistic fermions
    fermion_nonrel nrf;
    
    /// Thermodynamics of non-relativistic fermions with derivatives
    fermion_deriv_nr nrfd;
    
#ifndef DOXYGEN_NO_O2NS
    
  protected:

    /// \name Functions and parameters for calpar()
    //@{
    int calparfun(size_t nv, const ubvector &x, ubvector &y);
    int calparfun2(size_t nv, const ubvector &x, ubvector &y);
    double fixn0, fixeoa, fixesym, fixcomp, fixmsom;
    //@}

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
