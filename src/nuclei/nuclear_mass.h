/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCLEAR_MASS_H
#define O2SCL_NUCLEAR_MASS_H

#include <cmath>
#include <string>
#include <map>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/nucleus.h>
#include <o2scl/constants.h>
#include <o2scl/table.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/root_cern.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nuclear mass information

      This class exists to provide some basic information on 
      nuclei to nuclear mass classes which are children of
      \ref nuclear_mass. 
      
      Elements 113, 115, 117 and 118 are named "Uut", "Uup", "Uus",
      and "Uuo", respectively, until the IUPAC decides on their names.
      Note that some of the nuclear mass tables use older or
      alternative names for the heavier elements, so \ref Ztoel() may
      return something different than is stored in, e.g., \ref
      ame_mass::entry::el.

      \future Add the full names for each element.
  */
  class nuclear_mass_info {
  public:
    
    nuclear_mass_info();
    
    /** \brief Parse a string representing an element

	Accepts strings of one of the following forms:
	- <tt>Pb208</tt>
	- <tt>pb208</tt>
	- <tt>Pb 208</tt>
	- <tt>Pb-208</tt>
	- <tt>pb 208</tt>
	- <tt>pb-208</tt>
	or one of the special strings <tt>n</tt>, <tt>p</tt>, <tt>d</tt>
	or <tt>t</tt> for the neutron, proton, deuteron, and triton, 
	respectively. This function also allows the value of A
	to precede the element symbol.

	\note At present, this allows nuclei which don't make sense
	because A<Z, such as Carbon-5.

	\future Warn about malformed combinations like Carbon-5
	\future Right now, <tt>n4</tt> is interpreted incorrectly
	as Nitrogen-4, rather than the tetraneutron.
	\future Interpret strings with the full name rather
	than just the abbreviation.
    */
    int parse_elstring(std::string ela, int &Z, int &N, int &A);    

    /** \brief Return Z given the element name abbreviation
       
        If the string parameter \c el is invalid, the error handler is
        called and the value -1 is returned.
    */
    int eltoZ(std::string el);
    
    /** \brief Return the element name abbreviation given Z
       
        \note This function returns \c "n" indicating the neutron for
        Z=0, and if the argument \c Z is greater than 118, an empty
        string is returned after calling the error handler.
    */
    std::string Ztoel(size_t Z);

    /** \brief Return a string of the form "Pb208" for a given Z and N

        Note that if \c Z is zero, then and \c 'n' is 
        used to indicate the a nucleus composed entirely of neutrons 
        and if the argument \c Z is greater than 118, an
        empty string is returned (independ.
    */
    std::string tostring(size_t Z, size_t N);

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The number of elements (proton number)
    static const int nelements=119;
    
    /** \brief A map containing the proton numbers organized by 
	element abbreviation
    */
    std::map<std::string,int,string_comp> element_table;
    
    /// A convenient typedef for an iterator for element_table
    typedef std::map<std::string,int,string_comp>::iterator table_it;
    
    /// The list of elements organized by proton number
    std::string element_list[nelements];

#endif
    
  };
  
  /** \brief Nuclear mass formula base [abstract base]

      (See also the discussion in \ref nuclei_section.)

      This is abstract base class for the nuclear mass formulas. Some
      mass formulas are undefined for sufficiently exotic nuclei. You
      can use the function is_included() to find if a particular
      \nucleus is included or not in a particular mass formula.

      The quantities below are returned in units of MeV. The functions
      include a version which takes Z and N as integers and a version
      with a suffix <tt>"_d"</tt> which takes Z and N as
      double-precision numbers.

      The mass excess is given by \ref mass_excess() and \ref
      mass_excess_d() .

      Binding energies (\ref binding_energy() and \ref
      binding_energy_d() ) are determined from mass excesses by
      \f[
      \mathrm{binding~energy} = A u - Z \left(m_p + m_e\right) - N m_n +
      \mathrm{mass~excess}
      \f]
      The neutron, proton, and electron masses and atomic mass unit
      are stored in \ref m_prot, \ref m_neut, \ref m_elec, and \ref
      m_amu . By default, this are assigned to the values in \ref
      o2scl_mks times \ref o2scl_const::hc_mev_fm , but these default
      values are modified in the constructors of some children
      classes.

      Total masses, as returned by \ref total_mass() and \ref
      total_mass_d() , are the mass of the nuclide without the
      electron mass or binding energy contribution
      \f[
      \mathrm{total~mass} = \mathrm{mass~excess} + A u - Z m_e
      \f]

      Atomic masses are the total mass with the electron mass and
      binding energy contributions (see \ref atomic_mass() and \ref
      atomic_mass_d() ). Electron binding energies are computed
      in \ref electron_binding() and approximated
      with
      \f[
      14.4381 \times 10^{-6} Z^{2.39} + 1.55468 \times 10^{-12} 
      Z^{5.35}~\mathrm{MeV}
      \f]
      as in Eq. A4 of \ref Lunney03 . 

      Generally, descendants of this class only need to provide an
      implementation of \ref mass_excess() and \ref mass_excess_d()
      and possibly a new version of \ref is_included() to be fully
      functional.

      \comment
      \future It might be useful to consider a fudge factor 
      to ensure no problems with finite precision arithmetic
      when converting \c double to \c int. 
      11/22/09 - I waffle back and forth on this. I think
      maybe its best to let the user deal with this their own
      way.
      \endcomment

  */
  class nuclear_mass : public nuclear_mass_info {

  public:

    nuclear_mass();
    
    virtual ~nuclear_mass() {};

    /// Return the type, \c "nuclear_mass".
    virtual const char *type() { return "nuclear_mass"; }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N) {
      return true;
    }
    
    /** \brief Fill \c n with the information from nucleus with the given
	neutron and proton number
	
	All masses are given in \f$\mathrm{fm}^{-1}\f$. The total mass
	(withouth the electrons) is put in part::m and part::ms, the
	binding energy is placed in nucleus::be, the mass excess in
	nucleus::mex and the degeneracy (part::g) is arbitrarily set
	to 1 for even A nuclei and 2 for odd A nuclei.
    */
    virtual int get_nucleus(int Z, int N, nucleus &n);
    
    /// Given \c Z and \c N, return the mass excess in MeV [abstract]
    virtual double mass_excess(int Z, int N)=0;

    /// Given \c Z and \c N, return the mass excess in MeV [abstract]
    virtual double mass_excess_d(double Z, double N)=0;

    /** \brief Return the approximate electron binding energy in MeV
     */
    virtual double electron_binding(double Z) {
      return (14.4381*pow(Z,2.39)+1.55468e-6*pow(Z,5.35))*1.0e-6;
    }
    
    /** \brief Return the binding energy in MeV
	
	The binding energy is defined to be negative for bound 
	nuclei, thus the binding energy per baryon of Pb-208
	is about -8*208 = -1664 MeV.
    */
    virtual double binding_energy(int Z, int N) {
      return (mass_excess(Z,N)+((Z+N)*m_amu-Z*m_elec-N*m_neut-Z*m_prot));
    }

    /** \brief Return the binding energy in MeV
	
	The binding energy is defined to be negative for bound 
	nuclei, thus the binding energy per baryon of Pb-208
	is about -8*208 = -1664 MeV.
    */
    virtual double binding_energy_d(double Z, double N) {
      return (mass_excess_d(Z,N)+((Z+N)*m_amu-Z*m_elec-N*m_neut-Z*m_prot));
    }
    
    /** \brief Return the total mass of the nucleus (without the
	electrons) in MeV
    */
    virtual double total_mass(int Z, int N) {
      return (mass_excess(Z,N)+((Z+N)*m_amu-Z*m_elec));
    }
    
    /** \brief Return the total mass of the nucleus (without the electrons) 
	in MeV
    */
    virtual double total_mass_d(double Z, double N) {
      return (mass_excess_d(Z,N)+((Z+N)*m_amu-Z*m_elec));
    }

    /** \brief Return the atomic mass of the nucleus in MeV
	(includes electrons and their binding energy)
    */
    virtual double atomic_mass(int Z, int N) {
      return total_mass(Z,N)+Z*m_elec-electron_binding(Z);
    }
    
    /** \brief Return the atomic mass of the nucleus in MeV
	(includes electrons and their binding energy)
    */
    virtual double atomic_mass_d(double Z, double N) {
      return total_mass_d(Z,N)+Z*m_elec-electron_binding(Z);
    }

    /// \name Base masses
    //@{
    /** \brief Neutron mass in \f$ \mathrm{MeV} \f$ 
	(defaults to o2scl_mks::mass_neutron times o2scl_const::hc_mev_fm)
    */
    double m_neut;
    
    /** \brief Proton mass in \f$ \mathrm{MeV} \f$ 
	(defaults to o2scl_mks::mass_proton times o2scl_const::hc_mev_fm)
    */
    double m_prot;
    
    /** \brief Electron mass in \f$ \mathrm{MeV} \f$ 
	(defaults to o2scl_mks::mass_electron times o2scl_const::hc_mev_fm)
    */
    double m_elec;

    /** \brief Atomic mass unit in \f$ \mathrm{MeV} \f$ 
	(defaults to o2scl_mks::mass_amu times o2scl_const::hc_mev_fm)
    */
    double m_amu;
    //@}

  };

  /** \brief Tabulated nuclear masses [abstract base]

      This uses simple linear interpolation to obtain masses of nuclei
      with non-integer value of Z and N.

      Generally, descendants of this class only need to provide an
      implementation of \ref mass_excess() and possibly a version
      of \ref nuclear_mass::is_included()

  */
  class nuclear_mass_table : public nuclear_mass {

  public:
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess_d(double Z, double N);

  };

  /** \brief Fittable mass formula [abstract base]

      Nuclear mass formulas which are descendants of this class
      can be fit to experiment using \ref mass_fit.

      Within \o2p, this class has only two children,
      \ref frdm_mass and \ref semi_empirical_mass. There
      is also a child \ref ldrop_mass in \o2e.
  */
  class nuclear_mass_fit : public nuclear_mass {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Return the type, \c "nuclear_mass_fit".
    virtual const char *type() { return "nuclear_mass_fit"; }

    /// Number of fitting parameters
    size_t nfit;

    /// Fix parameters from an array for fitting [abstract]
    virtual int fit_fun(size_t nv, const ubvector &x)=0;
    
    /// Fill array with guess from present values for fitting [abstract]
    virtual int guess_fun(size_t nv, ubvector &x)=0;

  };

  /** \brief Semi-empirical mass formula

      A simple semi-empirical mass formula of the form
      \f[
      E/A = B + S_s \frac{1}{A^{1/3}}+E_c \frac{Z^2}{A^{4/3}}
      + S_v \left(1-\frac{2Z}{A}\right)^2+E_{\mathrm{pair}}(Z,N)
      \f]
      where the pairing energy is given by
      \f[
      E_{\mathrm{pair}}(Z,N) = - \frac{E_{\mathrm{pair}}}{2 A^{3/2}}
      \left[ \cos \left( \pi Z \right)+\cos \left( \pi N \right) \right]
      \f]
      which is equivalent to the traditional prescription
      \f[
      E_{\mathrm{pair}}(Z,N) = \frac{E_{\mathrm{pair}}}{A^{3/2}}
      \times
      \left\{ 
      \begin{array}{rl}
      -1 & \mathrm{N~and~Z~even} \\
      +1 & \mathrm{N~and~Z~odd} \\
      0 & \mathrm{otherwise} 
      \end{array}
      \right.
      \f]
      when \f$ Z \f$ and \f$ N \f$ and integers.

      \note The default parameters are arbitrary, and are not
      determined from a fit.

      There is an example of the usage of this class given in 
      \ref ex_mass_fit_sect.
  */
  class semi_empirical_mass : public nuclear_mass_fit {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Binding energy (negative and in MeV, default -16)
    double B;

    /// Symmetry energy (in MeV, default 23.7)
    double Sv;

    /// Surface energy (in MeV, default 18)
    double Ss;

    /// Coulomb energy (in MeV, default 0.7)
    double Ec;

    /// Pairing energy (MeV, default 13.0)
    double Epair;
    
    /// Return the type, \c "semi_empirical_mass".
    virtual const char *type() { return "semi_empirical_mass"; }

    semi_empirical_mass();

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess_d(double Z, double N);
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      return mass_excess_d(Z,N);
    }

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);

    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);
    
  };

  /** \brief An approximation of shell effects in nuclei based on
      the interacting boson model

      Shell effects from \ref Dieperink09 based on the interacting
      boson model, with corrections as suggested by \ref Duflo95. 
  */
  class ibm_shell_energy {

  public:

    ibm_shell_energy();

    virtual ~ibm_shell_energy() {}

    /** \name Shell correction coefficients in MeV (defaults are 
	-1.39, 0.02, 0.03, and 0.075, respectively)
    */
    //@{
    double s_a1;
    double s_a2;
    double s_a3;
    double s_anp;
    //@}
    
    /// Number of magic numbers
    static const size_t nshells=11;
    
    /// Magic numbers
    int shells[nshells];
    
    /// Most recently computed shell energy
    double shell;

    /// Compute the shell energy for nucleus Z and N
    virtual double shell_energy(int Z, int N);

    /** \brief Compute the shell energy for specified values of Z and N 
	using bilinear interpolation
    */
    virtual double shell_energy_interp(double Z, double N);

  };

  /** \brief Nuclear mass formula from Dieperink and van Isacker (2009)
   */
  class dvi_mass : public nuclear_mass_fit, public ibm_shell_energy {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Volume energy coefficient
    double av;
    /// Surface energy coefficient
    double as;
    /// Symmetry energy coefficient
    double sv;
    /// Coulomb energy coefficient
    double ac;
    /// Pairing energy coefficient
    double ap;
    /// Surface symmetry energy coefficient
    double y;
    
    /// Return the type, \c "dvi_mass".
    virtual const char *type() { return "dvi_mass"; }

    dvi_mass();

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess_d(double Z, double N);
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      return mass_excess_d(Z,N);
    }

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);

    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);
    
  };
  
  /** \brief Compute the RMS radius of a Fermi-Dirac density distribution
      with fixed diffusiveness
      
      This class computes the RMS radius given either the central density
      or the radius specified in the Fermi function. This class assumes
      the density distribution function is of the form
      \f[
      N = 4 \pi \rho_0 \int r^2~dr~\left\{1+\exp
      \left[\left(r-R_{\mathrm{fermi}}\right)\right]\right\}^{-1}
      \f]
      where \f$ N \f$ is the total number of particles, and 
      \f$ \rho_0 \f$ is the central density.
      
      The radius assuming constant density,
      \f[
      R_{\mathrm{cd}} = \left(\frac{3 N}{4 \pi \rho_0}\right)^3 \, ,
      \f]
      is also given.
  */
  class rms_radius {

  protected:

    /// The central denstiy
    double urho0;
    /// The diffusiveness
    double ud;
    /** \brief Store the user-specified value of the radius in the
	Fermi distribution
	
	This is used in the integrands \ref iand() and \ref iand2().
    */
    double uRfermi;
    /// The total number of particles
    double uN;

    /// The integrator
    inte_qagiu_gsl<funct> it;
    /// The solver
    root_cern<funct> cr;
    
    /// The function \f$ 4 \pi r^4 \rho(r) \f$
    double iand(double r);

    /// The function \f$ 4 \pi r^2 \rho(r) \f$
    double iand2(double r);

    /// The function to fix the total number of particles
    double solve(double x);

  public:

    rms_radius();
    
    /** \brief Compute the RMS radius from the central density

	Computes the RMS radius \c Rrms from the central density \c
	rho0, the number of particles \c N, and the diffusiveness \c
	d.  This function also computes the radius in the Fermi
	distribution function, \c Rfermi and the radius assuming
	constant density, \c Rcd.
    */
    int eval_rms_rho(double rho0, double N, double d,
		     double &Rcd, double &Rfermi, double &Rrms);
    
    /** \brief Compute the RMS radius from the Fermi distribution radius

	Computes the RMS radius \c Rrms from the radius \c Rfermi in
	the Fermi distribution assuming a total number of particles \c
	N, a diffusiveness paramter \c d. This function also produces
	the central density \c rho0, and the radius assuming constant
	density, \c Rcd.
    */
    int eval_rms_rsq(double Rfermi, double N, double d,
		     double &rho0, double &Rcd, double &Rrms);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
