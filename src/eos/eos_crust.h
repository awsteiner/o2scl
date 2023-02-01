/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
/** \file eos_crust.h
    \brief File defining \ref o2scl::eos_crust
*/
#ifndef O2SCL_BPS_EOS_H
#define O2SCL_BPS_EOS_H

#include <o2scl/part.h>
#include <o2scl/mroot.h>
#include <o2scl/eos_base.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/fermion.h>
#include <o2scl/nucmass.h>

namespace o2scl {

  /** \brief Baym-Pethick-Sutherland equation of state
    
      This calculates the equation of state of electrons and nuclei
      using the approach of Baym et al. (1971; based on the discussion
      in Shapiro and Teukolsky (1983)) between about \f$8 \times
      10^{6} ~\mathrm{g}/\mathrm{cm}^3\f$ and \f$4.3 \times 10^{11}
      ~\mathrm{g}/\mathrm{cm}^3\f$. Below these densities, more
      complex Coulomb corrections need to be considered, and above
      these densities, neutron drip is important.

      The default mass formula is semi-empirical
      \f{eqnarray*}
      M(A,Z)&=&(A-Z) m_n+Z (m_p+m_e)-
      15.76 A-17.81 A^{2/3} \\
      && -0.71 Z^2 /A^{1/3}-
      94.8/A \left(A/2-Z\right)^2+E_{\mathrm{pair}}
      \f}
      where 
      \f[
      E_{\mathrm{pair}} = \pm 39/A^{3/4}
      \f]
      if the nucleus is odd-odd (plus sign) or even-even (minus sign)
      and \f$E_{\mathrm{pair}}\f$ is zero for odd-even and even-odd
      nuclei. The nuclei are assumed not to contribute to the
      pressure. The electronic contribution to the pressure is assumed
      to be equal to the Fermi gas contribution plus a "lattice"
      contribution
      \f[
      \varepsilon_L = -1.444 Z^{2/3} e^2 n_e^{4/3}
      \f]
      This is Eq. 2.7.2 in Shapiro and Teukolsky (1983). The rest mass
      energy of the nucleons is included in the energy density.

      \verbatim embed:rst
      See [Shapiro83]_.
      \endverbatim
      
      The original results from Baym et al. (1971) are stored as a
      \ref table in file <tt>data/o2scl/bps.eos</tt>. The testing code
      for this class compares the calculations to the table and
      matches to within .2 percent for the energy density and 9
      percent for the pressure (for a fixed baryon number density).

      \future Can the pressure be made to match more closely?
      \future Convert to a \ref o2scl::eos_had_base object and offer an
      associated interface?

  */
  class eos_crust : public eos_base {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    eos_crust();
    
    virtual ~eos_crust() {};
    
    /** \brief Calculate the equation of state as a function of
	the baryon number density \c barn
	
	This calculates the equation of state as a function of the
	baryon number density in \f$\mathrm{fm}^{-3}\f$, returning the
	representative nucleus with proton number \c Z and atomic
	number \c A.  The pressure and energy density are returned 
	in \c th in \f$\mathrm{fm}^{-4}\f$.
    */
    virtual int calc_density(double barn, thermo &th, int &Z, int &A);

    /** \brief Calculate the equation of state as a function of
	the pressure

	This calculates the equation of state as a function of the
	pressure, returning the representative nucleus with proton
	number \c Z and atomic number \c A and the baryon number
	density \c barn in \f$\mathrm{fm}^{-3}\f$. The energy density
	is also returned in \f$\mathrm{fm}^{-4}\f$ in \c th.
    */
    virtual int calc_pressure(thermo &th, double &barn, int &Z, int &A);

    /** \brief The electron lattice energy */
    virtual double lattice_energy(int Z);

    /** \brief Get a pointer to the electron
     */
    virtual const fermion &get_electron() { return e; }

    /** \brief The mass formula

	The nuclear mass without the contribution of the rest mass
	of the electrons. The electron rest mass energy is included
	in the electron thermodynamics elsewhere.
     */
    virtual double mass_formula(int Z, int A);
    
    /// Return string denoting type ("eos_crust")
    virtual const char *type() { return "eos_crust"; }
    
    /// Default mass formula
    nucmass_semi_empirical def_mass;

    /// Set the nuclear mass formula to be used
    int set_mass_formula(nucmass &nm) {
      nmp=&nm;
      return 0;
    }

    /// Compute the ground state assuming a fixed atomic number
    int calc_density_fixedA(double barn, thermo &th, int &Z, int A);
    
    /** \brief The electron thermodynamics
	
	\note The electron rest mass is included by default in 
	the energy density and the chemical potential
    */
    fermion e;

  protected:

#ifndef DOXYGEN_INTERNAL

    /// Zero-temperature thermodynamics for the electrons
    fermion_zerot fzt;

    /// Solve Equation 2.7.4 for a given pressure
    virtual int eq274(size_t nv, const ubvector &nx, ubvector &ny, 
		      int Zt);

    /// The Gibbs free energy
    double gibbs(int Z, int A);

    /// The energy density
    double energy(double barn, int Z, int A);

    /// A solver to solve Eq. 2.7.4
    mroot_hybrids<mm_funct> gs;

    /// The nuclear mass formula
    nucmass *nmp;
    
#endif

  };

}

#endif
