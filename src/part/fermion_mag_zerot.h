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
#ifndef O2SCL_MAG_FERMION_ZEROT_H
#define O2SCL_MAG_FERMION_ZEROT_H

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/fermion.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/root_brent_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Thermodynamics of fermions in a magnetic field at 
      zero temperature

      Using the notation of \ref Broderick00, the
      effective mass of a fermion in a magnetic field is 
      \f[
      \tilde{m}^2 = m^2 + 2 \left( n + \frac{1}{2} - 
      \frac{q \sigma}{2 |q|} \right) |q| B
      \f]
      where \f$ m \f$ is the bare mass, \f$ q \f$ is the charge,
      and \f$ \sigma \f$ is the z-component of the spin along
      the axis of the magnetic field (for spin 1/2 fermions, for example,
      either +1 or -1). 

      \note This only works for spin 1/2 particles at the moment,
      so the discussion below assumes this is the case.

      \note The function calc_density_zerot_mag() 
      will fail if the density is small enough. 

      The Fermi momentum and the chemical potential are related
      in the usual way
      \f[
      k_{F,\sigma}^2 = \mu^2 - \tilde{m}_{\sigma}^2
      \f]

      The density is then given by 
      \f[
      n = \frac{|q| B}{2 \pi^2} \sum_{\sigma} 
      \sum_{n=0}^{n_{\mathrm{max},\sigma}} k_{F,\sigma}
      \f]
      where \f$ n_{\mathrm{max},\sigma} \f$ is the integer value of
      \f$ n \f$ for which the next largest integer makes \f$ \mu <
      \tilde{m} \f$. The value of \f$ n_{\mathrm{max}} \f$ is stored
      internally in this class as an integer (\ref nmax_dn and \ref
      nmax_up), and if there are no integers for which \f$ \mu <
      \tilde{m} \f$, then the corresponding integer will be -1, to
      indicate that no terms in the sum are present. For any fermion,
      at least one Landau level is always filled for at least one of
      the spin states at every density.
      
      When the number of Landau levels for \f$ \sigma=1 \f$ is larger
      than \ref sum_limit, the B=0 result is given.

      <b>Units:</b>

      It is useful to think of the magnetic field as typically being
      multiplied by the electron charge, so if magnetic field 
      is measured in Gauss, then for a 1 Gauss field in 
      "Gaussian" units, 
      \f[
      e B = 1.5~\times~10^{-19}~\mathrm{fm}^{-2}
      \f]
      This conversion factor is given in \ref o2scl_const::ec_gauss_fm2 .

      \todo Comment on Gaussian vs. Heaviside-Lorentz units.

  */
  class fermion_mag_zerot : public fermion_zerot {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:

    /// Solver to compute chemical potential from density
    mroot<mm_funct11,boost::numeric::ublas::vector<double>, 
      jac_funct<> > *density_root;
    
    /// The charge times the magnetic field in \f$ \mathrm{fm}^{-2} \f$
    double qBt;
    
    /// The anomalous magnetic moment
    double kt;

    /// The target density, set by calc_density_zerot_mag()
    double dent;
    
    /// Function to compute chemical potential from density
    int solve_fun(size_t nv, const ubvector &x,
		  ubvector &y);

    /// Pointer to the data object
    fermion *fp;

  public:

    /// Maximum Landau level for spin up particles (default 0)
    int nmax_up;

    /// Maximum Landau level for spin down particles (default 0)
    int nmax_dn;

    /// Limit on the sum (default \f$ 10^{6} \f$)
    int sum_limit;
    
    /// Create a fermion with mass \c mass and degeneracy \c dof.
    fermion_mag_zerot() {
      nmax_up=0;
      nmax_dn=0;
      sum_limit=1000000;
      density_root=&def_density_root;
      def_density_root.tol_abs/=1.0e4;
    }
    
    virtual ~fermion_mag_zerot() {
    }
    
    /**	\brief Thermodynamics in a magnetic field using the chemical 
	potential
	
	The parameter \c qB is the charge (in units of the positron
	charge) times the magnetic field strength. Thus, for example,
	\c qB should be negative for electrons.
    */
    virtual void calc_mu_zerot_mag(fermion &f, double qB, double kappa=0.0);

    /** \brief Thermodynamics in a magnetic field using the density
    */
    virtual void calc_density_zerot_mag(fermion &f, double qB, 
					double kappa=0.0);

    /** \brief Set the solver for use in calculating the chemical
        potential from the density 
    */
    int set_density_root
      (mroot<mm_funct11,boost::numeric::ublas::vector<double>,
       jac_funct<> > &rp) {
      density_root=&rp;
      return 0;
    }

    /** \brief The default solver for calc_density().
     */
    mroot_hybrids<mm_funct11,
      boost::numeric::ublas::vector<double>, 
      boost::numeric::ublas::matrix<double>,
      jac_funct<> > def_density_root;

    /// Return string denoting type ("fermion_mag_zerot")
    virtual const char *type() { return "fermion_mag_zerot"; }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
