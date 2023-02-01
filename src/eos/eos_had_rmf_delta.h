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
/** \file eos_had_rmf_delta.h
    \brief File defining \ref o2scl::eos_had_rmf_delta
*/
#ifndef O2SCL_RMF_DELTA_EOS_H
#define O2SCL_RMF_DELTA_EOS_H

#include <o2scl/eos_had_rmf.h>

namespace o2scl {
  
  /** \brief Field-theoretical EOS with scalar-isovector meson, 
      \f$ \delta \f$.

      \verbatim embed:rst
      See also [Kubis97]_ and [Gaitanos04]_.
      \endverbatim

      This essentially follows the notation in Kubis et al. (1997),
      except that our definitions of \c b and \c c follow their \f$
      \bar{b} \f$ and \f$ \bar{c} \f$, respectively.

      Also discussed in Gaitanos et al. (2004), where they take \f$
      m_{\delta}=980 \f$ MeV.

      The full Lagragian is:
      \f[
      {\cal L} = {\cal L}_{Dirac} + {\cal L}_{\sigma} + 
      {\cal L}_{\omega} + {\cal L}_{\rho} + {\cal L}_{\delta}
      \f]
    
      \f{eqnarray*}
      {\cal L}_{Dirac} &=& 
      \bar{\Psi} \left[ i {{\partial}\!\!\!{\backslash}} - 
      g_{\omega} {{\omega}\!\!\!{\backslash}} - \frac{g_{\rho}}{2} 
      {{\vec{\rho}}\!\!\!{\backslash}}~
      \vec{\tau} - M + g_{\sigma} \sigma - \frac{e}{2} 
      \left( 1 + \tau_3 \right) A_{\mu} \right] \Psi \nonumber \\
      {\cal L}_{\sigma} &=& 
      {\textstyle \frac{1}{2}} \left( \partial_{\mu} \sigma \right)^2 
      - {\textstyle \frac{1}{2}} m^2_{\sigma} \sigma^2 
      - \frac{b M}{3} \left( g_{\sigma} \sigma\right)^3 
      - \frac{c}{4} \left( g_{\sigma} \sigma\right)^4  \nonumber \\
      {\cal L}_{\omega} &=& 
      - {\textstyle \frac{1}{4}} f_{\mu \nu} f^{\mu \nu} 
      + {\textstyle \frac{1}{2}} m^2_{\omega}\omega^{\mu}\omega_{\mu} 
      + \frac{\zeta}{24} g_{\omega}^4 \left(\omega^\mu \omega_\mu\right)^2
      \nonumber \\
      {\cal L}_{\rho} &=& 
      - {\textstyle \frac{1}{4}} \vec{B}_{\mu \nu} \cdot \vec{B}^{\mu \nu}
      + {\textstyle \frac{1}{2}} m^2_{\rho} \vec{\rho}^{~\mu} \cdot 
      \vec{\rho}_{~\mu} 
      + \frac{\xi}{24} g_{\rho}^4 \left(\vec{\rho}^{~\mu}\right) \cdot 
      \vec{\rho}_{~\mu} 
      + g_{\rho}^2 f (\sigma, \omega) \vec{\rho}^{~\mu} \cdot 
      \vec{\rho}_{~\mu} \nonumber \\
      \f}
      where the additional terms are

      \f[
      {\cal L}_{\delta} = \bar{\Psi} \left( g_{\delta} \vec{\delta} \cdot 
      \vec{\tau} \right) \Psi 
      + \frac{1}{2} (\partial_{\mu} \vec{\delta})^2 - 
      \frac{1}{2} m_{\delta}^2 \vec{\delta}^{~2}
      \f]

      The new field equation for the delta meson is
      \f[
      m_{\delta}^2 \delta = g_{\delta} (n_{s,p} - n_{s,n})
      \f]

      \future Finish the finite temperature EOS 

   */
  class eos_had_rmf_delta : public eos_had_rmf {
  public:

    /// The mass of the scalar-isovector field
    double md;

    /// The coupling of the scalar-isovector field to the nucleons
    double cd;

    /// The value of the scalar-isovector field
    double del;
    
    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lth);

    /** \brief Equation of state as a function of chemical potentials
    */
    virtual int calc_eqd_p(fermion &neu, fermion &p, 
                           double sig, double ome, double rho, double delta,
                           double &f1, double &f2, double &f3, double &f4,
                           thermo& th);
    
    /** \brief Finite temperature (unfinished)
     */
    int calc_temp_eqd_p(fermion &ne, fermion &pr, double temper,
                        double sig, double ome, double lrho, 
                        double delta, double &f1, double &f2, 
                        double &f3, double &f4, thermo& lth);
      
    /** \brief Set a guess for the fields for the next call to calc_e(), 
	calc_p(), or saturation()
    */
    virtual int set_fields(double sig, double ome, double lrho, 
			   double delta) {
      sigma=sig;
      omega=ome;
      rho=lrho;
      del=delta;
      guess_set=true;
      return 0;
    }
    
    /** \brief Calculate saturation properties for nuclear matter 
	at the saturation density
	
	This requires initial guesses to the chemical 
	potentials, etc.
    */
    virtual int saturation();

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The function for calc_e()
    virtual int calc_e_solve_fun(size_t nv, const ubvector &ex, 
				 ubvector &ey);
    
    /// Compute matter at zero pressure (for saturation())
    virtual int zero_pressure(size_t nv, const ubvector &ex, 
			      ubvector &ey);
    
    
  private:
    
    /** \brief Forbid setting the guesses to the fields unless all four
	fields are specified
    */
    virtual int set_fields(double sig, double ome, double lrho) {
      return 0;
    }

#endif
    
  };
  
}
  
#endif
