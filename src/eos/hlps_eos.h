#ifndef MODELS_HLPS_H
#define MODELS_HLPS_H

#include <iostream>

#include <mpi.h>

#include <o2scl/poly.h>

#include "misc.h"
#include "entry.h"
#include "models.h"

namespace o2scl {

  /** \brief Schematic EOS from Hebeler et al.

      From an energy per baryon of
      \f[
      E/A = \left( 3 \pi^2 n_0/2 \right)^{2/3} \frac{1}{2 M}
      \left\{ \frac{3}{5} \left[ x^{5/3} + (1-x)^5/3 \right] (2 u)^{2/3}-
      [(2 \alpha - 4 \alpha_L) x (1-x)+\alpha_L] +
      \left[ (2 \eta - 4 \eta_L) x (1-x) + \eta_L \right] u^{\gamma} \right\}
      \f]

      One can fix the values of \f$ \alpha, \eta, \f$ and \f$ \gamma \f$
      by the requirement that the pressure is zero at saturation and
      by fixing the binding energy and incompressibility.

      Note that the original reference appears to have a typo in 
      the pressure. The \f$ 2/5 \f$ factor in front should be 
      \f$ 1/5 \f$ .
  */
  class hlps_eos : public o2scl::hadronic_eos_eden {

  protected:

    /// To solve quadratic equation for 'gamma'
    o2scl::quadratic_real_coeff_gsl quad;

  public:

    /// \name Constants (all unitless)
    //@{
    double gamma;
    double alpha;
    double eta;
    double alphaL;
    double etaL;
    //@}

    hlps_eos();

    virtual ~hlps_eos() {};

    /** \brief Fix 'alpha', 'eta' and 'gamma' from saturation properties

	All inputs must be in \f$ \mathrm{fm}^{-1} \f$.
    */
    void fix_coeffs(double M, double B, double K);

    /** \brief Fix 'alphaL' and 'etaL' from neutron matter EOS and its
	derivative

	The parameters \c M and \c Eneut must be in
	\f$ \mathrm{fm}^{-1} \f$ and \c dEneut must be in 
	\f$ \mathrm{fm}^{2} \f$
    */
    void fix_neutron_matter(double M, double Eneut, double dEneut);

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(o2scl::fermion &ln, o2scl::fermion &lp, 
		       o2scl::thermo &lth);

    /// Return string denoting type ("hlps_eos")
    virtual const char *type() { return "hlps_eos"; }

  };

}

#endif
