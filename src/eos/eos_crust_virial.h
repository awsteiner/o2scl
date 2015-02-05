/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2015, Andrew W. Steiner
  
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
/** \file eos_crust_virial.h
    \brief File defining \ref o2scl::eos_crust_virial
*/
#ifndef O2SCL_VIRIAL_EOS_H
#define O2SCL_VIRIAL_EOS_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/lib_settings.h>
#include <o2scl/interp.h>
#include <o2scl/boson.h>
#include <o2scl/eos_had_base.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Virial EOS for neutrons, protons, deuterons, and alpha 
      particles
      
      Virial EOS from \ref Horowitz06 and \ref Horowitz06b .
      
      \warning This class is implemented as a eos_had_base object
      because it might be helpful to be able to use \ref
      o2scl::eos_had_temp_base::calc_temp_e(), but because of the
      alpha particles and deuterons, some of the other \ref
      o2scl::eos_had_base methods don't have the correct
      interpretation.
  */
  class eos_crust_virial : public eos_had_temp_pres_base {

  protected:

    /// \name Interpolation for virial coefficients
    //@{
    std::vector<double> Tv, bnv, Tbnpv, bpnv, Tbpnpv;
    std::vector<double> banv, Tbanpv, bav, Tbapv;
    interp_vec<std::vector<double> > ibn, iTbnp, ibpn, iTbpnp, 
      iban, iTbanp, iba, iTbap;
    //@}

  public:

    eos_crust_virial();

    virtual ~eos_crust_virial() {
    }

    /// Internal alpha particle
    boson alpha;

    /// Internal deuteron
    boson deuteron;

    /** \brief Equation of state as a function of the chemical potentials
    */
    virtual int calc_p(fermion &ne, fermion &pr, thermo &th) {
      O2SCL_ERR("Virial EOS does not work at T=0",exc_efailed);
      return 0;
    }

    /** \name The virial coefficients and their temperature derivatives

	These functions assume that the temperature is specified in MeV
     */
    //@{
    virtual double bn(double T);
    virtual double ban(double T);
    virtual double ba(double T);
    virtual double bpn(double T);
    virtual double Tbn_prime(double T);
    virtual double Tban_prime(double T);
    virtual double Tba_prime(double T);
    virtual double Tbpn_prime(double T);
    //@}

    /** \brief Equation of state as a function of the chemical potentials
	at finite temperature
    */
    virtual int calc_temp_p(fermion &n, fermion &p, double T, 
			    thermo &th);

    /** \brief Equation of state as a function of the chemical
	potentials at finite temperature with alpha particles and
	deuterons
    */
    virtual int calc_temp_p_alpha(fermion &n, fermion &p, boson &d, boson &a, 
				  double T, thermo &th);

    /// Desc
    void fit();
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
