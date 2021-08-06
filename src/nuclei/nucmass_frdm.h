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
#ifndef NUCMASS_FRDM_H
#define NUCMASS_FRDM_H

/** \file nucmass_frdm.h \brief File defining \ref
    o2scl::nucmass_frdm, \ref o2scl::nucmass_mnmsk, and \ref
    o2scl::nucmass_mnmsk_exp
*/

#include <cmath>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/constants.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief FRDM semi-empirical mass formula (macroscopic part only
      with no deformation)

      \verbatim embed:rst
      The spherically-symmetric, macroscopic part of the finite-range
      droplet model from [Moller95]_.
      \endverbatim
  
      Using the relations
      \f[
      \bar{\delta} = (n_n - n_p)/n
      \f]
      and
      \f[
      \bar{\epsilon} = - (n-n_0)/3/n_0
      \f]
      we get
      \f[
      n_n = \frac{1}{2} (1+\bar{\delta}) (1-3 \bar{\epsilon}) n_0
      \f]
      and
      \f[
      n_p = \frac{1}{2} (1-\bar{\delta}) (1-3 \bar{\epsilon}) n_0
      \f]
      Assuming that 
      \f[
      \frac{4 \pi}{3} R_n^3 n_n = N
      \f]
      and 
      \f[
      \frac{4 \pi}{3} R_p^3 n_p = Z
      \f]
      we get
      \f[
      R_n^3 = 3 N / \alpha_n
      \f]
      \f[
      R_p^3 = 3 Z / \alpha_p
      \f]
      where \f$ \alpha \f$'s are
      \f[
      \alpha_n = 2 \pi (1+ \bar{\delta})(1 - 3 \bar{\epsilon}) n_0
      \f]
      \f[
      \alpha_p = 2 \pi (1- \bar{\delta})(1 - 3 \bar{\epsilon}) n_0
      \f]
      Note that the above relations are somehow self-consistent
      because they imply
      \f[
      R^3 n = R_n^3 n_n + R_p^3 n_p
      \f]

      Since we're using (is there a better way?)
      \f[
      R = r_0 A^{1/3}
      \f]
      with \f$ r_0 = 1.16 \f$ fm, then 
      \f$ n_0 = 0.152946 \mathrm{fm}^{-3} \f$.

      \verbatim embed:rst
      .. todo:: 

         In class nucmass_frdm:

         - Fix pairing energy and double vs. int
         - Document drip_binding_energy(), etc.
         - Decide on number of fit parameters (10 or 12?) or
           let the user decide
         - Document the protected variables
         - Set the neutron and proton masses and hbarc to Moller et al.'s 
           values

      \endverbatim

      \future Add microscopic part.
  */
  class nucmass_frdm : public nucmass_fit_base {

  public:

    nucmass_frdm();

    /// Volume-energy constant in MeV (default 16.247)
    double a1;
    /// Symmetry-energy constant in MeV (default 32.73)
    double J;
    /// Nuclear compressibility constant in MeV (default 240)
    double K;
    /// Surface-energy constant in MeV (default 22.92)
    double a2;
    /// Effective surface-stiffness constant in MeV (default 29.21)
    double Q;
    /// Curvature-energy constant in MeV (default 0)
    double a3;
    /// Charge-asymmetry constant in MeV (default 0.436)
    double ca;
    /// Wigner constant in MeV (default 30)
    double W;
    /** \brief electronic-binding constant in MeV 
        (default \f$ 1.433 \times 10^{-5} \f$ ).
    */
    double ael;
    /// Proton root-mean-square radius in fm (default 0.80)
    double rp;
    /// Nuclear-radius constant in fm (default 1.16)
    double r0;
    /// Hydrogen atom mass excess, 7.289034 MeV
    double MH;
    /// Neutron mass excess, 8.071431 MeV
    double Mn;
    /// Electronic charge squared, 1.4399764 MeV fm
    double e2;
    /// Range of Yukawa-plus-exponential potential, 0.68 fm
    double a;
    /** \brief Range of Yukawa function used to generate nuclear 
        charge distribution, 0.70 fm
    */
    double aden;
    /// Average pairing-gap constant, 4.80 MeV
    double rmac;
    /// Neutron-proton interaction constant, 6.6 MeV
    double h;
    /// Density-symmetry constant, 0 MeV
    double L;
    /// Pre-exponential compressibility-term constant, 60 MeV
    double C;
    /// Exponential compressibility-term range constant, 0.831
    double gamma;
    /// Atomic mass unit, 931.5014 MeV
    double amu;
    
    /// Internal average neutron density
    double nn;
    
    /// Internal average proton density
    double np;

    /// Neutron radius
    double Rn;

    /// Proton radius
    double Rp;

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

    /** \brief Return the binding energy in MeV
     */
    virtual double drip_binding_energy_d(double Z, double N,
                                         double npout, double nnout,
                                         double chi);

    /** \brief Given \c Z and \c N, return the mass excess in MeV
        in a many-body environment

        This is an experimental version of mass_excess_d which removes
        pairing, computes nn, np, Rn, and Rp, and attempts to correct
        the surface. This function probably doesn't work at the
        moment. It's not currently used by \ref
        drip_binding_energy_d().
    */
    virtual double drip_mass_excess_d(double Z, double N,
                                      double np_out, double nn_out,
                                      double chi);

  protected:

    /// Conversion from kg to inverse fm
    double kg_to_invfm;

    /// Proton pairing coefficient
    double Deltap;
    /// Neutron pairing coefficient
    double Deltan;
    /// Isubvector pairing coefficient
    double deltanp;

    /// Average bulk nuclear asymmetry
    double deltabar;
    /// Average relative deviation of bulk density
    double epsbar;

    /// Desc
    double Bs;
    /// Desc
    double Bk;
    /// Desc
    double Br;
    /// Desc
    double Bw;
    /// Desc
    double Bv;

    /// Coulomb energy coefficient
    double c1;
    /// Volume redistribution energy coefficient
    double c2;
    /// Coulomb exchange correction coefficient
    double c4;
    /// Surface redistribution energy coefficient
    double c5;

    /** \brief Coefficient for the proton form-factor correction to the 
        Coulomb energy
    */
    double f0;

    /// Desc
    double a0;

    /// Desc
    double B1;
    /// Desc
    double B2;
    /// Desc
    double B3;
    /// Desc
    double B4;

  };
  
  /** \brief Nuclear masses from Moller, et al.

      \verbatim embed:rst
      This is based on the tables given in [Moller95]_, 
      [Moller97]_, and [Moller16ng]_.
      \endverbatim

      In order to allow easier coordination of file I/O across
      multiple MPI tasks the constructor does not automatically load
      nuclear mass data. To load data from the \o2 HDF5 data files,
      use <tt>o2scl_hdf::mnmsk_load()</tt>. If no data is loaded, 
      then \ref o2scl::nucmass_table::is_loaded() will return
      <tt>false</tt> and calls to get_ZN() will call the error 
      handler.

      There are several entries in the original table which are
      blank because they are in some way not known, measured, or
      computable. These entries are filled with a positive number 
      larger than 1.0e90, given by the functions \ref blank(),
      \ref neither(), \ref beta_stable(), \ref beta_plus_and_minus(),
      \ref greater_100(), or \ref very_large() .

      \note This class requires data stored in an HDF file and
      thus requires HDF support for normal usage.
  */
  class nucmass_mnmsk : public nucmass_table {

  public:
    
    nucmass_mnmsk();

    virtual ~nucmass_mnmsk();
    
    /** \brief Entry structure for Moller, et al. masses
     */
    struct entry {
      
      /// Neutron number
      int N;
    
      /// Proton number
      int Z;
    
      /// Atomic number
      int A;
    
      /** \name Ground state deformations (perturbed-spheroid parameterization)
       */
      //@{
      /// Quadrupole
      double eps2;
      /// Octupole
      double eps3;
      /// Hexadecapole
      double eps4;
      /// Hexacontatetrapole
      double eps6;
      /// Hexacontatetrapole without mass asymmetry
      double eps6sym;
      //@}
    
      /** \name Ground state deformations in the spherical-harmonics expansion
       */
      //@{
      /// Quadrupole
      double beta2;
      /// Octupole
      double beta3;
      /// Hexadecapole
      double beta4;
      /// Hexacontatetrapole
      double beta6;
      //@}
    
      /// The ground-state microscopic energy
      double Emic;
    
      /// The theoretical mass excess (in MeV)
      double Mth;
    
      /// The experimental mass excess (in MeV)
      double Mexp;
    
      /// Experimental mass excess error
      double sigmaexp;
    
      /// The ground-state microscopic energy in the FRLDM
      double EmicFL;
    
      /// The theoretical mass excess in the FRLDM
      double MthFL;
    
      /// Spin and pairity of odd proton 
      char spinp[6];
    
      /// Spin and pairity of odd neutron
      char spinn[6];

      /// Lipkin-Nogami proton gap
      double gapp;
    
      /// Lipkin-Nogami neutron gap
      double gapn;

      /// Total binding energy
      double be;

      /// One neutron separation energy
      double S1n;
    
      /// Two neutron separation energy
      double S2n;

      /** \brief Percentage of daughters generated in beta decay after
          beta-delayed neutron emission
      */
      double PA;

      /// Desc
      double PAm1;

      /// Desc
      double PAm2;

      /// Energy released in beta-decay
      double Qbeta;

      /// Half-life w.r.t. GT beta-decay
      double Tbeta;

      /// One proton separation energy
      double S1p;

      /// Two proton separation energy
      double S2p;

      /// Energy released in alpha-decay
      double Qalpha;

      /// Half-life w.r.t. alpha-decay
      double Talpha;
    
    };
  
    /** \brief Return false if the mass formula does not include 
        specified nucleus
    */
    virtual bool is_included(int Z, int N);
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
    /** \brief Get the entry for the specified proton and neutron number
        
        This method searches the table using a cached binary search
        algorithm. It is assumed that the table is sorted first by
        proton number and then by neutron number.
    */
    nucmass_mnmsk::entry get_ZN(int l_Z, int l_N);
    
    /// The value which corresponds to a blank entry
    double blank() { return 1.0e99; };

    /// Neither beta+ or beta- is possible
    double neither() { return 1.0e98; };

    /// The value which corresponds to a blank entry
    double beta_stable() { return 1.0e97; };
    
    /// Both beta+ and beta- are possible
    double beta_plus_and_minus() { return 1.0e96; };
    
    /// The value is greater than 100
    double greater_100() { return 1.0e95; };

    /// The value is greater than \f$ 10^{20} \f$
    double very_large() { return 1.0e94; };

    /// Return the type, \c "nucmass_mnmsk".
    virtual const char *type() { return "nucmass_mnmsk"; }
    
    /** \brief Set data

        This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass, nucmass_mnmsk::entry *m, std::string ref);

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The reference for the original data
    std::string reference;
    
    /// The array containing the mass data of length ame::n
    nucmass_mnmsk::entry *mass;
    
    /// The last table index for caching
    int last;
    
#endif
    
  };

  /** \brief The experimental values from the Moller et al.
      mass tables

      \verbatim embed:rst 
      This mass formula only includes the experimental mass excesses
      tabulated in [Moller95]_, [Moller97]_, and [Moller16ng]_.
      \endverbatim
      
      \note This class requires data stored in an HDF file and
      thus requires HDF support for normal usage.
  */
  class nucmass_mnmsk_exp : public nucmass_mnmsk {

  public:

    /** \brief Return false if the mass formula does not include 
        specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);

  };
  
  /** \brief Nuclear masses from two tables patched together
      
      Experimental.
   */
  class nucmass_patch : public nucmass_table {
    
  protected:
    
    nucmass_ame ame;
    
    nucmass_table *nt;
    
    nucmass_fit_base *nf;

    bool inc_fit;
    
  public:

    nucmass_patch() {
      nt=0;
      nf=0;
      nt=&def_table;
      nf=&def_fit;
    }
    
    void load(bool include_fit=true);
    
    nucmass_mnmsk def_table;
    
    nucmass_frdm def_fit;
    
    /// Return the type, \c "nucmass_patch".
    virtual const char *type() { return "nucmass_patch"; }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N) {
      if (ame.is_included(Z,N) || nt->is_included(Z,N) ||
          nf->is_included(Z,N)) {
        return true;
      }
      return false;
    }

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      if (ame.is_included(Z,N)) {
        return ame.mass_excess(Z,N);
      }
      if (nt->is_included(Z,N)) {
        return nt->mass_excess(Z,N);
      }
      if (inc_fit && nf->is_included(Z,N)) {
        return nf->mass_excess(Z,N);
      }
      O2SCL_ERR("Failed to find nucleus in nucmass_patch::mass_excess().",
                o2scl::exc_einval);
      return 0.0;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

