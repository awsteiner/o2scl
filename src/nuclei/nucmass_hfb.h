/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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
#ifndef HFB_MASS_H
#define HFB_MASS_H

/** \file nucmass_hfb.h
    \brief File defining \ref o2scl::nucmass_hfb
*/

#include <cmath>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/constants.h>

namespace o2scl {
    
  /** \brief HFB Mass formula 

      \verbatim embed:rst
      .. todo:: 

         In class nucmass_hfb:

         - Mg40 is present in some tables but not others. Compare
           hfb14-plain with hfb14-plain_v0. This may be related to the
           fact that the mass excess of Mg40 differs significantly between
           the 2003 and 2013 Audi et al. tables?
         - Update to include hfb17. 

      \endverbatim
  */
  class nucmass_hfb : public nucmass_table {
    
  public:

    /** \brief Entry structure for HFB mass formula
     */
    struct entry {
    
      /// Neutron number
      int N;
    
      /// Proton number
      int Z;
    
      /// Atomic number
      int A;
    
      /// Beta 2 deformation
      double bet2;

      /// Beta 4 deformation
      double bet4;

      /// RMS charge radius
      double Rch;

      /// Deformation and Wigner energies
      double def_wig;

      /// Neutron separation energy
      double Sn;

      /// Proton separation energy
      double Sp;

      /// Beta-decay energy
      double Qbet;

      /// Calculated mass excess
      double Mcal;

      /// Error between experimental and calculated mass excess
      double Err;
    
    };
    
    /** \brief Create a new mass formula object 
     */
    nucmass_hfb();

    virtual ~nucmass_hfb();

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
    nucmass_hfb::entry get_ZN(int l_Z, int l_N);
    
    /// The value which corresponds to a blank entry
    double blank() { return 1.0e99; };

    /// Return the type, \c "nucmass_hfb".
    virtual const char *type() const { return "nucmass_hfb"; }

    /** \brief Set data
        
        This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass,
                 std::vector<nucmass_hfb::entry> &m, std::string ref);

    /// Desc
    virtual void clear() {
      mass.clear();
      this->n=0;
      last=0;
      return;
    }
    
  protected:
    
    /// The array containing the mass data of length ame::n
    std::vector<nucmass_hfb::entry> mass;
    
    /// The last table index for caching
    int last;
    
  };

  /** \brief HFB Mass formula with spin and parity information

      This class stores the HFB17 through HFB27 tables, loaded by
      \ref o2scl_hdf::hfb_sp_load() , and also the newer
      Brussels-Skyrme-on-a-Grid (BSkG1 through BSkG4) tables,
      loaded by \ref o2scl_hdf::bskg_load() . The BSkG tables
      provide several additional deformation, separation-energy,
      odd-even staggering, rotational-correction, pairing-gap, and
      moment-of-inertia fields beyond what the earlier HFB tables
      provide; see \ref nucmass_hfb_sp::entry for details. Not all
      of these additional fields are provided by every BSkG model:
      \ref nucmass_hfb_sp::entry::S2n through \ref
      nucmass_hfb_sp::entry::I3 are only filled in for BSkG3, while
      \ref nucmass_hfb_sp::entry::Erot through \ref
      nucmass_hfb_sp::entry::par_n are only filled in for BSkG1,
      BSkG2, and BSkG4.
  */
  class nucmass_hfb_sp : public nucmass_table {
    
  public:
    
    /** \brief Create a new mass formula object
     */
    nucmass_hfb_sp();

    virtual ~nucmass_hfb_sp();

    /** \brief Version of \ref nucmass_hfb::entry with spin and parity

        \note This cannot be a child of nucmass_hfb::entry in order
        for the HDF I/O preprocessor macros, like HOFFSET, to work

        \note The fields after \ref Pth (from \ref gamma to \ref
        I3) are only filled in by tables which provide the extra
        deformation, separation-energy, odd-even staggering, and
        moment of inertia information, e.g. the Brussels-Skyrme-
        on-a-Grid (BSkG) tables read by \ref
        o2scl_hdf::bskg_load() . For tables which do not provide
        this information, e.g. HFB17 through HFB27, these fields
        are unused and are left at zero. Similarly, \ref def_wig
        is not provided by the BSkG tables and is left at zero
        for entries obtained from \ref o2scl_hdf::bskg_load() .
    */
    struct entry {

      /// Neutron number
      int N;

      /// Proton number
      int Z;

      /// Atomic number
      int A;

      /// Beta 2 deformation
      double bet2;

      /// Beta 4 deformation
      double bet4;

      /// RMS charge radius
      double Rch;

      /// Deformation and Wigner energies
      double def_wig;

      /// Neutron separation energy
      double Sn;

      /// Proton separation energy
      double Sp;

      /// Beta-decay energy
      double Qbet;

      /// Calculated mass excess
      double Mcal;

      /// Error between experimental and calculated mass excess
      double Err;

      /// Experimental spin
      double Jexp;

      /// Theoretical spin
      double Jth;

      /// Experimental parity
      int Pexp;

      /// Theoretical parity
      int Pth;

      /** \brief Triaxial deformation angle gamma, in degrees
          (BSkG tables only)
      */
      double gamma;

      /// Axial quadrupole deformation, beta_20 (BSkG tables only)
      double beta20;

      /// Non-axial quadrupole deformation, beta_22 (BSkG tables only)
      double beta22;

      /// Axial octupole deformation, beta_30 (BSkG tables only)
      double beta30;

      /// Non-axial octupole deformation, beta_32 (BSkG tables only)
      double beta32;

      /// Two-neutron separation energy (BSkG tables only)
      double S2n;

      /// Two-proton separation energy (BSkG tables only)
      double S2p;

      /// Three-point neutron odd-even mass staggering (BSkG tables only)
      double delta3n;

      /// Three-point proton odd-even mass staggering (BSkG tables only)
      double delta3p;

      /// Five-point neutron odd-even mass staggering (BSkG tables only)
      double delta5n;

      /// Five-point proton odd-even mass staggering (BSkG tables only)
      double delta5p;

      /** \brief Fourth radial moment of the charge density, to the
          one-fourth power, i.e. <r_c^4>^(1/4) (BSkG tables only)
      */
      double rc4;

      /// Moment of inertia about the first axis (BSkG tables only)
      double I1;

      /// Moment of inertia about the second axis (BSkG tables only)
      double I2;

      /// Moment of inertia about the third axis (BSkG tables only)
      double I3;

      /** \brief Rotational correction energy (BSkG1, BSkG2, and
          BSkG4 tables only)
      */
      double Erot;

      /// Average neutron pairing gap (BSkG1, BSkG2, and BSkG4 tables only)
      double avgap_n;

      /// Average proton pairing gap (BSkG1, BSkG2, and BSkG4 tables only)
      double avgap_p;

      /** \brief Experimental RMS charge radius (BSkG1, BSkG2, and
          BSkG4 tables only)
      */
      double rc_exp;

      /** \brief Error between \ref Rch and \ref rc_exp (BSkG1,
          BSkG2, and BSkG4 tables only)
      */
      double rc_err;

      /** \brief Moment of inertia (BSkG1, BSkG2, and BSkG4 tables
          only, distinct from \ref I1, \ref I2, and \ref I3 which
          are only given for BSkG3)
      */
      double MOI;

      /// Parity of the proton subsystem (BSkG1, BSkG2, and BSkG4 tables only)
      int par_p;

      /// Parity of the neutron subsystem (BSkG1, BSkG2, and BSkG4 tables only)
      int par_n;

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
    nucmass_hfb_sp::entry get_ZN(int l_Z, int l_N);
    
    /// Return the type, \c "nucmass_hfb".
    virtual const char *type() const { return "nucmass_hfb_sp"; }

    /** \brief Set data
        
        This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass,
                 std::vector<nucmass_hfb_sp::entry> &m, std::string ref);

    /// Desc
    virtual void clear() {
      mass.clear();
      this->n=0;
      last=0;
      return;
    }
    
  protected:
    
    /// The array containing the mass data of length ame::n
    std::vector<nucmass_hfb_sp::entry> mass;
    
    /// The last table index for caching
    int last;
    
  };
  
}

#endif
