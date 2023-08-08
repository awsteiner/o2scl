/*
  ───────────────────────────────────────────────────────────────────
  
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
    virtual const char *type() { return "nucmass_hfb"; }

    /** \brief Set data
        
        This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass, nucmass_hfb::entry *m, std::string ref);

  protected:
    
    /// The array containing the mass data of length ame::n
    nucmass_hfb::entry *mass;
    
    /// The last table index for caching
    int last;
    
  };

  /** \brief HFB Mass formula with spin and parity information
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
    virtual const char *type() { return "nucmass_hfb_sp"; }

    /** \brief Set data
        
        This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass, nucmass_hfb_sp::entry *m, std::string ref);

  protected:
    
    /// The array containing the mass data of length ame::n
    nucmass_hfb_sp::entry *mass;

    /// The last table index for caching
    int last;
    
  };
  
}

#endif
