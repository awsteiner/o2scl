/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef HFB_MASS_H
#define HFB_MASS_H

#include <cmath>

#include <o2scl/nucleus.h>
#include <o2scl/nuclear_mass.h>
#include <o2scl/constants.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Entry structure for HFB mass formula
   */
  typedef struct {
    
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
    
  } hfb_mass_entry;

  /** \brief Version of \ref hfb_mass_entry with spin and parity

      \note This cannot be a child of hfb_mass_entry in order
      for the HDF I/O preprocessor macros, like HOFFSET, to work
  */
  typedef struct {
    
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

  } hfb_sp_mass_entry;
  
  /** \brief HFB Mass formula 

      \todo Mg40 is present in some tables but not others. Compare
      hfb14-plain with hfb14-plain_v0. This may be related to the
      fact that the mass excess of Mg40 differs significantly between
      the 2003 and 2013 Audi et al. tables?

      \todo Update to include hfb17. 
   */
  class hfb_mass : public nuclear_mass_table {
    
  public:
    
    /** \brief Create a new mass formula object 
     */
    hfb_mass();

    virtual ~hfb_mass();

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
    hfb_mass_entry get_ZN(int l_Z, int l_N);
    
    /// Verify that the constructor properly loaded the table
    bool is_loaded() { return (n>0); }
    
    /// The value which corresponds to a blank entry
    double blank() { return 1.0e99; };

    /// Return the type, \c "hfb_mass".
    virtual const char *type() { return "hfb_mass"; }

    /** \brief Set data
	
	This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass, hfb_mass_entry *m, std::string ref);

    /// Return number of entries
    int get_nentries() { return n; }
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The number of entries (about 3000).
    int n;
    
    /// The reference for the original data
    std::string reference;
    
    /// The array containing the mass data of length ame::n
    hfb_mass_entry *mass;
    
    /// The last table index for caching
    int last;
    
#endif
    
  };

  /** \brief HFB Mass formula with spin and parity information
   */
  class hfb_sp_mass : public nuclear_mass_table {
    
  public:
    
    /** \brief Create a new mass formula object
     */
    hfb_sp_mass();

    virtual ~hfb_sp_mass();

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
    hfb_sp_mass_entry get_ZN(int l_Z, int l_N);
    
    /// Return the type, \c "hfb_mass".
    virtual const char *type() { return "hfb_sp_mass"; }

    /** \brief Set data
	
	This function is used by the HDF I/O routines.
    */
    int set_data(int n_mass, hfb_sp_mass_entry *m, std::string ref);

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array containing the mass data of length ame::n
    hfb_sp_mass_entry *mass;

    /// The number of entries (about 3000).
    int n;
    
    /// The reference for the original data
    std::string reference;
    
    /// The last table index for caching
    int last;
    
#endif
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
