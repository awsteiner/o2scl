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
#ifndef O2SCL_KTUY_MASS_H
#define O2SCL_KTUY_MASS_H

#include <o2scl/nuclear_mass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Mass formula entry structure for KTUY mass formula

      Nuclear masses from \ref Koura00 and \ref Koura05
      as originally specified in the files <tt>KTUY04_m246.dat</tt>
      and <tt>KTUY05_m246.dat</tt> obtained from
      http://wwwndc.jaea.go.jp/nucldata/mass/KTUY04_E.html
   */
  typedef struct {
    
    /// Neutron number
    int N;
    
    /// Proton number
    int Z;
    
    /// Atomic number
    int A;
    
    /// Calculated mass excess
    double Mcal;

    /// Shell energy
    double Esh;

    /// Alpha 2 deformation
    double alpha2;

    /// Alpha 4 deformation
    double alpha4;

    /// Alpha 6 deformation
    double alpha6;

  } ktuy_mass_entry;

  /** \brief KTUY Mass formula 
   */
  class ktuy_mass : public nuclear_mass_table {
    
  public:
    
    /** \brief Create a new mass formula object using the specified model
	number
    */
    ktuy_mass(std::string model="05", bool external=false);

    virtual ~ktuy_mass();

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
    ktuy_mass_entry get_ZN(int l_Z, int l_N);
    
    /// Verify that the constructor properly loaded the table
    bool is_loaded() { return (n>0); }
    
    /// Return the type, \c "ktuy_mass".
    virtual const char *type() { return "ktuy_mass"; }

    /// Return number of entries
    int get_nentries() { return n; }
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The number of entries (about 3000).
    int n;
    
    /// The reference for the original data
    std::string reference;
    
    /// The array containing the mass data of length ame::n
    ktuy_mass_entry *mass;
    
    /// The last table index for caching
    int last;
    
#endif
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
