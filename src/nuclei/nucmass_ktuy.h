/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_KTUY_MASS_H
#define O2SCL_KTUY_MASS_H

/** \file nucmass_ktuy.h
    \brief File defining \ref o2scl::nucmass_ktuy
*/

#include <o2scl/nucmass.h>

namespace o2scl {

  /** \brief KTUY Mass formula 
   */
  class nucmass_ktuy : public nucmass_table {
    
  public:
    
    nucmass_ktuy();
    
    virtual ~nucmass_ktuy();

    /** \brief Load masses using the specified model number
     */
    int load(std::string model="05", bool external=false);
    
    /** \brief Mass formula entry structure for KTUY mass formula
        
        Nuclear masses from Koura et al. (2000) and Koura et al.
        (2005) as originally specified in the files
        <tt>KTUY04_m246.dat</tt> and <tt>KTUY05_m246.dat</tt> obtained
        from http://wwwndc.jaea.go.jp/nucldata/mass/KTUY04_E.html
        
        \verbatim embed:rst
        See [Koura00]_ and [Koura05]_.
        \endverbatim
    */
    struct entry {
    
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
    nucmass_ktuy::entry get_ZN(int l_Z, int l_N);
    
    /// Verify that the constructor properly loaded the table
    bool is_loaded() { return (n>0); }
    
    /// Return the type, \c "nucmass_ktuy".
    virtual const char *type() { return "nucmass_ktuy"; }

    /// Return number of entries
    virtual size_t get_nentries() { return n; }
    
    /// Desc
    virtual void clear() {
      mass.clear();
      this->n=0;
      last=0;
      return;
    }
    
  protected:
    
    /// The array containing the mass data of length ame::n
    std::vector<entry> mass;
    
    /// The last table index for caching
    int last;
    
  };
  
}

#endif
