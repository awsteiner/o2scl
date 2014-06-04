/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCMASS_WLW_H
#define O2SCL_NUCMASS_WLW_H

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief Nuclear structure from Wang et al.
  */
  class nucmass_wlw : public nucmass_table {
    
  public:
    
    nucmass_wlw(std::string model="", bool external=false);

    ~nucmass_wlw();
    
    /** \brief Entry structure
	
    */
    class entry {

    public:

      /// Proton number
      int Z;
      /// Neutron number
      int N;
      /// Desc
      double Mth;
    };
  
    /// Return the type, \c "nucmass_wlw".
    virtual const char *type() { return "nucmass_wlw"; }

    /// Returns true if data has been loaded
    bool is_loaded() { return (n>0); }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Return number of entries
    int get_nentries() { return n; }
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The number of entries (about 3000).
    int n;
    
    /// The reference for the original data
    std::string reference;
    
    /// The array containing the mass data of length n
    entry *mass;

    /// The last table index for caching
    int last;
    
#endif

  };
  
#ifndef DOXYGENP
}
#endif

#endif
