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

/** \file nucmass_wlw.h
    \brief File defining \ref o2scl::nucmass_wlw
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nuclear structure from Wang et al.

      \todo Unfinished.

      Models
      - "WS3.2" in file "wlw10.o2" from \ref Wang10 
      - "WS3.3" in file "wllw10.o2" from \ref Wang10b
      - "WS3.6" in file "lwdw11.o2" from \ref Liu11
      - "WS3_RBF" in file "wl11.o2" from \ref Wang11
      - "WS4_RBF" in file "wlwm14.o2" from \ref Wang14 
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
      /*
	/// Deformation
	double Beta2;
	/// Deformation
	double Beta4;
	/// Deformation
	double Beta6;
	/// Shell correction energy
	double Esh;
	double Dres;
	/// Experimental binding energy (MeV)
	double Eexp;
	/// Theoretical binding energy (MeV)
	double Eth;
	/// Experimental mass excess (MeV)
	double Mexp;
       */
      /// Theoretical mass excess (MeV)
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
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
