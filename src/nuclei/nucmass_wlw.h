/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2014-2024, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCMASS_WLW_H
#define O2SCL_NUCMASS_WLW_H

/** \file nucmass_wlw.h
    \brief File defining \ref o2scl::nucmass_wlw
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

namespace o2scl {

  /** \brief Nuclear structure from Wang et al.

      \verbatim embed:rst
      .. todo:: 

         Class nucmass_wlw is unfinished.

      Models
      - "WS3.2" in file "wlw10.o2" from [Wang10]_
      - "WS3.3" in file "wllw10.o2" from [Wang10b]_
      - "WS3.6" in file "lwdw11.o2" from [Liu11]_
      - "WS3_RBF" in file "wl11.o2" from [Wang11]_
      - "WS4_RBF" in file "wlwm14.o2" from [Wang14]_ 
      \endverbatim

  */
  class nucmass_wlw : public nucmass_table {
    
  public:
    
    nucmass_wlw();

    ~nucmass_wlw();

    int load(std::string model="", bool external=false);
    
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

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
  protected:

    /// The array containing the mass data of length n
    entry *mass;

    /// The last table index for caching
    int last;
    
  };
  
}

#endif
