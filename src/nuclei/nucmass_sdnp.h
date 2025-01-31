/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2014-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCMASS_SDNP_H
#define O2SCL_NUCMASS_SDNP_H

/** \file nucmass_sdnp.h
    \brief File defining \ref o2scl::nucmass_sdnp
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

namespace o2scl {

  /** \brief Nuclear structure from Stoitsov et al.
      
      \verbatim embed:rst
      .. todo:: 

         Class nucmass_sdnp is unfinished.

      \endverbatim
      
      Models
      - "sdnp03" - from Skyrme model SkM*
      - "sd_skp_04" - from Skyrme model SkP
      - "sd_sly4_04" - from Skyrme model SLy4

      \verbatim embed:rst
      See [Stoitsov03]_ and [Dobaczewski04]_ and
      http://www.fuw.edu.pl/~dobaczew/thodri/thodri.html .
      \endverbatim
  */
  class nucmass_sdnp : public nucmass_table {
    
  public:
    
    nucmass_sdnp();

    ~nucmass_sdnp();

    int load(std::string model="", bool external=false);
    
    /** \brief Entry structure
	
    */
    class entry {

    public:

      /// Proton number
      int Z;
      /// Neutron number
      int N;
      /// HFB energy minimum (MeV) 
      double ENERGY;
      /*
	///
	double ACCURACY;
	///
	double S_2N;
	///
	double S_2P;
	///
	double LAM_N;
	///
	double LAM_P;
	///
	double LA2_N;
	///
	double LA2_P;
	///
	double DEL_N;
	///
	double DEL_P;
	///
	double RAD_N;
	///
	double RAD_P;
	///
	double BETA;
	///
	double Q20_N;
	///
	double Q20_P;
       */

    };
  
    /// Return the type, \c "nucmass_sdnp".
    virtual const char *type() { return "nucmass_sdnp"; }

    /// Returns true if data has been loaded
    bool is_loaded() { return (n>0); }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Return number of entries
    virtual size_t get_nentries() { return n; }
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
    /// Desc
    void clear() {
      mass.clear();
      this->n=0;
      last=0;
      return;
    }
    
  protected:

    /// The number of entries (about 3000).
    size_t n;
    
    /// The array containing the mass data of length n
    std::vector<entry> mass;

    /// The last table index for caching
    int last;
    
  };
  
}

#endif
