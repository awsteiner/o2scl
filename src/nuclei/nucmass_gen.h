/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014-2020, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCMASS_GEN_H
#define O2SCL_NUCMASS_GEN_H

/** \file nucmass_gen.h
    \brief File defining \ref o2scl::nucmass_gen
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nuclear properties 
   */
  class nucmass_gen : public nucmass_table {
    
  public:
    
    nucmass_gen();

    ~nucmass_gen();

    /** \brief Load a file with binding energies
     */
    int load_be(std::string fname, std::string be_col,
		double be_units, bool external=false);
    
    /// Return the type, \c "nucmass_gen".
    virtual const char *type() { return "nucmass_gen"; }

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
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double get_string(int Z, int N, std::string column);
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The \ref o2scl::table object containing the data
    o2scl::table<> data;

    /// Column which refers to the mass excess
    size_t mex_col_ix;
    
    /// The last table index for caching
    int last;
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
