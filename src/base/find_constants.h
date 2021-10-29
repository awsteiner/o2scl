/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_FIND_CONSTANTS_H
#define O2SCL_FIND_CONSTANTS_H
#include <iostream>
#include <vector>

#include <o2scl/constants.h>

namespace o2scl {

  /** \brief A searchable database of constants with units
   */
  class find_constants {

  public:

    /// Type for constant database (also used for list of matches)
    typedef struct const_entry_s {
      /// List of names for the constant, with the preferred name first
      std::vector<std::string> names;
      /// Unit
      std::string unit;
      /// Flag (currently in the range 0 to 4)
      int unit_flag;
      /// Value
      double val;
      /// Source or reference for value
      std::string source;
      /// Power of length
      int m;
      /// Power of mass
      int k;
      /// Power of time
      int s;
      /// Power of temperature
      int K;
      /// Power of current
      int A;
      /// Power of moles
      int mol;
      /// Power of luminous intensity
      int cd;
    } const_entry;

    /// \name Return values for find_nothrow()
    //@{
    static const int one_exact_match_unit_match=0;
    static const int one_exact_match_unit_diff=1;
    static const int exact_matches_no_unit=2;
    static const int exact_matches_unit_match=3;
    static const int exact_matches_unit_diff=4;
    static const int one_pattern_match_unit_match=5;
    static const int one_pattern_match_unit_diff=6;
    static const int pattern_matches_no_unit=7;
    static const int pattern_matches_unit_match=8;
    static const int pattern_matches_unit_diff=9;
    static const int no_matches=10;
    //@}

  protected:

    /// \name List of constants and unit match function [protected]
    //@{
    /// Database of constant values
    std::vector<const_entry> list;

    /** \brief The function which decides if the requested unit matches
        the specified list entry

        Units match if 
        - the unit is unspecified (string of length zero) or 
          "none" and the flag is \ref fc_none
        - the unit is equal to "any" (case-insensitive comparison)
        - the unit is equal to the list unit (case-insensitive comparison)
        - the unit is "mks" (case-insensitive comparison) and the unit
          flag is either o2scl_mks or fc_none
        - the unit is "cgs" (case-insensitive comparison) and the unit
          flag is either o2scl_cgs or fc_none
    */
    bool unit_match_logic(std::string unit,
                          const const_entry &f) const;
    //@}
    
  public:
  
    find_constants();

    // FYI, from constants.h, we have:
    //
    // static const size_t o2scl_mks=1;
    // static const size_t o2scl_cgs=2;
    
    /// \name Other possible values of the unit flag
    //@{
    static const int fc_unknown=0;
    static const int fc_none=3;
    static const int fc_other=4;
    //@}

    /// \name Basic usage
    //@{
    /** \brief Search for constants matching \c name with unit
	\c unit (possibly empty) and store matches in \c indexes
    */
    int find_nothrow(std::string name, std::string unit,
		     std::vector<const_entry> &matches,
		     int verbose=0) const;
  
    /** \brief Search for constants matching \c name with unit \c unit
	and output result(s) with precision \c prec
    */
    void find_print(std::string name, std::string unit="", size_t prec=6,
		    int verbose=0) const;
  
    /** \brief Find a unique match and return the numerical value
     */
    double find_unique(std::string name, std::string unit="") const;
    //@}

    /// \name Functions to show or modify the constant list
    //@{
    /** \brief Output the full list of constants to \c os 
    */
    void output_list(std::ostream &os) const;

    /** \brief Output the full list of constants to \c os 
    */
    void output_list_full(std::ostream &os) const;

    /** \brief Output the full list of constants to 
        \c std::cout
    */
    void output_list_cout() const {
      output_list(std::cout);
      return;
    }

    /** \brief Output one entry from the constant database
        to \c os
    */
    void output(const find_constants::const_entry &c,
                std::ostream &os) const;
    
    /** \brief Add a constant
     */
    void add_constant(const const_entry &f, int verbose=0);
    
    /** \brief Remove a constant
     */
    void del_constant(std::string &name, int verbose=0);
    //@}
    
  };

}

#endif
