/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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

  /** \brief Find constant values which match a search term
   */
  class find_constants {

  protected:

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

    /// Type for constant database (also used for list of matches)
    typedef struct find_constants_list_s {
      /// List of names for the constant, with the preferred name first
      std::vector<std::string> names;
      /// Unit
      std::string unit;
      /// Flag (either 0, o2scl_mks, or o2scl_cgs)
      int unit_flag;
      /// Value
      double val;
    } find_constants_list;

    /// Database of constant values
    std::vector<find_constants_list> list;
  
  public:
  
    find_constants();

    /** \brief Search for constants matching \c name with unit
	\c unit (possibly empty) and store matches in \c indexes
    */
    int find_nothrow(std::string name, std::string unit,
		     std::vector<find_constants_list> &matches);
  
    /** \brief Search for constants matching \c name with unit \c unit
	and output result(s) with precision \c prec
    */
    void find_print(std::string name, std::string unit="", size_t prec=6);
  
    /** \brief Find a unique match and return the numerical value
     */
    double find_unique(std::string name, std::string unit="");
  
  };

}

#endif
