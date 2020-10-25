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
#ifndef O2SCL_AUTO_TABLE_H
#define O2SCL_AUTO_TABLE_H
/** \file auto_table.h
    \brief Desc
*/

#include <iostream>
#include <string>
#include <vector>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl_auto_table {
#endif

  /** \brief Lazy table formatting
   */
  class auto_table {
    
  protected:

    /// 
    std::vector<std::string> lines;

    ///
    size_t row_max;

    ///
    std::vector<int> aligns;

    ///
    bool inside_table;
    
  public:

    auto_table();
    
    /** \brief Desc
     */
    void add_string(std::string s, bool endl=false);

    /** \brief Desc
     */
    void done();

  };

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, double d);

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, float f);

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, int i);

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, char ch);

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, size_t s);

  /** \brief Desc
   */
  auto_table &operator<<(auto_table &c, std::string s);

  /** \brief Desc
   */
  //auto_table &operator<<(auto_table &c, o2scl_auto_table::auto_table &(*f)(o2scl_auto_table::auto_table&));

  /** \brief Desc
   */
  //auto_table &endl(auto_table &c);

  static const char endo='\n';
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

