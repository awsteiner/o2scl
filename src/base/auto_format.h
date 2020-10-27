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
#ifndef O2SCL_AUTO_FORMAT_H
#define O2SCL_AUTO_FORMAT_H
/** \file auto_format.h
    \brief Desc
*/

#include <iostream>
#include <string>
#include <vector>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl_auto_format {
#endif

  // Declarations for friendship
  class auto_format;
  auto_format &operator<<(auto_format &at, double d);

  /** \brief Automatically format output

      This class is a wrapper around output streams which performs
      automatic spacing and table formatting. Only scientific
      formatting for floating-point numbers is supported at present.

      \note Experimental.

      \note This class caches each line before sending to cout,
      so issuing <tt>cout << std::flush</tt> in the middle of 
      a line will not output the buffer to the screen.

      \todo Support char *'s, c-style arrays, matrix types, 
      ublas objects, multiprecision types.
      \todo Implement row_max
      \todo Allow output to files ("fout") in addition to cout.

      \future For automatic table detection: allow user to change the
      number of rows which must have the same number of 'words' to
      verify a table.
      \future Make internal algorithm more efficient.
   */
  class auto_format {
    
  protected:

    /** \brief If true, automatic formatting is enabled (default true)
     */
    bool enabled;
    
    /// \name Standard buffer
    //@{
    /// Output line buffer
    std::vector<std::string> lines;

    // Index of line for next output
    //size_t next_line;
    //@}
    size_t precision_;
    
    /// \name Table mode
    //@{
    /// If true, try to automatically detect tables (default true)
    bool auto_tables;
    
    /// The number of table header rows
    size_t n_headers;
    
    /// Headers for table mode
    std::vector<std::vector<std::string> > headers;

    /// If true, we are currently inside a table
    bool inside_table;

    /// Columns for table mode
    std::vector<std::vector<std::string> > columns;

    /// Index of next column
    size_t next_column;

    /// Maximum number of table rows (default 1000)
    size_t row_max;

    /// Alignment specifications for table columns
    std::vector<int> aligns;
    //@}

    // Ensure this operator is a friend to access precision_
    friend auto_format &o2scl_auto_format::operator<<(auto_format &at,
						      double d);
    
  public:

    auto_format();
    
    /** \brief Add a string to the output buffer
     */
    void add_string(std::string s);

    /** \brief Desc
     */
    void precision(size_t p);
    
    /** \brief Disable formatting and send all output 
	directly to \c cout
     */
    void off();
    
    /** \brief Turn on automatic formatting (on by default)
     */
    void on();
    
    /** \brief Add an endline
     */
    void endline();
    
    /** \brief Flush all buffered output to the screen
     */
    void done();

    /** \brief Start a table
     */
    void start_table();

    /** \brief Debug the table
     */
    void debug_table();

    /** \brief End a table
     */
    void end_table();

    /// Verbosity parameter (default 0)
    int verbose;

    /// Parameter for table line output
    int table_lines;

  };

  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &c, double d);

  /** \brief Output a single-precision number
   */
  auto_format &operator<<(auto_format &c, float f);

  /** \brief Output an integer
   */
  auto_format &operator<<(auto_format &c, int i);

  /** \brief Output a character
   */
  auto_format &operator<<(auto_format &c, char ch);

  /** \brief Output a \c size_t
   */
  auto_format &operator<<(auto_format &c, size_t s);

  /** \brief Output a string
   */
  auto_format &operator<<(auto_format &c, std::string s);
  
  /** \brief Output a vector of doubles
   */
  auto_format &operator<<(auto_format &c, const std::vector<double> &vd);
  
  /** \brief Output a vector of ints
   */
  auto_format &operator<<(auto_format &c, const std::vector<int> &vd);
  
  /** \brief Output a vector of size_ts
   */
  auto_format &operator<<(auto_format &c, const std::vector<size_t> &vd);
  
  /** \brief Output a vector of chars
   */
  auto_format &operator<<(auto_format &c, const std::vector<char> &vd);
  
  /** \brief Output a vector of std::strings
   */
  auto_format &operator<<(auto_format &c,
			  const std::vector<std::string> &vd);
  
  /** \brief Desc
   */
  //auto_format &operator<<(auto_format &c, o2scl_auto_format::auto_format &(*f)(o2scl_auto_format::auto_format&));

  /** \brief Desc
   */
  //auto_format &endl(auto_format &c);

  /// End the current line
  static const char endo='\n';
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

