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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

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

      \note The attach() function stores a pointer to the output
      stream, so the user must take care to make sure this pointer is
      valid.

      \todo Allow user-specified table alignments
      \todo switch to columnify::add_spaces() and add more complicated
      table line specifications
      \todo Support multiprecision types.
      \todo Implement row_max

      \future Create a replacement for std::flush
      \future Finish automatic table detection
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

    /// The output precision for floating-point numbers
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

    /// Pointer to the output stream
    std::ostream *outs;
    
    // Ensure these operators are friends
    friend auto_format &o2scl_auto_format::operator<<(auto_format &at,
						      double d);
    template<class data_t>
    friend auto_format &operator<<(auto_format &at,
				   const boost::numeric::ublas::matrix<data_t> &vu);
    template<class data_t>
    friend auto_format &operator<<(auto_format &at,
				   const std::vector<std::vector<data_t> > &vv);
    
  public:

    auto_format();

    /// If true, align the output of matrices (default true)
    bool align_matrices;
    
    /** \brief Desc
     */
    void attach(std::ostream &out);
    
    /** \brief Desc
     */
    void unattach();
    
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

  /// End the current line
  static const char endo='\n';
  
  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, double d);

  /** \brief Output a single-precision number
   */
  auto_format &operator<<(auto_format &at, float f);

  /** \brief Output an integer
   */
  auto_format &operator<<(auto_format &at, int i);

  /** \brief Output a character
   */
  auto_format &operator<<(auto_format &at, char ch);

  /** \brief Output a \c size_t
   */
  auto_format &operator<<(auto_format &at, size_t s);

  /** \brief Output a string
   */
  auto_format &operator<<(auto_format &at, std::string s);
  
  /** \brief Output a vector of doubles
   */
  auto_format &operator<<(auto_format &at, const std::vector<double> &vd);
  
  /** \brief Output a vector of ints
   */
  auto_format &operator<<(auto_format &at, const std::vector<int> &vd);
  
  /** \brief Output a vector of size_ts
   */
  auto_format &operator<<(auto_format &at, const std::vector<size_t> &vd);
  
  /** \brief Output a vector of chars
   */
  auto_format &operator<<(auto_format &at, const std::vector<char> &vd);
  
  /** \brief Output a vector of std::strings
   */
  auto_format &operator<<(auto_format &at,
			  const std::vector<std::string> &vd);

  /** \brief Output a ublas vector
   */
  template<class data_t>
  auto_format &operator<<(auto_format &at,
			  const boost::numeric::ublas::vector<data_t> &vu) {
    for(size_t i=0;i<vu.size();i++) {
      at << vu[i];
    }
    return at;
  }

  /** \brief Output a ublas matrix

      If \ref auto_format::align_matrices is true, then 
      the output is organized into a table.
   */
  template<class data_t>
  auto_format &operator<<(auto_format &at,
			  const boost::numeric::ublas::matrix<data_t> &vu) {

    bool table_started=false;
    if (at.align_matrices && !at.inside_table) {
      at.done();
      at.start_table();
      table_started=true;
    }
    for(size_t i=0;i<vu.size1();i++) {
      for(size_t j=0;j<vu.size2();j++) {
	at << vu(i,j);
      }
      at << endo;
    }
    if (table_started) {
      at.end_table();
    }
    return at;
  }

  /** \brief Output a vector of vectors

      If \ref auto_format::align_matrices is true and all of
      the vectors in the list have the same length, then 
      the output is organized into a table.
   */
  template<class data_t>
  auto_format &operator<<(auto_format &at,
			  const std::vector<std::vector<data_t> > &vv) {
    bool table_started=false;
    if (at.align_matrices && !at.inside_table && vv.size()>0) {
      size_t nc=vv[0].size();
      bool cols_match=true;
      for(size_t i=0;i<vv.size();i++) {
	if (vv[i].size()!=nc) cols_match=false;
      }
      if (cols_match) {
	at.done();
	at.start_table();
	table_started=true;
      }
    }
    for(size_t i=0;i<vv.size();i++) {
      for(size_t j=0;j<vv[i].size();j++) {
	at << vv[i][j];
      }
      at << endo;
    }
    if (table_started) {
      at.end_table();
    }
    return at;
  }

  
  /** \brief Desc
   */
  //auto_format &operator<<(auto_format &c, o2scl_auto_format::auto_format &(*f)(o2scl_auto_format::auto_format&));

  /** \brief Desc
   */
  //auto_format &endl(auto_format &c);

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

