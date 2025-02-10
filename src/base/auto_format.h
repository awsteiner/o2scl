/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2020-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_AUTO_FORMAT_H
#define O2SCL_AUTO_FORMAT_H
/** \file auto_format.h
    \brief File describing the class \ref auto_format
*/

#include <iostream>
#include <string>
#include <vector>

#include <o2scl/err_hnd.h>
#include <o2scl/columnify.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#ifdef O2SCL_MULTIP
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

namespace o2scl_auto_format {

#ifdef O2SCL_MULTIP
  typedef
  boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
  cpp_dec_float_35;
  typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
#endif
  
  // Declarations for friendship
  class auto_format;
  auto_format &operator<<(auto_format &at, double d);
  auto_format &operator<<(auto_format &at, long double d);
#ifdef O2SCL_MULTIP
  auto_format &operator<<(auto_format &at, const cpp_dec_float_35 &d);
  auto_format &operator<<(auto_format &at, const cpp_dec_float_50 &d);
  auto_format &operator<<(auto_format &at, const cpp_dec_float_100 &d);
#endif

  /** \brief Automatically format output

      This class is a wrapper around output streams which performs
      automatic spacing and table formatting. Only scientific
      formatting for floating-point numbers is supported at present.

      Experimental.

      This class caches each line before sending to cout, so issuing
      ``cout << std::flush`` in the middle of a line will not
      output the buffer to the screen.

      The attach() function stores a pointer to the output
      stream, so the user must take care to make sure this pointer is
      valid.

      For now, this class is in its own namespace,
      <tt>o2scl_auto_format</tt>.

      \verbatim embed:rst
      .. note::

      .. todo:: 

         In class auto_format:

         - Allow user-specified table alignments

         Future:

         - Create a replacement for std::flush
         - Finish automatic table detection
         - For automatic table detection: allow user to change the
           number of rows which must have the same number of 'words' to
           verify a table.
         - Make internal algorithm more efficient.

      \endverbatim
   */
  class auto_format {
    
  protected:

    /// \name Other variables (protected)
    //@{
    /** \brief If true, automatic formatting is enabled (default true)
     */
    bool enabled;
    
    /// Output line buffer
    std::vector<std::string> lines;

    /// The output precision for floating-point numbers
    size_t precision_;
    //@}
    
    /// \name Table mode (protected)
    //@{
    /// If true, try to automatically detect tables (default true)
    bool auto_tables;
    
    /// The number of table header rows (default 0)
    size_t n_headers;
    
    /// Pointer to the output stream
    std::ostream *outs;
    
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

    // Ensure these operators are friends
    /// \name Friend operators
    //@{
    friend auto_format &o2scl_auto_format::operator<<
    (auto_format &at, double d);
     
    friend auto_format &o2scl_auto_format::operator<<
    (auto_format &at, long double d);
    
#if defined (O2SCL_MULTIP) || defined (DOXYGEN)
    friend auto_format &o2scl_auto_format::operator<<
    (auto_format &at, const cpp_dec_float_35 &d);
     
    friend auto_format &o2scl_auto_format::operator<<
    (auto_format &at, const cpp_dec_float_50 &d);
     
    friend auto_format &o2scl_auto_format::operator<<
    (auto_format &at, const cpp_dec_float_100 &d);
#endif
    
    template<class data_t>
    friend auto_format &operator<<
    (auto_format &at,
     const boost::numeric::ublas::matrix<data_t> &vu);
    
    template<class data_t>
    friend auto_format &operator<<
    (auto_format &at,
     const std::vector<std::vector<data_t> > &vv);
    //@}
    
  public:

    /// \name Basic usage
    //@{
    /// Constructor
    auto_format();

    /** \brief Attach an ostream, so that all future output goes
        there
     */
    void attach(std::ostream &out);
    
    /** \brief Clear the buffer and put all future output to 
        \c std::cout 
     */
    void unattach();
    
    /** \brief Disable formatting and send all output 
        directly to \c cout
     */
    void off();
    
    /** \brief Turn on automatic formatting (on by default)
     */
    void on();
    
    /** \brief Start a table
     */
    void start_table();

    /** \brief End a table
     */
    void end_table();

    /** \brief Flush all buffered output to the screen
     */
    void done();
    //@}

    /// \name Other functions
    //@{
    /** \brief Add an endline
     */
    void endline();
    
    /** \brief Debug the table by outputting detailed table information
        to the attached stream.

        This function is useful just before end_table() is called.
     */
    void debug_table();

    /** \brief Add a string to the output buffer
     */
    void add_string(std::string s);
    //@}

    /// \name Settings
    //@{
    /// The columnify object for table formatting
    o2scl::columnify col;
    
    /// If true, align the output of matrices (default true)
    bool align_matrices;
    
    /// Verbosity parameter (default 0)
    int verbose;

    /** \brief Parameter for table line output

        This parameter is used for \ref columnify::table_lines.
     */
    int table_lines;

    /** \brief Set the precision of floating-point values
     */
    void precision(size_t p);
    //@}

  };

  /// End the current line
  static const char endo='\n';
  
  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, double d);

  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, long double d);

#if defined (O2SCL_MULTIP) || defined (DOXYGEN)
  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, const cpp_dec_float_25 &d);

  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, const cpp_dec_float_35 &d);

  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, const cpp_dec_float_50 &d);

  /** \brief Output a double-precision number
   */
  auto_format &operator<<(auto_format &at, const cpp_dec_float_100 &d);
#endif
  
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

      If \ref o2scl_auto_format::auto_format::align_matrices is true, then 
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

      If \ref o2scl_auto_format::auto_format::align_matrices is true and all of
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

  // Attempts to create an 'endl' similar to that for cout
  
  /* \brief Desc
   */
  //auto_format &operator<<(auto_format &c,
  //o2scl_auto_format::auto_format &(*f)(o2scl_auto_format::auto_format&));

  /* \brief Desc
   */
  //auto_format &endl(auto_format &c);

}

#endif

